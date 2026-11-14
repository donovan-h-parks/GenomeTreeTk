###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

import os
import sys
import gzip
import logging
import multiprocessing as mp
from pathlib import Path

import numpy as np
from tqdm import tqdm

from gtdblib.tree.bootstrap_merge import bootstrap_merge_replicates
from gtdblib.util.shell.gtdbshutil import make_sure_path_exists

from genometreetk.fasttree import FastTreeError, infer_trees, is_complete_tree


# Alignment shared with forked worker processes. It is stored as a module
# level global so workers inherit it through copy-on-write memory instead of
# receiving a pickled copy. Since the numpy buffer is never written to and
# has no per-row Python objects, no pages are duplicated in the workers.
_SHARED = {}

# Number of alignment rows gathered and written at a time by a worker.
ROW_BLOCK_SIZE = 4096

TREE_EXTENSIONS = ('.tree', '.tre', '.treefile')
MSA_EXTENSIONS = ('.fna', '.faa', '.fasta', '.fa')
NT_CHARS = b'ACGTUNacgtun'


class BootstrapError(Exception):
    pass


def read_msa(msa_file):
    """Read MSA into a 2D uint8 matrix.

    Parameters
    ----------
    msa_file : str
      FASTA file with aligned sequences (optionally gzipped).

    Returns
    -------
    list of bytes
      Sequence IDs in file order.
    np.ndarray
      Matrix of shape (num_seqs, alignment_len) holding one byte per base.
    """

    open_file = gzip.open if msa_file.endswith('.gz') else open

    # first pass: count sequences and determine alignment length so the
    # matrix can be allocated once (avoids holding two copies in memory)
    num_seqs = 0
    aln_len = 0
    with open_file(msa_file, 'rb') as f:
        for line in f:
            if line[:1] == b'>':
                num_seqs += 1
            elif num_seqs == 1:
                aln_len += len(line.strip())

    if num_seqs == 0:
        raise BootstrapError(f'No sequences found in {msa_file}.')
    if aln_len == 0:
        raise BootstrapError(f'First sequence in {msa_file} is empty.')

    msa = np.empty((num_seqs, aln_len), dtype=np.uint8)
    seq_ids = []

    def store(row_idx, pieces):
        seq = b''.join(pieces)
        if len(seq) != aln_len:
            raise BootstrapError(
                'Sequences in MSA have different lengths: {} has length {:,} but expected {:,}.'.format(
                    seq_ids[row_idx].decode(), len(seq), aln_len))
        msa[row_idx] = np.frombuffer(seq, dtype=np.uint8)

    # second pass: fill matrix
    row_idx = -1
    pieces = []
    with open_file(msa_file, 'rb') as f:
        for line in f:
            if line[:1] == b'>':
                if row_idx >= 0:
                    store(row_idx, pieces)
                row_idx += 1
                fields = line[1:].split(None, 1)
                if not fields:
                    raise BootstrapError(f'Sequence {row_idx + 1:,} in {msa_file} has an empty header.')
                seq_ids.append(fields[0])
                pieces = []
            else:
                line = line.strip()
                if line:
                    pieces.append(line)
        store(row_idx, pieces)

    if len(set(seq_ids)) != len(seq_ids):
        raise BootstrapError(f'MSA contains duplicate sequence IDs: {msa_file}')

    return seq_ids, msa


def nucleotide_fraction(msa, max_rows=1000):
    """Estimate fraction of non-gap characters that are nucleotides."""

    counts = np.bincount(msa[:max_rows].ravel(), minlength=256)
    nt = sum(counts[c] for c in NT_CHARS)
    non_gap = counts.sum() - counts[ord('-')] - counts[ord('.')]
    return nt / non_gap if non_gap else 0.0


def replicate_seed(base_seed, rep_idx):
    """Seed sequence for a replicate (independent of worker scheduling)."""

    return np.random.SeedSequence(entropy=base_seed, spawn_key=(rep_idx,))


def expected_replicate_size(headers, sample_len):
    """Exact size in bytes of a complete replicate MSA file."""

    return sum(len(h) for h in headers) + len(headers) * (sample_len + 1)


def write_replicate(rep_idx, output_msa):
    """Write a single bootstrap replicate of the shared MSA.

    Columns are sampled with replacement. The file is written to a temporary
    path and atomically renamed, so an interrupted run never leaves a file
    that looks complete.
    """

    msa = _SHARED['msa']
    headers = _SHARED['headers']
    sample_len = _SHARED['sample_len']

    rng = np.random.default_rng(replicate_seed(_SHARED['seed'], rep_idx))
    cols = rng.integers(0, msa.shape[1], size=sample_len)

    num_seqs = msa.shape[0]
    block = np.empty((min(ROW_BLOCK_SIZE, num_seqs), sample_len + 1), dtype=np.uint8)
    block[:, -1] = ord('\n')

    tmp_file = output_msa + '.tmp'
    with open(tmp_file, 'wb', buffering=8 * 1024 * 1024) as fout:
        for start in range(0, num_seqs, ROW_BLOCK_SIZE):
            end = min(start + ROW_BLOCK_SIZE, num_seqs)
            out = block[:end - start]
            np.take(msa[start:end], cols, axis=1, out=out[:, :-1], mode='clip')
            parts = []
            for i in range(end - start):
                parts.append(headers[start + i])
                parts.append(memoryview(out[i]))
            fout.writelines(parts)
    os.replace(tmp_file, output_msa)

    return rep_idx


class Bootstrap(object):
    """Assess robustness of genome tree by bootstrapping multiple sequence alignment."""

    def __init__(self, cpus, multiple_processors=False):
        """Initialization.

        Parameters
        ----------
        cpus : int
          Number of cpus to use.
        multiple_processors : bool
          Use FastTreeMP instead of FastTree.
        """

        self.logger = logging.getLogger()

        self.cpus = cpus
        self.mp = multiple_processors

    def _infer_trees(self, msa_files, tree_files):
        """Infer trees for replicate MSAs in parallel."""

        self.logger.info(f'Inferring {len(msa_files):,} bootstrap trees using {self.cpus:,} CPUs:')
        infer_trees(msa_files, tree_files,
                    cpus=self.cpus,
                    base_type=self.base_type,
                    model=self.model,
                    gamma=self.gamma,
                    multithreaded=self.mp)

    def _create_msa_replicates(self, msa_file, num_replicates, frac, seed):
        """Create bootstrap replicates of MSA in parallel."""

        self.logger.info('Reading MSA:')
        seq_ids, msa = read_msa(msa_file)
        self.logger.info(' - read {:,} sequences with an alignment length of {:,}'.format(
            msa.shape[0], msa.shape[1]))

        if self.base_type == 'prot' and nucleotide_fraction(msa) > 0.9:
            self.logger.warning('MSA appears to contain nucleotides, but bases are being treated '
                                'as amino acids. Use --base_type nt for nucleotide alignments.')

        sample_len = int(msa.shape[1] * frac)
        if sample_len < 1:
            raise BootstrapError(f'Fraction {frac} of alignment results in no columns being sampled.')

        headers = [b'>' + sid + b'\n' for sid in seq_ids]
        expected_size = expected_replicate_size(headers, sample_len)

        msa_files = []
        todo = []
        for rep_idx in range(num_replicates):
            output_msa = os.path.join(self.replicate_dir, f'bootstrap_msa.r_{rep_idx}.fna')
            msa_files.append(output_msa)
            if os.path.exists(output_msa) and os.path.getsize(output_msa) == expected_size:
                continue
            elif os.path.exists(output_msa):
                self.logger.warning(f'Recreating incomplete or mismatched replicate: {output_msa}')
            todo.append((rep_idx, output_msa))

        if len(todo) < num_replicates:
            self.logger.info(f'Skipping {num_replicates - len(todo):,} MSA replicates which already exist.')

        if todo:
            self.logger.info('Creating {:,} MSA replicates of length {:,} using {:,} CPUs (seed = {}):'.format(
                len(todo), sample_len, self.cpus, seed))

            _SHARED.update(msa=msa, headers=headers, sample_len=sample_len, seed=seed)
            try:
                if self.cpus > 1 and len(todo) > 1:
                    # fork is used explicitly so workers share the MSA in memory
                    # (spawn/forkserver would pickle a full copy to each worker)
                    ctx = mp.get_context('fork')
                    with ctx.Pool(min(self.cpus, len(todo))) as pool:
                        for _ in tqdm(pool.imap_unordered(_write_replicate_task, todo),
                                      total=len(todo), unit=' replicates'):
                            pass
                else:
                    for task in tqdm(todo, unit=' replicates'):
                        _write_replicate_task(task)
            finally:
                _SHARED.clear()

        return msa_files

    def run(self,
            input_tree,
            msa_file,
            num_replicates,
            model,
            gamma,
            base_type,
            frac,
            boot_dir,
            msa_only,
            output_dir,
            seed=None):
        """Bootstrap multiple sequence alignment.

        Parameters
        ----------
        input_tree : str
          Tree inferred with all data, or 'none'.
        msa_file : str
          File containing multiple sequence alignment for all taxa.
        num_replicates : int
          Number of replicates to perform.
        model : str
          Desired model of evolution.
        gamma : bool
          Use GAMMA model.
        base_type : str
          Indicates if bases are nucleotides or amino acids.
        frac : float
          Fraction of alignment to subsample.
        boot_dir : str
          Directory with pre-computed bootstrap trees or MSA replicates.
        msa_only : bool
          Only create MSA replicates.
        output_dir : str
          Directory for bootstrap trees.
        seed : int
          Seed for random number generator.

        Returns
        -------
        str
          Tree decorated with bootstrap support, or None.
        """

        assert(model in ['wag', 'lg', 'jtt'])
        assert(base_type in ['nt', 'prot'])

        self.model = model
        self.gamma = gamma
        self.base_type = base_type

        self.replicate_dir = os.path.join(output_dir, 'replicates')
        make_sure_path_exists(self.replicate_dir)

        rep_tree_files = []
        try:
            if not boot_dir:
                if seed is None:
                    seed = np.random.SeedSequence().entropy
                msa_files = self._create_msa_replicates(msa_file, num_replicates, frac, seed)

                if msa_only:
                    self.logger.info(f'MSA replicates written to: {self.replicate_dir}')
                    return None

                rep_tree_files = [os.path.join(self.replicate_dir, f'bootstrap_tree.r_{i}.tree')
                                  for i in range(num_replicates)]
                self._infer_trees(msa_files, rep_tree_files)
            else:
                files = sorted(os.listdir(boot_dir))
                rep_tree_files = [os.path.join(boot_dir, f) for f in files if f.endswith(TREE_EXTENSIONS)]
                msa_files = [os.path.join(boot_dir, f) for f in files if f.endswith(MSA_EXTENSIONS)]

                if rep_tree_files:
                    self.logger.info(f'Read {len(rep_tree_files):,} bootstrap trees from: {boot_dir}')
                    incomplete = [t for t in rep_tree_files if not is_complete_tree(t)]
                    if incomplete:
                        raise BootstrapError('Found {:,} empty or incomplete bootstrap trees, e.g.: {}'.format(
                            len(incomplete), incomplete[0]))
                elif msa_files:
                    self.logger.info(f'Read {len(msa_files):,} MSA replicates from: {boot_dir}')
                    if msa_only:
                        self.logger.info('No trees inferred as --msa_only was specified.')
                        return None

                    rep_tree_files = [os.path.join(self.replicate_dir, Path(m).stem + '.tree')
                                      for m in msa_files]
                    self._infer_trees(msa_files, rep_tree_files)
                else:
                    raise BootstrapError(f'No bootstrap trees or MSA replicates found in: {boot_dir}')
        except (BootstrapError, FastTreeError) as e:
            self.logger.error(str(e))
            sys.exit(1)

        # calculate support values
        output_tree = None
        if input_tree.lower() != 'none' and not msa_only:
            self.logger.info(f'Calculating bootstrap support values from {len(rep_tree_files):,} trees.')
            output_tree = os.path.join(output_dir, Path(input_tree).stem + '.bootstrap.tree')
            bootstrap_merge_replicates(input_tree, output_tree, rep_tree_files, self.cpus)

        return output_tree


def _write_replicate_task(task):
    """Pool entry point for writing a replicate."""

    rep_idx, output_msa = task
    return write_replicate(rep_idx, output_msa)
