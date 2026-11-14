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


import logging
import multiprocessing as mp
import os
import random
import sys
from math import floor
from pathlib import Path

from gtdblib.tree.bootstrap_merge import bootstrap_merge_replicates
from gtdblib.util.bio import seq_io
from gtdblib.util.shell.gtdbshutil import make_sure_path_exists
from tqdm import tqdm

from genometreetk.fasttree import FastTreeError, infer_trees

# State shared with forked worker processes (see JackknifeMarkers._create_replicates).
_SHARED = {}


def _write_replicate_task(task):
    """Pool entry point for writing a jackknifed alignment."""

    output_msa = task
    JackknifeMarkers.jackknife_alignment(_SHARED['msa'],
                                         _SHARED['perc_markers_to_keep'],
                                         _SHARED['marker_lengths'],
                                         output_msa)
    return output_msa


class JackknifeMarkers(object):
    """Assess robustness by jackkifing genes in alignment."""

    def __init__(self, cpus):
        """Initialization.

        Parameters
        ----------
        cpus : int
          Number of cpus to use.
        """

        self.logger = logging.getLogger()

        self.cpus = cpus

    @staticmethod
    def jackknife_alignment(msa, perc_markers_to_keep, marker_lengths, output_file):
        """Jackknife alignment to a subset of marker genes.

        The marker_lengths must be specified in the order
        in which genes were concatenated.

        Parameters
        ----------
        msa : d[seq_id] -> seq
          Full multiple sequence alignment.
        perc_markers_to_keep : float
          Percentage of marker genes to keep in each replicate [0, 1].
        marker_lengths : list
          Length of each marker gene.
        output_file : str
          File to write bootstrapped alignment.
        """

        markers_to_keep = random.sample(list(range(0, len(marker_lengths))),
                                        int(floor(perc_markers_to_keep * len(marker_lengths))))

        start_pos = [0]
        for index, ml in enumerate(marker_lengths):
            start_pos.append(start_pos[index] + ml)

        # alignment columns of the kept markers, in their original order
        segments = sorted((start_pos[i], start_pos[i] + marker_lengths[i]) for i in markers_to_keep)

        tmp_file = output_file + '.tmp'
        with open(tmp_file, 'w') as fout:
            for seq_id, seq in msa.items():
                fout.write('>' + seq_id + '\n')
                fout.write(''.join(seq[start:end] for start, end in segments) + '\n')
        os.replace(tmp_file, output_file)

    def _create_replicates(self, msa_files):
        """Create jackknifed alignments in parallel."""

        _SHARED.update(msa=self.msa,
                       perc_markers_to_keep=self.perc_markers_to_keep,
                       marker_lengths=self.marker_lengths)
        try:
            if self.cpus > 1 and len(msa_files) > 1:
                # fork is used so workers share the MSA in memory
                ctx = mp.get_context('fork')
                with ctx.Pool(min(self.cpus, len(msa_files))) as pool:
                    for _ in tqdm(pool.imap_unordered(_write_replicate_task, msa_files),
                                  total=len(msa_files), unit=' replicates'):
                        pass
            else:
                for msa_file in tqdm(msa_files, unit=' replicates'):
                    _write_replicate_task(msa_file)
        finally:
            _SHARED.clear()

    def run(self, input_tree,
            msa_file,
            marker_info_file,
            mask_file,
            perc_markers_to_keep,
            num_replicates,
            model,
            jk_dir,
            output_dir):
        """Jackknife marker genes.

        Marker file should have the format:
          <marker id>\t<marker name>\t<marker desc>\t<length>\n

        Parameters
        ----------
        input_tree : str
          Tree inferred with all data.
        msa_file : str
          File containing multiple sequence alignment for all taxa.
        marker_info_file : str
          File indicating database id, HMM name, description and length of each marker in the alignment.
        mask_file : str
          File indicating masking of multiple sequence alignment.
        perc_markers_to_keep : float [0, 1]
          Percentage of marker genes to keep in each replicate.
        num_replicates : int
          Number of replicates to perform.
        model : str
          Desired model of evolution.
        jk_dir : str
          Directory with pre-computed jackknife trees.
        output_dir : str
          Output directory for jackkife trees.
        """

        assert(model in ['wag', 'jtt'])

        self.model = model
        self.perc_markers_to_keep = perc_markers_to_keep

        rep_tree_files = []
        if not jk_dir:
            self.replicate_dir = os.path.join(output_dir, 'replicates')
            make_sure_path_exists(self.replicate_dir)

            # determine length of each marker gene in alignment
            marker_lengths = []
            total_len = 0
            with open(marker_info_file) as f:
                f.readline()
                for line in f:
                    line_split = line.split('\t')
                    ml = int(line_split[3])
                    marker_lengths.append(ml)
                    total_len += ml

            self.logger.info('Concatenated length of markers: %d' % total_len)

            # read mask
            with open(mask_file) as f:
                mask = f.readline().strip()
            start = 0
            self.marker_lengths = []
            total_mask_len = 0
            for ml in marker_lengths:
                end = start + ml
                zeros = mask[start:end].count('0')
                start = end

                self.marker_lengths.append(ml - zeros)
                total_mask_len += ml - zeros

            self.logger.info('Concatenated length of filtered MSA: %d' % total_mask_len)

            # read full multiple sequence alignment
            self.msa = seq_io.read(msa_file)

            if len(list(self.msa.values())[0]) != total_mask_len:
                self.logger.error('Length of MSA does not meet length of mask.')
                sys.exit(1)

            # create jackknifed alignments and infer trees
            msa_files = [os.path.join(self.replicate_dir, f'jk_markers.msa.{i}.faa') for i in range(num_replicates)]
            rep_tree_files = [os.path.join(self.replicate_dir, f'jk_markers.tree.{i}.tre') for i in range(num_replicates)]

            self.logger.info('Calculating jackknife marker replicates:')
            self._create_replicates(msa_files)
            self.msa = None

            self.logger.info(f'Inferring {num_replicates:,} jackknife trees using {self.cpus:,} CPUs:')
            try:
                infer_trees(msa_files, rep_tree_files, cpus=self.cpus, base_type='prot', model=self.model)
            except FastTreeError as e:
                self.logger.error(str(e))
                sys.exit(1)
        else:
            for f in sorted(os.listdir(jk_dir)):
                if f.endswith('.tree') or f.endswith('.tre'):
                    rep_tree_files.append(os.path.join(jk_dir, f))

            self.logger.info('Read %d jackknife replicates.' % len(rep_tree_files))
            if not rep_tree_files:
                self.logger.error(f'No jackknife trees found in: {jk_dir}')
                sys.exit(1)

        # calculate support
        self.logger.info(f'Calculating support from {len(rep_tree_files):,} replicates.')
        output_tree = os.path.join(output_dir, Path(input_tree).stem + '.jk_markers.tree')
        bootstrap_merge_replicates(input_tree, output_tree, rep_tree_files, self.cpus)

        return output_tree
