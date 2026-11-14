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

"""Run FastTree on replicate alignments."""

import logging
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from gtdblib.util.shell.execute import check_dependencies
from tqdm import tqdm


class FastTreeError(Exception):
    pass


def fasttree_cmd(msa_file, log_file, base_type='prot', model='wag', gamma=False, multithreaded=False):
    """Create FastTree command line.

    Nucleotide trees are inferred under the GTR model. Protein trees are
    inferred under the WAG, LG, or JTT (FastTree default) model.
    """

    if base_type not in ('nt', 'prot'):
        raise ValueError(f'Unknown base type: {base_type}')

    cmd = ['FastTreeMP' if multithreaded else 'FastTree', '-quiet', '-nosupport']
    if gamma:
        cmd.append('-gamma')

    if base_type == 'nt':
        cmd += ['-nt', '-gtr']
    elif model == 'wag':
        cmd.append('-wag')
    elif model == 'lg':
        cmd.append('-lg')
    elif model != 'jtt':
        raise ValueError(f'Unknown model: {model}')

    cmd += ['-log', log_file, msa_file]

    return cmd


def is_complete_tree(tree_file):
    """Check if file contains a complete Newick tree."""

    if not os.path.exists(tree_file) or os.path.getsize(tree_file) == 0:
        return False

    with open(tree_file, 'rb') as f:
        f.seek(max(0, os.path.getsize(tree_file) - 64))
        return f.read().rstrip().endswith(b';')


def run_fasttree(msa_file, output_tree, base_type='prot', model='wag', gamma=False, multithreaded=False):
    """Infer a tree with FastTree.

    The tree is written to a temporary file and renamed once FastTree has
    finished successfully, so an interrupted run never leaves a file that looks
    complete. FastTree's log and stderr are written next to the output tree as
    <prefix>.fasttree.log and <prefix>.fasttree.stderr.

    Returns
    -------
    bool
      True if the tree was inferred, False if a complete tree already existed.
    """

    if is_complete_tree(output_tree):
        return False

    out_dir = os.path.dirname(output_tree)
    prefix = Path(output_tree).stem
    log_file = os.path.join(out_dir, prefix + '.fasttree.log')
    stderr_file = os.path.join(out_dir, prefix + '.fasttree.stderr')
    tmp_tree = output_tree + '.tmp'

    cmd = fasttree_cmd(msa_file, log_file, base_type, model, gamma, multithreaded)
    with open(tmp_tree, 'wb') as fout, open(stderr_file, 'wb') as ferr:
        proc = subprocess.run(cmd, stdout=fout, stderr=ferr)

    if proc.returncode != 0 or not is_complete_tree(tmp_tree):
        os.remove(tmp_tree)
        with open(stderr_file, errors='replace') as f:
            stderr_tail = ''.join(f.readlines()[-5:]).strip()
        raise FastTreeError('FastTree failed on {} (exit code {}): {}'.format(
            msa_file, proc.returncode, stderr_tail))

    os.replace(tmp_tree, output_tree)

    return True


def infer_trees(msa_files, tree_files, cpus=1, base_type='prot', model='wag', gamma=False, multithreaded=False):
    """Infer a tree for each alignment in parallel.

    Trees which already exist and are complete are not inferred again.

    Returns
    -------
    int
      Number of trees which already existed.

    Raises
    ------
    FastTreeError
      If FastTree fails for any alignment.
    """

    logger = logging.getLogger()
    check_dependencies(['FastTreeMP' if multithreaded else 'FastTree'])

    if multithreaded and cpus > 1 and 'OMP_NUM_THREADS' not in os.environ:
        logger.warning('Running {} FastTreeMP processes without OMP_NUM_THREADS set; '
                       'each process may use all available cores.'.format(cpus))

    # threads are sufficient as the work is done by FastTree subprocesses
    failures = []
    num_skipped = 0
    with ThreadPoolExecutor(max_workers=cpus) as executor:
        futures = [executor.submit(run_fasttree, m, t, base_type, model, gamma, multithreaded)
                   for m, t in zip(msa_files, tree_files)]
        for future in tqdm(as_completed(futures), total=len(futures), unit=' trees'):
            try:
                if not future.result():
                    num_skipped += 1
            except FastTreeError as e:
                failures.append(str(e))

    if num_skipped:
        logger.info(f'Skipped {num_skipped:,} trees which were already inferred.')

    if failures:
        for msg in failures:
            logger.error(msg)
        raise FastTreeError(f'Failed to infer {len(failures):,} of {len(msa_files):,} trees.')

    return num_skipped
