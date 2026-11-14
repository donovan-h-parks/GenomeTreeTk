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
import logging
import random
from math import floor
from pathlib import Path

from gtdblib.util.bio import seq_io
from gtdblib.util.shell.gtdbshutil import make_sure_path_exists

from genometreetk.fasttree import infer_trees

from genometreetk.tree_support import TreeSupport


class JackknifeTaxa(object):
    """Assess robustness by jackknifing taxa in alignment."""

    def __init__(self, cpus):
        """Initialization.

        Parameters
        ----------
        cpus : int
          Number of cpus to use.
        """

        self.logger = logging.getLogger()

        self.cpus = cpus

    def _create_replicates(self, num_replicates):
        """Create jackknifed alignments and infer a tree for each."""

        msa_files = [os.path.join(self.replicate_dir, f'jk_taxa.msa.{i}.fna') for i in range(num_replicates)]
        tree_files = [os.path.join(self.replicate_dir, f'jk_taxa.tree.{i}.tre') for i in range(num_replicates)]
        for msa_file in msa_files:
            self.jackknife_taxa(self.msa, self.perc_taxa_to_keep, self.outgroup_ids, msa_file)

        infer_trees(msa_files, tree_files, cpus=self.cpus, base_type='prot', model=self.model)

    def jackknife_taxa(self, msa, perc_taxa_to_keep, outgroup_ids, output_file):
        """Jackknife alignment to a subset of taxa.

        Parameters
        ----------
        msa : d[seq_id] -> seq
          Full multiple sequence alignment.
        perc_taxa_to_keep : float
          Percentage of marker genes to keep in each replicate.
        outgroup_ids : set
          Labels of outgroup taxa.
        output_file : str
          File to write bootstrapped alignment.
        """

        # randomly select ingroup taxa
        ingroup_taxa = set(msa.keys()) - outgroup_ids
        taxa_to_keep = random.sample(sorted(ingroup_taxa), int(floor(len(ingroup_taxa) * perc_taxa_to_keep)))

        taxa_to_keep = set(taxa_to_keep).union(outgroup_ids)

        fout = open(output_file, 'w')
        for seq_id, seq in msa.items():
            if seq_id in taxa_to_keep:
                fout.write('>' + seq_id + '\n')
                fout.write(seq + '\n')
        fout.close()

    def run(self, input_tree, msa_file, outgroup_file, perc_taxa_to_keep, num_replicates, model, output_dir):
        """Jackknife taxa.

        Parameters
        ----------
        input_tree : str
          Tree inferred with all data.
        msa_file : str
          File containing multiple sequence alignment for all taxa.
        outgroup_file : str
          File indicating labels of outgroup taxa.
        perc_taxa_to_keep : float
          Percentage of taxa to keep in each replicate.
        num_replicates : int
          Number of replicates to perform.
        model : str
          Desired model of evolution.
        output_dir : str
          input_tree directory for bootstrap trees.
        """

        assert(model in ['wag', 'jtt'])

        self.perc_taxa_to_keep = perc_taxa_to_keep
        self.model = model
        self.replicate_dir = os.path.join(output_dir, 'replicates')
        make_sure_path_exists(self.replicate_dir)

        # read outgroup taxa
        self.outgroup_ids = set()
        if outgroup_file:
            for line in open(outgroup_file):
                self.outgroup_ids.add(line.strip())

        # read full multiple sequence alignment
        self.msa = seq_io.read(msa_file)

        # calculate replicates
        #***self.logger.info('Calculating jackknife taxa replicates:')
        #***self._create_replicates(num_replicates)

        # calculate support
        rep_tree_files = []
        for rep_index in range(num_replicates):
            rep_tree_files.append(os.path.join(self.replicate_dir, 'jk_taxa.tree.' + str(rep_index) + '.tre'))

        tree_support = TreeSupport()
        output_tree = os.path.join(output_dir, Path(input_tree).stem + '.jk_taxa.tree')
        tree_support.subset_taxa(input_tree, rep_tree_files, output_tree)

        return output_tree