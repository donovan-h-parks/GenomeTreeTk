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
from collections import defaultdict

from math import floor

import dendropy

from biolib.common import is_float
from biolib.taxonomy import Taxonomy


class PhylogeneticDiversity():
    """Calculate phylogenetic diversity."""

    def __init__(self):
        """Initialize."""
        self.logger = logging.getLogger()
        
    def _read_taxa_list(self, taxa_list):
        """Read taxa from file."""
        
        taxa = set()
        for line in open(taxa_list):
            taxa.add(line.strip().split('\t')[0])
            
        return taxa
        
    def _total_pd(self, tree):
        """Calculate PD over entire tree."""
        
        total_pd = 0
        total_taxa = 0
        for node in tree.preorder_node_iter():
            if node.parent_node is not None:
                total_pd += node.edge.length
                if node.is_leaf():
                    total_taxa += 1
                    
        return total_pd, total_taxa

    def _taxon_pd(self, tree, ingroup, outgroup):
        """Calculate phylogenetic gain of each ingroup taxon relative to outgroup."""

        pg_taxon = {}
        for leaf in tree.leaf_node_iter():
            if leaf.taxon.label in ingroup:
                # find first internal node containing an outgroup taxon
                pg = 0
                parent = leaf
                outgroup_taxon = 'None'
                while parent:
                    # check for outgroup taxon
                    stop = False
                    for tip in parent.leaf_iter():
                        if tip.taxon.label in outgroup:
                            outgroup_taxon = tip.taxon.label
                            if outgroup_taxon in ingroup:
                                outgroup_taxon += ' (one or more outgroup taxa are assigned to this ingroup taxon)'
                            stop = True
                            
                    if stop:
                        break
                    
                    pg += parent.edge.length
                    parent = parent.parent_node
                    
                pg_taxon[leaf.taxon.label] = [pg, outgroup_taxon]

        return pg_taxon
        
    def _taxa_pd(self, tree, taxa):
        """Calculate phylogenetic diversity of a set of taxa."""
        
        pd = 0
        visited_edges = set()
        for cur_node in taxa:
            while cur_node.parent_node is not None:
                if cur_node.edge in visited_edges:
                    break

                pd += cur_node.edge.length
                visited_edges.add(cur_node.edge)
                cur_node = cur_node.parent_node

        return pd

    def pd(self, tree, taxa_list, per_taxa_pg_file):
        """Calculate phylogenetic diversity of extant taxa."""

        self.logger.info('Reading tree.')
        tree = dendropy.Tree.get_from_path(tree,
                                            schema='newick',
                                            rooting='force-rooted',
                                            preserve_underscores=True)
                                                                            
        # get total branch length of tree
        self.logger.info('Calculating total PD.')
        total_pd, total_taxa = self._total_pd(tree)
        self.logger.info('Total PD for %d taxa = %.2f' % (total_taxa, total_pd))

        # get PD of ingroup taxa
        self.logger.info('Calculating total PD for specified taxa.')
        ingroup = self._read_taxa_list(taxa_list)
        
        in_taxa = set()
        for leaf in tree.leaf_node_iter():
            if leaf.taxon.label in ingroup:
                in_taxa.add(leaf)
        
        self.logger.info('Specified ingroup taxa: %d' % len(ingroup))
        self.logger.info('Ingroup taxa found in tree: %d' % len(in_taxa))

        # calculate PD for ingroup
        in_pd = self._taxa_pd(tree, in_taxa)

        # get PD of outgroup taxa
        outgroup = set()
        out_taxa = set()
        for leaf in tree.leaf_node_iter():
            if leaf.taxon.label not in ingroup:
                out_taxa.add(leaf)
                outgroup.add(leaf.taxon.label)
                
        self.logger.info('Outgroup taxa found in tree: %d' % len(out_taxa))
        out_pd = self._taxa_pd(tree, out_taxa)

        # calculate PG of each ingroup taxon relative to outgroup in requested
        if per_taxa_pg_file:
            self.logger.info('Calculating PG of each ingroup taxon relative to outgroup.')
            pg_taxon = self._taxon_pd(tree, ingroup, outgroup)
            
            fout = open(per_taxa_pg_file, 'w')
            fout.write('Taxon\tPG\tPercent PG\tFirst outgroup taxon\n')
            for taxon, pg_stats in pg_taxon.iteritems():
                pg, outgroup_taxon = pg_stats
                fout.write('%s\t%f\t%f\t%s\n' % (taxon, pg, pg * 100.0 / total_pd, outgroup_taxon))
            fout.close()
                   
        return total_pd, len(in_taxa), in_pd, len(out_taxa), out_pd
        
    def _clade_pd(self, tree, ingroup, outgroup):
        """Calculate PD for named clades."""
        
        pd = {}
        for node in tree.preorder_node_iter():
            if not node.label:
                continue

            taxon = None
            if ':' in node.label:
                _support, taxon = node.label.split(':')
            else:
                if not is_float(node.label):
                    taxon = node.label
                    
            if taxon:
                taxon_pd = 0
                taxon_count = 0
                in_taxon_pd = 0
                in_taxon_count = 0
                out_taxon_pd = 0
                out_taxon_count = 0
                for nn in node.postorder_iter():
                    if nn == node:
                        continue
                        
                    # check if group contains taxa from
                    # the ingroup and/or outgroup
                    ingroup_leaves = False
                    outgroup_leaves = False
                    for leaf in nn.leaf_iter():
                        genome_id = leaf.taxon.label
                        if genome_id in ingroup:
                            ingroup_leaves = True
                        
                        if genome_id in outgroup:
                            outgroup_leaves = True
                            
                    if ingroup_leaves:
                        in_taxon_pd += nn.edge.length

                    if outgroup_leaves:
                        out_taxon_pd += nn.edge.length
                        
                    if nn.is_leaf():
                        genome_id = nn.taxon.label
                        
                        if genome_id in ingroup:
                            in_taxon_count += 1
                        else:
                            out_taxon_count += 1
                            
                    taxon_pd += nn.edge.length

                pd[taxon] = [taxon_pd, in_taxon_pd, in_taxon_count, out_taxon_pd, out_taxon_count]

        return pd
        
    def pd_clade(self, decorated_tree, taxa_list, output_file):
        """Calculate phylogenetic diversity of named groups."""
        
        # calculate PD for entire tree
        self.logger.info('Reading tree.')
        tree = dendropy.Tree.get_from_path(decorated_tree,
                                            schema='newick',
                                            rooting='force-rooted',
                                            preserve_underscores=True)
 
        self.logger.info('Calculating total PD.')
        total_pd, total_taxa = self._total_pd(tree)
        
        # get ingroup and outgroup taxa
        ingroup = self._read_taxa_list(taxa_list)
        outgroup = set()
        in_tree = set()
        for leaf in tree.leaf_node_iter():
            if leaf.taxon.label not in ingroup:
                outgroup.add(leaf.taxon.label)
            if leaf.taxon.label in ingroup:
                in_tree.add(leaf.taxon.label)
                
        self.logger.info('Specified %d ingroup taxa.' % len(ingroup))
        self.logger.info('Identified %d ingroup taxa in tree.' % len(in_tree))
        self.logger.info('Identified %d outgroup taxa in tree.' % len(outgroup))

        # PD for named groups
        self.logger.info('Calculating PD for named clades.')
        pd_clade = self._clade_pd(tree, ingroup, outgroup)
        self.logger.info(' ... identified %d named clades.' % len(pd_clade))

        # report results
        fout = open(output_file, 'w')
        fout.write('Clade')
        fout.write('\tTaxa\tPD\tPercent PD')
        fout.write('\tOut Taxa\tOut PD\tOut Percent PD')
        fout.write('\tIn Taxa\tIn PD\tIn Percent PD')
        fout.write('\tIn PG\tIn Percent PG\n')
        
        ordered_taxa = Taxonomy().sort_taxa(pd_clade.keys())
        for taxon in ordered_taxa:
            taxon_pd, in_taxon_pd, in_taxon_count, out_taxon_pd, out_taxon_count = pd_clade[taxon]
            taxon_count = in_taxon_count + out_taxon_count
            in_taxon_pg = taxon_pd - out_taxon_pd
            
            taxon_pd = max(taxon_pd, 1e-9) # make sure PD is never exactly zero to avoid division errors
            
            row = taxon
            row += '\t%d\t%.2f\t%.2f' % (taxon_count, taxon_pd, taxon_pd * 100 / total_pd)
            row += '\t%d\t%.2f\t%.2f' % (out_taxon_count, out_taxon_pd, out_taxon_pd * 100 / taxon_pd)
            row += '\t%d\t%.2f\t%.2f' % (in_taxon_count, in_taxon_pd, in_taxon_pd * 100 / taxon_pd)
            row += '\t%.2f\t%.2f' % (in_taxon_pg, in_taxon_pg * 100 / taxon_pd)
            fout.write(row + '\n')
        fout.close()
