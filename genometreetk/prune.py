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

import sys
import logging

import dendropy

from biolib.newick import parse_label


class Prune(object):
    """Prune tree."""

    def __init__(self):
        """Initialization."""
        
        self.logger = logging.getLogger('timestamp')

    def run(self, input_tree, taxa_to_retain_file, output_tree):
        """Prune tree.
        
        Parameters
        ----------
        input_tree : str
          File containing Newick tree to rerooted.
        taxa_to_retain_file : str
          File specifying taxa to retain.
        output_tree : str
          Name of file for rerooted tree.
        """
        
        self.logger.info('Reading input tree.')
        tree = dendropy.Tree.get_from_path(input_tree, 
                                            schema='newick', 
                                            rooting='force-rooted', 
                                            preserve_underscores=True)
                                            
        # read taxa to retain
        taxa_to_retain = set()
        for line in open(taxa_to_retain_file):
            if line[0] == '#' or not line.strip():
                continue
                
            line_split = line.strip().split('\t')
            taxa_to_retain.add(line_split[0])
  
        # find taxa to retain
        self.logger.info('Identifying taxa to retain.')
        taxa_in_tree = set()
        for node in tree.postorder_node_iter():
            if node.is_leaf():
                support, taxon, _auxiliary_info = parse_label(node.taxon.label)
                if taxon in taxa_to_retain:
                    taxa_in_tree.add(node.taxon)
                    taxa_to_retain.remove(taxon)
            else:
                support, taxon, _auxiliary_info = parse_label(node.label)
                if taxon in taxa_to_retain:
                    for leaf in node.leaf_iter():
                        taxa_in_tree.add(leaf.taxon)
                    taxa_to_retain.remove(taxon)
                        
            # check if all outgroup taxa have been identified          
            if not taxa_to_retain:
                break
                
        self.logger.info('Identified {:,} extant taxa to retain in tree.'.format(
                            len(taxa_in_tree)))
        
        if taxa_to_retain:
            self.logger.warning('Failed to identify {:,} taxa: {}'.format(
                                    len(taxa_to_retain), 
                                    ','.join(taxa_to_retain)))
                
        # prune tree
        self.logger.info('Pruning tree.')
        tree.retain_taxa(taxa_in_tree)
        
        # write out results
        self.logger.info('Writing output tree.')  
        tree.write_to_path(output_tree, 
                            schema='newick', 
                            suppress_rooting=True, 
                            unquoted_underscores=True)
          