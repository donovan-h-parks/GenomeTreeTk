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

import csv
import logging
import sys

import dendropy
from biolib.common import check_file_exists, make_sure_path_exists, canonical_gid
from biolib.external.execute import check_dependencies
from biolib.newick import parse_label
from biolib.taxonomy import Taxonomy

from genometreetk.arb import Arb
from genometreetk.bootstrap import Bootstrap
from genometreetk.combine_support import CombineSupport
from genometreetk.common import read_gtdb_metadata
from genometreetk.derep_tree import DereplicateTree
from genometreetk.jackknife_markers import JackknifeMarkers
from genometreetk.jackknife_taxa import JackknifeTaxa
from genometreetk.phylogenetic_diversity import PhylogeneticDiversity
from genometreetk.arb import Arb
from genometreetk.derep_tree import DereplicateTree
from genometreetk.prune import Prune
from genometreetk.reroot_tree import RerootTree
from genometreetk.rna_workflow import RNA_Workflow

csv.field_size_limit(sys.maxsize)


class OptionsParser():
    def __init__(self):
        """Initialization"""
        self.logger = logging.getLogger()

    def ssu_tree(self, options):
        """Infer 16S tree spanning GTDB genomes."""

        check_dependencies(['mothur', 'ssu-align', 'ssu-mask', 'FastTreeMP', 'blastn'])

        check_file_exists(options.gtdb_metadata_file)
        check_file_exists(options.gtdb_ssu_file)
        make_sure_path_exists(options.output_dir)

        rna_workflow = RNA_Workflow(options.cpus)
        rna_workflow.run('ssu',
                         options.gtdb_metadata_file,
                         options.gtdb_ssu_file,
                         options.min_ssu_length,
                         options.min_scaffold_length,
                         options.min_quality,
                         options.max_contigs,
                         options.min_N50,
                         not options.disable_tax_filter,
                         options.genome_list,
                         options.output_dir,
                         options.align_method)

        self.logger.info('Results written to: %s' % options.output_dir)

    def lsu_tree(self, options):
        """Infer 23S tree spanning GTDB genomes."""

        check_dependencies(['esl-sfetch', 'cmsearch', 'cmalign', 'esl-alimask', 'FastTreeMP', 'blastn'])

        check_file_exists(options.gtdb_metadata_file)
        check_file_exists(options.gtdb_lsu_file)
        make_sure_path_exists(options.output_dir)

        rna_workflow = RNA_Workflow(options.cpus)
        rna_workflow.run('lsu',
                         options.gtdb_metadata_file,
                         options.gtdb_lsu_file,
                         options.min_lsu_length,
                         options.min_scaffold_length,
                         options.min_quality,
                         options.max_contigs,
                         options.min_N50,
                         not options.disable_tax_filter,
                         # options.reps_only,
                         # options.user_genomes,
                         options.genome_list,
                         options.output_dir)

        self.logger.info('Results written to: %s' % options.output_dir)

    def rna_tree(self, options):
        """Infer 16S + 23S tree spanning GTDB genomes."""

        check_dependencies(['FastTreeMP'])

        check_file_exists(options.ssu_msa)
        check_file_exists(options.ssu_tree)
        check_file_exists(options.lsu_msa)
        check_file_exists(options.lsu_tree)
        make_sure_path_exists(options.output_dir)

        rna_workflow = RNA_Workflow(options.cpus)
        rna_workflow.combine(options.ssu_msa,
                             options.ssu_tree,
                             options.lsu_msa,
                             options.lsu_tree,
                             options.output_dir)

        self.logger.info('Results written to: %s' % options.output_dir)

    def rna_dump(self, options):
        """Dump all 5S, 16S, and 23S sequences to files."""

        check_file_exists(options.genomic_file)
        make_sure_path_exists(options.output_dir)

        rna_workflow = RNA_Workflow(1)
        rna_workflow.dump(options.genomic_file,
                          options.gtdb_taxonomy,
                          options.min_5S_len,
                          options.min_16S_ar_len,
                          options.min_16S_bac_len,
                          options.min_23S_len,
                          options.min_contig_len,
                          options.include_user,
                          options.genome_list,
                          options.output_dir)

        self.logger.info('Results written to: %s' % options.output_dir)

    def derep_tree(self, options):
        """Dereplicate tree."""

        check_file_exists(options.input_tree)
        check_file_exists(options.gtdb_metadata)
        check_file_exists(options.msa_file)
        make_sure_path_exists(options.output_dir)

        derep_tree = DereplicateTree()
        derep_tree.run(options.input_tree,
                       options.lineage_of_interest,
                       options.outgroup,
                       options.gtdb_metadata,
                       options.taxa_to_retain,
                       options.msa_file,
                       options.keep_unclassified,
                       options.output_dir)

    def bootstrap(self, options):
        """Bootstrap multiple sequence alignment."""

        check_dependencies(['FastTree'])

        if options.input_tree.lower() != 'none':
            check_file_exists(options.input_tree)
        if options.msa_file.lower() != 'none':
            check_file_exists(options.msa_file)
        make_sure_path_exists(options.output_dir)

        bootstrap = Bootstrap(options.cpus)
        output_tree = bootstrap.run(options.input_tree,
                                    options.msa_file,
                                    options.num_replicates,
                                    options.model,
                                    options.gamma,
                                    options.base_type,
                                    options.fraction,
                                    options.boot_dir,
                                    options.msa_only,
                                    options.output_dir)

        self.logger.info('Bootstrapped tree written to: %s' % output_tree)

    def jk_markers(self, options):
        """Jackknife marker genes."""

        check_file_exists(options.input_tree)
        if options.msa_file != 'NONE':
            check_file_exists(options.msa_file)
        make_sure_path_exists(options.output_dir)

        jackknife_markers = JackknifeMarkers(options.cpus)
        output_tree = jackknife_markers.run(options.input_tree,
                                            options.msa_file,
                                            options.marker_info_file,
                                            options.mask_file,
                                            options.perc_markers,
                                            options.num_replicates,
                                            options.model,
                                            options.jk_dir,
                                            options.output_dir)

        self.logger.info('Jackknifed marker tree written to: %s' % output_tree)

    def jk_taxa(self, options):
        """Jackknife taxa."""

        check_file_exists(options.input_tree)
        check_file_exists(options.msa_file)
        make_sure_path_exists(options.output_dir)

        jackknife_taxa = JackknifeTaxa(options.cpus)
        output_tree = jackknife_taxa.run(options.input_tree,
                                         options.msa_file,
                                         options.outgroup_ids,
                                         options.perc_taxa,
                                         options.num_replicates,
                                         options.model,
                                         options.output_dir)

        self.logger.info('Jackknifed taxa tree written to: %s' % output_tree)

    def combine(self, options):
        """Combine support values into a single tree."""

        combineSupport = CombineSupport()
        combineSupport.run(options.support_type,
                           options.bootstrap_tree,
                           options.jk_marker_tree,
                           options.jk_taxa_tree,
                           options.output_tree)

    def support_wf(self, options):
        """"Perform entire tree support workflow."""

        self.bootstrap(options)
        self.jk_markers(options)
        self.jk_taxa(options)
        self.combine(options)

    def midpoint(self, options):
        """"Midpoint root tree."""

        reroot = RerootTree()
        reroot.midpoint(options.input_tree, options.output_tree)

    def outgroup(self, options):
        """Reroot tree with outgroup."""

        check_file_exists(options.taxonomy_file)

        self.logger.info('Identifying genomes from the specified outgroup.')
        outgroup = set()
        for genome_id, taxa in Taxonomy().read(options.taxonomy_file).items():
            if options.outgroup_taxon in taxa:
                outgroup.add(genome_id)
        self.logger.info('Identifying %d genomes in the outgroup.' % len(outgroup))

        reroot = RerootTree()
        reroot.root_with_outgroup(options.input_tree, options.output_tree, outgroup)

    def fill_ranks(self, options):
        """Ensure taxonomy strings contain all 7 canonical ranks."""

        check_file_exists(options.input_taxonomy)

        fout = open(options.output_taxonomy, 'w')
        taxonomy = Taxonomy()
        t = taxonomy.read(options.input_taxonomy)

        for genome_id, taxon_list in t.items():
            full_taxon_list = taxonomy.fill_missing_ranks(taxon_list)

            taxonomy_str = ';'.join(full_taxon_list)
            if not taxonomy.check_full(taxonomy_str):
                sys.exit(-1)

            fout.write('%s\t%s\n' % (genome_id, taxonomy_str))

        fout.close()

        self.logger.info('Revised taxonomy written to: %s' % options.output_taxonomy)

    def propagate(self, options):
        """Propagate labels to all genomes in a cluster."""

        check_file_exists(options.input_taxonomy)
        check_file_exists(options.metadata_file)

        user_to_uba = {}
        if options.uba_mapping_file:
            self.logger.info('Parsing genome ID mapping file.')
            with open(options.uba_mapping_file) as f:
                for line in f:
                    tokens = line.strip().split('\t')
                    if len(tokens) == 2:
                        user_to_uba[tokens[0]] = tokens[1]
            self.logger.info(' - found mappings for {:,} genomes.'.format(len(user_to_uba)))

        # get representative genome information
        rep_metadata = read_gtdb_metadata(options.metadata_file, ['gtdb_representative',
                                                                  'gtdb_clustered_genomes'])

        rep_metadata = {canonical_gid(gid): values
                        for gid, values in rep_metadata.items()}

        rep_metadata = {user_to_uba.get(gid, gid): values
                        for gid, values in rep_metadata.items()}

        explict_tax = Taxonomy().read(options.input_taxonomy)

        self.logger.info(f' - identified {len(rep_metadata):,} genomes')

        # sanity check all representatives have a taxonomy string
        rep_count = 0
        for gid in rep_metadata:
            is_rep_genome, clustered_genomes = rep_metadata.get(gid, (None, None))
            if is_rep_genome:
                rep_count += 1
                if gid not in explict_tax:
                    self.logger.error(
                        'Expected to find {} in input taxonomy as it is a GTDB representative.'.format(gid))
                    sys.exit(-1)

        self.logger.info('Identified {:,} representatives in metadata file and {:,} genomes in input taxonomy file.'.format(
            rep_count,
            len(explict_tax)))

        # propagate taxonomy to genomes clustered with each representative
        fout = open(options.output_taxonomy, 'w')
        for rid, taxon_list in explict_tax.items():
            taxonomy_str = ';'.join(taxon_list)
            rid = canonical_gid(rid)
            rid = user_to_uba.get(rid, rid)

            is_rep_genome, clustered_genomes = rep_metadata[rid]
            if is_rep_genome:
                # assign taxonomy to representative and all genomes in the cluster
                fout.write('{}\t{}\n'.format(rid, taxonomy_str))
                for cid in [gid.strip() for gid in clustered_genomes.split(';')]:
                    cid = canonical_gid(cid)
                    cid = user_to_uba.get(cid, cid)
                    if cid != rid:
                        if cid in rep_metadata:
                            fout.write('{}\t{}\n'.format(cid, taxonomy_str))
                        else:
                            self.logger.warning('Skipping {} as it is not in GTDB metadata file.'.format(cid))
            else:
                self.logger.error(
                    'Did not expected to find {} in input taxonomy as it is not a GTDB representative.'.format(rid))
                sys.exit(-1)

        self.logger.info('Taxonomy written to: {}'.format(options.output_taxonomy))

    def strip(self, options):
        """Remove taxonomic labels from tree."""

        check_file_exists(options.input_tree)

        outgroup_in_tree = set()
        tree = dendropy.Tree.get_from_path(options.input_tree,
                                           schema='newick',
                                           rooting='force-rooted',
                                           preserve_underscores=True)

        for node in tree.internal_nodes():
            if node.label:
                if ':' in node.label:
                    support, _taxa = node.label.split(':')
                    node.label = support
                else:
                    try:
                        # if number if a float (or int) treat
                        # it as a support value
                        f = float(node.label)
                    except ValueError:
                        node.label = None

        tree.write_to_path(options.output_tree,
                           schema='newick',
                           suppress_rooting=True,
                           unquoted_underscores=True)

        self.logger.info('Stripped tree written to: %s' % options.output_tree)

    def rm_support(self, options):
        """Remove support values from tree."""

        check_file_exists(options.input_tree)

        outgroup_in_tree = set()
        tree = dendropy.Tree.get_from_path(options.input_tree,
                                           schema='newick',
                                           rooting='force-rooted',
                                           preserve_underscores=True)

        for node in tree.internal_nodes():
            if node.label:
                if ':' in node.label:
                    support, taxa = node.label.split(':')
                    node.label = taxa
                else:
                    try:
                        # if number if a float (or int) treat
                        # it as a support value
                        f = float(node.label)
                        node.label = None
                    except ValueError:
                        pass  # keep other labels

        tree.write_to_path(options.output_tree,
                           schema='newick',
                           suppress_rooting=True,
                           unquoted_underscores=True)

        self.logger.info('Stripped tree written to: %s' % options.output_tree)

    def pull(self, options):
        """Create taxonomy file from a decorated tree."""

        check_file_exists(options.input_tree)

        if options.no_validation:
            tree = dendropy.Tree.get_from_path(options.input_tree,
                                               schema='newick',
                                               rooting="force-rooted",
                                               preserve_underscores=True)

            taxonomy = {}
            for leaf in tree.leaf_node_iter():
                taxon_id = leaf.taxon.label

                node = leaf.parent_node
                taxa = []
                while node:
                    support, taxon, aux_info = parse_label(node.label)
                    if taxon:
                        for t in list(map(str.strip, taxon.split(';')))[::-1]:
                            taxa.append(t)
                    node = node.parent_node

                taxonomy[taxon_id] = taxa[::-1]
        else:
            taxonomy = Taxonomy().read_from_tree(options.input_tree)

        Taxonomy().write(taxonomy, options.output_taxonomy)

        self.logger.info('Stripped tree written to: %s' % options.output_taxonomy)

    def append(self, options):
        """Append command"""

        check_file_exists(options.input_tree)
        check_file_exists(options.input_taxonomy)

        taxonomy = Taxonomy().read(options.input_taxonomy)

        tree = dendropy.Tree.get_from_path(options.input_tree,
                                           schema='newick',
                                           rooting='force-rooted',
                                           preserve_underscores=True)

        for n in tree.leaf_node_iter():
            taxa_str = taxonomy.get(n.taxon.label, None)
            if taxa_str == None:
                self.logger.error('Taxonomy file does not contain an entry for %s.' % n.label)
                sys.exit(-1)
            n.taxon.label = n.taxon.label + '|' + '; '.join(taxonomy[n.taxon.label])

        tree.write_to_path(options.output_tree,
                           schema='newick',
                           suppress_rooting=True,
                           unquoted_underscores=True)

        self.logger.info('Decorated tree written to: %s' % options.output_tree)

    def prune(self, options):
        """Prune tree."""

        check_file_exists(options.input_tree)
        check_file_exists(options.taxa_to_retain)

        prune = Prune()
        prune.run(options.input_tree,
                  options.taxa_to_retain,
                  options.output_tree)

    def phylogenetic_diversity(self, options):
        """Calculate phylogenetic diversity of extant taxa."""

        check_file_exists(options.tree)
        check_file_exists(options.taxa_list)

        pd = PhylogeneticDiversity()
        rtn = pd.pd(options.tree, options.taxa_list, options.per_taxa_pg_file)
        total_pd, num_in_taxa, in_pd, num_out_taxa, out_pd = rtn
        total_taxa = num_in_taxa + num_out_taxa
        in_pg = total_pd - out_pd

        # report phylogenetic diversity (PD) and gain (PG)
        print('')
        print('\tNo. Taxa\tPD\tPercent PD')
        print('%s\t%d\t%.2f\t%.2f%%' % ('Full tree', total_taxa, total_pd, 100))

        print('%s\t%d\t%.2f\t%.3f%%' % ('Outgroup taxa (PD)',
                                        num_out_taxa,
                                        out_pd,
                                        out_pd * 100 / total_pd))

        print('%s\t%d\t%.2f\t%.3f%%' % ('Ingroup taxa (PD)',
                                        num_in_taxa,
                                        in_pd,
                                        (in_pd) * 100 / total_pd))

        print('%s\t%d\t%.2f\t%.3f%%' % ('Ingroup taxa (PG)',
                                        num_in_taxa,
                                        in_pg,
                                        in_pg * 100 / total_pd))

    def phylogenetic_diversity_clade(self, options):
        """Calculate phylogenetic diversity of named groups."""

        check_file_exists(options.decorated_tree)

        pd = PhylogeneticDiversity()
        pd.pd_clade(options.decorated_tree,
                    options.taxa_list,
                    options.output_file)

    def arb_records(self, options):
        """Create an ARB records file from GTDB metadata."""

        check_file_exists(options.metadata_file)

        arb = Arb()
        arb.create_records(options.metadata_file,
                           options.msa_file,
                           options.taxonomy_file,
                           options.genome_list,
                           options.output_file)

    def parse_options(self, options):
        """Parse user options and call the correct pipeline(s)"""

        logging.basicConfig(format='', level=logging.INFO)

        if options.subparser_name == 'ssu_tree':
            self.ssu_tree(options)
        elif options.subparser_name == 'lsu_tree':
            self.lsu_tree(options)
        elif options.subparser_name == 'rna_tree':
            self.rna_tree(options)
        elif options.subparser_name == 'rna_dump':
            self.rna_dump(options)
        elif options.subparser_name == 'derep_tree':
            self.derep_tree(options)
        elif options.subparser_name == 'bootstrap':
            self.bootstrap(options)
        elif options.subparser_name == 'jk_markers':
            self.jk_markers(options)
        elif options.subparser_name == 'jk_taxa':
            self.jk_taxa(options)
        elif options.subparser_name == 'combine':
            self.combine(options)
        elif options.subparser_name == 'midpoint':
            self.midpoint(options)
        elif options.subparser_name == 'outgroup':
            self.outgroup(options)
        elif options.subparser_name == 'propagate':
            self.propagate(options)
        elif options.subparser_name == 'fill_ranks':
            self.fill_ranks(options)
        elif options.subparser_name == 'strip':
            self.strip(options)
        elif options.subparser_name == 'rm_support':
            self.rm_support(options)
        elif options.subparser_name == 'pull':
            self.pull(options)
        elif options.subparser_name == 'append':
            self.append(options)
        elif options.subparser_name == 'prune':
            self.prune(options)
        elif options.subparser_name == 'pd':
            self.phylogenetic_diversity(options)
        elif options.subparser_name == 'pd_clade':
            self.phylogenetic_diversity_clade(options)
        elif options.subparser_name == 'arb_records':
            self.arb_records(options)
        else:
            self.logger.error('Unknown GenomeTreeTk command: ' + options.subparser_name + '\n')
            sys.exit(-1)

        return 0
