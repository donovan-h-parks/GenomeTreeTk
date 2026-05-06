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
import logging
import random
from typing import Tuple, List, Dict
from collections import defaultdict, Counter
from dataclasses import dataclass

from biolib.seq_io import read_fasta

from tqdm import tqdm


GAP_CHARS = {'-', '.', '_', '*'}
STANDARD_AMINO_ACIDS = set('ACDEFGHIKLMNPQRSTVWY')


@dataclass
class MarkerInfo:
    """Information about marker genes."""
    id: str
    name: str
    description: str
    length: int


@dataclass
class MarkerFilterInfo:
    id: str
    length: int
    valid_cols: List[int]
    gap_filtered: int
    aa_conserved_filtered: int
    aa_unconserved_filtered: int


class MSA_Filter(object):
    """Filter MSA based on gaps and conservation of bases."""

    def __init__(self,
                    msa_length: int = 5000, 
                    min_perc_aa: int = 50, 
                    min_consensus: int = 25, 
                    max_consensus: int = 95):
        """Initialization."""
        
        self.logger = logging.getLogger('timestamp')

        self.target_msa_length = msa_length
        self.max_gaps_perc = (100 - min_perc_aa) / 100.0
        self.min_consensus_perc = min_consensus / 100.0
        self.max_consensus_perc = max_consensus / 100.0

    def parse_marker_info(self, marker_info_file: str) -> List[MarkerInfo]:
        """Parse marker info file.
        
        Parameters
        ----------
        marker_info_file : str
          File indicate marker genes in concatenated alignment (e.g. gtdb_r<cur>_bac120_markers_info.tsv).
        
        Returns
        -------
        List[MarkerInfo]
          List of MarkerInfo objects.
        """
        
        marker_info_list = []
        with open(marker_info_file) as f:
            header = f.readline().strip().split('\t')

            marker_id_idx = header.index('Marker Id')
            marker_name_idx = header.index('Name')
            marker_desc_idx = header.index('Description')
            marker_len_idx = header.index('Length (bp)')

            for line in f:
                tokens = line.strip().split('\t')
                marker_info = MarkerInfo(id=tokens[marker_id_idx],
                                            name=tokens[marker_name_idx],
                                            description=tokens[marker_desc_idx],
                                            length=int(tokens[marker_len_idx]))
                marker_info_list.append(marker_info)
        
        return marker_info_list

    def find_marker(self, col_idx: int, marker_info: List[MarkerInfo]) -> int:
        """Find marker for given column index.
        
        Parameters
        ----------
        col_idx : int
          Column index in concatenated MSA.
        marker_info : List[MarkerInfo]
          List of MarkerInfo objects.

        Returns
        -------
        int
          Marker gene corresponding to column index.
        """

        marker_end = 0
        for cur_marker_info in marker_info:
            marker_end += cur_marker_info.length
            if col_idx <= marker_end:
                return cur_marker_info.id

        raise ValueError(f"Column index {col_idx} exceeds total MSA length.")

    def parse_witchi_info(self, marker_info: List[MarkerInfo], witchi_info_file: str) -> Tuple[Dict[str, int], List[int]]:
        """Parse WitChi info file.
        
        Parameters
        ----------
        witchi_info_file : str
          File indicating columns pruned by WitChi (e.g. gtdb_r<cur>_bac120_full_concatenated_quartic_s400_touchdown_pruned.tsv).

        Returns
        -------
        Dict[str, int]
          Number of columns trimmed from each marker.
        List[bool]
          Mask indicating columns retained (1) or pruned (0).
        """

        total_marker_len = sum([m.length for m in marker_info])

        num_trimmed_cols = defaultdict(int)
        mask = [1]*total_marker_len
        with open(witchi_info_file) as f:
            header = f.readline().strip().split('\t')

            orig_index_idx = header.index('Original Index')

            for line in f:
                tokens = line.strip().split('\t')

                col_idx = int(tokens[orig_index_idx])
                marker_id = self.find_marker(col_idx, marker_info)
                num_trimmed_cols[marker_id] += 1

                mask[col_idx] = 0

        return num_trimmed_cols, mask

    def write_witchi_stats(self, 
                          marker_info: List[MarkerInfo], 
                          witchi_num_trimmed_cols: Dict[str, int], 
                          output_dir: str):
        """Write WitChi trimming statistics."""

        fout = open(os.path.join(output_dir, 'witchi_gtdb_marker_trimming_info.tsv'), 'w')
        fout.write('marker_id\tdescription\tlength\ttrimmed_columns\ttrimmed_columns_prec\n')
        for cur_marker_info in marker_info:
            num_trimmed_cols = witchi_num_trimmed_cols[cur_marker_info.id]
            fout.write('{}\t{}\t{}\t{}\t{:.2f}\n'.format(
                cur_marker_info.id,
                cur_marker_info.description,
                cur_marker_info.length,
                num_trimmed_cols,
                100 * num_trimmed_cols / cur_marker_info.length
            ))

        fout.close()

    def combine_mask(self, original_mask: List[int], filtered_mask: List[int]) -> List[int]:
        """Create mask with additional filtered columns.
        
        original_mask : List[int]
            The original mask.

        filtered_mask : List[int]
            The filtered mask that indicates additional columns to filter
            in the original mask.
        """

        # create mapping from retained cols to index position
        # in original mask
        mask_idx_map = {}
        mask_idx = 0
        for idx, state in enumerate(original_mask):
            if state == 0:
                continue
            
            mask_idx_map[mask_idx] = idx
            mask_idx += 1

        new_mask = list(original_mask)
        for idx, state in enumerate(filtered_mask):
            if state == 1:
                continue
            
            orig_mask_idx = mask_idx_map[idx]
            assert new_mask[orig_mask_idx] == 1

            new_mask[orig_mask_idx] = 0

        return new_mask

    def identify_valid_columns(self, start: int, end: int, seqs: Dict[str, str]) -> Tuple[set, int, int, int]:
        """Identify columns meeting gap and amino acid ubiquity criteria (micro-optimized)."""

        seq_values = list(seqs.values())
        num_seqs = len(seq_values)
        if num_seqs == 0 or end <= start:
            self.logger.error('Invalid execution of identify_valid_columns: no sequences or invalid range.')
            sys.exit(1)

        # Pre-slice and uppercase once
        slices = [s[start:end].upper() for s in seq_values]

        # Iterate column-wise; zip(*slices) yields tuples of characters for each column
        gap_filtered = 0
        aa_conserved_filtered = 0
        aa_unconserved_filtered = 0
        valid_cols = set()

        for idx, col in enumerate(zip(*slices)):
            gap_cnt = 0
            counts = {}
            for ch in col:
                if ch in GAP_CHARS:
                    gap_cnt += 1
                else:
                    counts[ch] = counts.get(ch, 0) + 1

            if float(gap_cnt) / num_seqs <= self.max_gaps_perc:
                if not counts:
                    gap_filtered += 1
                    continue

                # Find most common amino acid
                letter, count = max(counts.items(), key=lambda t: t[1])
                if letter not in STANDARD_AMINO_ACIDS:
                    self.logger.warning(f'Most common amino acid was not in standard alphabet: {letter}')

                aa_ratio = float(count) / (num_seqs - gap_cnt)
                if aa_ratio < self.min_consensus_perc:
                    aa_unconserved_filtered += 1
                elif aa_ratio >= self.max_consensus_perc:
                    aa_conserved_filtered += 1
                elif self.min_consensus_perc <= aa_ratio < self.max_consensus_perc:
                    valid_cols.add(idx)
                else:
                    self.logger.error('Unexpected case. Should never occur.')
                    sys.exit(1)
            else:
                gap_filtered += 1

        return valid_cols, gap_filtered, aa_conserved_filtered, aa_unconserved_filtered

    def filter_columns(self,
                        witchi_msa: Dict[str, str], 
                        marker_info: List[MarkerInfo], 
                        witchi_num_trimmed_cols: Dict[str, int]) -> List[MarkerFilterInfo]:
        """Filter columns based on gap and amino acid conservation criteria."""

        start = 0
        filter_info = []
        for cur_marker_info in tqdm(marker_info, desc='Filtering markers', unit='marker'):
            marker_len = cur_marker_info.length - witchi_num_trimmed_cols[cur_marker_info.id]
            end = start + marker_len

            rtn = self.identify_valid_columns(start, end, witchi_msa)
            valid_cols, cur_gap_filtered, cur_aa_conserved_filtered, cur_aa_unconserved_filtered = rtn
            
            assert len(valid_cols) <= marker_len # sanity check

            offset_valid_cols = [i + start for i in valid_cols]

            filter_info.append(MarkerFilterInfo(
                id=cur_marker_info.id,
                length=marker_len,
                valid_cols=offset_valid_cols,
                gap_filtered=cur_gap_filtered,
                aa_conserved_filtered=cur_aa_conserved_filtered,
                aa_unconserved_filtered=cur_aa_unconserved_filtered
            ))

            start = end

        return filter_info

    def subsample_columns(self, filter_info: List[MarkerFilterInfo]) -> Tuple[int, List[int]]:
        """Subsample columns from filtered markers."""

        # determine number of columns that should be subsampled in order
        # to obtain MSA of desired length
        for cols_to_sample in range(1, 1000):
            msa_length = 0
            for cur_filter_info in filter_info:
                msa_length += min(len(cur_filter_info.valid_cols), cols_to_sample)

            if msa_length >= self.target_msa_length:
                break

        # randomly sample columns from each marker
        rnd_sampled_cols = []
        for cur_filter_info in filter_info:
            sel_cols = random.sample(cur_filter_info.valid_cols, min(cols_to_sample, len(cur_filter_info.valid_cols)))
            rnd_sampled_cols.extend(sel_cols)

        return cols_to_sample, rnd_sampled_cols

    def run(self, 
            marker_info_file: str, 
            witchi_info_file: str, 
            witch_msa_file: str, 
            output_dir: str):
        """Filter MSA based on gaps and conservation of bases.
        
        Parameters
        ----------
        marker_info_file : str
          File indicate marker genes in concatenated alignment (e.g. gtdb_r<cur>_bac120_full_markers_info.tsv).
        witchi_info_file : str
          File indicating columns pruned by WitChi (e.g. gtdb_r<cur>_bac120_full_concatenated_quartic_s400_touchdown_pruned.tsv).
        witch_msa_file : str
          File with WitChi pruned MSA (e.g. gtdb_r<cur>_bac120_full_concatenated_quartic_s400_touchdown_pruned.fasta).
        output_dir : str
          Output directory.
        """

        # determine number of columns for each marker gene
        self.logger.info('Parsing marker info file:')
        marker_info = self.parse_marker_info(marker_info_file)
        concatenated_marker_len = sum([mi.length for mi in marker_info])
        self.logger.info(f' - concatenated marker length: {concatenated_marker_len:,}')

        # determine columns trimmed by WitChi
        self.logger.info('Determining columns trimmed by WitChi:')
        witchi_num_trimmed_cols, witchi_mask = self.parse_witchi_info(marker_info, witchi_info_file)
        self.logger.info(f' - retained columns: {sum(witchi_mask):,}')

        if len(witchi_mask) != concatenated_marker_len:
            self.logger.error(f'Concatenated marker length and length of WitChi mask are different: {concatenated_marker_len:,} vs {len(witchi_mask):,}')
            sys.exit(1)
        
        # write out WitChi mask and filtering statistics
        with open(os.path.join(output_dir, 'witchi_mask.txt'), 'w') as mask_file:
            mask_file.write(''.join([str(n) for n in witchi_mask]))

        self.write_witchi_stats(marker_info, witchi_num_trimmed_cols, output_dir)

        # read WitChi filtered MSA
        self.logger.info('Reading WitChi filtered MSA:')
        witchi_msa = read_fasta(witch_msa_file)
        first_seq = next(iter(witchi_msa.values()))
        self.logger.info(' - MSA contains {:,} sequences and has {:,} columns'.format(len(witchi_msa), len(first_seq)))

        if sum(witchi_mask) != len(first_seq):
            self.logger.error(f'WitChi mask and length WitChi-filtered MSA are different: {sum(witchi_mask):,} vs {len(first_seq):,}')
            sys.exit(1)

        # filter columns in WitChi MSA
        self.logger.info('Identifying columns to retain based on gap and AA conservation criteria:')
        filter_info = self.filter_columns(witchi_msa, marker_info, witchi_num_trimmed_cols)

        gap_filtered = sum([f.gap_filtered for f in filter_info])
        aa_conserved_filtered = sum([f.aa_conserved_filtered for f in filter_info])
        aa_unconserved_filtered = sum([f.aa_unconserved_filtered for f in filter_info])
        
        self.logger.info(f' - gappy columns: {gap_filtered:,}')
        self.logger.info(f' - conserved columns: {aa_conserved_filtered:,}')
        self.logger.info(f' - unconserved columns: {aa_unconserved_filtered:,}')

        total_valid_cols = sum([len(f.valid_cols) for f in filter_info])
        min_valid_cols = min([len(f.valid_cols) for f in filter_info])
        max_valid_cols = max([len(f.valid_cols) for f in filter_info])
        self.logger.info(f' - retained {total_valid_cols:,} columns (min: {min_valid_cols:,}, max: {max_valid_cols:,}) after filtering')

        sampled_cols = set([i for f in filter_info for i in f.valid_cols])
        filtered_mask = [1 if i in sampled_cols else 0 for i in range(len(first_seq))]
        witchi_filtered_mask = self.combine_mask(witchi_mask, filtered_mask)
        with open(os.path.join(output_dir, 'filtered_mask.txt'), 'w') as mask_file:
            mask_file.write(''.join([str(n) for n in witchi_filtered_mask]))

        # randomly subsample columns to obtain MSA of desired length
        self.logger.info('Randomly subsampling columns to obtain MSA of ~{:,}:'.format(self.target_msa_length))
        cols_to_sample, rnd_sampled_cols = self.subsample_columns(filter_info)

        num_markers_below_cols_to_sample = sum([1 for f in filter_info if len(f.valid_cols) < cols_to_sample])
        self.logger.info(f' - sampling {cols_to_sample:,} columns from each marker')
        self.logger.info(f' - {num_markers_below_cols_to_sample:,} markers had fewer than {cols_to_sample:,} columns available for sampling')
        self.logger.info(f' - subsampled MSA has {len(rnd_sampled_cols):,} columns')

        subsampled_mask = [1 if i in rnd_sampled_cols else 0 for i in range(len(first_seq))]
        final_mask = self.combine_mask(witchi_mask, subsampled_mask)
        num_mask_cols = sum(final_mask)
        with open(os.path.join(output_dir, 'final_mask.txt'), 'w') as mask_file:
            mask_file.write(''.join([str(n) for n in final_mask]))

        if num_mask_cols != len(rnd_sampled_cols):
            self.logger.error(f"Final mask and MSA do not have the sample number of cols: {num_mask_cols} vs {len(rnd_sampled_cols)}")
            sys.exit(1)

        # create filtered and subsampled MSA
        final_msa = {}
        for seq_id, seq in witchi_msa.items():
            masked_seq = ''.join([seq[i] for i in sorted(rnd_sampled_cols)])
            final_msa[seq_id] = masked_seq

        # write out filtered MSA
        fout = open(os.path.join(output_dir, 'genome_msa_stats.tsv'), 'w')
        fout.write('Genome ID\tMSA length\tAmino acids\tAmino acids (%)\n')
        with open(os.path.join(output_dir, 'filtered_msa.faa'), 'w') as filter_file:
            for gid, seq in final_msa.items():
                fasta_outstr = ">%s\n%s\n" % (gid, seq)
                filter_file.write(fasta_outstr)

                aa_len = sum([1 for c in seq if c.isalpha()])
                fout.write('{}\t{}\t{}\t{:.2f}\n'.format(
                    gid,
                    len(seq),
                    aa_len,
                    100.0 * aa_len / len(seq)
                ))
        fout.close()

        # write out filtering statistics
        with open(os.path.join(output_dir, 'filtered_msa_stats.tsv'), 'w') as fout:
            fout.write('marker_id\tdescription\tlength\twitchi_filtered_columns\twitchi_filtered_columns_perc')
            fout.write('\tfiltered_columns\tfiltered_columns_prec')
            fout.write('\tgappy_columns\tconserved_columns\tunconserved_columns\n')

            for idx, cur_filter_info in enumerate(filter_info):
                cur_marker_info = marker_info[idx]
                witchi_filtered_columns = witchi_num_trimmed_cols[cur_filter_info.id]

                assert cur_filter_info.length + witchi_filtered_columns == cur_marker_info.length

                fout.write('{}\t{}\t{}\t{}\t{:.2f}\t{}\t{:.2f}\t{}\t{}\t{}\n'.format(
                    cur_filter_info.id,
                    cur_marker_info.description,
                    cur_marker_info.length,
                    witchi_filtered_columns,
                    100 * witchi_filtered_columns / cur_marker_info.length,
                    len(cur_filter_info.valid_cols),
                    100 * len(cur_filter_info.valid_cols) / cur_marker_info.length,
                    cur_filter_info.gap_filtered,
                    cur_filter_info.aa_conserved_filtered,
                    cur_filter_info.aa_unconserved_filtered,
                ))
