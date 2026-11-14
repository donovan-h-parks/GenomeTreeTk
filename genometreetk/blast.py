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

"""Minimal BLAST wrapper used to identify rRNA sequences with incongruent taxonomy."""

import subprocess
from collections import namedtuple

from gtdblib.util.shell.execute import check_dependencies

# fields written by blastn with BLAST_CUSTOM_FMT
BLAST_CUSTOM_FMT = '6 qseqid qlen sseqid stitle slen length pident evalue bitscore'
BlastHit = namedtuple('BlastHit', ['query_id',
                                   'query_len',
                                   'subject_id',
                                   'subject_annotation',
                                   'subject_len',
                                   'alignment_len',
                                   'perc_identity',
                                   'evalue',
                                   'bitscore'])


class BlastError(Exception):
    pass


def _run(cmd):
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    if proc.returncode != 0:
        raise BlastError('{} failed with exit code {}: {}'.format(
            cmd[0], proc.returncode, proc.stderr.strip()))


def create_blastn_db(fasta_file, db_file):
    """Create nucleotide BLAST database."""

    check_dependencies(['makeblastdb'])
    _run(['makeblastdb', '-dbtype', 'nucl', '-in', fasta_file, '-out', db_file])


def blastn(query_seqs, nucl_db, output_file, cpus=1, evalue=1e-3, max_matches=500, task='megablast'):
    """Search query sequences against a nucleotide database, writing hits in BLAST_CUSTOM_FMT."""

    check_dependencies(['blastn'])
    _run(['blastn',
          '-num_threads', str(cpus),
          '-query', query_seqs,
          '-db', nucl_db,
          '-out', output_file,
          '-evalue', str(evalue),
          '-max_target_seqs', str(max_matches),
          '-task', task,
          '-outfmt', BLAST_CUSTOM_FMT])


def read_hits(table):
    """Read hits from a table written by blastn()."""

    with open(table) as f:
        for line in f:
            tokens = line.rstrip('\n').split('\t')
            yield BlastHit(query_id=tokens[0],
                           query_len=int(tokens[1]),
                           subject_id=tokens[2],
                           subject_annotation=tokens[3],
                           subject_len=int(tokens[4]),
                           alignment_len=int(tokens[5]),
                           perc_identity=float(tokens[6]),
                           evalue=float(tokens[7]),
                           bitscore=float(tokens[8]))
