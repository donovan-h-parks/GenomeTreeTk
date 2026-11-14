"""End-to-end tests of GenomeTreeTk commands.

Each case runs bin/genometreetk on the small inputs in tests/data and compares
every file it writes (except the log) with tests/data/expected/<case>/.

To regenerate the expected outputs after an intended change in behaviour:

    GENOMETREETK_UPDATE_EXPECTED=1 python -m unittest tests.test_cli

Run from the repository root with biolib/gtdblib importable, e.g.:

    python -m unittest discover -s tests -t .
"""

import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
BIN = REPO / 'bin' / 'genometreetk'
DATA = REPO / 'tests' / 'data'
EXPECTED = DATA / 'expected'
UPDATE = os.environ.get('GENOMETREETK_UPDATE_EXPECTED') == '1'

# name -> command line; {data} is the input directory and {out} the output directory
CASES = {
    'strip': 'strip {data}/decorated.tree {out}/stripped.tree',
    'rm_support': 'rm_support {data}/decorated.tree {out}/no_support.tree',
    'pull': 'pull {data}/decorated.tree {out}/taxonomy.tsv',
    'pull_no_validation': 'pull --no_validation {data}/decorated.tree {out}/taxonomy.tsv',
    'pull_undecorated': 'pull {data}/support_a.tree {out}/taxonomy.tsv',
    'append': 'append {data}/support_a.tree {data}/taxonomy.tsv {out}/appended.tree',
    'fill_ranks': 'fill_ranks {data}/partial_taxonomy.tsv {out}/filled.tsv',
    'propagate': 'propagate {data}/rep_taxonomy.tsv {data}/metadata.tsv {out}/propagated.tsv',
    'midpoint': 'midpoint {data}/decorated.tree {out}/midpoint.tree',
    'outgroup': 'outgroup {data}/decorated.tree {data}/taxonomy.tsv p__Halobacteriota {out}/outgroup.tree',
    'prune': 'prune {data}/decorated.tree {data}/taxa_subset.txt {out}/pruned.tree',
    'pd': 'pd {data}/decorated.tree {data}/ingroup.txt --per_taxa_pg_file {out}/pg.tsv',
    'pd_clade': 'pd_clade {data}/decorated.tree {data}/ingroup.txt {out}/pd_clade.tsv',
    'combine': 'combine {data}/support_a.tree {data}/support_b.tree {data}/support_c.tree {out}/combined.tree',
    'combine_no_root_support': ('combine {data}/support_no_root.tree {data}/support_b.tree {data}/support_c.tree '
                                '{out}/combined.tree'),
    'arb_records': ('arb_records {data}/metadata.tsv {out}/arb_records.txt --taxonomy_file {data}/taxonomy.tsv '
                    '--msa_file {data}/derep_msa.faa'),
    'jk_markers_jk_dir': 'jk_markers {data}/support_a.tree NONE {data}/marker_info.tsv NONE {out} --jk_dir {data}/boot_trees',
    'sh_uf_hack': 'sh_uf_hack {data}/sh_uf.tree {out}/sh_uf_combined.tree',
    'bootstrap_boot_dir': 'bootstrap {data}/support_a.tree None {out} --boot_dir {data}/boot_trees',
    'derep_tree': ('derep_tree {data}/decorated.tree g__Bacillus d__Archaea {data}/metadata.tsv {out} '
                   '--taxa_to_retain 1 --msa_file {data}/derep_msa.faa'),
    'msa_filter': 'msa_filter {data}/marker_info.tsv None {data}/msa.faa {out} --msa_length 1000',
}


def run_genometreetk(args, cwd):
    env = dict(os.environ)
    env['PYTHONPATH'] = os.pathsep.join([str(REPO)] + [p for p in env.get('PYTHONPATH', '').split(os.pathsep) if p])
    return subprocess.run([sys.executable, str(BIN)] + args, cwd=cwd, env=env,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)


def program_output(stdout):
    """Lines printed by a command, excluding timestamped log lines and progress bars."""
    return ''.join(line + '\n' for line in stdout.splitlines()
                   if line.strip() and not line.startswith('[') and '\r' not in line and '━' not in line)


class TestCommands(unittest.TestCase):

    def check_case(self, name):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / 'out'
            out.mkdir()
            args = CASES[name].format(data=DATA, out=out).split()
            proc = run_genometreetk(args, cwd=tmp)
            self.assertEqual(proc.returncode, 0, f'{name} failed:\n{proc.stdout}\n{proc.stderr}')

            stdout = program_output(proc.stdout)
            if stdout:
                (out / 'stdout.txt').write_text(stdout)

            produced = sorted(p.relative_to(out) for p in out.rglob('*')
                              if p.is_file() and p.name != 'genometreetk.log')

            expected_dir = EXPECTED / name
            if UPDATE:
                if expected_dir.exists():
                    shutil.rmtree(expected_dir)
                for rel in produced:
                    (expected_dir / rel).parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy(out / rel, expected_dir / rel)
                return

            expected = sorted(p.relative_to(expected_dir) for p in expected_dir.rglob('*') if p.is_file())
            self.assertEqual([str(p) for p in produced], [str(p) for p in expected], f'{name}: output files differ')
            for rel in produced:
                self.assertEqual((out / rel).read_text(), (expected_dir / rel).read_text(),
                                 f'{name}: {rel} differs from expected output')


def _make_test(name):
    def test(self):
        self.check_case(name)
    return test


for _name in CASES:
    setattr(TestCommands, f'test_{_name}', _make_test(_name))


class TestErrors(unittest.TestCase):

    def test_missing_input_file_exits_with_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            proc = run_genometreetk(['strip', str(Path(tmp) / 'missing.tree'), str(Path(tmp) / 'out.tree')], cwd=tmp)
            self.assertNotEqual(proc.returncode, 0)

    def test_help(self):
        with tempfile.TemporaryDirectory() as tmp:
            proc = run_genometreetk(['bootstrap', '-h'], cwd=tmp)
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertIn('--boot_dir', proc.stdout)
            # CustomHelpFormatter appends defaults to help strings
            self.assertIn('(default: 100)', ' '.join(proc.stdout.split()))


if __name__ == '__main__':
    unittest.main()
