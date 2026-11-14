# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

GenomeTreeTk is a Python 3 command-line toolbox for working with genome trees (rRNA tree inference, MSA filtering, tree support/stability, rerooting, taxonomy decoration, pruning, phylogenetic diversity). It is primarily used by the GTDB curation team, so many commands assume GTDB-formatted inputs (GTDB metadata TSV, GTDB genome IDs, 7-rank `d__;p__;...;s__` taxonomy strings).

## Running and development

There is no build step or linter config. Run the CLI directly from a checkout:

```bash
python3 bin/genometreetk -h
python3 bin/genometreetk <command> -h
```

Install for use on `PATH` with `pip install .` (packaging uses `setup.py` with setuptools).

Python dependencies: `gtdblib` (≥ 1.12.1, the GTDB team's shared library at `~/git/gtdb-lib`; supplies file checks, FASTA I/O, Newick label parsing, taxonomy helpers, bootstrap support calculation, and the CLI help formatter), `dendropy` (all tree manipulation), `numpy` (≥ 1.17), and `tqdm`. Several commands also shell out to external tools, checked at runtime with `gtdblib.util.shell.execute.check_dependencies`: FastTree/FastTreeMP, makeblastdb/blastn, ssu-align/ssu-mask, mothur, and Infernal/Easel (`cmsearch`, `cmalign`, `esl-sfetch`, `esl-alimask`). The mothur path in `rna_workflow.py` hard-codes `/srv/sw/Mothur/...` and `/srv/db/mothur/silva_128/...`.

### Tests

`tests/test_cli.py` runs each command through `bin/genometreetk` on the small inputs in `tests/data/` and compares every file it writes with `tests/data/expected/<case>/`. The tests use `unittest`, so pytest is optional:

```bash
python -m unittest discover -s tests
python -m unittest discover -s tests -k test_bootstrap_boot_dir
```

After an intended change in output, regenerate the expected files with `GENOMETREETK_UPDATE_EXPECTED=1 python -m unittest discover -s tests -k <case>` and review the diff. The tests don't run FastTree, BLAST or the rRNA workflows (`ssu_tree`, `lsu_tree`, `rna_tree`, `rna_dump`).

Keep test trees free of ties between equally long paths. dendropy breaks such ties differently between runs, which makes `midpoint` output nondeterministic.

## Architecture

The CLI follows a three-layer pattern; adding or changing a command usually touches all three:

1. **`bin/genometreetk`** defines every subcommand's argparse parser and the hand-written `print_help()` menu (update both when adding a command). It sets up logging to stdout and, if the command has an `output_dir` argument, to `<output_dir>/genometreetk.log`. It then calls `OptionsParser().parse_options(args)`.
2. **`genometreetk/main.py`** — `OptionsParser` has one method per command, plus an `if/elif` chain in `parse_options` that dispatches on `options.subparser_name`. Note that some method names differ from the command name (`pd` → `phylogenetic_diversity`, `pd_clade` → `phylogenetic_diversity_clade`). These methods validate inputs (`check_file_exists`, `make_sure_path_exists`, `check_dependencies`) and then either call a worker class or, for the simpler taxonomy/tree commands (`fill_ranks`, `propagate`, `strip`, `rm_support`, `pull`, `append`), do the work inline with dendropy and gtdblib's taxonomy functions.
3. **Worker modules** in `genometreetk/` each expose one class with a `run(...)` method (e.g. `MSA_Filter`, `Bootstrap`, `JackknifeMarkers`, `JackknifeTaxa`, `CombineSupport`, `RerootTree`, `Prune`, `PhylogeneticDiversity`, `RNA_Workflow`, `Arb`, `DereplicateTree`).

Other notes:

- `common.py` holds the GTDB metadata readers `read_gtdb_metadata(metadata_file, fields)` (returns a dict keyed by genome ID of namedtuples with the requested columns) and `read_gtdb_taxonomy`.
- `fasttree.py` runs FastTree for `bootstrap`, `jk_markers` and `jk_taxa`. `run_fasttree` checks the exit code, writes the tree to `*.tmp` and renames it when complete, and keeps stderr in `<prefix>.fasttree.stderr`. `infer_trees` runs many trees in parallel and skips trees that are already complete. `blast.py` is a minimal BLAST wrapper for `rna_workflow.py`.
- Bootstrap and jackknife-marker support values are calculated with `gtdblib.tree.bootstrap_merge.bootstrap_merge_replicates`. It keeps the input tree's rooting and taxon labels and writes whole-number support, and it leaves the root label unchanged. `combine` therefore has to cope with nodes that have no support value. `jk_taxa` uses its own `tree_support.TreeSupport.subset_taxa`, and its replicate-inference step is deliberately commented out (`#***`).
- Internal-node labels in Newick trees often combine support and taxon names (e.g. `'95.3/80.2:INT1'` or `'100:p__Foo'`). Use `gtdblib.util.bio.newick.parse_label`/`create_label` to split and rebuild them instead of parsing by hand. Trees are read with `preserve_underscores=True` and written with `unquoted_underscores=True` so genome IDs and taxon names survive unchanged.
- `msa_filter` works on a concatenated marker MSA. It uses a marker-info TSV (marker column spans), an optional WitChi pruning TSV (pass `None` to skip it), and the WitChi-pruned MSA. It filters gappy and over- or under-conserved columns per marker, then randomly subsamples to `--msa_length` while preserving the column order within each marker. It writes `witchi_mask.txt`, `filtered_mask.txt`, `final_mask.txt`, `filtered_msa.faa`, and stats TSVs to the output directory.
- `bootstrap` has to handle production MSAs of around 200,000 × 5,000:
  - It loads the MSA into a single numpy `uint8` matrix and shares it with workers through a module-level global and a `fork` pool. Spawn or forkserver would copy the matrix to every worker, and a dict of `str` gets duplicated page by page by refcount writes.
  - Each replicate is seeded from `(seed, replicate index)`, so output does not depend on `-c`.
  - On a rerun, an MSA replicate is skipped only if its byte size matches exactly, and a tree only if it is complete.
- Logging: gtdblib logs to its own `rich` logger, which it gives a RichHandler and stops from propagating. `logger_setup` in `bin/genometreetk` removes that handler and turns propagation back on, so gtdblib's messages use GenomeTreeTk's format, go to `genometreetk.log`, and obey `--silent`. `bin/genometreetk` also reports a `GtdbLibException` (e.g. a taxon without a rank prefix) as a logged error with exit code 1 instead of a traceback.
- Errors should log the problem and call `sys.exit(1)`. `bin/genometreetk` passes the exit code through (as do gtdblib's `check_file_exists` etc.), so batch jobs can detect failures.
- `cm_files/` contains Infernal covariance models for LSU and ships as package data.

## Versioning and commits

- `genometreetk/VERSION` is both the version source and the changelog. Its first line is the version, read by `genometreetk/__init__.py` and `setup.py`. Below it, each release is listed newest first with `- ` bullet notes. To bump the version, add a new version line and its notes at the top.
- Commit messages use conventional prefixes: `feat:`, `fix:`, `chore:`.
