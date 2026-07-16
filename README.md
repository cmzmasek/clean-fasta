# clean-fasta

[![PyPI](https://img.shields.io/pypi/v/clean-fasta.svg)](https://pypi.org/project/clean-fasta/)
[![CI](https://github.com/cmzmasek/clean-fasta/actions/workflows/ci.yml/badge.svg)](https://github.com/cmzmasek/clean-fasta/actions/workflows/ci.yml)
[![Python](https://img.shields.io/pypi/pyversions/clean-fasta.svg)](https://pypi.org/project/clean-fasta/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

A small, dependency-free command-line tool (and Python library) for **cleaning
and filtering FASTA sequence files**. It strips gaps and whitespace, drops
sequences that are too short, too long, malformed, or below a valid-character
threshold, optionally removes duplicate identifiers, and writes tidy,
line-wrapped FASTA output — along with a **before/after statistics report** of
what happened (length distribution, N50/L50, GC content, extremes, and more).
It can also profile a file without cleaning it (`--stats-only`) or preview a run
without writing anything (`--dry-run`).

## What it does

For every sequence in the input, `clean-fasta`:

1. Removes gap characters (`-`) and all whitespace from the sequence body.
2. Collapses runs of whitespace in the identifier line.
3. Drops the sequence if it:
   - has an **empty identifier**;
   - **contains digits** (unless `--allow-digits`);
   - is **shorter** than `--min-length` (default `20`);
   - is **longer** than `--max-length` (if set);
   - has a **valid-character ratio** below `--min-ratio` (default `0.99`);
   - has a **duplicate identifier** (unless `--no-unique`).
4. Writes the surviving sequences, wrapped at `--wrap` characters (default `80`).

The **valid-character ratio** is:
- for `--type aa`: the fraction of characters that are *not* `_ - ? X x * .`;
- for `--type na`: the fraction of characters that are `A`, `C`, `G`, or `T`.

## Installation

Requires Python 3.9+.

```bash
pip install clean-fasta
```

This installs a `clean-fasta` command on your PATH. Verify it with:

```bash
clean-fasta --version
```

To install the latest development version from source instead:

```bash
git clone https://github.com/cmzmasek/clean-fasta
cd clean-fasta
pip install .

# or, for development (editable install + test/lint tools)
pip install -e ".[dev]"
```

## Usage

```
clean-fasta [OPTIONS] INPUT [OUTPUT]
```

Use `-` for `INPUT` or `OUTPUT` to read from stdin / write to stdout. `OUTPUT`
is required for a normal run, but not with `--dry-run` or `--stats-only`.

### Options

| Option | Default | Description |
| --- | --- | --- |
| `-t`, `--type {aa,na}` | `aa` | Sequence type: amino acid or nucleic acid. |
| `-m`, `--min-length N` | `20` | Minimum sequence length. |
| `-M`, `--max-length N` | none | Maximum sequence length. |
| `-r`, `--min-ratio R` | `0.99` | Minimum ratio of valid characters (0.0–1.0). |
| `-u`, `--unique` / `--no-unique` | `--unique` | Drop / keep sequences with duplicate ids. |
| `--allow-digits` | off | Keep sequences that contain digit characters. |
| `-w`, `--wrap N` | `80` | Wrap output lines at N chars (`0` = no wrapping). |
| `-f`, `--force` | off | Overwrite the output file if it exists. |
| `-q`, `--quiet` | off | Suppress the summary report. |
| `-n`, `--dry-run` | off | Analyze and report, but write no output. |
| `--stats-only` | off | Print statistics for `INPUT` only; write nothing (no `OUTPUT` needed). |
| `--version` | | Print the version and exit. |

The summary report is written to **stderr**, so piping via stdout stays clean.

### Statistics report

Every run (unless `--quiet`) prints a **before/after** report: descriptive
statistics over the input set next to the same statistics over the kept
(surviving) set. This includes sequence count, total residues, the length
distribution (min / Q1 / median / Q3 / max / mean / stddev), **N50** and
**L50**, **GC content** (nucleic-acid input only), the overall valid-character
ratio, and the named longest / shortest / most- and least-"bad"-character
sequences.

```bash
# Profile a file without cleaning it (no OUTPUT argument needed)
clean-fasta genes.fasta --stats-only -t na

# Preview what a cleaning run would do, writing nothing
clean-fasta genes.fasta genes.clean.fasta --dry-run -t na -m 200
```

### Examples

```bash
# Clean a protein FASTA, keeping sequences >= 50 residues
clean-fasta proteins.fasta proteins.clean.fasta -m 50

# Nucleotide sequences, allow up to 2% ambiguous bases, cap length at 5000
clean-fasta genes.fasta genes.clean.fasta -t na -r 0.98 -M 5000

# Use in a pipeline (read stdin, write stdout, no report)
gunzip -c raw.fasta.gz | clean-fasta - - -q > clean.fasta

# Overwrite an existing output file
clean-fasta in.fasta out.fasta --force
```

## Use as a library

The filtering logic is I/O-free and importable:

```python
from clean_fasta import stream_fasta, filter_sequences, CleanStats

stats = CleanStats()
with open("in.fasta") as handle:
    for seq in filter_sequences(stream_fasta(handle), min_length=50, stats=stats):
        print(seq.to_fasta_wrapped(80))

print("passed:", stats.passed, "of", stats.total)
```

Or clean a file in one call:

```python
from clean_fasta import clean_fasta_file

stats = clean_fasta_file("in.fasta", "out.fasta", min_length=50, is_aa=True, force=True)
```

## Development

```bash
pip install -e ".[dev]"
pytest          # run the test suite
ruff check .    # lint
```

### Releasing (maintainers)

Releases are published to PyPI automatically by
[`.github/workflows/release.yml`](.github/workflows/release.yml) using PyPI
[trusted publishing](https://docs.pypi.org/trusted-publishers/) — no API tokens
are stored anywhere.

To cut a release:

1. Bump `__version__` in `src/clean_fasta/__init__.py` and update `CHANGELOG.md`.
2. Commit and push to `main`.
3. Create a GitHub Release whose tag is `v<version>` (e.g. `v2.0.1`).

Publishing the release triggers the workflow, which builds the sdist and wheel
and uploads them to PyPI. The workflow fails fast if the release tag does not
match the packaged version, so a forgotten version bump can't ship.

## License

[MIT](LICENSE) © Christian M. Zmasek
