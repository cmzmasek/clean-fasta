# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0.0] - 2026-07-08

First release of `clean-fasta` as its own installable package. This version
grew out of the original `clean_fasta.py` / `molseq.py` scripts.

### Added
- Installable package with a `clean-fasta` console command
  (`pip install .` provides it on your PATH).
- Public Python API: `stream_fasta`, `filter_sequences`, `clean_fasta_file`,
  `format_stats`, `CleanStats`, and `MolSeq`.
- Read from stdin / write to stdout using `-` for the input/output path.
- `--force` to overwrite an existing output file (previously the tool exited
  silently with a success status when the output already existed).
- `--wrap N` to control the output line width (`0` disables wrapping).
- `--allow-digits` to keep sequences containing digit characters.
- `--no-unique` as the explicit counterpart to `--unique`.
- A full `pytest` test suite and GitHub Actions CI across Python 3.9-3.13.

### Changed
- **Breaking:** modernized the command-line interface.
  - `-t aa|na` is now validated via `--type {aa,na}`.
  - `-u t|f` is now the boolean flag `--unique` / `--no-unique`.
  - `-ml` / `-mal` are now `-m/--min-length` / `-M/--max-length`.
  - `-r` is now `-r/--min-ratio`.
  - The summary report is written to **stderr** so that `-` output can be piped.
- Split the code into `clean_fasta.molseq`, `clean_fasta.cleaner`, and
  `clean_fasta.cli` so the filtering logic can be tested without file I/O.

### Fixed
- A bare `>` header (no identifier) no longer raises `AttributeError`; the
  record is counted as "empty id" and skipped.
- Zero-length sequences no longer cause a `ZeroDivisionError` when
  `--min-length 0` is used.
- Files are now opened with explicit UTF-8 encoding and always closed, even on
  error.
- Overwriting an existing output file now reports an error on stderr and exits
  with a non-zero status instead of exiting silently with status 0.

### Notes on filtering behavior (unchanged from the original tool)
- Gaps (`-`) and whitespace are removed from every sequence body.
- Sequences containing digits are dropped (unless `--allow-digits`).
- For amino acids, the "valid" ratio is the fraction of characters that are not
  in `_ - ? X x * .`; for nucleic acids it is the fraction of `A/C/G/T`.
