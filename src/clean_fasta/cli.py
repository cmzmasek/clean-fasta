# Copyright (c) 2025 Christian M. Zmasek
#
# Permission is hereby granted, free of charge, to any person
# obtaining a copy of this software and associated
# documentation files (the "Software"), to deal in the Software
# without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to
# whom the Software is furnished to do so, subject to the
# following conditions:
#
# The above copyright notice and this permission notice shall
# be included in all copies or substantial portions of the
# Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY
# KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
# WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
# PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
# OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

"""Command-line interface for ``clean-fasta``."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence

from clean_fasta import __version__
from clean_fasta.cleaner import (
    DEFAULT_MIN_LENGTH,
    DEFAULT_MIN_RATIO,
    DEFAULT_WRAP,
    CleanStats,
    clean_fasta_file,
    format_stats,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="clean-fasta",
        description=(
            "Clean and filter a FASTA file: strip gaps/whitespace, filter by "
            "length and valid-character ratio, and optionally deduplicate ids."
        ),
        epilog="Use '-' as INPUT or OUTPUT to read from stdin / write to stdout.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("input", help="input FASTA file (or '-' for stdin)")
    parser.add_argument(
        "output", nargs="?", default=None,
        help="output FASTA file (or '-' for stdout); "
        "not required with --dry-run or --stats-only",
    )

    parser.add_argument(
        "-t", "--type", choices=("aa", "na"), default="aa",
        help="sequence type: amino acid or nucleic acid",
    )
    parser.add_argument(
        "-m", "--min-length", type=int, default=DEFAULT_MIN_LENGTH, metavar="N",
        help="minimum sequence length",
    )
    parser.add_argument(
        "-M", "--max-length", type=int, default=None, metavar="N",
        help="maximum sequence length (default: no limit)",
    )
    parser.add_argument(
        "-r", "--min-ratio", type=float, default=DEFAULT_MIN_RATIO, metavar="R",
        help="minimum ratio of valid characters (0.0-1.0)",
    )
    parser.add_argument(
        "-u", "--unique", action=argparse.BooleanOptionalAction, default=True,
        help="drop sequences whose id was already seen",
    )
    parser.add_argument(
        "--allow-digits", action="store_true",
        help="keep sequences that contain digit characters (dropped by default)",
    )
    parser.add_argument(
        "-w", "--wrap", type=int, default=DEFAULT_WRAP, metavar="N",
        help="wrap output sequence lines at N characters (0 disables wrapping)",
    )
    parser.add_argument(
        "-f", "--force", action="store_true",
        help="overwrite the output file if it already exists",
    )
    parser.add_argument(
        "-q", "--quiet", action="store_true",
        help="do not print the summary report to stderr",
    )

    mode_group = parser.add_mutually_exclusive_group()
    mode_group.add_argument(
        "-n", "--dry-run", action="store_true",
        help="analyze and report, but do not write any output",
    )
    mode_group.add_argument(
        "--stats-only", action="store_true",
        help="print statistics for INPUT only; write no output (OUTPUT not needed)",
    )

    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}",
    )
    return parser


def _format_header(args: argparse.Namespace, mode: str) -> str:
    max_length = "no limit" if args.max_length is None else str(args.max_length)
    return "\n".join(
        (
            f"clean-fasta {__version__}",
            f"  mode            : {mode}",
            f"  type            : {'amino acid' if args.type == 'aa' else 'nucleic acid'}",
            f"  min length      : {args.min_length}",
            f"  max length      : {max_length}",
            f"  min valid ratio : {args.min_ratio}",
            f"  unique ids      : {'yes' if args.unique else 'no'}",
            f"  allow digits    : {'yes' if args.allow_digits else 'no'}",
        )
    )


def _status_line(args: argparse.Namespace, stats: CleanStats) -> str | None:
    """The trailing one-line summary of what was (or would have been) written."""
    if args.dry_run:
        if args.output not in (None, "-"):
            return (
                f"Dry run: {stats.passed} sequences would be written to "
                f"{args.output} (nothing written)."
            )
        return f"Dry run: {stats.passed} of {stats.total} sequences would pass (nothing written)."
    if args.stats_only:
        return (
            f"Statistics only: {stats.passed} of {stats.total} sequences "
            "pass the filters (nothing written)."
        )
    if args.output != "-":
        return f"Wrote {stats.passed} sequences to {args.output}"
    return None


def main(argv: Sequence[str] | None = None) -> int:
    """Run the command-line tool. Returns a process exit code."""
    parser = build_parser()
    args = parser.parse_args(argv)

    if not 0.0 <= args.min_ratio <= 1.0:
        parser.error("--min-ratio must be between 0.0 and 1.0")
    if args.min_length < 0:
        parser.error("--min-length must be >= 0")
    if args.max_length is not None and args.max_length <= args.min_length:
        parser.error("--max-length must be greater than --min-length")

    writing = not (args.dry_run or args.stats_only)
    if writing and args.output is None:
        parser.error(
            "the output argument is required "
            "(use '-' for stdout, or --dry-run / --stats-only to analyze without writing)"
        )
    if args.stats_only and args.quiet:
        parser.error("--quiet cannot be combined with --stats-only (there would be no output)")

    if args.dry_run:
        mode = "dry run (no output written)"
    elif args.stats_only:
        mode = "statistics only (no output written)"
    else:
        mode = "clean"

    try:
        stats: CleanStats = clean_fasta_file(
            args.input,
            args.output if args.output is not None else "-",
            min_length=args.min_length,
            max_length=args.max_length,
            min_ratio=args.min_ratio,
            is_aa=(args.type == "aa"),
            unique_ids=args.unique,
            allow_digits=args.allow_digits,
            wrap=args.wrap,
            force=args.force,
            write_output=writing,
        )
    except FileExistsError as exc:
        print(f"clean-fasta: error: {exc}", file=sys.stderr)
        return 2
    except FileNotFoundError as exc:
        print(f"clean-fasta: error: input file not found: {exc.filename}", file=sys.stderr)
        return 2
    except OSError as exc:
        print(f"clean-fasta: error: {exc}", file=sys.stderr)
        return 1

    if not args.quiet:
        print(_format_header(args, mode), file=sys.stderr)
        print(format_stats(stats), file=sys.stderr)
        status = _status_line(args, stats)
        if status is not None:
            print(status, file=sys.stderr)

    return 0


def run() -> None:
    """Console-script entry point."""
    sys.exit(main())


if __name__ == "__main__":
    run()
