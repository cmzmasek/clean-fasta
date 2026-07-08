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

"""Core FASTA parsing and filtering logic, kept free of command-line concerns.

The public surface is:

* :func:`stream_fasta`     -- parse a stream of lines into :class:`MolSeq` records.
* :func:`filter_sequences` -- filter an iterable of records, collecting statistics.
* :func:`clean_fasta_file` -- wire the two together with file/stream I/O.
* :class:`CleanStats`      -- counters gathered during a run.
* :func:`format_stats`     -- render a :class:`CleanStats` as a text report.
"""

from __future__ import annotations

import contextlib
import re
import sys
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from pathlib import Path

from clean_fasta.molseq import MolSeq

DEFAULT_MIN_LENGTH = 20
DEFAULT_MIN_RATIO = 0.99
DEFAULT_WRAP = 80

_ID_RE = re.compile(r">\s*(.+)")
_GAP_RE = re.compile(r"[-\s]+")
_WS_RE = re.compile(r"\s+")
_DIGIT_RE = re.compile(r"\d")


@dataclass
class CleanStats:
    """Counters describing a single cleaning run."""

    total: int = 0
    passed: int = 0
    ignored_no_name: int = 0
    ignored_numbers: int = 0
    ignored_too_short: int = 0
    ignored_too_long: int = 0
    ignored_irregular: int = 0
    ignored_duplicate_id: int = 0
    min_input_length: int | None = None
    max_input_length: int | None = None
    sum_input_length: int = 0
    min_ratio_passing: float | None = None
    max_ratio_passing: float | None = None

    @property
    def average_input_length(self) -> float:
        """Mean length of all input sequences (0.0 when there were none)."""
        return self.sum_input_length / self.total if self.total else 0.0

    def _record_input_length(self, length: int) -> None:
        if self.min_input_length is None or length < self.min_input_length:
            self.min_input_length = length
        if self.max_input_length is None or length > self.max_input_length:
            self.max_input_length = length
        self.sum_input_length += length

    def _record_ratio(self, ratio: float) -> None:
        if self.min_ratio_passing is None or ratio < self.min_ratio_passing:
            self.min_ratio_passing = ratio
        if self.max_ratio_passing is None or ratio > self.max_ratio_passing:
            self.max_ratio_passing = ratio


def stream_fasta(lines: Iterable[str], remove_gaps: bool = True) -> Iterator[MolSeq]:
    """Parse ``lines`` (any iterable of strings) into :class:`MolSeq` records.

    Blank lines are skipped. Within a sequence body, whitespace is always
    removed; when ``remove_gaps`` is true, gap characters (``-``) are removed as
    well. A header line with no identifier (e.g. a bare ``>``) yields a record
    with an empty id rather than raising.
    """
    seq_id: str | None = None
    seq: list[str] = []
    for raw in lines:
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if seq_id is not None:
                yield MolSeq(seq_id, "".join(seq))
            match = _ID_RE.search(line)
            seq_id = match.group(1) if match else ""
            seq = []
        else:
            line = (_GAP_RE if remove_gaps else _WS_RE).sub("", line)
            seq.append(line)
    if seq_id is not None:
        yield MolSeq(seq_id, "".join(seq))


def filter_sequences(
    sequences: Iterable[MolSeq],
    *,
    min_length: int = DEFAULT_MIN_LENGTH,
    max_length: int | None = None,
    min_ratio: float = DEFAULT_MIN_RATIO,
    is_aa: bool = True,
    unique_ids: bool = True,
    allow_digits: bool = False,
    stats: CleanStats | None = None,
) -> Iterator[MolSeq]:
    """Yield the sequences that pass every filter, updating ``stats`` in place.

    A record is dropped (and the corresponding counter incremented) when it has
    an empty id, contains digit characters (unless ``allow_digits``), is shorter
    than ``min_length`` (or empty), is longer than ``max_length`` (when set), has
    a valid-character ratio below ``min_ratio``, or (when ``unique_ids``) repeats
    an id already seen. The passing id is whitespace-collapsed before output.

    The valid-character ratio is the fraction of unambiguous nucleotides
    (A/C/G/T) for ``is_aa=False``, or the fraction of non-gap/unknown characters
    for ``is_aa=True``.
    """
    if stats is None:
        stats = CleanStats()
    seen_ids: set[str] = set()

    for seq in sequences:
        stats.total += 1
        length = seq.get_length()
        stats._record_input_length(length)

        name = _WS_RE.sub(" ", seq.get_seq_id()).strip()
        if not name:
            stats.ignored_no_name += 1
            continue
        if not allow_digits and _DIGIT_RE.search(seq.get_seq()):
            stats.ignored_numbers += 1
            continue
        if length == 0 or length < min_length:
            stats.ignored_too_short += 1
            continue
        if max_length is not None and length > max_length:
            stats.ignored_too_long += 1
            continue

        if is_aa:
            regular = length - seq.count_irregular_chars_aa()
        else:
            regular = seq.count_regular_chars_na()
        ratio = regular / length
        stats._record_ratio(ratio)
        if ratio < min_ratio:
            stats.ignored_irregular += 1
            continue

        if unique_ids and name in seen_ids:
            stats.ignored_duplicate_id += 1
            continue

        seen_ids.add(name)
        seq.set_seq_id(name)
        stats.passed += 1
        yield seq


def clean_fasta_file(
    infile: str,
    outfile: str,
    *,
    min_length: int = DEFAULT_MIN_LENGTH,
    max_length: int | None = None,
    min_ratio: float = DEFAULT_MIN_RATIO,
    is_aa: bool = True,
    unique_ids: bool = True,
    allow_digits: bool = False,
    wrap: int = DEFAULT_WRAP,
    force: bool = False,
    encoding: str = "utf-8",
) -> CleanStats:
    """Clean ``infile`` into ``outfile`` and return the run statistics.

    ``"-"`` may be used for ``infile``/``outfile`` to read from stdin / write to
    stdout. Existing output files are not overwritten unless ``force`` is true
    (raising :class:`FileExistsError` otherwise). ``wrap <= 0`` disables output
    line wrapping.
    """
    if outfile != "-" and not force and Path(outfile).exists():
        raise FileExistsError(
            f"output file already exists: {outfile!r} (use force=True / --force to overwrite)"
        )

    stats = CleanStats()
    with contextlib.ExitStack() as stack:
        in_stream = sys.stdin if infile == "-" else stack.enter_context(
            open(infile, encoding=encoding)
        )
        out_stream = sys.stdout if outfile == "-" else stack.enter_context(
            open(outfile, "w", encoding=encoding)
        )

        records = stream_fasta(in_stream, remove_gaps=True)
        for seq in filter_sequences(
            records,
            min_length=min_length,
            max_length=max_length,
            min_ratio=min_ratio,
            is_aa=is_aa,
            unique_ids=unique_ids,
            allow_digits=allow_digits,
            stats=stats,
        ):
            record = seq.to_fasta_wrapped(wrap) if wrap and wrap > 0 else seq.to_fasta()
            out_stream.write(record)
            out_stream.write("\n")

    return stats


def format_stats(stats: CleanStats) -> str:
    """Render ``stats`` as an aligned, human-readable multi-line report."""

    def num(value: int | None) -> str:
        return "n/a" if value is None else str(value)

    def ratio(value: float | None) -> str:
        return "n/a" if value is None else f"{value:.4f}"

    mean = f"{stats.average_input_length:.1f}" if stats.total else "n/a"

    rows = [
        ("Input sequences", str(stats.total)),
        ("  min length", num(stats.min_input_length)),
        ("  max length", num(stats.max_input_length)),
        ("  mean length", mean),
        ("Ignored: empty id", str(stats.ignored_no_name)),
        ("Ignored: digits in sequence", str(stats.ignored_numbers)),
        ("Ignored: too short", str(stats.ignored_too_short)),
        ("Ignored: too long", str(stats.ignored_too_long)),
        ("Ignored: too many irregular chars", str(stats.ignored_irregular)),
        ("Ignored: duplicate id", str(stats.ignored_duplicate_id)),
        ("  min passing ratio", ratio(stats.min_ratio_passing)),
        ("  max passing ratio", ratio(stats.max_ratio_passing)),
        ("Passed", str(stats.passed)),
    ]
    label_width = max(len(label) for label, _ in rows)
    return "\n".join(f"{label.ljust(label_width)} : {value}" for label, value in rows)
