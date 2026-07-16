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
import statistics
import sys
from collections.abc import Iterable, Iterator
from dataclasses import dataclass, field
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
class SeqSummary:
    """Descriptive statistics accumulated over one *set* of sequences.

    Two of these are gathered per run: one over every input record (the "before"
    picture) and one over just the records that passed every filter (the "after"
    picture). Only integer lengths are retained per record, so memory stays
    modest even for very large files; the sequences themselves are not held.
    """

    is_aa: bool = True
    count: int = 0
    total_residues: int = 0
    total_bad: int = 0
    total_gc: int = 0
    longest_len: int | None = None
    longest_name: str = ""
    shortest_len: int | None = None
    shortest_name: str = ""
    most_bad_count: int | None = None
    most_bad_name: str = ""
    least_bad_count: int | None = None
    least_bad_name: str = ""
    lengths: list[int] = field(default_factory=list)

    def add(self, name: str, length: int, bad: int, gc: int) -> None:
        """Fold one record's length / bad-char / GC counts into the summary."""
        self.count += 1
        self.total_residues += length
        self.total_bad += bad
        self.total_gc += gc
        self.lengths.append(length)
        if self.longest_len is None or length > self.longest_len:
            self.longest_len, self.longest_name = length, name
        if self.shortest_len is None or length < self.shortest_len:
            self.shortest_len, self.shortest_name = length, name
        if self.most_bad_count is None or bad > self.most_bad_count:
            self.most_bad_count, self.most_bad_name = bad, name
        if self.least_bad_count is None or bad < self.least_bad_count:
            self.least_bad_count, self.least_bad_name = bad, name

    @property
    def mean_length(self) -> float | None:
        return self.total_residues / self.count if self.count else None

    @property
    def median_length(self) -> float | None:
        return statistics.median(self.lengths) if self.lengths else None

    @property
    def length_stddev(self) -> float | None:
        if not self.lengths:
            return None
        return statistics.pstdev(self.lengths)

    @property
    def quartiles(self) -> tuple[float | None, float | None]:
        """(Q1, Q3) length quartiles, or (None, None) with fewer than 2 records."""
        if len(self.lengths) < 2:
            return (None, None)
        q1, _q2, q3 = statistics.quantiles(self.lengths, n=4, method="inclusive")
        return (q1, q3)

    def _n50_l50(self) -> tuple[int | None, int | None]:
        if not self.lengths:
            return (None, None)
        half = self.total_residues / 2
        cumulative = 0
        for count, length in enumerate(sorted(self.lengths, reverse=True), start=1):
            cumulative += length
            if cumulative >= half:
                return (length, count)
        return (self.lengths[-1], len(self.lengths))  # pragma: no cover

    @property
    def n50(self) -> int | None:
        return self._n50_l50()[0]

    @property
    def l50(self) -> int | None:
        return self._n50_l50()[1]

    @property
    def gc_content(self) -> float | None:
        """G+C fraction (0-1) for nucleic-acid input; None for amino acids."""
        if self.is_aa or not self.total_residues:
            return None
        return self.total_gc / self.total_residues

    @property
    def valid_ratio(self) -> float | None:
        """Overall fraction of valid (non-"bad") characters across the set."""
        if not self.total_residues:
            return None
        return (self.total_residues - self.total_bad) / self.total_residues


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
    is_aa: bool = True
    min_ratio_passing: float | None = None
    max_ratio_passing: float | None = None
    input_summary: SeqSummary = field(default_factory=SeqSummary)
    kept_summary: SeqSummary = field(default_factory=SeqSummary)

    @property
    def min_input_length(self) -> int | None:
        return self.input_summary.shortest_len

    @property
    def max_input_length(self) -> int | None:
        return self.input_summary.longest_len

    @property
    def sum_input_length(self) -> int:
        return self.input_summary.total_residues

    @property
    def average_input_length(self) -> float:
        """Mean length of all input sequences (0.0 when there were none)."""
        return self.input_summary.mean_length or 0.0

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
    stats.is_aa = is_aa
    stats.input_summary.is_aa = is_aa
    stats.kept_summary.is_aa = is_aa
    seen_ids: set[str] = set()

    for seq in sequences:
        stats.total += 1
        length = seq.get_length()

        name = _WS_RE.sub(" ", seq.get_seq_id()).strip()

        # Bad-character and GC counts, computed once for every input record so
        # they can feed both the input summary and (for survivors) the kept
        # summary. "Bad" is the complement of the valid-character ratio: gap /
        # unknown chars for amino acids, non-A/C/G/T for nucleic acids.
        if is_aa:
            bad = seq.count_irregular_chars_aa()
        else:
            bad = length - seq.count_regular_chars_na()
        gc = seq.count_gc()
        stats.input_summary.add(name, length, bad, gc)

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

        ratio = (length - bad) / length
        stats._record_ratio(ratio)
        if ratio < min_ratio:
            stats.ignored_irregular += 1
            continue

        if unique_ids and name in seen_ids:
            stats.ignored_duplicate_id += 1
            continue

        seen_ids.add(name)
        seq.set_seq_id(name)
        stats.kept_summary.add(name, length, bad, gc)
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
    write_output: bool = True,
) -> CleanStats:
    """Clean ``infile`` into ``outfile`` and return the run statistics.

    ``"-"`` may be used for ``infile``/``outfile`` to read from stdin / write to
    stdout. Existing output files are not overwritten unless ``force`` is true
    (raising :class:`FileExistsError` otherwise). ``wrap <= 0`` disables output
    line wrapping.

    When ``write_output`` is false the input is still fully parsed and filtered
    (so the returned statistics are complete) but nothing is written and
    ``outfile`` is not opened or checked -- this backs the ``--dry-run`` and
    ``--stats-only`` command-line modes.
    """
    if write_output and outfile != "-" and not force and Path(outfile).exists():
        raise FileExistsError(
            f"output file already exists: {outfile!r} (use force=True / --force to overwrite)"
        )

    stats = CleanStats()
    with contextlib.ExitStack() as stack:
        in_stream = sys.stdin if infile == "-" else stack.enter_context(
            open(infile, encoding=encoding)
        )
        out_stream = None
        if write_output:
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
            if out_stream is None:
                continue
            record = seq.to_fasta_wrapped(wrap) if wrap and wrap > 0 else seq.to_fasta()
            out_stream.write(record)
            out_stream.write("\n")

    return stats


def _format_extremes(title: str, summary: SeqSummary) -> list[str]:
    """Render the named-extreme lines (longest / shortest / bad-char) for one set."""

    def named(name: str, value: int | None) -> str:
        if value is None:
            return "(none)"
        return f"{name or '(unnamed)'} ({value})"

    rows = [
        ("longest", named(summary.longest_name, summary.longest_len)),
        ("shortest", named(summary.shortest_name, summary.shortest_len)),
        ("most bad chars", named(summary.most_bad_name, summary.most_bad_count)),
        ("least bad chars", named(summary.least_bad_name, summary.least_bad_count)),
    ]
    label_w = max(len(label) for label, _ in rows)
    return [title] + [f"  {label:<{label_w}} : {value}" for label, value in rows]


def format_stats(stats: CleanStats) -> str:
    """Render ``stats`` as a human-readable before/after statistics report."""
    inp = stats.input_summary
    kept = stats.kept_summary

    def as_int(value: int | None) -> str:
        return "n/a" if value is None else str(value)

    def as_float(value: float | None, places: int = 1) -> str:
        return "n/a" if value is None else f"{value:.{places}f}"

    def as_pct(value: float | None) -> str:
        return "n/a" if value is None else f"{value * 100:.1f}%"

    def as_ratio(value: float | None) -> str:
        return "n/a" if value is None else f"{value:.4f}"

    inp_q1, inp_q3 = inp.quartiles
    kept_q1, kept_q3 = kept.quartiles

    # (label, input value, kept value)
    rows = [
        ("sequences", as_int(inp.count), as_int(kept.count)),
        ("total residues", as_int(inp.total_residues), as_int(kept.total_residues)),
        ("length min", as_int(inp.shortest_len), as_int(kept.shortest_len)),
        ("length Q1", as_float(inp_q1), as_float(kept_q1)),
        ("length median", as_float(inp.median_length), as_float(kept.median_length)),
        ("length Q3", as_float(inp_q3), as_float(kept_q3)),
        ("length max", as_int(inp.longest_len), as_int(kept.longest_len)),
        ("length mean", as_float(inp.mean_length), as_float(kept.mean_length)),
        ("length stddev", as_float(inp.length_stddev), as_float(kept.length_stddev)),
        ("N50", as_int(inp.n50), as_int(kept.n50)),
        ("L50", as_int(inp.l50), as_int(kept.l50)),
        ("GC content", as_pct(inp.gc_content), as_pct(kept.gc_content)),
        ("valid-char ratio", as_ratio(inp.valid_ratio), as_ratio(kept.valid_ratio)),
    ]

    label_w = max(len(label) for label, _, _ in rows)
    in_w = max(len("Input"), *(len(v) for _, v, _ in rows))
    kept_w = max(len("Kept"), *(len(v) for _, _, v in rows))

    lines = [
        "Statistics (before / after filtering)",
        f"  {'':<{label_w}}   {'Input':>{in_w}}   {'Kept':>{kept_w}}",
    ]
    for label, in_v, kept_v in rows:
        lines.append(f"  {label:<{label_w}}   {in_v:>{in_w}}   {kept_v:>{kept_w}}")

    lines.append("")
    lines.extend(_format_extremes("Extremes (input)", inp))
    lines.append("")
    lines.extend(_format_extremes("Extremes (kept)", kept))

    filtered = [
        ("input", stats.total),
        ("ignored: empty id", stats.ignored_no_name),
        ("ignored: digits in sequence", stats.ignored_numbers),
        ("ignored: too short", stats.ignored_too_short),
        ("ignored: too long", stats.ignored_too_long),
        ("ignored: too many irregular chars", stats.ignored_irregular),
        ("ignored: duplicate id", stats.ignored_duplicate_id),
        ("passed", stats.passed),
    ]
    filt_w = max(len(label) for label, _ in filtered)
    lines.append("")
    lines.append("Filtering")
    lines.extend(f"  {label:<{filt_w}} : {value}" for label, value in filtered)

    return "\n".join(lines)
