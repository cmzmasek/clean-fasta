"""Tests for the core parsing and filtering logic."""

import io

import pytest

from clean_fasta.cleaner import (
    CleanStats,
    SeqSummary,
    clean_fasta_file,
    filter_sequences,
    format_stats,
    stream_fasta,
)
from clean_fasta.molseq import MolSeq


def parse(text, remove_gaps=True):
    return list(stream_fasta(io.StringIO(text), remove_gaps=remove_gaps))


# --------------------------------------------------------------------------- #
# stream_fasta
# --------------------------------------------------------------------------- #

def test_stream_basic_two_records():
    records = parse(">a\nACGT\n>b\nTTTT\n")
    assert records == [MolSeq("a", "ACGT"), MolSeq("b", "TTTT")]


def test_stream_joins_multiline_body_and_skips_blank_lines():
    records = parse(">a\nACGT\n\nTTTT\n\n")
    assert records == [MolSeq("a", "ACGTTTTT")]


def test_stream_removes_gaps_by_default():
    (record,) = parse(">a\nAC--GT NN\n")
    assert record.get_seq() == "ACGTNN"


def test_stream_keeps_dashes_when_gaps_not_removed():
    (record,) = parse(">a\nAC--GT NN\n", remove_gaps=False)
    assert record.get_seq() == "AC--GTNN"


def test_stream_collapses_header_whitespace_via_molseq():
    (record,) = parse(">   my id  \nACGT\n")
    assert record.get_seq_id() == "my id"


def test_stream_bare_header_yields_empty_id_without_error():
    records = parse(">\nACGT\n>b\nTTTT\n")
    assert records[0].get_seq_id() == ""
    assert records[0].get_seq() == "ACGT"
    assert records[1] == MolSeq("b", "TTTT")


def test_stream_no_trailing_newline():
    records = parse(">a\nACGT")
    assert records == [MolSeq("a", "ACGT")]


def test_stream_empty_input():
    assert parse("") == []


# --------------------------------------------------------------------------- #
# filter_sequences
# --------------------------------------------------------------------------- #

def test_filter_passes_clean_sequence():
    stats = CleanStats()
    out = list(filter_sequences([MolSeq("a", "ACGT" * 10)], min_length=5, stats=stats))
    assert [s.get_seq_id() for s in out] == ["a"]
    assert stats.passed == 1
    assert stats.total == 1


def test_filter_drops_too_short():
    stats = CleanStats()
    out = list(filter_sequences([MolSeq("a", "ACGT")], min_length=10, stats=stats))
    assert out == []
    assert stats.ignored_too_short == 1


def test_filter_drops_too_long():
    stats = CleanStats()
    out = list(
        filter_sequences([MolSeq("a", "A" * 50)], min_length=1, max_length=20, stats=stats)
    )
    assert out == []
    assert stats.ignored_too_long == 1


def test_filter_drops_empty_id():
    stats = CleanStats()
    out = list(filter_sequences([MolSeq("", "ACGT" * 10)], min_length=1, stats=stats))
    assert out == []
    assert stats.ignored_no_name == 1


def test_filter_drops_sequences_with_digits_by_default():
    stats = CleanStats()
    out = list(filter_sequences([MolSeq("a", "ACGT123")], min_length=1, stats=stats))
    assert out == []
    assert stats.ignored_numbers == 1


def test_filter_allow_digits_keeps_them():
    # Digits count as irregular for NA (8 of 10 chars valid -> 0.8), so use a low
    # ratio to isolate the effect of --allow-digits.
    stats = CleanStats()
    out = list(
        filter_sequences(
            [MolSeq("a", "ACGTACGT12")],
            min_length=1,
            is_aa=False,
            allow_digits=True,
            min_ratio=0.5,
            stats=stats,
        )
    )
    assert [s.get_seq_id() for s in out] == ["a"]
    assert stats.ignored_numbers == 0


def test_filter_aa_ratio_rejects_too_many_irregular():
    # 2 irregular ('X','X') out of 10 -> ratio 0.8 < 0.99
    stats = CleanStats()
    out = list(filter_sequences([MolSeq("a", "MKLMKLMKXX")], min_length=1, stats=stats))
    assert out == []
    assert stats.ignored_irregular == 1
    assert stats.min_ratio_passing == pytest.approx(0.8)


def test_filter_na_ratio_counts_acgt():
    # 8 of 10 are ACGT -> 0.8
    stats = CleanStats()
    out = list(
        filter_sequences(
            [MolSeq("a", "ACGTACGTNN")], min_length=1, is_aa=False, min_ratio=0.75, stats=stats
        )
    )
    assert [s.get_seq_id() for s in out] == ["a"]
    assert stats.max_ratio_passing == pytest.approx(0.8)


def test_filter_deduplicates_ids():
    stats = CleanStats()
    seqs = [MolSeq("dup", "ACGTACGT"), MolSeq("dup", "TTTTTTTT")]
    out = list(filter_sequences(seqs, min_length=1, stats=stats))
    assert len(out) == 1
    assert stats.ignored_duplicate_id == 1


def test_filter_allows_duplicate_ids_when_unique_off():
    stats = CleanStats()
    seqs = [MolSeq("dup", "ACGTACGT"), MolSeq("dup", "TTTTTTTT")]
    out = list(filter_sequences(seqs, min_length=1, unique_ids=False, stats=stats))
    assert len(out) == 2
    assert stats.ignored_duplicate_id == 0


def test_filter_input_length_stats():
    stats = CleanStats()
    list(
        filter_sequences(
            [MolSeq("a", "AA"), MolSeq("b", "AAAAA"), MolSeq("c", "AAA")],
            min_length=1,
            stats=stats,
        )
    )
    assert stats.min_input_length == 2
    assert stats.max_input_length == 5
    assert stats.average_input_length == pytest.approx(10 / 3)


def test_filter_collapses_whitespace_in_output_id():
    (out,) = list(filter_sequences([MolSeq("my    long   id", "ACGTACGT")], min_length=1))
    assert out.get_seq_id() == "my long id"


# --------------------------------------------------------------------------- #
# SeqSummary
# --------------------------------------------------------------------------- #

def _summary(lengths, is_aa=True):
    s = SeqSummary(is_aa=is_aa)
    for i, length in enumerate(lengths):
        s.add(f"seq{i}", length, bad=0, gc=0)
    return s


def test_summary_empty_metrics_are_none():
    s = SeqSummary()
    assert s.count == 0
    assert s.median_length is None
    assert s.mean_length is None
    assert s.length_stddev is None
    assert s.n50 is None
    assert s.l50 is None
    assert s.quartiles == (None, None)
    assert s.valid_ratio is None


def test_summary_basic_length_metrics():
    s = _summary([2, 3, 4, 5, 6])
    assert s.count == 5
    assert s.total_residues == 20
    assert s.shortest_len == 2
    assert s.longest_len == 6
    assert s.mean_length == pytest.approx(4.0)
    assert s.median_length == pytest.approx(4.0)
    assert s.length_stddev == pytest.approx(2**0.5)
    assert s.quartiles == (pytest.approx(3.0), pytest.approx(5.0))


def test_summary_n50_l50():
    # lengths 6,5,4,3,2 (total 20, half 10): cumulative reaches 10 at the 2nd seq.
    s = _summary([2, 3, 4, 5, 6])
    assert s.n50 == 5
    assert s.l50 == 2


def test_summary_gc_content_only_for_na():
    na = SeqSummary(is_aa=False)
    na.add("x", length=10, bad=0, gc=5)
    assert na.gc_content == pytest.approx(0.5)

    aa = SeqSummary(is_aa=True)
    aa.add("x", length=10, bad=0, gc=5)
    assert aa.gc_content is None


def test_summary_valid_ratio_and_bad_extremes():
    s = SeqSummary()
    s.add("clean", length=10, bad=0, gc=0)
    s.add("dirty", length=10, bad=4, gc=0)
    assert s.valid_ratio == pytest.approx((20 - 4) / 20)
    assert (s.most_bad_name, s.most_bad_count) == ("dirty", 4)
    assert (s.least_bad_name, s.least_bad_count) == ("clean", 0)


def test_filter_populates_input_and_kept_summaries():
    stats = CleanStats()
    seqs = [MolSeq("keep", "ACGT" * 10), MolSeq("drop", "AC")]  # 40 and 2
    list(filter_sequences(seqs, min_length=10, stats=stats))

    assert stats.input_summary.count == 2
    assert stats.input_summary.shortest_len == 2
    assert stats.input_summary.longest_len == 40
    # Only the 40-mer survives, so the kept summary sees just that one.
    assert stats.kept_summary.count == 1
    assert stats.kept_summary.shortest_len == 40
    assert stats.kept_summary.longest_len == 40


def test_filter_summary_bad_chars_track_type():
    # NA: 'NN' are the two bad chars in a 10-mer -> most bad = 2.
    stats = CleanStats()
    list(
        filter_sequences(
            [MolSeq("a", "ACGTACGTNN")],
            min_length=1,
            is_aa=False,
            min_ratio=0.0,
            stats=stats,
        )
    )
    assert stats.input_summary.most_bad_count == 2
    assert stats.input_summary.is_aa is False


# --------------------------------------------------------------------------- #
# clean_fasta_file
# --------------------------------------------------------------------------- #

def test_clean_fasta_file_end_to_end(tmp_path):
    infile = tmp_path / "in.fasta"
    outfile = tmp_path / "out.fasta"
    infile.write_text(">a\nAC--GT ACGTACGT\n>short\nAC\n>a\nTTTTTTTTTTTT\n")

    stats = clean_fasta_file(str(infile), str(outfile), min_length=5)

    assert stats.passed == 1
    assert stats.ignored_too_short == 1
    assert stats.ignored_duplicate_id == 1
    assert outfile.read_text() == ">a\nACGTACGTACGT\n"


def test_clean_fasta_file_refuses_to_overwrite(tmp_path):
    infile = tmp_path / "in.fasta"
    outfile = tmp_path / "out.fasta"
    infile.write_text(">a\nACGTACGT\n")
    outfile.write_text("do not clobber")

    with pytest.raises(FileExistsError):
        clean_fasta_file(str(infile), str(outfile), min_length=1)
    assert outfile.read_text() == "do not clobber"


def test_clean_fasta_file_force_overwrites(tmp_path):
    infile = tmp_path / "in.fasta"
    outfile = tmp_path / "out.fasta"
    infile.write_text(">a\nACGTACGT\n")
    outfile.write_text("old")

    clean_fasta_file(str(infile), str(outfile), min_length=1, force=True)
    assert outfile.read_text() == ">a\nACGTACGT\n"


def test_clean_fasta_file_wrap_zero_disables_wrapping(tmp_path):
    infile = tmp_path / "in.fasta"
    outfile = tmp_path / "out.fasta"
    infile.write_text(">a\n" + "A" * 200 + "\n")

    clean_fasta_file(str(infile), str(outfile), min_length=1, wrap=0)
    assert outfile.read_text() == ">a\n" + "A" * 200 + "\n"


def test_clean_fasta_file_stdout(tmp_path, capsys):
    infile = tmp_path / "in.fasta"
    infile.write_text(">a\nACGTACGT\n")

    clean_fasta_file(str(infile), "-", min_length=1)
    assert capsys.readouterr().out == ">a\nACGTACGT\n"


def test_clean_fasta_file_write_output_false_gathers_stats_without_writing(tmp_path):
    infile = tmp_path / "in.fasta"
    outfile = tmp_path / "out.fasta"
    infile.write_text(">a\nACGTACGT\n>short\nAC\n")

    stats = clean_fasta_file(str(infile), str(outfile), min_length=5, write_output=False)

    assert stats.total == 2
    assert stats.passed == 1
    assert not outfile.exists()  # nothing written


def test_clean_fasta_file_write_output_false_ignores_existing_output(tmp_path):
    infile = tmp_path / "in.fasta"
    outfile = tmp_path / "out.fasta"
    infile.write_text(">a\nACGTACGT\n")
    outfile.write_text("keep me")

    # No FileExistsError, and the existing file is left untouched.
    clean_fasta_file(str(infile), str(outfile), min_length=1, write_output=False)
    assert outfile.read_text() == "keep me"


# --------------------------------------------------------------------------- #
# format_stats
# --------------------------------------------------------------------------- #

def test_format_stats_handles_empty_run():
    report = format_stats(CleanStats())
    assert "sequences" in report
    assert "length median" in report
    assert "N50" in report
    assert "(none)" in report  # extremes have no records to report
    assert "n/a" in report  # numeric metrics are n/a with no input


def test_format_stats_reports_passed():
    stats = CleanStats(total=3, passed=2)
    report = format_stats(stats)
    assert "Filtering" in report
    assert "passed" in report


def test_format_stats_has_two_columns_and_extremes():
    stats = CleanStats()
    list(
        filter_sequences(
            [MolSeq("long", "ACGT" * 25), MolSeq("shorty", "ACGTAC")],
            min_length=10,
            stats=stats,
        )
    )
    report = format_stats(stats)
    assert "Input" in report and "Kept" in report
    assert "Extremes (input)" in report
    assert "Extremes (kept)" in report
    # 'long' is the longest input and the only survivor; 'shorty' the shortest input.
    assert "long (100)" in report
    assert "shorty (6)" in report
