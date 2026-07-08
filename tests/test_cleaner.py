"""Tests for the core parsing and filtering logic."""

import io

import pytest

from clean_fasta.cleaner import (
    CleanStats,
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


# --------------------------------------------------------------------------- #
# format_stats
# --------------------------------------------------------------------------- #

def test_format_stats_handles_empty_run():
    report = format_stats(CleanStats())
    lines = report.splitlines()
    assert any(line.startswith("Input sequences") and line.endswith(": 0") for line in lines)
    assert "mean length" in report
    assert "n/a" in report  # min/max length and ratios are n/a


def test_format_stats_reports_passed():
    stats = CleanStats(total=3, passed=2)
    assert "Passed" in format_stats(stats)
