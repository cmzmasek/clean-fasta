"""Tests for the command-line interface."""

import pytest

from clean_fasta import __version__
from clean_fasta.cli import main


def write_fasta(tmp_path, text):
    path = tmp_path / "in.fasta"
    path.write_text(text)
    return path


def test_main_success_writes_output(tmp_path):
    infile = write_fasta(tmp_path, ">a\nAC--GT ACGTACGT\n")
    outfile = tmp_path / "out.fasta"

    code = main([str(infile), str(outfile), "-m", "5"])

    assert code == 0
    assert outfile.read_text() == ">a\nACGTACGTACGT\n"


def test_main_reports_to_stderr_not_stdout(tmp_path, capsys):
    infile = write_fasta(tmp_path, ">a\nACGTACGTACGT\n")
    outfile = tmp_path / "out.fasta"

    main([str(infile), str(outfile), "-m", "1"])

    captured = capsys.readouterr()
    assert captured.out == ""
    assert "clean-fasta" in captured.err
    assert "Passed" in captured.err


def test_main_quiet_suppresses_report(tmp_path, capsys):
    infile = write_fasta(tmp_path, ">a\nACGTACGTACGT\n")
    outfile = tmp_path / "out.fasta"

    main([str(infile), str(outfile), "-m", "1", "-q"])
    assert capsys.readouterr().err == ""


def test_main_stdout_output(tmp_path, capsys):
    infile = write_fasta(tmp_path, ">a\nACGTACGTACGT\n")

    code = main([str(infile), "-", "-m", "1", "-q"])

    assert code == 0
    assert capsys.readouterr().out == ">a\nACGTACGTACGT\n"


def test_main_existing_output_returns_2(tmp_path, capsys):
    infile = write_fasta(tmp_path, ">a\nACGTACGT\n")
    outfile = tmp_path / "out.fasta"
    outfile.write_text("keep me")

    code = main([str(infile), str(outfile), "-m", "1"])

    assert code == 2
    assert outfile.read_text() == "keep me"
    assert "already exists" in capsys.readouterr().err


def test_main_force_overwrites(tmp_path):
    infile = write_fasta(tmp_path, ">a\nACGTACGT\n")
    outfile = tmp_path / "out.fasta"
    outfile.write_text("old")

    code = main([str(infile), str(outfile), "-m", "1", "--force"])
    assert code == 0
    assert outfile.read_text() == ">a\nACGTACGT\n"


def test_main_missing_input_returns_2(tmp_path, capsys):
    outfile = tmp_path / "out.fasta"
    code = main([str(tmp_path / "nope.fasta"), str(outfile), "-m", "1"])
    assert code == 2
    assert "not found" in capsys.readouterr().err


def test_main_na_type(tmp_path):
    infile = write_fasta(tmp_path, ">a\nACGTACGTNN\n")
    outfile = tmp_path / "out.fasta"

    # 8/10 valid; default ratio 0.99 would drop it, low ratio keeps it.
    code = main([str(infile), str(outfile), "-t", "na", "-m", "1", "-r", "0.5"])
    assert code == 0
    assert outfile.read_text() == ">a\nACGTACGTNN\n"


def test_main_no_unique_flag(tmp_path):
    infile = write_fasta(tmp_path, ">dup\nACGTACGT\n>dup\nTTTTTTTT\n")
    outfile = tmp_path / "out.fasta"

    code = main([str(infile), str(outfile), "-m", "1", "--no-unique"])
    assert code == 0
    assert outfile.read_text().count(">dup") == 2


@pytest.mark.parametrize(
    "extra",
    [
        ["-r", "1.5"],        # ratio out of range
        ["-m", "-1"],         # negative min length
        ["-M", "5", "-m", "10"],  # max <= min
    ],
)
def test_main_invalid_arguments_exit_2(tmp_path, extra):
    infile = write_fasta(tmp_path, ">a\nACGT\n")
    outfile = tmp_path / "out.fasta"
    with pytest.raises(SystemExit) as exc:
        main([str(infile), str(outfile), *extra])
    assert exc.value.code == 2


def test_version(capsys):
    with pytest.raises(SystemExit) as exc:
        main(["--version"])
    assert exc.value.code == 0
    assert __version__ in capsys.readouterr().out
