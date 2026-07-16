"""Tests for the MolSeq value object."""

from clean_fasta.molseq import MolSeq


def test_id_and_seq_are_stripped():
    seq = MolSeq("  my id  ", "  ACGT\n")
    assert seq.get_seq_id() == "my id"
    assert seq.get_seq() == "ACGT"


def test_length_and_dunder_len():
    seq = MolSeq("x", "ACGTA")
    assert seq.get_length() == 5
    assert len(seq) == 5


def test_set_seq_id():
    seq = MolSeq("old", "ACGT")
    seq.set_seq_id("new")
    assert seq.get_seq_id() == "new"


def test_to_fasta():
    assert MolSeq("id", "ACGT").to_fasta() == ">id\nACGT"


def test_to_fasta_wrapped():
    seq = MolSeq("id", "AAAAAAAAAA")  # 10 chars
    assert seq.to_fasta_wrapped(4) == ">id\nAAAA\nAAAA\nAA"


def test_count_regular_chars_na_is_case_insensitive():
    # A C G T a c g t -> 8 regular; N and n are not regular
    assert MolSeq("x", "ACGTacgtNn").count_regular_chars_na() == 8


def test_count_irregular_chars_aa():
    # irregular set is: _ - ? X x * .
    assert MolSeq("x", "MK_-?Xx*.").count_irregular_chars_aa() == 7


def test_count_gc_is_case_insensitive():
    # C c G g -> 4 GC bases; A, T, N do not count.
    assert MolSeq("x", "ACGTacgtNn").count_gc() == 4
    assert MolSeq("x", "AATT").count_gc() == 0


def test_equality_and_repr():
    assert MolSeq("id", "ACGT") == MolSeq("id", "ACGT")
    assert MolSeq("id", "ACGT") != MolSeq("id", "ACGA")
    assert MolSeq("id", "ACGT") != "not a molseq"
    assert repr(MolSeq("id", "ACGT")) == "MolSeq('id', 'ACGT')"


def test_str_wraps_at_60():
    body = "A" * 130
    text = str(MolSeq("id", body))
    lines = text.splitlines()
    assert lines[0] == ">id"
    assert [len(line) for line in lines[1:]] == [60, 60, 10]
