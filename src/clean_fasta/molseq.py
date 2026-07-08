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

"""A minimal representation of a molecular (biological) sequence."""

from __future__ import annotations

import textwrap

#: Characters treated as "irregular" for amino-acid sequences: gaps, unknowns,
#: stop codons, and separators. Spelled out explicitly to preserve the exact
#: historical counting behaviour of the original tool.
_IRREGULAR_AA_CHARS = ("_", "-", "?", "X", "x", "*", ".")

#: Regular (unambiguous) nucleotide characters, both cases.
_REGULAR_NA_CHARS = ("a", "A", "c", "C", "g", "G", "t", "T")


class MolSeq:
    """A molecular sequence with an identifier.

    Both the identifier and the sequence are coerced to ``str`` and stripped of
    surrounding whitespace on construction.

    Attributes:
        seq_id: The sequence identifier or name.
        seq: The molecular sequence itself.
    """

    __slots__ = ("_seq_id", "_seq")

    def __init__(self, seq_id: object, seq: object) -> None:
        self._seq_id = str(seq_id).strip()
        self._seq = str(seq).strip()

    def get_seq_id(self) -> str:
        return self._seq_id

    def set_seq_id(self, seq_id: object) -> None:
        self._seq_id = str(seq_id)

    def get_seq(self) -> str:
        return self._seq

    def get_length(self) -> int:
        return len(self._seq)

    def to_fasta(self) -> str:
        """Return the sequence as a single-line-body FASTA record (no wrapping)."""
        return f">{self._seq_id}\n{self._seq}"

    def to_fasta_wrapped(self, width: int) -> str:
        """Return the sequence as a FASTA record with the body wrapped at ``width``."""
        wrapped = textwrap.fill(self._seq, width=width)
        return f">{self._seq_id}\n{wrapped}"

    def count_regular_chars_na(self) -> int:
        """Count unambiguous nucleotide characters (A, C, G, T; both cases)."""
        return sum(self._seq.count(ch) for ch in _REGULAR_NA_CHARS)

    def count_irregular_chars_aa(self) -> int:
        """Count gap / unknown / stop characters (``_ - ? X x * .``)."""
        return sum(self._seq.count(ch) for ch in _IRREGULAR_AA_CHARS)

    def __str__(self) -> str:
        return self.to_fasta_wrapped(60)

    def __len__(self) -> int:
        return self.get_length()

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, MolSeq):
            return NotImplemented
        return self._seq_id == other._seq_id and self._seq == other._seq

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self._seq_id!r}, {self._seq!r})"


if __name__ == "__main__":
    example = MolSeq(" abcd ", " acgttgtca")
    print(example.to_fasta())
    print(example.get_length())
    print(str(example))
    print(repr(example))
