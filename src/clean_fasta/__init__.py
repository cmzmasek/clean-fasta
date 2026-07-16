"""clean_fasta -- clean and filter FASTA sequence files.

Public API::

    from clean_fasta import (
        MolSeq,
        CleanStats,
        stream_fasta,
        filter_sequences,
        clean_fasta_file,
    )
"""

from __future__ import annotations

__version__ = "2.1.0"

from clean_fasta.cleaner import (
    CleanStats,
    SeqSummary,
    clean_fasta_file,
    filter_sequences,
    format_stats,
    stream_fasta,
)
from clean_fasta.molseq import MolSeq

__all__ = [
    "MolSeq",
    "CleanStats",
    "SeqSummary",
    "stream_fasta",
    "filter_sequences",
    "clean_fasta_file",
    "format_stats",
    "__version__",
]
