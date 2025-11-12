"""
I/O helpers for intron-aware design workflows.
"""

from pathlib import Path
from typing import Iterable, Tuple


def _wrap_sequence(sequence: str, width: int = 60) -> Iterable[str]:
    for i in range(0, len(sequence), width):
        yield sequence[i:i + width]


def write_intron_multifasta(
    output_path: str,
    utr5: str,
    main_sequence: str,
    utr3: str,
    headers: Tuple[str, str, str] = ("5utr", "main", "3utr")
) -> Path:
    """
    Write a multi-FASTA file containing 5'UTR, main (exon/intron mix), and 3'UTR sequences.
    """
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w') as handle:
        for header, sequence in zip(headers, (utr5, main_sequence, utr3)):
            handle.write(f">{header}\n")
            for chunk in _wrap_sequence(sequence):
                handle.write(chunk + "\n")
    return path
