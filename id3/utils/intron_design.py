"""
Utilities for the IntronAwaredExonDesigner workflow.

This module parses multi-FASTA inputs (5'UTR / main / 3'UTR),
tracks exon/intron layouts, and provides helpers for rebuilding
full RNA sequences from optimized exon segments.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence, Tuple
import logging

from id3.utils.sequence_utils import rna_string_to_amino_acids

logger = logging.getLogger(__name__)


def _normalize_rna_with_case(sequence: str) -> str:
    """Convert DNA alphabets to RNA while preserving exon/intron casing."""
    mapping = {
        'A': 'A', 'a': 'a',
        'C': 'C', 'c': 'c',
        'G': 'G', 'g': 'g',
        'T': 'U', 't': 'u',
        'U': 'U', 'u': 'u'
    }
    cleaned = []
    for char in sequence.strip():
        if char in mapping:
            cleaned.append(mapping[char])
        else:
            raise ValueError(f"Unsupported nucleotide '{char}' in sequence: {sequence}")
    return ''.join(cleaned)


def _parse_multi_fasta(path: Path) -> Dict[str, str]:
    """Parse multi-FASTA file into a dict keyed by lowercase header."""
    entries: Dict[str, List[str]] = {}
    current_header = None

    with open(path, 'r') as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith('>'):
                current_header = line[1:].strip().lower()
                entries.setdefault(current_header, [])
            else:
                if current_header is None:
                    raise ValueError(f"FASTA sequence encountered before header in {path}")
                entries[current_header].append(line)

    joined = {header: ''.join(seq_parts) for header, seq_parts in entries.items()}
    required = ['5utr', 'main', '3utr']
    missing = [key for key in required if key not in joined]
    if missing:
        raise ValueError(f"Missing entries {missing} in multi-FASTA {path}")

    return joined


@dataclass(frozen=True)
class Segment:
    """Represents a contiguous exon or intron region within the main transcript."""
    start: int
    end: int   # exclusive
    is_exon: bool


class IntronAwaredExonDesignerContext:
    """
    Context for the IntronAwaredExonDesigner: map uppercase exons back into the full transcript.

    Attributes:
        utr5 / utr3: Fixed UTR sequences.
        main_sequence: Full main transcript (mixed case).
        exon_positions: Indices of designable nucleotides within main_sequence.
        intron_segments: List of contiguous intron spans.
    """

    def __init__(self, fasta_path: str, amino_acid_sequence: str):
        fasta = _parse_multi_fasta(Path(fasta_path))
        self.utr5 = _normalize_rna_with_case(fasta['5utr']).upper()
        self.utr3 = _normalize_rna_with_case(fasta['3utr']).upper()
        self.main_sequence = _normalize_rna_with_case(fasta['main'])
        self.amino_acid_sequence = amino_acid_sequence

        self.segments: List[Segment] = self._extract_segments()
        self.exon_positions = [idx for idx, base in enumerate(self.main_sequence) if base.isupper()]
        if not self.exon_positions:
            raise ValueError("No designable exon nucleotides (uppercase) were found in 'main' entry.")

        self.design_length = len(self.exon_positions)
        self.intron_segments = [seg for seg in self.segments if not seg.is_exon]

        self._warn_on_amino_acid_mismatch()

    def _extract_segments(self) -> List[Segment]:
        """Split main sequence into contiguous exon/intron segments."""
        segments: List[Segment] = []
        idx = 0
        while idx < len(self.main_sequence):
            is_exon = self.main_sequence[idx].isupper()
            start = idx
            while idx < len(self.main_sequence) and self.main_sequence[idx].isupper() == is_exon:
                idx += 1
            segments.append(Segment(start=start, end=idx, is_exon=is_exon))
        return segments

    def _warn_on_amino_acid_mismatch(self) -> None:
        """Warn if the initial exon sequence does not encode the provided amino acids."""
        exon_sequence = self.get_exon_sequence()
        translated = rna_string_to_amino_acids(exon_sequence)
        if translated != self.amino_acid_sequence:
            logger.warning(
                "Initial exon sequence translates to %s, which mismatches provided amino-acid sequence %s",
                translated,
                self.amino_acid_sequence
            )

    def get_exon_sequence(self) -> str:
        """Return the concatenated exon sequence (design space) as RNA (uppercase)."""
        return ''.join(self.main_sequence[idx].upper() for idx in self.exon_positions)

    def rebuild_main_with_exons(self, exon_design: str) -> str:
        """
        Merge a design exon sequence back into the main transcript (mixed intron case retained).
        """
        if len(exon_design) != self.design_length:
            raise ValueError(
                f"Exon design length {len(exon_design)} does not match expected {self.design_length}"
            )
        exon_design = exon_design.upper().replace('T', 'U')
        main_chars = list(self.main_sequence)
        for idx, pos in enumerate(self.exon_positions):
            main_chars[pos] = exon_design[idx]
        return ''.join(main_chars)

    def build_full_sequence(self, exon_design: str, uppercase: bool = True) -> str:
        """
        Construct the complete RNA sequence (UTR + main + UTR) for ViennaRNA evaluation.
        """
        main_seq = self.rebuild_main_with_exons(exon_design)
        full = f"{self.utr5}{main_seq}{self.utr3}"
        return full.upper() if uppercase else full

    def get_intron_window_ranges(self, upstream: int, downstream: int) -> List[Tuple[int, int]]:
        """
        Return absolute indices (relative to full sequence) for each intron window.
        """
        full_offset = len(self.utr5)
        full_length = full_offset + len(self.main_sequence) + len(self.utr3)
        windows: List[Tuple[int, int]] = []
        for segment in self.intron_segments:
            start = max(0, full_offset + segment.start - upstream)
            end = min(full_length, full_offset + segment.end + downstream)
            windows.append((start, end))
        return windows

    def get_boundary_indices(self, flank: int = 3) -> List[int]:
        """
        Return absolute indices for the first/last `flank` nucleotides within each intron segment.
        """
        full_offset = len(self.utr5)
        indices: List[int] = []
        for segment in self.intron_segments:
            intron_length = segment.end - segment.start
            if intron_length < flank * 2:
                logger.warning(
                    "Intron of length %d is shorter than 2*flank=%d; using available positions.",
                    intron_length,
                    flank * 2
                )
            front = set(range(segment.start, min(segment.start + flank, segment.end)))
            back_start = max(segment.start, segment.end - flank)
            back = set(range(back_start, segment.end))
            boundary_positions = sorted(front.union(back))
            for rel_idx in boundary_positions:
                indices.append(full_offset + rel_idx)
        return indices

    def describe(self) -> Dict[str, Sequence[int]]:
        """Return metadata summary for debugging/logging."""
        return {
            'utr5_length': len(self.utr5),
            'utr3_length': len(self.utr3),
            'main_length': len(self.main_sequence),
            'design_length': self.design_length,
            'num_introns': len(self.intron_segments)
        }
