"""
Utility helpers for reading/writing sequences used by the IntronAwaredExonDesigner workflow.
"""

from pathlib import Path


def load_sequence_from_file(path: Path) -> str:
    """Read a FASTA-like file and concatenate non-header lines."""
    if not path.exists():
        raise FileNotFoundError(f"Sequence file not found: {path}")

    lines = path.read_text().splitlines()
    return ''.join(line.strip() for line in lines if line and not line.startswith('>'))


def load_protein_sequence(path: Path) -> str:
    """Load an amino-acid sequence from FASTA or plain text."""
    return load_sequence_from_file(path)


def load_utr_sequence(path: Path) -> str:
    """Load a UTR sequence and normalize to uppercase RNA."""
    raw = load_sequence_from_file(path)
    return raw.upper().replace('T', 'U')


def rna_to_dna(seq: str) -> str:
    """Convert RNA tokens (U) to DNA tokens (T) while preserving case."""
    return seq.replace('U', 'T').replace('u', 't')


def get_default_utrs() -> tuple[str, str]:
    """Return default 5'/3' UTR sequences from the data templates, falling back when missing."""
    project_root = Path(__file__).resolve().parents[2]
    utr5_path = project_root / "data" / "utr_templates" / "5utr_templates.txt"
    utr3_path = project_root / "data" / "utr_templates" / "3utr_templates.txt"

    if utr5_path.exists() and utr3_path.exists():
        return load_utr_sequence(utr5_path), load_utr_sequence(utr3_path)

    # Fallback minimal templates
    return (
        "GGGAAAUAAGAGAGAAAAGAAGAGUAAGAAGAAAUAUAAGAGCCACC",
        "UGAA"
    )
