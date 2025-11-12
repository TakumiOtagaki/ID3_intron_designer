"""
ViennaRNA partition function helpers.

Provides wrappers around the ViennaRNA Python bindings to obtain
ensemble free energy values and base-pair probability matrices.
"""

from typing import List, Tuple
import numpy as np


class ViennaRNAPartition:
    """Thin wrapper around ViennaRNA fold compound utilities."""

    def __init__(self):
        try:
            import RNA  # type: ignore
        except ImportError as exc:  # pragma: no cover - depends on local install
            raise RuntimeError(
                "viennaRNA module not found. Install with `uv pip install ViennaRNA`."
            ) from exc
        self.RNA = RNA

    @staticmethod
    def _sanitize_sequence(sequence: str) -> str:
        return sequence.replace('T', 'U').upper()

    def _fold_compound(self, sequence: str):
        seq = self._sanitize_sequence(sequence)
        fc = self.RNA.fold_compound(seq)
        ss, mfe = fc.mfe()
        fc.exp_params_rescale(mfe)
        pf_result = fc.pf()
        if isinstance(pf_result, (list, tuple)):
            # ViennaRNA returns (structure, energy); keep the energy term
            energy = pf_result[-1]
        else:
            energy = pf_result
        return fc, float(energy)

    def ensemble_free_energy(self, sequence: str) -> float:
        """
        Compute ensemble free energy (kcal/mol) for the provided sequence.
        """
        _, energy = self._fold_compound(sequence)
        return energy  # kcal/mol

    def base_pair_probability_matrix(self, sequence: str) -> np.ndarray:
        """
        Compute the full base-pair probability matrix for the sequence.
        """
        fc, _ = self._fold_compound(sequence)
        raw = fc.bpp()  # tuple-of-tuples with 1-based indexing padding
        matrix = np.array(raw, dtype=np.float32)
        if matrix.ndim != 2:
            raise RuntimeError("Unexpected BPP matrix shape from ViennaRNA")
        # Drop the 0-index padding row/column
        matrix = matrix[1:, 1:]
        return matrix


def compute_intron_losses(
    vienna: ViennaRNAPartition,
    full_sequence: str,
    window_ranges: List[Tuple[int, int]],
    boundary_indices: List[int]
) -> Tuple[float, float]:
    """
    Compute structural losses: (-EFE) for intron windows and summed boundary BPP.

    Returns:
        (window_efe_loss, boundary_probability_loss)
    """
    efe_loss = 0.0
    sanitized = vienna._sanitize_sequence(full_sequence)
    for start, end in window_ranges:
        window_seq = sanitized[start:end]
        if not window_seq:
            continue
        efe = vienna.ensemble_free_energy(window_seq)
        efe_loss += -efe  # maximize EFE → minimize -EFE

    if boundary_indices:
        bpp_matrix = vienna.base_pair_probability_matrix(sanitized)
        boundary_loss = float(sum(bpp_matrix[idx, :].sum() for idx in boundary_indices))
    else:
        boundary_loss = 0.0

    return efe_loss, boundary_loss
