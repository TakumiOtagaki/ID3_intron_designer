## Intron-aware exon design requirements

### Goal
Design the exon segments of a pre-mRNA sequence such that:
- Introns remain fixed (provided in lowercase) and exons are designable (uppercase).
- 5' and 3' UTRs are fixed and prepended/appended to any output sequence.
- Structural penalties derived from ViennaRNA discourage structure formation in and around the intron.

### Inputs
- Multi-FASTA file with three entries in order: `>5utr`, `>main`, `>3utr`.
  - `main` contains uppercase exons (designable) and lowercase introns (fixed).
  - Provided amino-acid sequence must match the concatenated exon coding region. Raise a warning if the initial exon RNA violates the target AA sequence.
- Initial exon RNA sequence (acts as starting point for optimization).

### Optimization requirements
1. **Intron design context**
   - Parse multi-FASTA into 5'UTR, main (exon/intron), 3'UTR.
   - Track exon segments, intron segments, and mapping needed to merge optimized exon output back into the full pre-mRNA sequence.
   - Provide helpers to:
     - Extract designable exon-only sequence.
     - Reconstruct the full RNA (UTR + exon/intron) from any candidate exon design.
     - Validate exon sequence against provided amino-acid sequence.

2. **ViennaRNA wrapper**
   - Use `viennaRNA` Python module (already installed) to expose:
     - Ensemble free energy (EFE) for the intron window with flanking exons: `[intron_start-60, intron_end+30)`.
     - Base-pair probabilities (BPP) for the full concatenated sequence (UTR + main).
   - Guard against missing ViennaRNA import and provide informative errors.

3. **Loss terms**
   - Only two loss terms are needed; DeepRaccess/accessibility is **not** part of this flow.
   - **Window EFE loss**: `loss_efe = -ensemble_free_energy(window)` to maximize EFE and discourage structure.
   - **Boundary BPP loss**: sum of pairing probabilities for the three nucleotides at each intron boundary (5' and 3' ends), computed using the full-sequence BPP matrix.
   - Parameters/weights must be configurable via CLI/config (e.g., `--efe-weight`, `--boundary-weight`).

4. **Integration**
   - Extend `demo.py` (and the unified experiment path) to:
  - Load the multi-FASTA input and initialize `IntronAwaredExonDesignerContext`.
     - Run constraint.forward with `beta=1` when computing ViennaRNA-based losses so that discrete sequences feed ViennaRNA, while STE keeps gradients flowing.
     - Combine the weighted EFE and boundary losses into the objective optimized by the existing loop.
   - Ensure the loss hooks are available to config-based workflows so experiments can enable/disable structural penalties consistently.

### Non-functional
- Keep new utilities under `id3/utils`.
- Maintain ASCII files and concise comments; warn (but do not crash) when starter exon RNA encodes the wrong amino-acid sequence.
- Provide extensible interfaces so future structural objectives can reuse the same context/wrapper.
