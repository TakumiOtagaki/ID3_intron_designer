#!/usr/bin/env python3
"""
Lightweight CLI for the ID3 IntronAwaredExonDesigner demo workflows.

This script dispatches either the DeepRaccess accessibility demo or the intron-aware
structural sampler depending on the provided arguments.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

import argparse
import torch
import yaml

from id3.apps.accessibility_runner import run_accessibility_optimization
from id3.apps.exon_driven_structural import run_intron_awared_exon_structural_optimization


def apply_config_overrides(args):
    """Override argparse namespace with values from YAML config, if provided."""
    config_path = getattr(args, 'config', None)
    if not config_path:
        return args

    path = Path(config_path)
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {path}")

    with open(path, 'r') as handle:
        data = yaml.safe_load(handle) or {}

    if not isinstance(data, dict):
        raise ValueError(f"Config file {path} must contain a key-value mapping.")

    for key, value in data.items():
        attr = key.replace('-', '_')
        if hasattr(args, attr):
            setattr(args, attr, value)
        else:
            raise ValueError(f"Unknown config key '{key}' in {path}")

    return args


def main():
    parser = argparse.ArgumentParser(
        description='ID3 IntronAwaredExonDesigner demo driver',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python demo.py --protein MSKGEELFTGVVPILVELDGDVNGHKFSVSGEG
  python demo.py --structure-fasta data/intron_examples/egfp_with_intron.fasta --window-upstream 60
  python demo.py --config id3/config/intron_design_example.yaml
        """
    )

    parser.add_argument(
        '--protein-seq', '--protein',
        type=str,
        default='MSKGEELFTGVVPILVELDGDVNGHKFSVSGEG',
        help='Protein sequence (amino acids)'
    )
    parser.add_argument('--protein-file', type=str, help='Load protein from FASTA file')
    parser.add_argument('--config', type=str, help='YAML config file whose keys override CLI arguments')
    parser.add_argument(
        '--constraint',
        type=str,
        choices=['lagrangian', 'amino_matching', 'codon_profile'],
        default='lagrangian',
        help='Constraint mechanism'
    )
    parser.add_argument(
        '--structure-fasta',
        type=str,
        help="Multi-FASTA file with '>5utr', '>main', '>3utr' for the IntronAwaredExonDesigner structural demo"
    )
    parser.add_argument('--utr5-file', type=str, help='5\' UTR sequence file (default: data/utr_templates/5utr_templates.txt)')
    parser.add_argument('--utr3-file', type=str, help='3\' UTR sequence file (default: data/utr_templates/3utr_templates.txt)')
    parser.add_argument(
        '--mode',
        type=str,
        choices=['det.soft', 'det.hard', 'sto.soft', 'sto.hard'],
        default='sto.soft',
        help='Operational mode (see paper for details)'
    )
    parser.add_argument('--iterations', type=int, default=20, help='Optimization iterations')
    parser.add_argument('--learning-rate', type=float, default=0.01, help='Learning rate for optimizer')
    parser.add_argument('--cai-target', type=float, default=0.8, help='Target CAI value')
    parser.add_argument('--cai-weight', type=float, default=0.1, help='CAI weight in loss')
    parser.add_argument('--efe-weight', type=float, default=1.0, help='Weight for intron window -EFE loss')
    parser.add_argument('--boundary-weight', type=float, default=1.0, help='Weight for intron boundary pairing loss')
    parser.add_argument('--window-upstream', type=int, default=60, help='Upstream window size for EFE')
    parser.add_argument('--window-downstream', type=int, default=30, help='Downstream window size for EFE')
    parser.add_argument('--boundary-flank', type=int, default=3, help='Flank size for boundary BPP loss')
    parser.add_argument('--structure-output', type=str, help='Path to save optimized UTR/main multi-FASTA')
    parser.add_argument('--sample-count', type=int, default=0, help='Number of sampled exon sequences to save')
    parser.add_argument('--device', type=str, default='cuda' if torch.cuda.is_available() else 'cpu', help='Computation device')
    parser.add_argument('--deepraccess-model', type=str, help='Path to DeepRaccess model (optional)')
    parser.add_argument('--output', '-o', type=str, help='Output FASTA for optimized RNA sequence')
    parser.add_argument('--save-result', type=str, help='Save detailed optimization log to JSON')
    parser.add_argument('--verbose', '-v', action='store_true', help='Verbose output')

    args = parser.parse_args()
    args = apply_config_overrides(args)

    try:
        if args.structure_fasta:
            run_intron_awared_exon_structural_optimization(args)
        else:
            run_accessibility_optimization(args)
    except Exception as exc:
        print(f"\n❌ Error: {exc}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
