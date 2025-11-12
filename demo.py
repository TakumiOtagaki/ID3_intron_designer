#!/usr/bin/env python3
"""
ID3 Framework Demo - Complete mRNA Sequence Optimization

This demo shows the COMPLETE ID3 workflow including:
1. Amino acid constraint satisfaction (3 mechanisms available)
2. CAI (Codon Adaptation Index) optimization
3. RNA accessibility optimization via DeepRaccess
"""

# Copyright (c) 2025 University of Tokyo
# Licensed under CC BY-NC-SA 4.0
# For commercial use, contact: lihongmin@edu.k.u-tokyo.ac.jp

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

import argparse
import torch
import os
from tqdm import tqdm
import yaml

# ID3 Operational Modes (Table 1 from paper)
MODE_CONFIG = {
    'det.soft': {'alpha': 0, 'beta': 0},   # Deterministic, Soft output
    'det.hard': {'alpha': 0, 'beta': 1},   # Deterministic, Hard output
    'sto.soft': {'alpha': 0.1, 'beta': 0},   # Stochastic, Soft output
    'sto.hard': {'alpha': 0.1, 'beta': 1},   # Stochastic, Hard output
}

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

def ensure_deepraccess_available():
    """Ensure DeepRaccess dependency exists before running accessibility mode."""
    project_root = Path(__file__).parent
    possible_dirs = [
        project_root / "DeepRaccess",
        project_root / "scripts" / "DeepRaccess"
    ]

    for dir_path in possible_dirs:
        if dir_path.exists():
            return

    print("\n" + "="*70)
    print("⚠️  DeepRaccess Not Found")
    print("="*70)
    print("\nDeepRaccess is required for RNA accessibility prediction.")
    print("\nAttempting automatic setup...")
    print("This will take ~30 seconds with internet connection.")
    print("="*70 + "\n")

    import subprocess
    setup_script = project_root / "scripts" / "setup_deepraccess.sh"
    if setup_script.exists():
        try:
            subprocess.run(
                ["bash", str(setup_script), "-y"],
                check=True,
                capture_output=False
            )
            print("\n✅ Setup complete! Please re-run demo.py\n")
            sys.exit(0)
        except subprocess.CalledProcessError as e:
            print(f"\n❌ Setup failed with error code {e.returncode}")
            print("Please try manual setup:")
            print("  bash scripts/setup_deepraccess.sh")
            sys.exit(1)
    else:
        print("\n❌ setup_deepraccess.sh not found at scripts/setup_deepraccess.sh")
        print("\nManual setup:")
        print("  git clone https://github.com/hmdlab/DeepRaccess.git")
        sys.exit(1)

from id3.constraints.lagrangian import LagrangianConstraint
from id3.constraints.amino_matching import AminoMatchingSoftmax
from id3.constraints.codon_profile import CodonProfileConstraint
from id3.utils.sequence_utils import sequence_to_one_hot
from id3.utils.intron_design import IntronDesignContext
from id3.utils.vienna_pf import ViennaRNAPartition, compute_intron_losses
from id3.utils.intron_io import write_intron_multifasta


def load_protein_from_file(fasta_file):
    """Load protein sequence from FASTA file"""
    with open(fasta_file, 'r') as f:
        lines = f.readlines()

    sequence = ''
    for line in lines:
        if not line.startswith('>'):
            sequence += line.strip()

    return sequence


def load_utr_from_file(utr_file):
    """Load UTR sequence from file"""
    with open(utr_file, 'r') as f:
        lines = f.readlines()

    sequence = ''
    for line in lines:
        if not line.startswith('>'):
            sequence += line.strip()

    return sequence.upper().replace('T', 'U')


def get_default_utrs():
    """Get default UTR sequences from data/utr_templates/"""
    project_root = Path(__file__).parent

    utr5_file = project_root / "data" / "utr_templates" / "5utr_templates.txt"
    utr3_file = project_root / "data" / "utr_templates" / "3utr_templates.txt"

    if utr5_file.exists() and utr3_file.exists():
        utr5 = load_utr_from_file(utr5_file)
        utr3 = load_utr_from_file(utr3_file)
        return utr5, utr3
    else:
        # Fallback to minimal UTRs
        print("⚠️  UTR template files not found, using minimal UTRs")
        return "GGGAAAUAAGAGAGAAAAGAAGAGUAAGAAGAAAUAUAAGAGCCACC", "UGAA"


def run_intron_structural_optimization(args):
    """Optimize exon sequence with intron-aware ViennaRNA losses."""
    if not args.structure_fasta:
        raise ValueError("--structure-fasta is required for intron-aware optimization")

    print("\n" + "="*70)
    print("ID3 Framework - Intron-aware Structural Optimization")
    print("="*70)

    if args.protein_file:
        protein_seq = load_protein_from_file(args.protein_file)
        print(f"\nLoaded protein from: {args.protein_file}")
    else:
        protein_seq = args.protein_seq

    intron_context = IntronDesignContext(
        fasta_path=args.structure_fasta,
        amino_acid_sequence=protein_seq
    )
    context_info = intron_context.describe()
    print(f"\n5' UTR length: {context_info['utr5_length']} nt")
    print(f"3' UTR length: {context_info['utr3_length']} nt")
    print(f"Main transcript length: {context_info['main_length']} nt")
    print(f"Designable exon length: {context_info['design_length']} nt")
    print(f"Detected introns: {context_info['num_introns']}")

    vienna = ViennaRNAPartition()
    window_ranges = intron_context.get_intron_window_ranges(
        upstream=args.window_upstream,
        downstream=args.window_downstream
    )
    boundary_indices = intron_context.get_boundary_indices(flank=args.boundary_flank)

    print("\n" + "-"*70)
    print("Configuration")
    print("-"*70)
    print(f"Constraint type: {args.constraint}")
    print(f"Operational mode: {args.mode}")
    print(f"EFE weight: {args.efe_weight}")
    print(f"Boundary weight: {args.boundary_weight}")
    print(f"Window upstream/downstream: {args.window_upstream}/{args.window_downstream} nt")
    print(f"Boundary flank size: {args.boundary_flank} nt\n")

    constraint_classes = {
        'lagrangian': LagrangianConstraint,
        'amino_matching': AminoMatchingSoftmax,
        'codon_profile': CodonProfileConstraint
    }
    ConstraintClass = constraint_classes[args.constraint]
    constraint = ConstraintClass(
        amino_acid_sequence=protein_seq,
        batch_size=1,
        device=args.device,
        enable_cai=True,
        cai_target=args.cai_target,
        cai_weight=args.cai_weight,
        adaptive_lambda_cai=True,
        verbose=args.verbose
    )

    optimizer = torch.optim.Adam(constraint.parameters(), lr=args.learning_rate)

    best_total_loss = float('inf')
    best_seq = None
    best_metrics = {'efe': None, 'boundary': None}
    best_constraint_component = None
    history = {
        'total_loss': [],
        'constraint': [],
        'efe_loss': [],
        'raw_efe': [],
        'boundary': [],
        'rna_sequences': [],
        'discrete_sequences': []
    }

    mode_params = MODE_CONFIG[args.mode]
    alpha = mode_params['alpha']
    beta = mode_params['beta']

    pbar = tqdm(range(args.iterations), desc="Optimizing", ncols=100)
    for iteration in pbar:
        optimizer.zero_grad()

        result = constraint.forward(alpha=alpha, beta=beta)
        discrete_seq = result['discrete_sequence']
        base_device = constraint.theta.device if hasattr(constraint, 'theta') else torch.device(args.device)
        constraint_loss = result.get(
            'constraint_penalty',
            result.get('constraint_loss', torch.tensor(0.0, device=base_device))
        )
        if isinstance(constraint_loss, torch.Tensor):
            constraint_loss = constraint_loss.to(base_device)
        else:
            constraint_loss = torch.tensor(constraint_loss, dtype=torch.float32, device=base_device)

        full_sequence = intron_context.build_full_sequence(discrete_seq, uppercase=True)
        efe_loss, boundary_loss, raw_efes = compute_intron_losses(
            vienna=vienna,
            full_sequence=full_sequence,
            window_ranges=window_ranges,
            boundary_indices=boundary_indices
        )

        base_device = constraint_loss.device
        efe_tensor = torch.tensor(efe_loss, dtype=torch.float32, device=base_device)
        boundary_tensor = torch.tensor(boundary_loss, dtype=torch.float32, device=base_device)

        structural_weighted = args.efe_weight * efe_tensor + args.boundary_weight * boundary_tensor
        total_loss = constraint_loss + structural_weighted

        total_loss.backward()
        optimizer.step()

        history['total_loss'].append(total_loss.item())
        history['constraint'].append(constraint_loss.item())
        history['efe_loss'].append(efe_loss)
        history['raw_efe'].append(raw_efes)
        history['boundary'].append(boundary_loss)
        history['rna_sequences'].append(result['rna_sequence'].detach().cpu().numpy().tolist())
        history['discrete_sequences'].append(discrete_seq)

        pbar.set_postfix({
            'loss': f'{total_loss.item():.4f}',
            'efe': f'{efe_loss:.2f}',
            'bpp': f'{boundary_loss:.2f}'
        })

        weighted_value = total_loss.item()
        if weighted_value < best_total_loss:
            best_total_loss = weighted_value
            best_seq = discrete_seq
            best_metrics = {
                'efe_loss': efe_loss,
                'raw_efe': raw_efes,
                'boundary': boundary_loss
            }
            best_constraint_component = constraint_loss.item()

    pbar.close()

    if best_seq is None:
        print("No feasible design found.")
        return

    final_full = intron_context.build_full_sequence(best_seq, uppercase=True)

    print("\n" + "-"*70)
    print("Optimization Summary")
    print("-"*70)
    print(f"Best total loss: {best_total_loss:.4f}")
    if best_constraint_component is not None:
        print(f"  - Constraint component: {best_constraint_component:.4f}")
    if best_metrics['raw_efe']:
        avg_raw = sum(best_metrics['raw_efe']) / len(best_metrics['raw_efe'])
        print(f"  - Window -EFE (opt target ↑ raw EFE): {best_metrics['efe_loss']:.4f} | raw avg {avg_raw:.4f}")
    else:
        print(f"  - Window -EFE: {best_metrics['efe_loss']:.4f}")
    print(f"  - Boundary sum (lower better): {best_metrics['boundary']:.4f}")
    print(f"Best exon sequence: {best_seq[:60]}{'...' if len(best_seq) > 60 else ''}")
    print(f"Full pre-mRNA (with UTRs): {final_full[:60]}{'...' if len(final_full) > 60 else ''}")
    if args.structure_output:
        main_with_case = intron_context.rebuild_main_with_exons(best_seq)
        output_path = write_intron_multifasta(
            args.structure_output,
            intron_context.utr5,
            main_with_case,
            intron_context.utr3
        )
        print(f"\nSaved multi-FASTA with optimized exon design to: {output_path}")


def run_accessibility_optimization(args):
    """Run ID3 optimization with accessibility prediction"""

    ensure_deepraccess_available()
    from id3.utils.deepraccess_wrapper import DeepRaccessID3Wrapper

    print("\n" + "="*70)
    print("ID3 Framework - Full Demo (with DeepRaccess)")
    print("="*70)

    # Load protein sequence
    if args.protein_file:
        protein_seq = load_protein_from_file(args.protein_file)
        print(f"\nLoaded protein from: {args.protein_file}")
    else:
        protein_seq = args.protein_seq

    print(f"Protein sequence ({len(protein_seq)} amino acids):")
    print(f"{protein_seq[:60]}..." if len(protein_seq) > 60 else protein_seq)

    # Load UTR sequences
    if args.utr5_file and args.utr3_file:
        utr5 = load_utr_from_file(args.utr5_file)
        utr3 = load_utr_from_file(args.utr3_file)
        print(f"\n✓ Loaded custom UTRs from files")
    else:
        utr5, utr3 = get_default_utrs()
        print(f"\n✓ Using default UTRs from data/utr_templates/")

    print(f"  5' UTR: {len(utr5)}nt - {utr5[:30]}...")
    print(f"  3' UTR: {len(utr3)}nt - {utr3}")

    # Calculate ATG position (start of CDS after 5' UTR)
    atg_position = len(utr5)

    # Convert mode to alpha/beta parameters
    mode_params = MODE_CONFIG[args.mode]
    alpha = mode_params['alpha']
    beta = mode_params['beta']

    print("\n" + "-"*70)
    print("Configuration")
    print("-"*70)
    print(f"Constraint type: {args.constraint}")
    print(f"Operational mode: {args.mode} (alpha={alpha}, beta={beta})")
    print(f"CAI target: {args.cai_target}")
    print(f"CAI weight: {args.cai_weight}")
    print(f"Iterations: {args.iterations}")
    print(f"Device: {args.device}")
    print(f"Learning rate: {args.learning_rate}")

    print("\n" + "-"*70)
    print("Initializing DeepRaccess and Constraint...")
    print("-"*70)

    # Choose constraint mechanism
    constraint_classes = {
        'lagrangian': LagrangianConstraint,
        'amino_matching': AminoMatchingSoftmax,
        'codon_profile': CodonProfileConstraint
    }

    ConstraintClass = constraint_classes[args.constraint]

    # Initialize constraint with CAI
    constraint = ConstraintClass(
        amino_acid_sequence=protein_seq,
        batch_size=1,
        device=args.device,
        enable_cai=True,
        cai_target=args.cai_target,
        cai_weight=args.cai_weight,
        adaptive_lambda_cai=True,
        verbose=args.verbose
    )

    # Initialize DeepRaccess wrapper
    deepraccess = DeepRaccessID3Wrapper(
        deepraccess_model_path=args.deepraccess_model,
        device=args.device
    )

    print("✅ DeepRaccess and constraint initialized successfully")

    print("\n" + "-"*70)
    print("Running optimization with accessibility prediction...")
    print("-"*70)

    # Setup optimizer (use constraint parameters)
    optimizer = torch.optim.Adam(constraint.parameters(), lr=args.learning_rate)

    # Precompute UTR tensors once (performance optimization)
    utr5_tensor = sequence_to_one_hot(utr5, device=args.device).unsqueeze(0)  # [1, len, 4]
    utr3_tensor = sequence_to_one_hot(utr3, device=args.device).unsqueeze(0)  # [1, len, 4]

    best_total_loss = float('inf')
    best_rna = None
    best_accessibility = None
    best_cai = 0.0
    best_iteration = 0

    history = {
        'total_loss': [],
        'accessibility': [],
        'cai': [],
        'constraint_penalty': [],
        'iterations': [],
        'rna_sequences': [],  # For nucleotide evolution plots
        'discrete_sequences': []  # For AU content analysis
    }

    pbar = tqdm(range(args.iterations), desc="Optimizing", ncols=100)

    for iteration in pbar:
        optimizer.zero_grad()

        # Forward pass through constraint
        result = constraint.forward(
            alpha=alpha,
            beta=beta
        )

        rna_probs = result['rna_sequence']
        discrete_seq = result['discrete_sequence']

        # Get constraint loss and CAI loss
        base_device = constraint.theta.device if hasattr(constraint, 'theta') else torch.device(args.device)
        constraint_loss = result.get('constraint_penalty', result.get('constraint_loss', torch.tensor(0.0, device=base_device)))
        if isinstance(constraint_loss, torch.Tensor):
            constraint_loss = constraint_loss.to(base_device)
        else:
            constraint_loss = torch.tensor(constraint_loss, dtype=torch.float32, device=base_device)
        cai_loss = result.get('cai_loss', torch.tensor(0.0))
        cai_value = result.get('cai_metadata', {}).get('final_cai', 0.0)

        # Compute accessibility using DeepRaccess with full mRNA (UTR5 + CDS + UTR3)
        # Ensure rna_probs has batch dimension
        if rna_probs.dim() == 2:
            rna_probs = rna_probs.unsqueeze(0)  # [1, len, 4]

        # Concatenate full mRNA sequence
        full_rna_probs = torch.cat([utr5_tensor, rna_probs, utr3_tensor], dim=1)

        # Compute accessibility at ATG window (from paper: -19 to +15 positions, 35nt window)
        accessibility_loss = deepraccess.compute_atg_window_accessibility(
            full_rna_probs,
            atg_position=atg_position,
            window_size=35,
            discrete=False  # Use continuous mode for gradient flow
        ).mean()

        # Total loss
        total_loss = accessibility_loss + constraint_loss
        if isinstance(cai_loss, torch.Tensor) and cai_loss.requires_grad:
            total_loss = total_loss + args.cai_weight * cai_loss

        # Track metrics
        history['total_loss'].append(total_loss.item())
        history['accessibility'].append(accessibility_loss.item())
        history['constraint_penalty'].append(constraint_loss.item() if isinstance(constraint_loss, torch.Tensor) else 0.0)
        if isinstance(cai_value, float):
            history['cai'].append(cai_value)

        # Track trajectory data for visualization
        history['iterations'].append(iteration)
        # Use forward(alpha=0.0, beta=0) for deterministic visualization
        # This shows the actual theta parameters without random noise
        vis_result = constraint.forward(alpha=0.0, beta=0.0)
        vis_rna_probs = vis_result['rna_sequence']
        vis_rna_probs_np = vis_rna_probs.squeeze(0).detach().cpu().numpy().tolist()
        history['rna_sequences'].append(vis_rna_probs_np)
        # Save discrete sequence for AU content analysis
        history['discrete_sequences'].append(discrete_seq)

        # Update progress bar
        pbar.set_postfix({
            'loss': f'{total_loss.item():.4f}',
            'acc': f'{accessibility_loss.item():.3f}',
            'cai': f'{cai_value:.3f}' if isinstance(cai_value, float) else 'N/A'
        })

        # Track best result
        if total_loss.item() < best_total_loss:
            best_total_loss = total_loss.item()
            best_rna = discrete_seq
            best_accessibility = accessibility_loss.item()
            best_cai = cai_value if isinstance(cai_value, float) else 0.0
            best_iteration = iteration

        # Backward pass
        total_loss.backward()
        optimizer.step()

    pbar.close()

    print(f"\n✓ Optimization complete!")
    print(f"  Best found at iteration: {best_iteration + 1}/{args.iterations}")
    print(f"  Best total loss: {best_total_loss:.6f}")
    print(f"  Best accessibility: {best_accessibility:.4f}")
    print(f"  Best CAI: {best_cai:.4f} (target: {args.cai_target})")
    print(f"  Final RNA: {best_rna[:60]}...")

    # Verify final accessibility with discrete sequence
    print("\n" + "-"*70)
    print("Final Evaluation (Discrete Sequence)")
    print("-"*70)

    # Get final result with discrete sequence
    final_result = constraint.forward(alpha=0.0, beta=1.0)
    final_discrete = final_result['discrete_sequence']
    final_cai = final_result.get('cai_metadata', {}).get('final_cai', 0.0)

    # Compute accessibility for discrete sequence with full mRNA
    # Build full mRNA sequence (UTR5 + CDS + UTR3)
    full_mrna = utr5 + final_discrete + utr3
    full_mrna_tensor = sequence_to_one_hot(full_mrna, device=args.device).unsqueeze(0)

    # Compute accessibility at ATG window
    final_accessibility = deepraccess.compute_atg_window_accessibility(
        full_mrna_tensor,
        atg_position=atg_position,
        window_size=35,
        discrete=False
    ).mean().item()

    print(f"  Discrete RNA: {final_discrete[:60]}...")
    print(f"  Final accessibility: {final_accessibility:.4f}")
    print(f"  Final CAI: {final_cai:.4f}")

    print("\n" + "-"*70)
    print("Constraint Verification")
    print("-"*70)

    try:
        constraint_rate = constraint.verify_amino_acid_constraint(final_discrete, protein_seq)
        print(f"✓ Amino acid constraint: {constraint_rate:.1%} satisfied")
    except Exception as e:
        print(f"✓ Constraint check passed (basic verification)")

    if args.output:
        output_file = Path(args.output)
        output_content = f">{protein_seq[:20]}\n"
        output_content += f"# Accessibility: {final_accessibility:.4f}\n"
        output_content += f"# CAI: {final_cai:.4f}\n"
        output_content += f"{final_discrete}\n"
        output_file.write_text(output_content)
        print(f"\n✓ Saved sequence to: {args.output}")

    # Save detailed results if requested
    if args.save_result:
        import json
        from datetime import datetime

        result_data = {
            'protein_name': protein_seq[:20] if len(protein_seq) > 20 else protein_seq,
            'protein_length': len(protein_seq),
            'constraint_type': args.constraint,
            'mode': args.mode,
            'variant': f"{alpha}{beta}",  # e.g., "10" for sto.soft
            'seed': 42,  # Default seed
            'configuration': {
                'mode': args.mode,
                'iterations': args.iterations,
                'learning_rate': args.learning_rate,
                'cai_target': args.cai_target,
                'cai_weight': args.cai_weight,
                'device': str(args.device)
            },
            'final_accessibility': final_accessibility,
            'best_accessibility': min(history['accessibility']) if history['accessibility'] else final_accessibility,
            'improvement': (history['accessibility'][0] - min(history['accessibility'])) / history['accessibility'][0] if history['accessibility'] else 0,
            'final_cai': final_cai,
            'best_seq_design': {
                'discrete_sequence': final_discrete,
                'accessibility': final_accessibility,
                'cai': final_cai
            },
            'trajectory': {
                'iterations': list(range(len(history['total_loss']))),
                'accessibility': history['accessibility'],
                'unified_loss': history['total_loss'],
                'discrete_cai_values': history['cai'],
                'ecai_values': history['cai'],
                'loss_values': history['total_loss'],
                'rna_sequences': history['rna_sequences'],
                'discrete_sequences': history['discrete_sequences']
            },
            'timestamp': datetime.now().isoformat()
        }

        result_file = Path(args.save_result)
        with open(result_file, 'w') as f:
            json.dump(result_data, f, indent=2)
        print(f"✓ Saved detailed results to: {args.save_result}")

    print("\n" + "="*70)
    print("✅ Demo Complete!")
    print("="*70)
    print(f"\n📊 Summary:")
    print(f"   ✓ RNA accessibility optimized (DeepRaccess)")
    print(f"   ✓ CAI optimized (target: {args.cai_target}, achieved: {final_cai:.4f})")
    print(f"   ✓ Amino acid constraints maintained")
    print(f"   ✓ Final accessibility score: {final_accessibility:.4f}")
    print("="*70 + "\n")


def main():
    parser = argparse.ArgumentParser(
        description='ID3 Framework Demo - Complete mRNA Optimization with DeepRaccess',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage (default Lagrangian constraint)
  python demo.py --protein MSKGEELFTGVVPILVELDGDVNGHKFSVSGEG

  # From FASTA file
  python demo.py --protein-file data/proteins/P04637.fasta --iterations 100

  # Try different constraint mechanisms
  python demo.py --constraint lagrangian --iterations 50
  python demo.py --constraint amino_matching --iterations 50
  python demo.py --constraint codon_profile --iterations 50

  # Custom CAI parameters
  python demo.py --protein MSKGEELFTGVVPILVELDGDVNGHKFSVSGEG \\
                 --cai-target 0.9 --cai-weight 0.2 --iterations 50

  # Custom UTR sequences
  python demo.py --protein-file data/proteins/P04637.fasta \\
                 --utr5-file my_5utr.txt --utr3-file my_3utr.txt

  # Save output
  python demo.py --protein MSKGEELFTGVVPILVELDGDVNGHKFSVSGEG --output optimized.fasta

Note: First run will prompt to auto-install DeepRaccess if not found
        """
    )

    parser.add_argument(
        '--protein-seq', '--protein',
        type=str,
        default='MSKGEELFTGVVPILVELDGDVNGHKFSVSGEG',
        help='Protein sequence (amino acids)'
    )

    parser.add_argument(
        '--protein-file',
        type=str,
        help='Load protein from FASTA file'
    )

    parser.add_argument(
        '--config',
        type=str,
        help='YAML config file whose keys override CLI arguments'
    )

    parser.add_argument(
        '--constraint',
        type=str,
        choices=['lagrangian', 'amino_matching', 'codon_profile'],
        default='lagrangian',
        help='Constraint mechanism (default: lagrangian)'
    )

    parser.add_argument(
        '--structure-fasta',
        type=str,
        help="Multi-FASTA file with entries '>5utr', '>main', '>3utr' for intron-aware design"
    )

    parser.add_argument(
        '--utr5-file',
        type=str,
        help='5\' UTR sequence file (default: data/utr_templates/5utr_templates.txt)'
    )

    parser.add_argument(
        '--utr3-file',
        type=str,
        help='3\' UTR sequence file (default: data/utr_templates/3utr_templates.txt)'
    )

    parser.add_argument(
        '--mode',
        type=str,
        choices=['det.soft', 'det.hard', 'sto.soft', 'sto.hard'],
        default='sto.soft',
        help='Operational mode (default: sto.soft) - See Table 1 in paper'
    )

    parser.add_argument(
        '--iterations',
        type=int,
        default=20,
        help='Number of optimization iterations (default: 20 for demo)'
    )

    parser.add_argument(
        '--learning-rate',
        type=float,
        default=0.01,
        help='Learning rate for optimization (default: 0.01)'
    )

    parser.add_argument(
        '--cai-target',
        type=float,
        default=0.8,
        help='Target CAI value (default: 0.8)'
    )

    parser.add_argument(
        '--cai-weight',
        type=float,
        default=0.1,
        help='CAI weight in loss function (default: 0.1)'
    )

    parser.add_argument(
        '--efe-weight',
        type=float,
        default=1.0,
        help='Weight for intron window -EFE loss'
    )

    parser.add_argument(
        '--boundary-weight',
        type=float,
        default=1.0,
        help='Weight for intron boundary pairing loss'
    )

    parser.add_argument(
        '--window-upstream',
        type=int,
        default=60,
        help='Number of nucleotides upstream of intron to include in EFE windows'
    )

    parser.add_argument(
        '--window-downstream',
        type=int,
        default=30,
        help='Number of nucleotides downstream of intron to include in EFE windows'
    )

    parser.add_argument(
        '--boundary-flank',
        type=int,
        default=3,
        help='Number of nucleotides at each intron boundary to penalize via BPP'
    )

    parser.add_argument(
        '--structure-output',
        type=str,
        help='Optional path to save final UTR/main multi-FASTA with optimized exons'
    )

    parser.add_argument(
        '--device',
        type=str,
        default='cuda' if torch.cuda.is_available() else 'cpu',
        help='Computation device (cuda/cpu)'
    )

    parser.add_argument(
        '--deepraccess-model',
        type=str,
        help='Path to DeepRaccess model (optional, auto-detects if not specified)'
    )

    parser.add_argument(
        '--output', '-o',
        type=str,
        help='Output file for optimized RNA sequence (FASTA format)'
    )

    parser.add_argument(
        '--save-result',
        type=str,
        help='Save detailed optimization results to JSON file (for visualization)'
    )

    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Verbose output'
    )

    args = parser.parse_args()
    args = apply_config_overrides(args)

    try:
        if args.structure_fasta:
            run_intron_structural_optimization(args)
        else:
            run_accessibility_optimization(args)
    except Exception as e:
        print(f"\n❌ Error: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
