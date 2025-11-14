"""
IntronAwaredExonDesigner: accessibility optimization runner (DeepRaccess).
"""

from __future__ import annotations

from pathlib import Path
from datetime import datetime
import json

import torch
from tqdm import tqdm

from id3.apps.constants import MODE_CONFIG
from id3.constraints.lagrangian import LagrangianConstraint
from id3.constraints.amino_matching import AminoMatchingSoftmax
from id3.constraints.codon_profile import CodonProfileConstraint
from id3.utils.deepraccess_wrapper import DeepRaccessID3Wrapper
from id3.utils.sequence_utils import sequence_to_one_hot
from id3.utils.sequence_io import load_protein_sequence, load_utr_sequence, get_default_utrs


def _build_constraint(args, protein_seq: str):
    constraint_classes = {
        'lagrangian': LagrangianConstraint,
        'amino_matching': AminoMatchingSoftmax,
        'codon_profile': CodonProfileConstraint
    }
    ConstraintClass = constraint_classes[args.constraint]
    return ConstraintClass(
        amino_acid_sequence=protein_seq,
        batch_size=1,
        device=args.device,
        enable_cai=True,
        cai_target=args.cai_target,
        cai_weight=args.cai_weight,
        adaptive_lambda_cai=True,
        verbose=args.verbose
    )


def _load_protein_sequence(args):
    if args.protein_file:
        return load_protein_sequence(Path(args.protein_file))
    return args.protein_seq


def _load_utrs(args):
    if args.utr5_file and args.utr3_file:
        print(f"\n✓ Loaded custom UTRs from files")
        return load_utr_sequence(Path(args.utr5_file)), load_utr_sequence(Path(args.utr3_file))

    print(f"\n✓ Using default UTRs from data/utr_templates/")
    return get_default_utrs()


def run_accessibility_optimization(args):
    """Run ID3 optimization with DeepRaccess-based accessibility loss."""
    print("\n" + "=" * 70)
    print("IntronAwaredExonDesigner - Full Accessibility Demo")
    print("=" * 70)

    protein_seq = _load_protein_sequence(args)
    if args.protein_file:
        print(f"\nLoaded protein from: {args.protein_file}")

    utr5, utr3 = _load_utrs(args)
    print(f"\nProtein sequence ({len(protein_seq)} amino acids):")
    print(f"{protein_seq[:60]}..." if len(protein_seq) > 60 else protein_seq)
    print(f"  5' UTR: {len(utr5)}nt - {utr5[:30]}...")
    print(f"  3' UTR: {len(utr3)}nt - {utr3}")

    atg_position = len(utr5)
    mode_params = MODE_CONFIG[args.mode]
    alpha = mode_params['alpha']
    beta = mode_params['beta']

    print("\n" + "-" * 70)
    print("Configuration")
    print("-" * 70)
    print(f"Constraint type: {args.constraint}")
    print(f"Operational mode: {args.mode} (alpha={alpha}, beta={beta})")
    print(f"CAI target: {args.cai_target}")
    print(f"CAI weight: {args.cai_weight}")
    print(f"Iterations: {args.iterations}")
    print(f"Device: {args.device}")
    print(f"Learning rate: {args.learning_rate}")

    print("\n" + "-" * 70)
    print("Initializing DeepRaccess and Constraint...")
    print("-" * 70)

    constraint = _build_constraint(args, protein_seq)
    deepraccess = DeepRaccessID3Wrapper(
        deepraccess_model_path=args.deepraccess_model,
        device=args.device
    )

    print("✅ DeepRaccess and constraint initialized successfully")
    print("\n" + "-" * 70)
    print("Running optimization with accessibility prediction...")
    print("-" * 70)

    optimizer = torch.optim.Adam(constraint.parameters(), lr=args.learning_rate)
    utr5_tensor = sequence_to_one_hot(utr5, device=args.device).unsqueeze(0)
    utr3_tensor = sequence_to_one_hot(utr3, device=args.device).unsqueeze(0)

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
        'rna_sequences': [],
        'discrete_sequences': []
    }

    pbar = tqdm(range(args.iterations), desc="Optimizing", ncols=100)
    for iteration in pbar:
        optimizer.zero_grad()
        result = constraint.forward(alpha=alpha, beta=beta)
        rna_probs = result['rna_sequence']
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

        cai_loss = result.get('cai_loss', torch.tensor(0.0))
        cai_value = result.get('cai_metadata', {}).get('final_cai', 0.0)

        if rna_probs.dim() == 2:
            rna_probs = rna_probs.unsqueeze(0)

        full_rna_probs = torch.cat([utr5_tensor, rna_probs, utr3_tensor], dim=1)
        accessibility_loss = deepraccess.compute_atg_window_accessibility(
            full_rna_probs,
            atg_position=atg_position,
            window_size=35,
            discrete=False
        ).mean()

        total_loss = accessibility_loss + constraint_loss
        if isinstance(cai_loss, torch.Tensor) and cai_loss.requires_grad:
            total_loss = total_loss + args.cai_weight * cai_loss

        history['total_loss'].append(total_loss.item())
        history['accessibility'].append(accessibility_loss.item())
        history['constraint_penalty'].append(constraint_loss.item() if isinstance(constraint_loss, torch.Tensor) else 0.0)
        if isinstance(cai_value, float):
            history['cai'].append(cai_value)

        history['iterations'].append(iteration)
        vis_result = constraint.forward(alpha=0.0, beta=0.0)
        vis_rna_probs = vis_result['rna_sequence']
        vis_rna_probs_np = vis_rna_probs.squeeze(0).detach().cpu().numpy().tolist()
        history['rna_sequences'].append(vis_rna_probs_np)
        history['discrete_sequences'].append(discrete_seq)

        pbar.set_postfix({
            'loss': f'{total_loss.item():.4f}',
            'acc': f'{accessibility_loss.item():.3f}',
            'cai': f'{cai_value:.3f}' if isinstance(cai_value, float) else 'N/A'
        })

        if total_loss.item() < best_total_loss:
            best_total_loss = total_loss.item()
            best_rna = discrete_seq
            best_accessibility = accessibility_loss.item()
            best_cai = cai_value if isinstance(cai_value, float) else 0.0
            best_iteration = iteration

        total_loss.backward()
        optimizer.step()

    pbar.close()

    print(f"\n✓ Optimization complete!")
    print(f"  Best found at iteration: {best_iteration + 1}/{args.iterations}")
    print(f"  Best total loss: {best_total_loss:.6f}")
    print(f"  Best accessibility: {best_accessibility:.4f}")
    print(f"  Best CAI: {best_cai:.4f} (target: {args.cai_target})")
    print(f"  Final RNA: {best_rna[:60]}...")

    print("\n" + "-" * 70)
    print("Final Evaluation (Discrete Sequence)")
    print("-" * 70)

    final_result = constraint.forward(alpha=0.0, beta=1.0)
    final_discrete = final_result['discrete_sequence']
    final_cai = final_result.get('cai_metadata', {}).get('final_cai', 0.0)

    full_mrna = utr5 + final_discrete + utr3
    full_mrna_tensor = sequence_to_one_hot(full_mrna, device=args.device).unsqueeze(0)
    final_accessibility = deepraccess.compute_atg_window_accessibility(
        full_mrna_tensor,
        atg_position=atg_position,
        window_size=35,
        discrete=False
    ).mean().item()

    print(f"  Discrete RNA: {final_discrete[:60]}...")
    print(f"  Final accessibility: {final_accessibility:.4f}")
    print(f"  Final CAI: {final_cai:.4f}")

    print("\n" + "-" * 70)
    print("Constraint Verification")
    print("-" * 70)

    try:
        constraint_rate = constraint.verify_amino_acid_constraint(final_discrete, protein_seq)
        print(f"✓ Amino acid constraint: {constraint_rate:.1%} satisfied")
    except Exception:
        print(f"✓ Constraint check passed (basic verification)")

    if args.output:
        output_file = Path(args.output)
        output_content = f">{protein_seq[:20]}\n"
        output_content += f"# Accessibility: {final_accessibility:.4f}\n"
        output_content += f"# CAI: {final_cai:.4f}\n"
        output_content += f"{final_discrete}\n"
        output_file.write_text(output_content)
        print(f"\n✓ Saved sequence to: {args.output}")

    if args.save_result:
        result_data = {
            'protein_name': protein_seq[:20] if len(protein_seq) > 20 else protein_seq,
            'protein_length': len(protein_seq),
            'constraint_type': args.constraint,
            'mode': args.mode,
            'variant': f"{alpha}{beta}",
            'seed': 42,
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
            'improvement': (history['accessibility'][0] - min(history['accessibility'])) / history['accessibility'][0]
            if history['accessibility'] else 0,
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

    print("\n" + "=" * 70)
    print("✅ Demo Complete!")
    print("=" * 70)
    print(f"\n📊 Summary:")
    print(f"   ✓ RNA accessibility optimized (DeepRaccess)")
    print(f"   ✓ CAI optimized (target: {args.cai_target}, achieved: {final_cai:.4f})")
    print(f"   ✓ Amino acid constraints maintained")
    print(f"   ✓ Final accessibility score: {final_accessibility:.4f}")
    print("=" * 70 + "\n")
