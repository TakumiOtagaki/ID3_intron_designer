"""
IntronAwaredExonDesigner: intron-aware structural optimization runner.
"""

from __future__ import annotations

from pathlib import Path
import json
from typing import Dict, List

import torch
from tqdm import tqdm

from id3.apps.constants import MODE_CONFIG
from id3.constraints.lagrangian import LagrangianConstraint
from id3.constraints.amino_matching import AminoMatchingSoftmax
from id3.constraints.codon_profile import CodonProfileConstraint
from id3.utils.intron_design import IntronAwaredExonDesignerContext
from id3.utils.sequence_io import load_protein_sequence, rna_to_dna
from id3.utils.vienna_pf import ViennaRNAPartition, compute_intron_losses


def _compute_structural_metrics(full_sequence: str, vienna, window_ranges, boundary_indices) -> Dict[str, List[float] | float]:
    efe_loss, boundary_loss, raw_efes = compute_intron_losses(
        vienna=vienna,
        full_sequence=full_sequence,
        window_ranges=window_ranges,
        boundary_indices=boundary_indices
    )
    raw_avg = sum(raw_efes) / len(raw_efes) if raw_efes else 0.0
    return {
        'efe_loss': efe_loss,
        'boundary': boundary_loss,
        'raw_efe': raw_efes,
        'raw_efe_avg': raw_avg
    }


def _build_candidate_entry(label: str, exon_sequence: str, context: IntronAwaredExonDesignerContext, vienna, window_ranges, boundary_indices) -> Dict:
    main_rna = context.rebuild_main_with_exons(exon_sequence)
    full_sequence = context.build_full_sequence(exon_sequence, uppercase=True)
    metrics = _compute_structural_metrics(full_sequence, vienna, window_ranges, boundary_indices)
    return {
        'label': label,
        'exon_rna': exon_sequence,
        'main_rna': main_rna,
        'full_sequence': full_sequence,
        'metrics': metrics
    }


def _format_candidate_header(prefix: str, label: str, metrics: Dict[str, float]) -> str:
    return (
        f"{prefix} {label} | efe={metrics['efe_loss']:.2f} | "
        f"boundary={metrics['boundary']:.2f} | raw_avg={metrics['raw_efe_avg']:.2f}"
    )


def _write_design_multifasta(path: Path, utr5_dna: str, utr3_dna: str, candidates: List[Dict]) -> None:
    with open(path, 'w') as handle:
        for candidate in candidates:
            main_dna = rna_to_dna(candidate['main_rna'])
            for prefix, seq in (('5utr', utr5_dna), ('main', main_dna), ('3utr', utr3_dna)):
                header = _format_candidate_header(prefix, candidate['label'], candidate['metrics'])
                handle.write(f">{header}\n")
                for idx in range(0, len(seq), 60):
                    handle.write(seq[idx:idx + 60] + "\n")


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


def run_intron_awared_exon_structural_optimization(args):
    """Optimizes exon codons so that intron windows are destabilized."""
    if not args.structure_fasta:
        raise ValueError("--structure-fasta is required for intron-aware optimization")

    print("\n" + "=" * 70)
    print("IntronAwaredExonDesigner - Structural Optimization")
    print("=" * 70)

    protein_seq = _load_protein_sequence(args)
    if args.protein_file:
        print(f"\nLoaded protein from: {args.protein_file}")

    context = IntronAwaredExonDesignerContext(
        fasta_path=args.structure_fasta,
        amino_acid_sequence=protein_seq
    )
    context_info = context.describe()
    print(f"\n5' UTR length: {context_info['utr5_length']} nt")
    print(f"3' UTR length: {context_info['utr3_length']} nt")
    print(f"Main transcript length: {context_info['main_length']} nt")
    print(f"Designable exon length: {context_info['design_length']} nt")
    print(f"Detected introns: {context_info['num_introns']}")

    vienna = ViennaRNAPartition()
    window_ranges = context.get_intron_window_ranges(
        upstream=args.window_upstream,
        downstream=args.window_downstream
    )
    boundary_indices = context.get_boundary_indices(flank=args.boundary_flank)

    # Baseline metrics on the provided sequence
    original_exon = context.get_exon_sequence()
    initial_full_sequence = context.build_full_sequence(original_exon, uppercase=True)
    baseline_efe_loss, baseline_boundary_loss, baseline_raw_efes = compute_intron_losses(
        vienna=vienna,
        full_sequence=initial_full_sequence,
        window_ranges=window_ranges,
        boundary_indices=boundary_indices
    )

    print("\n" + "-" * 70)
    print("Configuration")
    print("-" * 70)
    print(f"Constraint type: {args.constraint}")
    print(f"Operational mode: {args.mode}")
    print(f"EFE weight: {args.efe_weight}")
    print(f"Boundary weight: {args.boundary_weight}")
    print(f"Window upstream/downstream: {args.window_upstream}/{args.window_downstream} nt")
    print(f"Boundary flank size: {args.boundary_flank} nt\n")

    constraint = _build_constraint(args, protein_seq)
    optimizer = torch.optim.Adam(constraint.parameters(), lr=args.learning_rate)

    best_total_loss = float('inf')
    best_seq = None
    best_constraint_component = None
    history: Dict[str, List] = {
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
        rna_probs_current = result['rna_sequence']
        base_device = constraint.theta.device if hasattr(constraint, 'theta') else torch.device(args.device)
        constraint_loss = result.get(
            'constraint_penalty',
            result.get('constraint_loss', torch.tensor(0.0, device=base_device))
        )
        if isinstance(constraint_loss, torch.Tensor):
            constraint_loss = constraint_loss.to(base_device)
        else:
            constraint_loss = torch.tensor(constraint_loss, dtype=torch.float32, device=base_device)

        full_sequence = context.build_full_sequence(discrete_seq, uppercase=True)
        efe_loss, boundary_loss, raw_efes = compute_intron_losses(
            vienna=vienna,
            full_sequence=full_sequence,
            window_ranges=window_ranges,
            boundary_indices=boundary_indices
        )

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
        if rna_probs_current.dim() == 2:
            rna_probs_store = rna_probs_current.detach().cpu().numpy().tolist()
        else:
            rna_probs_store = rna_probs_current.squeeze(0).detach().cpu().numpy().tolist()
        history['rna_sequences'].append(rna_probs_store)
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
            best_constraint_component = constraint_loss.item()

    pbar.close()

    if best_seq is None:
        print("No feasible design found.")
        return

    best_candidate = _build_candidate_entry(
        'best',
        best_seq,
        context,
        vienna,
        window_ranges,
        boundary_indices
    )
    best_metrics = best_candidate['metrics']
    final_full = best_candidate['full_sequence']
    final_vis = constraint.forward(alpha=0.0, beta=0.0)
    final_probs = final_vis['rna_sequence']
    if final_probs.dim() == 2:
        final_probs = final_probs.unsqueeze(0)

    mutations = []
    if len(original_exon) == len(best_seq):
        utr5_len = len(context.utr5)
        for i, (ref_base, alt_base) in enumerate(zip(original_exon, best_seq)):
            if ref_base != alt_base:
                main_index = context.exon_positions[i]
                full_index = utr5_len + main_index
                codon_index = i // 3 + 1
                pos_in_codon = i % 3
                mutations.append({
                    'design_index': i + 1,
                    'main_index': main_index,
                    'full_index': full_index,
                    'codon_index': codon_index,
                    'pos_in_codon': pos_in_codon,
                    'ref_rna': ref_base,
                    'alt_rna': alt_base,
                })
    else:
        print("Warning: exon length mismatch when computing mutations; skipping list.")

    output_file_path = args.output_file
    loss_png_path = None
    try:
        import matplotlib.pyplot as _plt

        base_dir = Path(output_file_path).parent if output_file_path else Path("outputs")
        base_dir.mkdir(parents=True, exist_ok=True)
        loss_png = base_dir / "intron_loss_curve.png"
        xs = list(range(len(history['total_loss'])))
        _plt.figure(figsize=(8, 5))
        _plt.plot(xs, history['total_loss'], label='total')
        _plt.plot(xs, history['efe_loss'], label='window_-EFE')
        _plt.plot(xs, history['boundary'], label='boundary_sum')
        _plt.plot(xs, history['constraint'], label='constraint')
        _plt.xlabel('iteration')
        _plt.ylabel('loss')
        _plt.title('Intron structural optimization - loss trajectory')
        _plt.legend()
        _plt.tight_layout()
        _plt.savefig(loss_png)
        _plt.close()
        loss_png_path = str(loss_png)
        print(f"Saved intron loss curve to: {loss_png}")
    except Exception as exc:
        print(f"Warning: failed to save intron loss curve: {exc}")

    print("\n" + "-" * 70)
    print("Optimization Summary")
    print("-" * 70)
    print(f"Best total loss: {best_total_loss:.4f}")
    if best_constraint_component is not None:
        print(f"  - Constraint component: {best_constraint_component:.4f}")
    avg_raw = best_metrics['raw_efe_avg']
    if best_metrics['raw_efe']:
        print(f"  - Window -EFE (opt target ↑ raw EFE): {best_metrics['efe_loss']:.4f} | raw avg {avg_raw:.4f}")
    else:
        print(f"  - Window -EFE: {best_metrics['efe_loss']:.4f}")
    print(f"  - Boundary sum (lower better): {best_metrics['boundary']:.4f}")
    init_raw_avg = (sum(baseline_raw_efes) / len(baseline_raw_efes)) if baseline_raw_efes else 0.0
    print(f"  - Initial window -EFE: {baseline_efe_loss:.4f} | raw avg {init_raw_avg:.4f}")
    print(f"  - Δ window -EFE (best - init, lower better): {best_metrics['efe_loss'] - baseline_efe_loss:+.4f}")
    print(f"  - Δ raw EFE avg (best - init, higher better): {(avg_raw - init_raw_avg):+.4f}")
    print(f"  - Initial boundary sum: {baseline_boundary_loss:.4f}")
    print(f"  - Δ boundary sum (best - init, lower better): {best_metrics['boundary'] - baseline_boundary_loss:+.4f}")
    best_exon_dna = rna_to_dna(best_seq)
    final_full_dna = rna_to_dna(final_full)
    print(f"Best exon sequence: {best_exon_dna[:60]}{'...' if len(best_exon_dna) > 60 else ''}")
    print(f"Full pre-mRNA (with UTRs): {final_full_dna[:60]}{'...' if len(final_full_dna) > 60 else ''}")

    if mutations:
        print("\n" + "-" * 70)
        print("Mutations (exon design, DNA view)")
        print("-" * 70)
        print(f"Total changes: {len(mutations)} (of {len(original_exon)})")
        preview = mutations[:20]
        for m in preview:
            ref_dna = rna_to_dna(m['ref_rna'])
            alt_dna = rna_to_dna(m['alt_rna'])
            print(
                f"pos(exon) {m['design_index']:>4d} | codon {m['codon_index']:>4d}.{m['pos_in_codon']} | "
                f"main_idx {m['main_index']:>4d} | full_idx {m['full_index']:>4d} : {ref_dna} -> {alt_dna}"
            )
        if len(mutations) > len(preview):
            print(f"... and {len(mutations) - len(preview)} more changes")
    else:
        print("\nNo base changes in exon design (identical to original).")

    if output_file_path:
        output_path = Path(output_file_path)
        sample_n = int(getattr(args, 'sample_count', 0) or 0)
        sampled_candidates = []
        if sample_n > 0:
            try:
                import numpy as _np

                final_probs_np = final_probs.squeeze(0).detach().cpu().numpy()
                for idx in range(sample_n):
                    bases = []
                    for i in range(len(context.exon_positions)):
                        p = final_probs_np[i]
                        p = _np.maximum(p, 1e-9)
                        p = p / p.sum()
                        base_idx = _np.random.choice(4, p=p)
                        bases.append('ACGU'[base_idx])
                    exon_rna = ''.join(bases)
                    sampled_candidates.append(
                        _build_candidate_entry(
                            f'sample_{idx + 1}',
                            exon_rna,
                            context,
                            vienna,
                            window_ranges,
                            boundary_indices
                        )
                    )
            except Exception as exc:
                print(f"Warning: failed to sample sequences: {exc}")

        candidates = [best_candidate] + sampled_candidates
        utr5_dna = rna_to_dna(context.utr5)
        utr3_dna = rna_to_dna(context.utr3)
        base_dir = output_path.parent
        base_dir.mkdir(parents=True, exist_ok=True)
        _write_design_multifasta(output_path, utr5_dna, utr3_dna, candidates)
        print(f"\nSaved multi-FASTA with optimized and sampled sequences to: {output_path}")

        try:
            mut_path = Path(str(output_path)).with_suffix(".mutations.tsv")
            with open(mut_path, 'w') as mf:
                mf.write("design_index\tmain_index\tfull_index\tcodon_index\tpos_in_codon\tref_dna\talt_dna\n")
                for m in mutations:
                    ref_dna = rna_to_dna(m['ref_rna'])
                    alt_dna = rna_to_dna(m['alt_rna'])
                    mf.write(
                        f"{m['design_index']}\t{m['main_index']}\t{m['full_index']}\t{m['codon_index']}\t"
                        f"{m['pos_in_codon']}\t{ref_dna}\t{alt_dna}\n"
                    )
            print(f"Saved mutation list to: {mut_path}")
        except Exception as exc:
            print(f"Warning: failed to save mutation list TSV: {exc}")

        try:
            summary = {
                'baseline': {
                    'efe_loss': baseline_efe_loss,
                    'boundary_sum': baseline_boundary_loss,
                    'raw_efe_values': baseline_raw_efes,
                },
                'best': {
                    'efe_loss': best_metrics['efe_loss'],
                    'boundary_sum': best_metrics['boundary'],
                    'raw_efe_values': best_metrics['raw_efe'],
                },
                'delta': {
                    'efe_loss': best_metrics['efe_loss'] - baseline_efe_loss,
                    'boundary_sum': best_metrics['boundary'] - baseline_boundary_loss,
                    'raw_efe_avg': best_metrics['raw_efe_avg'] - (
                        sum(baseline_raw_efes) / len(baseline_raw_efes)
                        if baseline_raw_efes else 0.0
                    )
                },
                'mutations': [
                    {
                        'design_index': m['design_index'],
                        'main_index': m['main_index'],
                        'full_index': m['full_index'],
                        'codon_index': m['codon_index'],
                        'pos_in_codon': m['pos_in_codon'],
                        'ref_dna': rna_to_dna(m['ref_rna']),
                        'alt_dna': rna_to_dna(m['alt_rna'])
                    } for m in mutations
                ],
                'samples': [
                    {
                        'label': candidate['label'],
                        'main_dna': rna_to_dna(candidate['main_rna']),
                        'efe_loss': candidate['metrics']['efe_loss'],
                        'boundary_sum': candidate['metrics']['boundary'],
                        'raw_efe_avg': candidate['metrics']['raw_efe_avg']
                    } for candidate in sampled_candidates
                ],
                'loss_plot': loss_png_path,
                'multi_fasta': str(output_path)
            }

            json_path = Path(str(output_path)).with_suffix('.summary.json')
            with open(json_path, 'w') as jf:
                json.dump(summary, jf, indent=2)
            print(f"Saved structural summary to: {json_path}")
        except Exception as exc:
            print(f"Warning: failed to save structural summary: {exc}")
