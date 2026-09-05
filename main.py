"""
main.py — FoldX double-mutant scanning pipeline for AF3 PPI structures.

Usage
-----
python main.py \
    --cif model.cif \
    --chain-a A \
    --chain-b B \
    --mode alanine \
    --confidence-json confidences.json \
    --foldx-bin /path/to/foldx \
    --workers 8

For full amino acid scanning (expensive):
    --mode full

Run python main.py --help for all options.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from config import PipelineConfig
from logger import get_logger
from structure_prep import cif_to_pdb, iterative_repair, load_plddt, filter_by_plddt
from interface import (
    get_interface_residues,
    generate_single_mutations,
    generate_double_mutations,
    chunk_mutation_file,
)
from runner import run_buildmodel_parallel
from analysis import (
    compute_epistasis,
    save_single_results,
    save_double_results,
    print_summary,
)


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="FoldX double-mutant scanning on AlphaFold 3 PPI structures."
    )
    p.add_argument("--cif", type=Path, help="AF3 .cif input file")
    p.add_argument("--pdb", type=Path, default=None,
                   help="Pre-converted .pdb (skips CIF conversion)")
    p.add_argument("--chain-a", default="A", help="Receptor chain ID (default: A)")
    p.add_argument("--chain-b", default="B", help="Ligand chain ID (default: B)")
    p.add_argument("--mode", choices=["alanine", "full"], default="alanine",
                   help="Scan mode (default: alanine)")
    p.add_argument("--confidence-json", type=Path, default=None,
                   help="AF3 confidence JSON for pLDDT filtering")
    p.add_argument("--plddt-min", type=float, default=70.0,
                   help="Minimum pLDDT to include a residue (default: 70)")
    p.add_argument("--cutoff", type=float, default=5.0,
                   help="Interface distance cutoff in Å (default: 5.0)")
    p.add_argument("--foldx-bin", type=Path, default=Path("foldx"),
                   help="Path to FoldX binary (default: foldx in PATH)")
    p.add_argument("--work-dir", type=Path, default=Path("foldx_run"),
                   help="Working directory for intermediate files")
    p.add_argument("--results-dir", type=Path, default=Path("results"),
                   help="Output directory for CSV results")
    p.add_argument("--repair-iters", type=int, default=3,
                   help="RepairPDB iterations (default: 3)")
    p.add_argument("--foldx-runs", type=int, default=5,
                   help="FoldX runs per mutant for averaging (default: 5)")
    p.add_argument("--workers", type=int, default=4,
                   help="Parallel FoldX worker processes (default: 4)")
    p.add_argument("--chunk-size", type=int, default=200,
                   help="Mutations per FoldX call (default: 200)")
    p.add_argument("--epistasis-threshold", type=float, default=0.5,
                   help="Min |epistasis| to flag coupling (default: 0.5 kcal/mol)")
    p.add_argument("--no-double", action="store_true",
                   help="Only run single mutant scan")
    p.add_argument("--singles-only-repair", action="store_true",
                   help="Re-use existing repaired PDB (skip repair step)")
    p.add_argument("--log-file", type=Path, default=None,
                   help="Optional log file path")
    return p.parse_args()


# ── Validation ────────────────────────────────────────────────────────────────

def validate_inputs(cfg: PipelineConfig, args: argparse.Namespace) -> None:
    errors = []

    if cfg.input_pdb is None and not cfg.input_cif.exists():
        errors.append(f"CIF not found: {cfg.input_cif}")
    if cfg.input_pdb is not None and not cfg.input_pdb.exists():
        errors.append(f"PDB not found: {cfg.input_pdb}")

    import shutil as _sh
    if _sh.which(str(cfg.foldx_bin)) is None and not cfg.foldx_bin.exists():
        errors.append(
            f"FoldX binary not found: {cfg.foldx_bin}. "
            "Set --foldx-bin or add it to PATH."
        )

    if not (1 <= cfg.repair_iterations <= 10):
        errors.append("--repair-iters must be between 1 and 10.")

    if cfg.chain_A == cfg.chain_B:
        errors.append("--chain-a and --chain-b must be different.")

    if errors:
        for e in errors:
            print(f"[ERROR] {e}", file=sys.stderr)
        sys.exit(1)


# ── Main pipeline ─────────────────────────────────────────────────────────────

def run(cfg: PipelineConfig, args: argparse.Namespace) -> None:
    log = get_logger("main", log_file=args.log_file)
    log.info("=" * 60)
    log.info("FoldX PPI mutagenesis pipeline starting")
    log.info("  Mode          : %s", cfg.scan_mode)
    log.info("  Double muts   : %s", cfg.double_mutations)
    log.info("  Chains        : %s / %s", cfg.chain_A, cfg.chain_B)
    log.info("  Workers       : %d", cfg.n_workers)
    log.info("=" * 60)

    cfg.work_dir.mkdir(parents=True, exist_ok=True)
    cfg.results_dir.mkdir(parents=True, exist_ok=True)

    # ── Step 1: Get a PDB ────────────────────────────────────────────────────
    if cfg.input_pdb:
        raw_pdb = cfg.input_pdb
        log.info("Using provided PDB: %s", raw_pdb)
    else:
        log.info("Converting CIF → PDB …")
        raw_pdb = cif_to_pdb(cfg.input_cif, cfg.work_dir / "converted")

    # ── Step 2: Repair ───────────────────────────────────────────────────────
    repaired_pdb_marker = cfg.work_dir / "repaired_final.pdb"
    if args.singles_only_repair and repaired_pdb_marker.exists():
        log.info("Reusing existing repaired PDB: %s", repaired_pdb_marker)
        repaired_pdb = repaired_pdb_marker
    else:
        log.info("Running %d RepairPDB iterations …", cfg.repair_iterations)
        repaired_pdb = iterative_repair(cfg, raw_pdb)
        import shutil
        shutil.copy2(repaired_pdb, repaired_pdb_marker)

    # ── Step 3: pLDDT filter ─────────────────────────────────────────────────
    plddt_map = load_plddt(args.confidence_json)

    # ── Step 4: Interface residues ───────────────────────────────────────────
    log.info("Detecting interface residues (cutoff=%.1f Å) …", cfg.interface_cutoff_angstrom)
    iface_a, iface_b = get_interface_residues(
        repaired_pdb, cfg.chain_A, cfg.chain_B, cfg.interface_cutoff_angstrom
    )
    all_interface = iface_a + iface_b

    if not all_interface:
        log.error(
            "No interface residues found between chains %s and %s. "
            "Check your chain IDs and cutoff.",
            cfg.chain_A, cfg.chain_B,
        )
        sys.exit(1)

    all_interface = filter_by_plddt(all_interface, plddt_map, cfg.plddt_min)

    if not all_interface:
        log.error("No residues survived the pLDDT filter (threshold=%.1f).", cfg.plddt_min)
        sys.exit(1)

    # ── Step 5: Generate single mutations ────────────────────────────────────
    log.info("Generating single mutation list …")
    single_muts = generate_single_mutations(
        all_interface, cfg.scan_mode, cfg.skip_ala_equivalents
    )

    if not single_muts:
        log.error("No single mutations generated — check scan mode and residue list.")
        sys.exit(1)

    single_chunks = chunk_mutation_file(
        single_muts,
        cfg.chunk_size,
        cfg.work_dir / "single_chunks",
        prefix="single",
    )

    # ── Step 6: Run single scan ──────────────────────────────────────────────
    log.info("Running single mutant scan (%d mutations) …", len(single_muts))
    single_results = run_buildmodel_parallel(cfg, repaired_pdb, single_chunks)

    if not single_results:
        log.error(
            "No single mutant results returned. "
            "Check FoldX output in %s/buildmodel_runs/", cfg.work_dir
        )
        sys.exit(1)

    save_single_results(single_results, cfg.results_dir / "single_mutants.csv")

    # ── Step 7: Double mutations (optional) ──────────────────────────────────
    epistasis_results = []

    if cfg.double_mutations:
        log.info("Generating double mutation list …")
        double_muts = generate_double_mutations(
            all_interface, cfg.scan_mode, cfg.skip_ala_equivalents
        )

        if double_muts:
            double_chunks = chunk_mutation_file(
                double_muts,
                cfg.chunk_size,
                cfg.work_dir / "double_chunks",
                prefix="double",
            )

            log.info("Running double mutant scan (%d mutations) …", len(double_muts))
            double_results = run_buildmodel_parallel(cfg, repaired_pdb, double_chunks)

            if double_results:
                epistasis_results = compute_epistasis(
                    single_results, double_results, cfg.epistasis_threshold
                )
                save_double_results(
                    double_results,
                    epistasis_results,
                    cfg.results_dir / "double_mutants.csv",
                )
            else:
                log.warning("No double mutant results returned.")
        else:
            log.warning("No double mutations generated (too few interface residues?).")

    # ── Step 8: Summary ──────────────────────────────────────────────────────
    print_summary(single_results, epistasis_results)
    log.info("Results written to: %s", cfg.results_dir)
    log.info("Pipeline complete.")


# ── Entry point ───────────────────────────────────────────────────────────────

def main() -> None:
    args = parse_args()

    cfg = PipelineConfig(
        input_cif=args.cif or Path("model.cif"),
        input_pdb=args.pdb,
        foldx_bin=args.foldx_bin,
        work_dir=args.work_dir,
        results_dir=args.results_dir,
        repair_iterations=args.repair_iters,
        chain_A=args.chain_a,
        chain_B=args.chain_b,
        interface_cutoff_angstrom=args.cutoff,
        plddt_min=args.plddt_min,
        scan_mode=args.mode,
        double_mutations=not args.no_double,
        foldx_runs=args.foldx_runs,
        n_workers=args.workers,
        chunk_size=args.chunk_size,
        epistasis_threshold=args.epistasis_threshold,
    )

    validate_inputs(cfg, args)
    run(cfg, args)


if __name__ == "__main__":
    main()
