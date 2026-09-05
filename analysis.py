"""
analysis.py — Post-processing: epistasis calculation, ranking, CSV export.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from logger import get_logger
from runner import MutantResult

log = get_logger(__name__)


# ── Epistasis ─────────────────────────────────────────────────────────────────

@dataclass
class EpistasisResult:
    mutation_pair: str           # "WA100A,YB50A"
    ddG_double: float            # ΔΔG of double mutant
    ddG_single_1: Optional[float]
    ddG_single_2: Optional[float]
    epistasis: Optional[float]   # ΔΔG_AB − ΔΔG_A − ΔΔG_B
    coupled: bool = False        # |epistasis| > threshold


def compute_epistasis(
    single_results: list[MutantResult],
    double_results: list[MutantResult],
    threshold: float = 0.5,
) -> list[EpistasisResult]:
    """
    Compute coupling energy for every double mutant.

    Epistasis = ΔΔG_double − ΔΔG_A − ΔΔG_B
      > 0 : synergistic destabilisation (both residues matter together)
      < 0 : one mutation partially rescues the other
      ≈ 0 : additive (independent contributions)
    """
    single_map: dict[str, float] = {r.mutation: r.ddG for r in single_results}

    epistasis_results: list[EpistasisResult] = []
    missing_singles = 0

    for dr in double_results:
        parts = dr.mutation.split(",")
        if len(parts) != 2:
            log.debug("Skipping non-double mutation: %s", dr.mutation)
            continue

        m1, m2 = parts[0].strip(), parts[1].strip()
        ddg1 = single_map.get(m1)
        ddg2 = single_map.get(m2)

        if ddg1 is None or ddg2 is None:
            missing_singles += 1
            epi = None
            coupled = False
        else:
            epi = dr.ddG - ddg1 - ddg2
            coupled = abs(epi) >= threshold

        epistasis_results.append(EpistasisResult(
            mutation_pair=dr.mutation,
            ddG_double=dr.ddG,
            ddG_single_1=ddg1,
            ddG_single_2=ddg2,
            epistasis=epi,
            coupled=coupled,
        ))

    if missing_singles:
        log.warning(
            "%d double mutants had no matching single mutant data "
            "(epistasis = None for those).",
            missing_singles,
        )

    coupled_count = sum(1 for r in epistasis_results if r.coupled)
    log.info(
        "Epistasis computed: %d pairs, %d coupled (|ε| ≥ %.2f kcal/mol)",
        len(epistasis_results), coupled_count, threshold,
    )
    return epistasis_results


# ── Hot spot identification ───────────────────────────────────────────────────

def rank_single_hotspots(
    results: list[MutantResult],
    ddg_threshold: float = 1.0,
) -> list[MutantResult]:
    """
    Sort single mutants by ΔΔG descending.
    Positive ΔΔG = destabilising = hot spot candidate.
    """
    ranked = sorted(results, key=lambda r: r.ddG, reverse=True)
    hotspots = [r for r in ranked if r.ddG >= ddg_threshold]
    log.info(
        "Single mutant hot spots (ΔΔG ≥ %.1f kcal/mol): %d / %d",
        ddg_threshold, len(hotspots), len(results),
    )
    return ranked


# ── CSV export ────────────────────────────────────────────────────────────────

def save_single_results(results: list[MutantResult], out_path: Path) -> None:
    """Write single mutant results to CSV."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["mutation", "ddG_kcal_mol", "sd"])
        for r in results:
            writer.writerow([r.mutation, f"{r.ddG:.4f}", f"{r.sd:.4f}"])
    log.info("Single mutant results → %s (%d rows)", out_path, len(results))


def save_double_results(
    results: list[MutantResult],
    epistasis: list[EpistasisResult],
    out_path: Path,
) -> None:
    """Write double mutant results + epistasis to CSV."""
    out_path.parent.mkdir(parents=True, exist_ok=True)

    epi_map = {e.mutation_pair: e for e in epistasis}

    with open(out_path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow([
            "mutation_pair",
            "ddG_double",
            "sd_double",
            "ddG_single_1",
            "ddG_single_2",
            "epistasis_kcal_mol",
            "coupled",
        ])
        for r in sorted(results, key=lambda x: x.ddG, reverse=True):
            e = epi_map.get(r.mutation)
            writer.writerow([
                r.mutation,
                f"{r.ddG:.4f}",
                f"{r.sd:.4f}",
                f"{e.ddG_single_1:.4f}" if e and e.ddG_single_1 is not None else "NA",
                f"{e.ddG_single_2:.4f}" if e and e.ddG_single_2 is not None else "NA",
                f"{e.epistasis:.4f}" if e and e.epistasis is not None else "NA",
                str(e.coupled) if e else "NA",
            ])
    log.info("Double mutant results → %s (%d rows)", out_path, len(results))


def print_summary(
    single_results: list[MutantResult],
    epistasis_results: list[EpistasisResult],
    top_n: int = 10,
) -> None:
    """Print a brief human-readable summary to stdout."""
    print("\n" + "=" * 60)
    print(f"  TOP {top_n} SINGLE MUTANT HOT SPOTS")
    print("=" * 60)
    print(f"  {'Mutation':<20} {'ΔΔG (kcal/mol)':>16} {'SD':>8}")
    print("-" * 60)
    for r in sorted(single_results, key=lambda x: x.ddG, reverse=True)[:top_n]:
        print(f"  {r.mutation:<20} {r.ddG:>16.3f} {r.sd:>8.3f}")

    coupled = [e for e in epistasis_results if e.coupled and e.epistasis is not None]
    coupled.sort(key=lambda e: abs(e.epistasis), reverse=True)

    print(f"\n{'=' * 60}")
    print(f"  TOP {top_n} COUPLED DOUBLE MUTANTS (by |epistasis|)")
    print("=" * 60)
    print(f"  {'Pair':<35} {'ΔΔG_AB':>8} {'ΔΔG_A':>7} {'ΔΔG_B':>7} {'ε':>8}")
    print("-" * 60)
    for e in coupled[:top_n]:
        print(
            f"  {e.mutation_pair:<35} "
            f"{e.ddG_double:>8.3f} "
            f"{e.ddG_single_1:>7.3f} "
            f"{e.ddG_single_2:>7.3f} "
            f"{e.epistasis:>8.3f}"
        )
    print("=" * 60 + "\n")
