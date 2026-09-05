"""
interface.py — Interface residue detection and mutation list generation.

Detects interface residues between two chains using a distance cutoff,
then generates single or double mutation lists for FoldX BuildModel.
"""

from __future__ import annotations

from itertools import combinations, product
from pathlib import Path
from typing import Iterator

from config import ALL_AAS, ALA_EQUIVALENTS, PipelineConfig
from logger import get_logger

log = get_logger(__name__)


# ── Interface detection ───────────────────────────────────────────────────────

def get_interface_residues(
    pdb_path: Path,
    chain_a: str,
    chain_b: str,
    cutoff: float,
) -> tuple[list[tuple[str, int, str]], list[tuple[str, int, str]]]:
    """
    Return interface residues on chain_a and chain_b within `cutoff` Å.

    Returns
    -------
    (residues_on_A, residues_on_B) — each item is (chain, resnum, one_letter_aa)
    """
    try:
        from Bio.PDB import PDBParser, NeighborSearch, is_aa  # type: ignore
        from Bio.PDB.Polypeptide import protein_letters_3to1  # type: ignore
    except ImportError as exc:
        raise ImportError(
            "BioPython is required for interface detection. "
            "Install with: pip install biopython"
        ) from exc

    if not pdb_path.exists():
        raise FileNotFoundError(f"PDB not found: {pdb_path}")

    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("ppi", str(pdb_path))
    except Exception as exc:
        raise RuntimeError(f"Failed to parse PDB {pdb_path}: {exc}") from exc

    model = structure[0]

    # Validate chains exist
    available_chains = [c.id for c in model.get_chains()]
    for ch in (chain_a, chain_b):
        if ch not in available_chains:
            raise ValueError(
                f"Chain '{ch}' not found in {pdb_path}. "
                f"Available chains: {available_chains}"
            )

    def _get_residues(chain_id: str) -> list:
        return [
            r for r in model[chain_id].get_residues()
            if is_aa(r, standard=True)
        ]

    res_a = _get_residues(chain_a)
    res_b = _get_residues(chain_b)

    # Build neighbour search on chain B atoms
    atoms_b = [a for r in res_b for a in r.get_atoms()]
    ns = NeighborSearch(atoms_b)

    def _one_letter(residue) -> str:
        name3 = residue.get_resname().strip()
        return protein_letters_3to1.get(name3, "X")

    def _resnum(residue) -> int:
        return residue.get_id()[1]

    # Find chain A residues with any atom within cutoff of chain B
    iface_a: list[tuple[str, int, str]] = []
    for res in res_a:
        if any(ns.search(atom.coord, cutoff) for atom in res.get_atoms()):
            aa = _one_letter(res)
            if aa != "X":
                iface_a.append((chain_a, _resnum(res), aa))

    # Symmetric: find chain B residues near chain A
    atoms_a = [a for r in res_a for a in r.get_atoms()]
    ns_a = NeighborSearch(atoms_a)
    iface_b: list[tuple[str, int, str]] = []
    for res in res_b:
        if any(ns_a.search(atom.coord, cutoff) for atom in res.get_atoms()):
            aa = _one_letter(res)
            if aa != "X":
                iface_b.append((chain_b, _resnum(res), aa))

    log.info(
        "Interface: %d residues on chain %s, %d on chain %s (cutoff %.1f Å)",
        len(iface_a), chain_a, len(iface_b), chain_b, cutoff,
    )
    return iface_a, iface_b


# ── Mutation string helpers ───────────────────────────────────────────────────

def _ala_scan_target(wt_aa: str, skip_equivalents: bool) -> str:
    """Return the target amino acid for alanine scanning."""
    if skip_equivalents and wt_aa in ALA_EQUIVALENTS:
        return "V"   # valine for Ala/Gly
    return "A"


def _mut_str(chain: str, resnum: int, wt_aa: str, target_aa: str) -> str:
    """FoldX mutation string: e.g. WA100A"""
    return f"{wt_aa}{chain}{resnum}{target_aa}"


# ── Single mutation list ──────────────────────────────────────────────────────

def generate_single_mutations(
    residues: list[tuple[str, int, str]],
    mode: str,
    skip_ala_equivalents: bool,
) -> list[str]:
    """
    Generate FoldX single-mutation strings.

    mode: 'alanine' → mutate each residue to Ala (or Val if Gly/Ala)
          'full'    → mutate each residue to all 19 other amino acids
    """
    mutations: list[str] = []
    for chain, resnum, wt_aa in residues:
        if mode == "alanine":
            target = _ala_scan_target(wt_aa, skip_ala_equivalents)
            if target != wt_aa:
                mutations.append(_mut_str(chain, resnum, wt_aa, target))
        elif mode == "full":
            for target in ALL_AAS:
                if target != wt_aa:
                    mutations.append(_mut_str(chain, resnum, wt_aa, target))
        else:
            raise ValueError(f"Unknown scan mode: '{mode}'. Use 'alanine' or 'full'.")

    log.info("Generated %d single mutations (mode=%s)", len(mutations), mode)
    return mutations


# ── Double mutation list ──────────────────────────────────────────────────────

def generate_double_mutations(
    residues: list[tuple[str, int, str]],
    mode: str,
    skip_ala_equivalents: bool,
) -> list[str]:
    """
    Generate FoldX double-mutation strings (comma-separated pairs).

    For 'alanine' mode: each residue → Ala/Val; enumerate all pairwise combos.
    For 'full' mode:    each pair × all (19 × 19) AA combinations.

    Warning: full mode scales as O(N² × 361) — use interface residues only.
    """
    def _targets(wt_aa: str) -> list[str]:
        if mode == "alanine":
            t = _ala_scan_target(wt_aa, skip_ala_equivalents)
            return [t] if t != wt_aa else []
        else:
            return [aa for aa in ALL_AAS if aa != wt_aa]

    mutations: list[str] = []

    for (c1, r1, wt1), (c2, r2, wt2) in combinations(residues, 2):
        targets1 = _targets(wt1)
        targets2 = _targets(wt2)
        if not targets1 or not targets2:
            continue
        for t1, t2 in product(targets1, targets2):
            m1 = _mut_str(c1, r1, wt1, t1)
            m2 = _mut_str(c2, r2, wt2, t2)
            mutations.append(f"{m1},{m2}")

    n_pairs = len(list(combinations(residues, 2)))
    log.info(
        "Generated %d double mutations from %d residues (%d pairs, mode=%s)",
        len(mutations), len(residues), n_pairs, mode,
    )
    if len(mutations) > 10_000:
        log.warning(
            "Large double-mutant set (%d). This will take significant compute. "
            "Consider reducing interface_cutoff_angstrom or using 'alanine' mode.",
            len(mutations),
        )
    return mutations


# ── Write individual_list.txt ─────────────────────────────────────────────────

def write_mutation_file(mutations: list[str], out_path: Path) -> Path:
    """Write FoldX individual_list.txt format (one mutation set per line, trailing ;)."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as fh:
        for mut in mutations:
            fh.write(mut + ";\n")
    log.info("Wrote %d mutations → %s", len(mutations), out_path)
    return out_path


def chunk_mutation_file(
    all_mutations: list[str],
    chunk_size: int,
    out_dir: Path,
    prefix: str = "chunk",
) -> list[tuple[Path, list[str]]]:
    """
    Split mutation list into chunks and write each to its own file.

    Returns list of (chunk_file_path, mutations_in_chunk).
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    chunks: list[tuple[Path, list[str]]] = []
    for i, start in enumerate(range(0, len(all_mutations), chunk_size)):
        chunk = all_mutations[start : start + chunk_size]
        path = out_dir / f"{prefix}_{i:04d}.txt"
        write_mutation_file(chunk, path)
        chunks.append((path, chunk))
    log.info("Split into %d chunks of ≤%d mutations", len(chunks), chunk_size)
    return chunks
