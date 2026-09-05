"""
structure_prep.py — AF3 CIF → repaired PDB, with pLDDT filtering.

Steps
-----
1. Convert .cif → .pdb  (via gemmi, with biopython fallback)
2. Iterative FoldX RepairPDB  (config.repair_iterations passes)
3. Load pLDDT scores from AF3 JSON confidence file
4. Return the repaired PDB path + high-confidence residue set
"""

from __future__ import annotations

import json
import shutil
import subprocess
import time
from pathlib import Path
from typing import Optional

from config import PipelineConfig
from logger import get_logger

log = get_logger(__name__)


# ── 1. CIF → PDB ─────────────────────────────────────────────────────────────

def _convert_cif_gemmi(cif_path: Path, out_pdb: Path) -> bool:
    """Try gemmi first — fastest and most faithful."""
    try:
        import gemmi  # noqa: F401
        result = subprocess.run(
            ["gemmi", "convert", "--to", "pdb", str(cif_path), str(out_pdb)],
            capture_output=True, text=True,
        )
        if result.returncode == 0 and out_pdb.exists():
            log.debug("gemmi conversion succeeded")
            return True
        log.warning("gemmi returned non-zero: %s", result.stderr.strip())
    except (ImportError, FileNotFoundError):
        log.debug("gemmi not available, trying BioPython")
    return False


def _convert_cif_biopython(cif_path: Path, out_pdb: Path) -> bool:
    """BioPython fallback."""
    try:
        from Bio.PDB import MMCIFParser, PDBIO  # type: ignore
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("af3", str(cif_path))
        io = PDBIO()
        io.set_structure(structure)
        io.save(str(out_pdb))
        log.debug("BioPython conversion succeeded")
        return True
    except Exception as exc:
        log.error("BioPython conversion failed: %s", exc)
    return False


def cif_to_pdb(cif_path: Path, out_dir: Path) -> Path:
    """Convert AF3 CIF to PDB; raise RuntimeError on failure."""
    if not cif_path.exists():
        raise FileNotFoundError(f"CIF not found: {cif_path}")
    out_dir.mkdir(parents=True, exist_ok=True)
    out_pdb = out_dir / (cif_path.stem + ".pdb")

    if out_pdb.exists():
        log.info("PDB already exists, skipping conversion: %s", out_pdb)
        return out_pdb

    if _convert_cif_gemmi(cif_path, out_pdb):
        return out_pdb
    if _convert_cif_biopython(cif_path, out_pdb):
        return out_pdb

    raise RuntimeError(
        f"Could not convert {cif_path} to PDB. "
        "Install gemmi (`pip install gemmi`) or biopython."
    )


# ── 2. Iterative RepairPDB ────────────────────────────────────────────────────

def _run_foldx_repair(foldx_bin: Path, pdb_path: Path, work_dir: Path,
                       ph: float, water: str, max_retries: int,
                       retry_delay: float) -> Path:
    """Run a single RepairPDB pass; return path to repaired PDB."""
    pdb_name = pdb_path.name
    repaired_name = pdb_path.stem + "_Repair.pdb"

    # FoldX writes output next to itself unless --output-dir is given
    run_dir = work_dir / f"repair_{pdb_path.stem}"
    run_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(pdb_path, run_dir / pdb_name)

    cmd = [
        str(foldx_bin),
        "--command=RepairPDB",
        f"--pdb={pdb_name}",
        f"--pH={ph}",
        f"--water={water}",
        "--ionStrength=0.05",
        "--pdb-dir=.",
        "--output-dir=.",
    ]

    for attempt in range(1, max_retries + 1):
        log.debug("RepairPDB attempt %d/%d: %s", attempt, max_retries, pdb_name)
        result = subprocess.run(
            cmd, cwd=run_dir,
            capture_output=True, text=True,
        )
        repaired = run_dir / repaired_name
        if result.returncode == 0 and repaired.exists():
            log.debug("RepairPDB succeeded → %s", repaired)
            return repaired
        log.warning(
            "RepairPDB attempt %d failed (rc=%d): %s",
            attempt, result.returncode,
            result.stderr.strip()[-300:],
        )
        if attempt < max_retries:
            time.sleep(retry_delay)

    raise RuntimeError(
        f"RepairPDB failed after {max_retries} attempts for {pdb_path}. "
        f"Last stderr:\n{result.stderr.strip()}"
    )


def iterative_repair(cfg: PipelineConfig, pdb_path: Path) -> Path:
    """Run RepairPDB cfg.repair_iterations times, chaining outputs."""
    current = pdb_path
    for i in range(1, cfg.repair_iterations + 1):
        log.info("RepairPDB pass %d/%d …", i, cfg.repair_iterations)
        out_dir = cfg.work_dir / f"repair_pass_{i}"
        repaired = _run_foldx_repair(
            foldx_bin=cfg.foldx_bin,
            pdb_path=current,
            work_dir=out_dir,
            ph=cfg.ph,
            water=cfg.water_mode,
            max_retries=cfg.max_retries,
            retry_delay=cfg.retry_delay_s,
        )
        # Copy to work_dir root for the next iteration
        next_pdb = cfg.work_dir / f"repaired_pass{i}.pdb"
        shutil.copy2(repaired, next_pdb)
        current = next_pdb
        log.info("Pass %d complete → %s", i, current)

    log.info("Final repaired structure: %s", current)
    return current


# ── 3. pLDDT confidence filtering ────────────────────────────────────────────

def load_plddt(confidence_json: Optional[Path]) -> dict[tuple[str, int], float]:
    """
    Parse AF3 confidences JSON → {(chain_id, res_seq): pLDDT}.

    AF3 JSON structure (subject to change between AF3 versions):
    {
      "atom_chain_ids": ["A", "A", ...],
      "atom_res_ids":   [1, 1, ...],
      "atom_plddts":    [85.2, 85.2, ...],
      ...
    }
    Returns per-residue *mean* pLDDT.
    """
    if confidence_json is None or not confidence_json.exists():
        log.warning(
            "No confidence JSON found — pLDDT filtering disabled. "
            "Pass --confidence-json to enable."
        )
        return {}

    try:
        with open(confidence_json) as fh:
            data = json.load(fh)
    except (json.JSONDecodeError, OSError) as exc:
        log.error("Could not read confidence JSON: %s", exc)
        return {}

    required_keys = {"atom_chain_ids", "atom_res_ids", "atom_plddts"}
    if not required_keys.issubset(data):
        log.warning(
            "Confidence JSON missing keys %s — pLDDT filtering disabled.",
            required_keys - set(data),
        )
        return {}

    chains = data["atom_chain_ids"]
    res_ids = data["atom_res_ids"]
    plddts = data["atom_plddts"]

    if not (len(chains) == len(res_ids) == len(plddts)):
        log.error("Confidence JSON arrays have mismatched lengths — skipping.")
        return {}

    from collections import defaultdict
    bucket: dict[tuple[str, int], list[float]] = defaultdict(list)
    for chain, rid, score in zip(chains, res_ids, plddts):
        bucket[(str(chain), int(rid))].append(float(score))

    return {key: sum(vals) / len(vals) for key, vals in bucket.items()}


def filter_by_plddt(
    residues: list[tuple[str, int, str]],   # [(chain, resnum, wt_aa), ...]
    plddt_map: dict[tuple[str, int], float],
    threshold: float,
) -> list[tuple[str, int, str]]:
    """Remove residues below the pLDDT threshold."""
    if not plddt_map:
        return residues  # no filtering data — return all

    filtered, dropped = [], []
    for chain, resnum, aa in residues:
        score = plddt_map.get((chain, resnum), 100.0)
        if score >= threshold:
            filtered.append((chain, resnum, aa))
        else:
            dropped.append((chain, resnum, aa, score))

    if dropped:
        log.info(
            "Dropped %d low-confidence residues (pLDDT < %.1f): %s",
            len(dropped),
            threshold,
            [(c, r, a) for c, r, a, _ in dropped[:10]],
        )
    log.info("%d residues pass pLDDT filter (≥ %.1f)", len(filtered), threshold)
    return filtered
