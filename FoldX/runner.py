"""
runner.py — Runs FoldX BuildModel jobs with retry logic and parallelism.

Each chunk of mutations is run as an independent FoldX process.
Failed chunks are retried up to cfg.max_retries times.
Results from all chunks are collected and returned as a unified list.
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from config import PipelineConfig
from logger import get_logger

log = get_logger(__name__)


# ── Result dataclass ──────────────────────────────────────────────────────────

@dataclass
class MutantResult:
    mutation: str          # e.g. "WA100A" or "WA100A,YB50A"
    ddG: float             # kcal/mol  (ΔΔG vs wildtype)
    sd: float              # standard deviation across FoldX runs
    raw_energies: list[float] = None   # per-run energies if available

    def __post_init__(self):
        if self.raw_energies is None:
            self.raw_energies = []


# ── Parsing FoldX output ──────────────────────────────────────────────────────

def _parse_dif_fxout(fxout_path: Path, mutation_file: Path) -> list[MutantResult]:
    """
    Parse Dif_*.fxout and map row indices back to mutation strings.
    FoldX names rows as <pdb_stem>_<mut_idx>_<run_idx>.pdb
    """
    if not fxout_path.exists():
        log.warning("Dif fxout not found: %s", fxout_path)
        return []

    # Load mutation strings from individual_list.txt
    mutations: list[str] = []
    if mutation_file.exists():
        with open(mutation_file) as fh:
            for line in fh:
                line = line.strip().rstrip(";")
                if line:
                    mutations.append(line)

    try:
        with open(fxout_path) as fh:
            lines = fh.readlines()
    except OSError as exc:
        log.error("Cannot read fxout %s: %s", fxout_path, exc)
        return []

    # Find data header
    data_start = None
    for idx, line in enumerate(lines):
        if line.strip().startswith("Pdb"):
            data_start = idx + 1
            break
    if data_start is None:
        log.warning("Could not find data header in %s", fxout_path)
        return []

    # Group rows by mutation index, average ddG across runs
    import re
    from collections import defaultdict
    mut_ddgs: dict[int, list[float]] = defaultdict(list)

    for line in lines[data_start:]:
        line = line.strip()
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        try:
            pdb_col = parts[0].strip()
            ddg = float(parts[1].strip())
            # Extract mutation index: pdb_stem_MUTIDX_RUNIDX.pdb
            m = re.search(r'_(\d+)_\d+\.pdb$', pdb_col)
            if m:
                mut_idx = int(m.group(1))
                mut_ddgs[mut_idx].append(ddg)
        except (ValueError, IndexError):
            continue

    # Build results
    results: list[MutantResult] = []
    for mut_idx, ddgs in sorted(mut_ddgs.items()):
        label_idx = mut_idx - 1
        label = mutations[label_idx] if label_idx < len(mutations) else f"mut_{mut_idx}"
        avg_ddg = sum(ddgs) / len(ddgs)
        sd = (sum((x - avg_ddg) ** 2 for x in ddgs) / max(len(ddgs) - 1, 1)) ** 0.5
        results.append(MutantResult(mutation=label, ddG=avg_ddg, sd=sd, raw_energies=ddgs))
    log.debug("Parsed %d mutations from %s (groups: %s)", len(results), fxout_path.name, list(mut_ddgs.keys())[:5])
    return results


def _parse_average_fxout(fxout_path: Path, dif_results: list[MutantResult]) -> None:
    """
    Supplement dif_results with SD from Average_*.fxout (in place).
    The Average file has the same row order as Dif.
    """
    if not fxout_path.exists():
        return
    try:
        with open(fxout_path) as fh:
            lines = fh.readlines()
    except OSError:
        return

    data_start = None
    for idx, line in enumerate(lines):
        if line.strip().startswith("Pdb"):
            data_start = idx + 1
            break
    if data_start is None:
        return

    data_lines = [l.strip() for l in lines[data_start:] if l.strip()]
    for i, (line, result) in enumerate(zip(data_lines, dif_results)):
        parts = line.split("\t")
        if len(parts) >= 3:
            try:
                result.sd = float(parts[2].strip())
            except ValueError:
                pass


# ── Single chunk runner ───────────────────────────────────────────────────────

def _run_chunk(
    foldx_bin: Path,
    repaired_pdb: Path,
    mutation_file: Path,
    run_dir: Path,
    n_runs: int,
    chains: tuple[str, str],
    water: str,
    ion_strength: float,
    write_pdbs: bool,
    max_retries: int,
    retry_delay: float,
) -> list[MutantResult]:
    """Run FoldX BuildModel on one mutation chunk; return parsed results."""
# Skip if already completed with matching mutations
    if run_dir.exists():
        existing_dif = list(run_dir.glob("Dif_*.fxout"))
        existing_mut_file = run_dir / "individual_list.txt"
        if existing_dif and existing_mut_file.exists():
            # Check mutation file matches what we planned
            planned = [l.strip() for l in mutation_file.read_text().splitlines() if l.strip()]
            existing = [l.strip() for l in existing_mut_file.read_text().splitlines() if l.strip()]
            if planned == existing:
                results = _parse_dif_fxout(existing_dif[0], existing_mut_file)
                if results:
                    log.info("Chunk already complete and verified, skipping: %s", run_dir.name)
                    return results
                    
    # Clean any previous run artifacts before starting
    if run_dir.exists():
        shutil.rmtree(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)

    pdb_name = repaired_pdb.name
    shutil.copy2(repaired_pdb, run_dir / pdb_name)
    shutil.copy2(mutation_file, run_dir / "individual_list.txt")

    rotabase = foldx_bin.parent / "rotabase.txt"
    if rotabase.exists():
        shutil.copy2(rotabase, run_dir / "rotabase.txt")

    analyse_complex_str = ",".join(chains)

    cmd = [
        str(foldx_bin.resolve()),
        "--command=BuildModel",
        f"--pdb={pdb_name}",
        f"--mutant-file={str((run_dir / 'individual_list.txt').resolve())}",
        f"--numberOfRuns={n_runs}",
        f"--water={water}",
        f"--ionStrength={ion_strength}",
        f"--out-pdb={'true' if write_pdbs else 'false'}",
        f"--pdb-dir={str(run_dir.resolve())}",
        f"--output-dir={str(run_dir.resolve())}",
    ]

    last_stderr = ""
    for attempt in range(1, max_retries + 1):
        log.debug(
            "Chunk %s — attempt %d/%d", mutation_file.stem, attempt, max_retries
        )
        log.info("Running command: %s", " ".join(cmd))
        proc = subprocess.run(
            cmd, cwd=run_dir, text=True,
            stdout=None, stderr=subprocess.PIPE,
        )
        
        last_stderr = proc.stderr

        if proc.returncode != 0:
            log.warning(
                "FoldX non-zero exit (rc=%d) for chunk %s, attempt %d: %s",
                proc.returncode, mutation_file.stem, attempt,
                last_stderr.strip()[-300:],
            )
            if attempt < max_retries:
                time.sleep(retry_delay * attempt)   # exponential back-off
            continue

        # Locate output files — FoldX names them Dif_<pdb_stem>.fxout
        dif_files = list(run_dir.glob("Dif_*.fxout"))
        avg_files = list(run_dir.glob("Average_*.fxout"))
        dif_file = dif_files[0] if dif_files else run_dir / f"Dif_{repaired_pdb.stem}.fxout"
        avg_file = avg_files[0] if avg_files else run_dir / f"Average_{repaired_pdb.stem}.fxout"

        if not dif_file.exists():
            log.warning("Expected Dif file not found: %s", dif_file)
            if attempt < max_retries:
                time.sleep(retry_delay * attempt)
            continue

        results = _parse_dif_fxout(dif_file, run_dir / "individual_list.txt")
        _parse_average_fxout(avg_file, results)

        if not results:
            log.warning("Parsed 0 results from %s — may be empty chunk", dif_file)

        log.debug("Chunk %s → %d results", mutation_file.stem, len(results))
        return results

    log.error(
        "Chunk %s failed after %d attempts. Last stderr:\n%s",
        mutation_file.stem, max_retries, last_stderr.strip()[-500:],
    )
    return []   # Return empty rather than crashing the whole pipeline


# ── Parallel chunk runner ─────────────────────────────────────────────────────

def _run_chunk_worker(args: dict) -> list[MutantResult]:
    """Top-level worker for ProcessPoolExecutor (must be picklable)."""
    try:
        return _run_chunk(**args)
    except Exception as exc:
        # Log inside worker; return empty so parent can continue
        import traceback
        print(f"[worker error] {exc}\n{traceback.format_exc()}")
        return []


def run_buildmodel_parallel(
    cfg: PipelineConfig,
    repaired_pdb: Path,
    chunks: list[tuple[Path, list[str]]],
) -> list[MutantResult]:
    """
    Run BuildModel in parallel across chunks.

    Returns aggregated list of MutantResult objects.
    """
    if not chunks:
        log.warning("No mutation chunks provided — nothing to run.")
        return []

    jobs = []
    for i, (chunk_file, _mutations) in enumerate(chunks):
        pdb_stem = repaired_pdb.stem
        run_dir = cfg.work_dir / "buildmodel_runs" / pdb_stem / f"chunk_{i:04d}"
        jobs.append({
            "foldx_bin": cfg.foldx_bin,
            "repaired_pdb": repaired_pdb,
            "mutation_file": chunk_file,
            "run_dir": run_dir,
            "n_runs": cfg.foldx_runs,
            "chains": (cfg.chain_A, cfg.chain_B),
            "water": cfg.water_mode,
            "ion_strength": cfg.ion_strength,
            "write_pdbs": cfg.write_output_pdbs,
            "max_retries": cfg.max_retries,
            "retry_delay": cfg.retry_delay_s,
        })

    all_results: list[MutantResult] = []
    n_workers = min(cfg.n_workers, len(jobs))

    log.info("Running %d chunks on %d workers …", len(jobs), n_workers)

    if n_workers <= 1:
        # Serial fallback — easier to debug
        for job in jobs:
            all_results.extend(_run_chunk(**job))
    else:
        with ProcessPoolExecutor(max_workers=n_workers) as pool:
            futures = {pool.submit(_run_chunk_worker, job): i
                       for i, job in enumerate(jobs)}
            completed = 0
            for future in as_completed(futures):
                chunk_idx = futures[future]
                completed += 1
                try:
                    results = future.result()
                    all_results.extend(results)
                    log.info(
                        "Chunk %d/%d done (%d results)",
                        completed, len(jobs), len(results),
                    )
                except Exception as exc:
                    log.error("Chunk %d raised an exception: %s", chunk_idx, exc)

    log.info("Total results collected: %d", len(all_results))
    return all_results
