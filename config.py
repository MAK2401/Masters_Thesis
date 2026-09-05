"""
config.py — Central configuration for the FoldX PPI mutagenesis pipeline.
Edit this file to match your environment before running.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


# ── Amino acid constants ──────────────────────────────────────────────────────

ALL_AAS = "ACDEFGHIKLMNPQRSTVWY"

# Residues that are already effectively alanine in a scan — use valine instead
ALA_EQUIVALENTS = {"A", "G"}

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

ONE_TO_THREE = {v: k for k, v in THREE_TO_ONE.items()}


# ── Pipeline configuration ────────────────────────────────────────────────────

@dataclass
class PipelineConfig:
    # Paths
    input_cif: Path = Path("model.cif")
    input_pdb: Optional[Path] = None          # skip conversion if already PDB
    foldx_bin: Path = Path("foldx")
    work_dir: Path = Path("foldx_run")
    results_dir: Path = Path("results")

    # Structure prep
    repair_iterations: int = 3                # RepairPDB passes
    ph: float = 7.0
    ion_strength: float = 0.05               # mol/L
    water_mode: str = "IGNORE"               # IGNORE | CRYSTAL | PREDICT

    # Interface detection
    chain_A: str = "A"                        # receptor chain
    chain_B: str = "B"                        # ligand chain
    interface_cutoff_angstrom: float = 5.0
    plddt_min: float = 70.0                  # ignore low-confidence residues

    # Mutagenesis
    scan_mode: str = "alanine"               # alanine | full
    double_mutations: bool = True
    skip_ala_equivalents: bool = True        # skip A/G in alanine scan
    foldx_runs: int = 5                      # repeats per mutant (average)
    write_output_pdbs: bool = False          # False = much faster

    # Parallelism
    n_workers: int = 4
    chunk_size: int = 200                    # mutations per FoldX call

    # Retry / robustness
    max_retries: int = 3
    retry_delay_s: float = 2.0

    # Epistasis threshold for reporting
    epistasis_threshold: float = 0.5        # kcal/mol

    def __post_init__(self):
        self.work_dir = Path(self.work_dir)
        self.results_dir = Path(self.results_dir)
        self.foldx_bin = Path(self.foldx_bin)
        if self.input_pdb:
            self.input_pdb = Path(self.input_pdb)
        self.input_cif = Path(self.input_cif)
