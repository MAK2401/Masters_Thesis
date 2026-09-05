# FoldX PPI Double-Mutant Scanning Pipeline

Stable, parallelised pipeline for running FoldX mutagenesis scanning
(single and double mutants) on AlphaFold 3 protein–protein interaction structures.

---

## File overview

```
foldx_ppi/
├── config.py          # All tunable parameters in one place
├── logger.py          # Centralised logging (console + file)
├── structure_prep.py  # CIF→PDB conversion, RepairPDB, pLDDT filter
├── interface.py       # Interface detection, mutation list generation
├── runner.py          # Parallel FoldX execution with retry logic
├── analysis.py        # Epistasis calculation, CSV export, summary
└── main.py            # CLI entry point
```

---

## Installation

```bash
pip install biopython gemmi
# FoldX must be separately downloaded from https://foldxsuite.crg.eu/
```

---

## Quick start

### Alanine scanning (singles + doubles)
```bash
python main.py \
  --cif model.cif \
  --chain-a A \
  --chain-b B \
  --mode alanine \
  --confidence-json confidences.json \
  --foldx-bin /path/to/foldx \
  --workers 8
```

### Skip CIF conversion (already have a PDB)
```bash
python main.py \
  --pdb model.pdb \
  --chain-a A --chain-b B \
  --foldx-bin ./foldx
```

### Full amino acid scan (expensive — interface residues only)
```bash
python main.py \
  --cif model.cif \
  --mode full \
  --workers 16 \
  --chunk-size 100
```

### Single mutants only
```bash
python main.py --cif model.cif --no-double
```

---

## Output files

```
results/
├── single_mutants.csv    # mutation, ddG, sd
└── double_mutants.csv    # pair, ddG_double, ddG_A, ddG_B, epistasis, coupled
```

**Interpreting results:**

| Column | Meaning |
|---|---|
| `ddG` | ΔΔG in kcal/mol. Positive = destabilising = hot spot |
| `epistasis` | ΔΔG_AB − ΔΔG_A − ΔΔG_B. Non-zero = coupling |
| `coupled` | True if \|epistasis\| ≥ threshold (default 0.5 kcal/mol) |

---

## Stability features

| Feature | Where |
|---|---|
| Multi-pass RepairPDB | `structure_prep.iterative_repair` |
| Retry with exponential back-off | `runner._run_chunk` |
| pLDDT confidence filter | `structure_prep.filter_by_plddt` |
| Graceful empty-chunk handling | `runner.run_buildmodel_parallel` |
| Ala/Gly → Val in alanine scan | `interface._ala_scan_target` |
| Input validation before any work | `main.validate_inputs` |
| Centralised logging to file | `logger.get_logger` |
| Serial fallback when workers=1 | `runner.run_buildmodel_parallel` |

---

## Caveats

- FoldX ΔΔG values on AF3 models carry ~1–2 kcal/mol uncertainty.
  Use results for **ranking**, not as absolute energetics.
- Low-pLDDT regions (< 70) are excluded by default — adjust `--plddt-min`.
- Full double scanning scales as O(N² × 361) — use with care.
  For > 20 interface residues, alanine mode is strongly recommended first.
- FoldX must be licensed separately from https://foldxsuite.crg.eu/
