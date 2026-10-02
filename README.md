# Contact-graph definitions shape what a model describes but not what it predicts: a controlled comparison across five encoder capacities

This repository contains the complete, reproducible analysis pipeline for a
controlled benchmark of **residue-level contact-graph representations** used to
predict protein stability changes (ΔΔG) upon single-point mutation.

## What this study asks

The residue contact graph is a standard ingredient of structure-based ΔΔG
predictors. It is defined by choosing **one representative atom per residue** and
a **distance cutoff** — a choice that determines whether mutation-induced
side-chain changes are visible to the model. This repository evaluates that
choice explicitly.

The analysis separates **three questions that are usually conflated**:

| Question | Depends on model capacity? | Depends on FoldX? |
|---|---|---|
| **(A)** What structural change does a definition *report*? | No — a structural fact | Yes |
| **(B)** Does the definition change *prediction*? | Yes — tested across a capacity ladder | No |
| **(C)** Does structure add increment over a sequence baseline? | Yes | No |

Separating these matters: a weak encoder can only threaten claim (B), never (A).

## Design

- **Four graph definitions** — Cα, Cβ, side-chain centroid (heavy-atom centroid
  excluding hydrogens), all-atom minimum inter-residue distance
- **A five-rung capacity ladder** — from a 5.8 K-parameter mean-pooled GCN to a
  509 K-parameter E(3)-equivariant network:

  | Rung | Encoder | Parameters |
  |---|---|---|
  | 1 | GCN + global mean pooling | 5.8 K |
  | 2 | GCN + mutation-site local pooling | 5.9 K |
  | 3 | Shallow edge-aware GINEConv | 50 K |
  | 4 | Deep GINE + attention pooling + local readout | 460 K |
  | 5 | E(3)-equivariant EGNN (3D coordinates) | 509 K |

- **Two benchmarks** — S669 (common intersection of 511 mutations / 88 proteins)
  and ssym (342 mutations / 15 proteins, not used during development)
- **Two modelling engines cross-checked** — FoldX BuildModel and SCWRL4
- **Sequence baselines** — ESM-2 zero-shot and ESM-2 embeddings with a supervised head

## Repository layout

```
.
├── src/                        # analysis code
│   ├── paths.py                # relocatable path handling (GED_ROOT override)
│   ├── seed_utils.py           # full determinism control
│   ├── contact_graph_defs.py   # the four residue-level graph definitions
│   ├── ladder_common.py        # shared data loading for the capacity ladder
│   ├── strong_backbones.py     # DeepGINE and EGNN encoders
│   ├── run_ladder.py           # capacity-ladder experiment
│   ├── esm2_zeroshot.py        # ESM-2 masked-marginal scoring
│   ├── esm2_supervised.py      # ESM-2 embeddings + regression head
│   ├── submission_gate.py      # 55-check pre-submission quality gate
│   └── test_*.py               # test suite
├── results/                    # 86 derived result tables + README.md describing each
├── figures/publication/        # figures plus per-figure source-data CSV
├── docs/                       # audit trail and internal notes
├── requirements.lock           # pinned environment
└── ruff.toml                   # lint policy with documented exemptions
```

**What this repository does and does not contain.** It ships the analysis code
and every derived result table, so all reported numbers can be recomputed
directly:

```bash
python src/reproduce_headline.py     # recomputes 29 reported values from results/
```

It does **not** ship the raw inputs — wild-type PDB structures, FoldX mutant
structures, contact-graph `.npz` files or ESM embedding caches — which together
total roughly 21 GB. Those are rebuilt from public sources by the scripts in
`src/` (see "Reproducing the analysis" below for the limitations and required
inputs). Scripts that read only `results/` run immediately; scripts that
rebuild from raw inputs need those inputs first.

## Quick start

```bash
# 1) environment (Python 3.11)
python -m venv .venv
.venv\Scripts\activate            # Windows
# source .venv/bin/activate       # Linux / macOS

# 2) PyTorch first, from the wheel index (otherwise you get the CPU build)
pip install torch==2.6.0 --index-url https://download.pytorch.org/whl/cu124
# ... or .../whl/cpu for a CPU-only install

# 3) remaining dependencies
pip install -r requirements.lock
```

`src/paths.py` resolves the project root from its own location. It sets `DATA`
to `data/` when that directory exists (the full working copy, where the large
inputs live) and otherwise to `results/` (this release). Set `GED_DATA` to point
at a data directory elsewhere, or `GED_ROOT` to relocate the project root
entirely. External tools are expected under `tools/` and can be redirected with
`FOLDX_BIN` and `SCWRL4_BIN`.

## Reproducing the analysis

The commands below are **not** executable end-to-end from this derived-results-only
archive. Obtain the public sequences/labels and structures named in the Methods,
install and configure FoldX separately, and build the `data/structures/`,
`data/mutant_structures_*`, `data/contact_graphs_*`, and `data/esm_emb_cache/`
intermediates before attempting training. No one-command downloader for the full
input bundle or trained model checkpoints are supplied. The released per-sample
predictions support numerical reproduction without retraining (see
`src/reproduce_headline.py`). For example, the raw-input helper requires explicit
arguments: `download_pdbs.py --ids <comma-separated-PDB-IDs> --out_dir
<data/structures>`, and FoldX requires `--label_csv`, `--struct_dir`, `--foldx_bin`,
and `--out_dir`; consult each script's `--help` before use.

The following commands describe analysis **after** those prerequisites are met:

```bash
# contact graphs for one definition (repeat for ca, cb, centroid, allatom)
python src/build_graphs_atomdef.py --atom_def centroid

# capacity ladder: three small encoders use four definitions
python src/run_ladder.py \
    --models gnn_global,gnn_local,gnn_edge \
    --defs ca,cb,centroid,allatom --seeds 42,123,2024 --epochs 20
# the two largest rungs omit all-atom graphs
python src/run_ladder.py \
    --models deep_gine,egnn \
    --defs ca,cb,centroid --seeds 42,123,2024 --epochs 20
# NOTE: inspect --tag / output destinations first; results/ contains the
# published, architecture-matched 508934-parameter EGNN run, and a separate
# 508938-parameter stability-control run. Never aggregate them together.

# sequence baselines
python src/esm2_zeroshot.py  --model esm2_650m
python src/esm2_supervised.py --stage all --model esm2_650m

# pre-submission quality gate
python src/submission_gate.py
```

## Reproducibility

- **Seeding.** `src/seed_utils.py` seeds `random`, `numpy`, `torch`, CUDA and
  cuDNN; it sets `CUBLAS_WORKSPACE_CONFIG` *before* the CUDA context is created
  and enables `torch.use_deterministic_algorithms`. `src/test_determinism.py`
  verifies empirically that two runs with the same seed produce bit-identical
  predictions, and reports the CPU-versus-GPU difference.
- **Quality checks.** `src/submission_gate.py` checks a defined set of manuscript
  assertions, runs selected tests and attempts a DOCX render; it does **not**
  prove the absence of every inconsistency. The local DOCX render requires a
  separately configured LibreOffice-kit executable (`DSH_LO_CLI`, `DSH_NODE`),
  and its absence may cause the gate to fail in other environments. The
  independent `src/reproduce_headline.py` recalculates 29 reported quantities
  from released tables. Neither check replaces full raw-input reconstruction.
- **Linting.** `ruff check` passes on the analysis pipeline. The handful of
  exemptions in `ruff.toml` are declared per file with a written reason rather
  than by disabling whole rule families.

## Data provenance, and a defect we found and fixed

All inputs are public: RCSB PDB structures, MegaScale (Tsuboyama et al., 2023),
ThermoMutDB (Xavier et al., 2021), and the S669 and ssym benchmark sets.
Training data were homology-filtered against both test sets by BLAST before use.

While preparing the structure-aware encoder we audited the training graphs and
found that the **MegaScale half of the training set had been built with a
hydrogen-including side-chain centroid**, whereas every test set used a
hydrogen-excluding centroid. AlphaFold models do contain hydrogens (for example,
1A32 has 566 of them among 1095 atoms), so this was a genuine train/test
representation mismatch affecting half the training data. The affected graphs are
preserved in the full working copy under
`data/contact_graphs_megascale_sc_Hincluded_bug/` for audit (not shipped here —
this release contains result tables only), and every result in
the manuscript was regenerated from the corrected graphs. The finding is written
up in `docs/` alongside the other audit records. We report it here because it is
exactly the kind of defect a reader should be able to check.

## License

Code is released under the MIT License (see `LICENSE`). Processed data derived
from public sources remain subject to the terms of the original databases.

## Contact

Yilin Huang — Yilin.Huang24@student.xjtlu.edu.cn
