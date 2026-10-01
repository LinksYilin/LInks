# Mutation-sensitive residue contact graphs: a controlled benchmark of graph representations for ΔΔG prediction

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
│   ├── submission_gate.py      # 51-check pre-submission quality gate
│   └── test_*.py               # test suite
├── data/                       # processed results and per-sample predictions
├── figures/publication/        # figures plus per-figure source-data CSV
├── docs/                       # audit trail and internal notes
├── requirements.lock           # pinned environment
└── ruff.toml                   # lint policy with documented exemptions
```

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

Paths resolve automatically from the repository location. Set the `GED_ROOT`
environment variable only if the data live elsewhere. Model checkpoints for the
external tools are expected at `tools/` by default and can be redirected with
`FOLDX_BIN` and `SCWRL4_BIN`.

## Reproducing the analysis

```bash
# contact graphs for one definition (repeat for ca, cb, centroid, allatom)
python src/build_graphs_atomdef.py --atom_def centroid

# capacity ladder: 5 encoders x 4 definitions x N seeds
python src/run_ladder.py \
    --models gnn_global,gnn_local,gnn_edge,deep_gine,egnn \
    --defs ca,cb,centroid,allatom --seeds 42,123,2024 --epochs 20

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
- **Quality gate.** `src/submission_gate.py` derives every expected number *from
  the result files* (nothing is hard-coded), searches the manuscript for those
  values, asserts that superseded values are absent, runs the test suite, and
  renders the document. It exits non-zero on any failure, so it can be used as a
  commit hook or CI step.
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
preserved under
`data/contact_graphs_megascale_sc_Hincluded_bug/` for audit, and every result in
the manuscript was regenerated from the corrected graphs. The finding is written
up in `docs/` alongside the other audit records. We report it here because it is
exactly the kind of defect a reader should be able to check.

## License

Code is released under the MIT License (see `LICENSE`). Processed data derived
from public sources remain subject to the terms of the original databases.

## Contact

Yilin Huang — Yilin.Huang24@student.xjtlu.edu.cn
