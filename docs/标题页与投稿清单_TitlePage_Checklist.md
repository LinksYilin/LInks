# Title Page

**Manuscript title**
Contact-graph definitions shape what a model describes but not what it predicts: a controlled comparison across five encoder capacities

**Running title** (≤ 50 characters)
Contact-graph definitions and ΔΔG prediction

**Article type**
Research article

**Target journal**
BMC Bioinformatics

---

## Authors

| Order | Full name | Affiliation | ORCID | Email |
|---|---|---|---|---|
| 1 | *to be supplied* | *to be supplied* | *to be supplied* | *to be supplied* |
| 2 | *to be supplied* | *to be supplied* | *to be supplied* | *to be supplied* |

**Corresponding author**
Yilin Huang — Yilin.Huang24@student.xjtlu.edu.cn

*Replace the placeholder rows above with the final author list, affiliations and ORCIDs before submission.*

---

## Abstract

*See the manuscript. Structured abstract, 345 words (journal limit 350).*

**Keywords**
protein stability; mutation effect prediction; side-chain-centroid contact graphs; graph neural networks; FoldX; protein-cluster bootstrap

---

## Declarations

**Ethics approval and consent to participate**
Not applicable. This study used only publicly available protein structures and experimental stability datasets; no human or animal subjects were involved.

**Consent for publication**
Not applicable.

**Availability of data and materials**
All datasets are publicly available: PDB structures (RCSB Protein Data Bank), MegaScale, ThermoMutDB, and the S669 and ssym benchmarks. Code, model configurations and the scripts that regenerate every table and figure are released under the MIT licence and archived at Zenodo.

- Repository: https://github.com/LinksYilin/LInks
- Archive DOI: DOI to be assigned at acceptance

The repository also contains exploratory graph-edit-distance code that is not differentiable and contributes to no result reported here.

**Competing interests**
The authors declare no competing interests.

**Funding**
*To be supplied. If the work received no specific funding, state: "This research received no specific grant from any funding agency in the public, commercial or not-for-profit sectors."*

**Authors' contributions**
*To be completed once the author list is final, using the CRediT taxonomy (conceptualisation, methodology, software, validation, formal analysis, investigation, data curation, writing — original draft, writing — review and editing, visualisation, supervision).*

**Acknowledgements**
*To be supplied (or omitted if none).*

---

## Figures and tables

| Item | Title | File |
|---|---|---|
| Figure 1 | Study design and evaluation protocol | `figures/publication/figure1_study_design.png` |
| Figure 2 | Side-chain-centroid graphs expose reproducible contact rewiring | `figures/publication/figure2_contact_rewiring.png` |
| Figure 3 | A capacity ladder: the representation effect does not grow with encoder size | `figures/publication/figure3_capacity_ladder.png` |
| Figure 4 | Feature-block accounting for ΔΔG prediction | `figures/publication/figure4_information_budget.png` |
| Figure 5 | Contact-edit diagnostics for the side-chain-centroid definition | `figures/publication/figure5_contact_diagnostics.png` |
| Table 1 | Predictive performance on the common S669 intersection and on ssym | in manuscript |

---

# Submission Checklist

## Before upload

- [ ] Author list, affiliations, ORCIDs and corresponding-author details inserted on this page and in the manuscript
- [ ] GitHub repository created and made public (MIT licence), exact URL captured
- [ ] Zenodo DOI minted from the release tag, both DOI and URL inserted in the Abstract Availability statement, Declarations, and Supplementary S9
- [ ] Funding statement finalised
- [ ] Authors' contributions written (CRediT)
- [ ] Acknowledgements finalised or removed

## Manuscript files

- [ ] Main manuscript (DOCX): 234 paragraphs, 5 figures, 1 table
- [ ] Structured abstract within 350 words (currently 343)
- [ ] Declarations complete (ethics, consent, availability, competing interests, funding, contributions)
- [ ] Supplementary material (10 sections: S1 audit trail, S2 provenance, S3 encoder definitions, S4 increment bounds, S5 ladder and multiplicity, S6 threshold sensitivity, S7 EGNN instability, S8 cross-engine, S9 reproducibility, S10 locality statistic)
- [ ] Cover letter
- [ ] Figures supplied at submission resolution; source data files present

## Verification already performed

- [ ] `submission_gate.py` — 55/55 checks pass (includes stale-value removal, test suite, LibreOffice render)
- [ ] `verify_key_numbers.py` — 75/75 checks pass, each value verified both against its source file and for presence in the manuscript text
- [ ] `verify_supplementary.py` — 29/29 supplementary numbers trace to a named results file
- [ ] `forensic_audit.py` — 0 confirmed internal errors
- [ ] `audit_docx_images.py` — 0 stale or duplicated images
- [ ] `audit_docx_health.py` — 0 structural defects
- [ ] Render check — 16 pages, no missing fonts
- [ ] Reference integrity — 17 references, all cited, none orphaned

## Statements to keep accurate at proof stage

- [ ] The abstract's comparison count (one of thirteen significant before correction, none after Holm correction)
- [ ] The structure-increment bounds are **directional** (+0.010 upper bound on S669; +0.070 on ssym), not a single symmetric bound
- [ ] The Cα zero-change control holds **at the 8 Å reference cutoff only**
- [ ] No published structure-based predictor was reproduced
- [ ] The training data contained a hydrogen-inconsistency defect that was found and corrected
- [ ] The equivariant encoder is reported as unstable across seeds, not as a representation effect
