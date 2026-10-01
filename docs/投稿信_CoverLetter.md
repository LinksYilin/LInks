# Cover Letter

**To:** The Editor, *BMC Bioinformatics*

**Subject:** Submission of a research article — "What does structure contribute to ΔΔG prediction? A controlled information-budget analysis across residue contact-graph representations"

---

Dear Editor,

We submit the enclosed manuscript for consideration as a Research Article in *BMC Bioinformatics*.

**What the paper asks.** Structure-based ΔΔG predictors all begin from a residue contact graph, and the definition of that graph — one representative atom per residue plus a distance cutoff — is normally treated as an implementation detail. We show that this assumption deserves separate scrutiny, because the definition controls two things that the literature routinely conflates: what a model can *describe* about a mutation, and what it can *predict*.

**What we found.**

1. **The description is highly definition-dependent.** With wild-type and modelled mutant structures processed identically, Cα graphs recorded no contact change in any of 505 quality-controlled pairs, whereas Cβ, side-chain-centroid and all-atom graphs changed in 9.9%, 90.7% and 81.6% of pairs. Broken contacts were strongly localised (mean 4.03 Å from the mutated residue, versus 15.70 Å for contacts that persisted). This pattern replicated on the development-independent ssym benchmark and across two side-chain modelling engines (FoldX and SCWRL4).

2. **The prediction is not.** We repeated the four-definition comparison at five encoder capacities spanning 5.8 k to 509 k parameters. Twelve of the thirteen paired definition comparisons were non-significant, and the mean absolute effect did not grow with capacity — the pattern that would be absent if limited encoder strength were the reason no effect was detected.

3. **Structure adds no reliable increment over a modern sequence baseline.** Adding a side-chain-centroid graph to an ESM-2 650M baseline produced Δr = −0.013 (95% CI −0.036 to +0.010) on S669 and +0.025 (95% CI −0.007 to +0.070) on ssym. Read as an equivalence statement, these intervals exclude structural increments above |Δr| ≈ 0.04. The graph branch alone reaches 0.478 on the independent ssym benchmark, above the physicochemical baseline (0.305) and the narrower sequence head (0.458), so the null increment is not an artefact of a weak graph encoder.

4. **The one structural quantity with univariate signal is redundant.** Broken hydrophobic contacts correlate with ΔΔG at r = 0.135 (95% CI 0.044 to 0.201), the only contact-derived feature to exclude zero, yet they add no incremental information over five physicochemical features (Δr = −0.009, 95% CI −0.047 to +0.037).

**Why this is useful.** The practical consequence is a quantitative guide for where representation-layer effort is likely to pay off, together with the audit protocol that produced it (capacity ladder, paired equivalence testing, information decomposition, cross-engine check). Practitioners can apply the same protocol to their own representations.

**Transparency.** During the study we identified and corrected two data-handling defects in our own pipeline: half the training graphs had been built from side-chain centroids that included hydrogens while the test graphs excluded them, and one exploratory script fitted its evaluation set in-sample. Both are documented, all affected results were regenerated from corrected inputs, and the audit trail is included in the released repository. We report the equivariant encoder as unstable across seeds rather than presenting its single best run.

**Data and code.** All data are public (PDB, MegaScale, ThermoMutDB, S669, ssym). Code, model configurations and scripts that regenerate every table and figure are released under the MIT licence. We have no competing interests, and the work involved no human or animal subjects.

We confirm this manuscript is original, is not under consideration elsewhere, and has been approved by all authors.

Yours sincerely,

**Yilin Huang** (corresponding author)
*on behalf of all authors*

Correspondence: Yilin.Huang24@student.xjtlu.edu.cn

---

## Suggested reviewers

*(to be completed by the authors)*

| Name | Affiliation | Expertise |
|---|---|---|
| | | structure-based ΔΔG prediction |
| | | protein language models |
| | | graph representation learning for proteins |

## Opposed reviewers

None.
