# Cover Letter

**To:** The Editor, *BMC Bioinformatics*

**Subject:** Submission of a research article — "Contact-graph definitions shape what a model describes more than what it predicts: a controlled comparison across five encoder capacities"

---

Dear Editor,

We submit the enclosed manuscript for consideration as a Research Article in *BMC Bioinformatics*.

**What the paper asks.** Structure-based ΔΔG predictors all begin from a residue contact graph, and the definition of that graph — one representative atom per residue plus a distance cutoff — is normally treated as an implementation detail. We show that this assumption deserves separate scrutiny, because the definition controls two things that the literature routinely conflates: what a model can *describe* about a mutation, and what it can *predict*.

**What we found.**

1. **What a definition reports is set by the definition, not by the mutation.** With wild-type and modelled mutant structures processed identically, Cα graphs recorded no contact change in any of 505 quality-controlled pairs at the reference 8 Å cutoff, whereas Cβ, side-chain-centroid and all-atom graphs changed in 9.9%, 90.7% and 81.6% of pairs. A representation can therefore be vacuous by construction. Contact-change distances are localised — a mean of 7.9 Å from the mutated residue versus 15.9 Å for contacts that persisted, excluding contacts incident on the mutated residue itself — but the absolute rate is engine-dependent: under SCWRL4 the Cα result reproduced while the Cβ contact-change rate moved from 9.9% to 97.8%.

2. **What a definition predicts does not follow from what it reports.** We repeated the comparison at five encoder capacities spanning 5.8 k to 509 k parameters. One of the thirteen paired definition comparisons was significant before correction and none survived Holm correction across the thirteen. The mean absolute effect over the definition pairs common to every stable rung did not grow with capacity (exact permutation test, Spearman ρ = −0.80, two-sided exact P = 0.333 over four rungs), which is the pattern that would be absent if limited encoder strength were the reason no effect was detected.

3. **Structure adds no reliable increment over a modern sequence baseline.** Adding a side-chain-centroid graph to an ESM-2 650M baseline produced Δr = −0.013 (95% CI −0.036 to +0.009) on S669 and +0.025 (95% CI −0.007 to +0.070) on ssym. The intervals are asymmetric, and we report them directionally as confidence limits rather than as a pre-specified equivalence test: on S669 the interval excludes positive increments above +0.009, on ssym above +0.070. The graph branch alone reaches 0.460 on the development-independent ssym benchmark, above the physicochemical baseline (0.305) and the narrower sequence head (0.458), so the null increment is not an artefact of a weak graph encoder.

4. **The one structural quantity with univariate signal is redundant.** Broken hydrophobic contacts correlate with ΔΔG at r = 0.135 (95% CI 0.044 to 0.201), the only contact-derived feature to survive Bonferroni correction across the six typed counts, yet they add no incremental information over five physicochemical features (Δr = −0.009, 95% CI −0.047 to +0.037).

**Why this is useful.** The practical consequence is a quantitative account of where representation-layer effort is and is not likely to pay off, together with the audit protocol that produced it (capacity ladder, paired testing with multiplicity control, feature-block accounting, cross-engine check). Practitioners can apply the same protocol to their own representations.

**Transparency.** During the study we identified and corrected two data-handling defects in our own pipeline: half the training graphs had been built from side-chain centroids that included hydrogens while the test graphs excluded them, and one exploratory script fitted its evaluation set in-sample. Both are documented, all affected results were regenerated from corrected inputs, and the audit trail is included in the released repository. We report the equivariant encoder as unstable across seeds rather than presenting its single best run. We also document that the distance-to-broken-contact statistic is sensitive to whether contacts incident on the mutated residue are counted, and report both values.

**Scope.** The conclusions are bounded by what was tested: five encoder families from 5.8 k to 509 k parameters, contact graphs built from FoldX-modelled mutants on fixed backbones, and an encoder family trained in this study rather than any published structure-based predictor. We state these boundaries explicitly in the Discussion rather than leaving them implicit.

**Data and code.** The input datasets are publicly available (PDB, MegaScale, ThermoMutDB, S669, ssym). Code, per-sample predictions, derived tables and figure source data are released under the MIT licence at https://github.com/LinksYilin/LInks and archived at https://doi.org/10.5281/zenodo.23086956. The archive omits the large raw structures, contact graphs, embedding caches and trained checkpoints; rebuilding them requires obtaining the cited public inputs. The work involved no human or animal subjects.

**Author confirmation required before submission:** confirm that the manuscript is original, is not under consideration elsewhere, the competing-interest declaration is accurate, and every named author has approved the final submitted version.

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
