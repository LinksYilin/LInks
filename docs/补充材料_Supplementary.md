# Supplementary Material

**Manuscript:** Contact-graph definitions shape what a model describes but not what it predicts: a controlled comparison across five encoder capacities

---

## S1. Data pipeline audit trail

During this study we identified and corrected two defects in our own processing pipeline. We report both here rather than only in the released code, because each affected a subset of the numerical results.

### S1.1 Hydrogen inconsistency in training-graph construction

**Defect.** The side-chain-centroid definition requires one representative point per residue. Half of the training set (the 232 MegaScale proteins that survive sampling, 239 graph files on disk) had been built from candidates that **included hydrogen atoms** in the centroid average, whereas the ThermoMutDB training graphs and all test graphs (S669, ssym) **excluded** hydrogens. AlphaFold models contain substantial hydrogen content (for 1A32, 566 of 1,095 atoms), so the two centroid definitions are numerically different quantities.

**Consequence.** Training and test inputs differed in definition for half the training set, and the resulting edges differed (1A32: 206 contacts with hydrogens versus 192 without).

**Correction.** All MegaScale training graphs were rebuilt with hydrogen exclusion, and the superseded directory is retained as `contact_graphs_megascale_sc_Hincluded_bug`. Every predictive number in this manuscript was regenerated from the corrected graphs.

**Before/after for the affected quantities.** The hydrogen rebuild changes graph-derived predictions only, because the physicochemical ridge baseline never reads the graphs.

| Quantity | Before the hydrogen rebuild | After the rebuild |
|---|---|---|
| Edge-aware GINE, range across definitions (S669) | 0.289–0.354 (three seeds, H-included training graphs) | 0.340–0.373 (means of three per-seed correlations) |
| Mutation-site GCN, side-chain-centroid definition (S669) | 0.365–0.382 (three seeds) | 0.378 (mean of three per-seed correlations) |
| EGNN, side-chain-centroid definition (S669) | 0.393 / −0.081 (two seeds) | 0.393 / −0.081 / 0.079 (three seeds, same architecture) |
| Physicochemical ridge baseline (S669) | 0.392 | 0.390 |
| Physicochemical ridge baseline (ssym) | 0.316 | 0.305 |

The values quoted before the rebuild are those in the superseded table, in which the GINE range spanned both definitions and seeds; the after column reports the same quantity recomputed on the corrected graphs. The S669 ridge baseline is unchanged to two decimals (0.392 to 0.390). The ssym entry moves from 0.316 to 0.305, but those two values come from different evaluation pipelines rather than from the rebuild: the hydrogen correction cannot affect a ridge model that never reads the graphs. The hydrogen audit counts are recorded in `data/hydrogen_bias_audit.csv`, and every per-stage exclusion is logged in `data/data_flow_skip_log.csv`.

### S1.2 In-sample fitting in an exploratory script

**Defect.** An early version of the information-decomposition analysis fitted its ridge models and evaluated them on the same mutations, which inflated the apparent contribution of high-dimensional feature blocks (the frozen ESM-2 block reached r = 0.912 under that error).

**Correction.** The analysis now fits on the training set and evaluates out of sample. The affected manuscript quantities are the linear-readout values in Section 4.4: the frozen ESM-2 block is reported at r = 0.238 (S669) and 0.345 (ssym) rather than 0.912, and Figure 4a plots the corrected values. No result outside that section used the in-sample fit.

---

## S2. Sample provenance

The S669 evaluation sets were built by successive filters from the 543-mutation benchmark. The filters are not strictly nested: the 505-pair contact-change set was constructed from the 538 modelled mutant structures independently of the 511- and 508-mutation predictive sets, so it contains four pairs that the predictive sets lack and omits seven that they contain.

| Stage | n | Excluded | Reason |
|---|---|---|---|
| S669 benchmark | 543 | — | — |
| Residue-index verification passed | 512 | 31 | inconsistent structure (single protein, 3DV0) |
| Contact graph available | 511 | 1 | graph could not be built (1G3P) |
| ESM-2 embedding cached | 508 | 3 | sequence not representable (1O6X ×2, 2HBB ×1) |
| FoldX output with matching residue count | 505 | 33 removed, 2 re-admitted | incomplete or mismatched FoldX output, measured against the 538 modelled pairs |

Consequently: **505** mutations for contact-change statistics, **511** for predictive comparisons without sequence features, **508** for sequence and fusion comparisons. The 505-pair set was built from the 538 modelled pairs (538 − 33 + 2 = 507 candidate pairs, of which 505 carry a complete graph). The ssym benchmark contributed **342** mutations with no losses at any stage.

Training data came from MegaScale plus ThermoMutDB after BLAST leakage filtering, sampled 1:1 to 7,905 examples drawn from **420 distinct proteins** (232 sampled from MegaScale and 191 from ThermoMutDB, with three proteins present in both sources). The leakage criterion was applied as implemented: sequence identity above 25% together with an E-value below 0.01 and query coverage above 0.5.

---

## S3. Encoder definitions

All encoders receive 24 node features (20-residue one-hot of the mutant identity at the mutated position, three normalised physicochemical properties, and a binary mutation-position flag) and 2 edge features (hydrophobic indicator, electrostatic indicator).

| Name | Layers | Hidden width | Readout | Parameters |
|---|---|---|---|---|
| GCN, global pooling | 2 × GCNConv, ReLU | 64 | global mean pool → linear | 5,825 |
| GCN, mutation-site pooling | 2 × GCNConv, ReLU | 64 | concat(global mean pool, mutated-residue representation) → linear | 5,889 |
| GINE, edge-aware | 2 × GINEConv (edge features mixed into messages), ReLU | 64 | attention pool over nodes → concat(global, k-hop local, k = 1) → linear | 50,050 |
| Deep GINE, attention pooling | 6 × GINEConv with LayerNorm residuals | 128 | jumping knowledge concat → attention pool → concat(global, k-hop local, k = 2) → 2-layer head | 460,034 |
| EGNN, E(3)-equivariant | 4 × EGNN layer (edge MLP on pair distances, coordinate update) | 128 | attention pool → concat(global, k-hop local, k = 2) → 2-layer head | 508,934 |

Training settings were identical for every encoder: 20 epochs, Adam with learning rate 1 × 10⁻³, batch size 64, mean-squared-error loss, and the same data, splits and graph definitions (Cα, Cβ, side-chain centroid with hydrogens excluded, all-atom minimum distance). No hyperparameter search was performed; the ladder varies capacity by architecture, not by tuning, so that any performance difference is attributable to the encoder family and not to per-rung optimisation. The equivariant encoder additionally receives per-residue coordinates. Its coordinate-update step is **unscaled** in the reported ladder, which yields 508,934 parameters; the stability-control variant described in S7 multiplies the coordinate update by a learnable scalar initialised at zero, which adds four parameters (508,938) and is not used in any result reported in the main text.

---

## S4. Structure increment and its directional bounds

| Comparison | Δr | 95% bootstrap interval | Excluded below | Excluded above |
|---|---|---|---|---|
| S669, fusion − ESM-2 only | −0.013 | [−0.036, +0.009] | −0.036 | +0.009 |
| ssym, fusion − ESM-2 only | +0.025 | [−0.007, +0.070] | −0.007 | +0.070 |

Intervals are 95% percentile intervals from a protein-cluster bootstrap (B = 2,000) in which whole proteins are resampled with replacement. The intervals are **asymmetric**, so a single symmetric bound cannot be quoted for either benchmark: on S669 the interval excludes positive increments above +0.009 and negative increments below −0.036, and on ssym the negative direction is constrained only above −0.007. No equivalence margin was pre-specified, and no formal equivalence test was performed; the intervals are reported as confidence limits rather than as a demonstrated equivalence claim.

---

## S5. Capacity ladder: paired comparisons and multiplicity

Thirteen paired definition comparisons (each non-reference definition against the side-chain centroid) were computed with a joint seed and protein bootstrap (B = 2,000, resampling both proteins and the seed ensemble).

| Encoder | Parameters | Comparison | Δr | Seed-aware 95% CI | Uncorrected P | Holm P |
|---|---|---|---|---|---|---|
| GCN, global pooling | 5,825 | all-atom vs centroid | −0.101 | [−0.186, −0.002] | 0.045 | 0.585 |
| GCN, global pooling | 5,825 | Cα vs centroid | −0.079 | [−0.164, +0.034] | 0.153 | 1.000 |
| GCN, global pooling | 5,825 | Cβ vs centroid | −0.018 | [−0.087, +0.060] | 0.622 | 1.000 |
| GCN, mutation-site pooling | 5,889 | Cα vs centroid | −0.064 | [−0.128, +0.016] | 0.107 | 1.000 |
| GCN, mutation-site pooling | 5,889 | all-atom vs centroid | −0.049 | [−0.107, +0.031] | 0.176 | 1.000 |
| GCN, mutation-site pooling | 5,889 | Cβ vs centroid | −0.031 | [−0.085, +0.029] | 0.283 | 1.000 |
| GINE, edge-aware | 50,050 | Cα vs centroid | −0.024 | [−0.105, +0.070] | 0.686 | 1.000 |
| GINE, edge-aware | 50,050 | all-atom vs centroid | −0.014 | [−0.109, +0.094] | 0.848 | 1.000 |
| GINE, edge-aware | 50,050 | Cβ vs centroid | +0.019 | [−0.090, +0.104] | 0.630 | 1.000 |
| Deep GINE | 460,034 | Cβ vs centroid | −0.042 | [−0.139, +0.070] | 0.436 | 1.000 |
| Deep GINE | 460,034 | Cα vs centroid | −0.044 | [−0.163, +0.067] | 0.380 | 1.000 |
| EGNN | 508,934 | Cα vs centroid | −0.052 | [−0.151, +0.086] | 0.439 | 1.000 |
| EGNN | 508,934 | Cβ vs centroid | −0.065 | [−0.549, +0.520] | 0.753 | 1.000 |

**One of the thirteen comparisons was significant before correction and none survived Holm correction** (smallest corrected P = 0.585). The uncorrected case is the all-atom definition at the smallest encoder, where the all-atom graph performs worse than the side-chain centroid.

**Capacity trend.** Restricting to the definition pairs available at every stable rung (Cα and Cβ against the centroid) gives mean absolute effects of 0.048 (5.8 k), 0.048 (5.9 k), 0.021 (50 k) and 0.043 (460 k). An exact permutation test over all 24 orderings of these four values gives Spearman ρ = −0.80 with a two-sided exact P = 0.333. With four rungs the test is underpowered, so this is a failure to detect a growth trend rather than evidence that none exists; the point estimate lies in the direction opposite to the one that limited encoder capacity would produce.

---

## S6. Threshold sensitivity

Predictions under the side-chain-centroid definition as a function of the contact cutoff (GCN with mutation-site pooling, seeds 42 and 123, predictions averaged over the two seeds before correlating; the ladder in S5 uses three seeds for the same encoder). The "mean edges" column is the mean number of contacts per mutation in the evaluated test graphs, not a graph size.

| Cutoff | S669 r | ssym r | S669 mean edges per mutation |
|---|---|---|---|
| 6 Å | 0.410 | 0.381 | 140 |
| 7 Å | 0.386 | 0.371 | 237 |
| 8 Å | 0.382 | 0.410 | 338 |
| 9 Å | 0.353 | 0.359 | 465 |
| 10 Å | 0.353 | 0.348 | 615 |

The 8 Å row reproduces the mutation-site GCN value used in the main ladder (0.378 mean of three per-seed correlations versus 0.382 two-seed mean). Shifts across cutoffs are small relative to the bootstrap uncertainty but are larger than the structure-increment interval, and the main text reports the increment as a bounded quantity rather than as a demonstrated equivalence for this reason.

---

## S7. EGNN seed instability

The E(3)-equivariant encoder was unstable across seeds under every definition tested. On the side-chain-centroid definition, using the same architecture (508,934 parameters) for all three seeds:

| Variant | Seed 42 | Seed 123 | Seed 2024 | Mean | Sample SD (n − 1) |
|---|---|---|---|---|---|
| Original | 0.393 | −0.081 | 0.079 | 0.131 | 0.241 |
| Coordinate update rescaled from zero at initialisation | −0.182 | 0.346 | 0.076 | 0.080 | 0.264 |

The same pattern appears on the other definitions: Cα gives 0.343, −0.053 and −0.016, and Cβ gives −0.072, 0.358 and −0.014. All three seeds share the same 508,934-parameter architecture. Rescaling the coordinate update did not stabilise the encoder. We therefore report the highest capacity rung using the 460 k-parameter attention-pooled GINE network, which was stable across seeds, and report the equivariant results as an optimisation instability rather than as a representation effect.

---

## S8. Cross-engine contact-change reproducibility

Repeated FoldX BuildModel runs (three per mutation) on three case-study mutations gave broken-contact agreement of 92.3%, 96.2% and 100.0% (mean 96.2%) and formed-contact agreement of 40.0%, 57.1% and 100.0% (mean 65.7%), measured as intersection over union of the changed-contact sets. The three-mutation basis is illustrative rather than an estimate of population reproducibility, and the formed-contact figure rests on the two non-saturated cases.

A matched-scope SCWRL4 comparison covered 507 of 543 S669 pairs. Cα conclusions were preserved (maximum displacement 0 Å) and the side-chain-centroid definition changed in 100% of pairs, whereas the Cβ contact-change rate moved from 9.9% (FoldX) to 97.8% (SCWRL4) because SCWRL4 rebuilds Cβ coordinates and FoldX BuildModel retains them. The absolute contact-change rate is therefore engine-dependent, and only the Cα comparison is engine-invariant.

---

## S9. Reproducibility

| Resource | Location |
|---|---|
| Code, model configurations, analysis scripts | https://github.com/LinksYilin/LInks (MIT licence), archived at Zenodo (DOI and URL to be inserted at submission) |
| Per-sample predictions for every comparison | `data/*_predictions.csv` |
| Ladder per-seed results and audited aggregates | `data/ladder_results.csv` (seeds 42, 123), `data/ladder_egnn_legacy_s2024_results.csv` (the architecture-matched 508,934-parameter EGNN seed 2024 run used in the audited aggregates), `data/ladder_seed_summary_audited.csv`. `data/ladder_s2024_results.csv` holds the rescaled 508,938-parameter variant and is excluded before aggregation. |
| Paired effects with multiplicity correction | `data/ladder_paired_effects_audited.csv`, `data/ladder_paired_effects_audited_holm.csv` |
| Audit records for S1 | `data/hydrogen_bias_audit.csv`, `data/data_flow_skip_log.csv` |
| Figure source data | `figures/publication/*_source_data.csv` |

Seed handling is stated per analysis. The capacity ladder used seeds 42, 123 and 2024 for the four stable encoders and the same three seeds for an architecture-matched equivariant encoder; the sequence-baseline and fusion models used seeds 42, 123 and 2024. Two aggregations are used and are distinguished throughout. The correlation reported for each encoder and graph definition is the unweighted mean of the per-seed correlations (`seed_mean_r`). The paired definition effects reported in S5 are computed from the joint seed and protein bootstrap, in which the seed-averaged prediction enters the resampling. Averaging per-seed correlations and correlating seed-averaged predictions give different values; for the four stable encoders the difference is at most 0.033, and for the equivariant network it reaches 0.216 because its per-seed correlations are bimodal. `set_seed` fixes Python, NumPy, PyTorch and CUDA state and enables deterministic algorithms. Segmented softmax pooling uses a deterministic clamp/exp/`index_add_` formulation rather than `index_reduce_('amax')`.

---

## S10. Locality statistic: definition and sensitivity

The distance from a mutation to the contacts it loses admits more than one definition, and the reported value depends on the choice. For each mutation we measured, for each contact present in the wild-type graph and absent from the modelled mutant graph, the distance from the mutated residue to the nearest endpoint of that contact.

A contact incident on the mutated residue itself is at distance zero under this rule, and such self-contacts account for **382 of the 1,068 broken contacts** and appear in **253 of the 358** mutation pairs with at least one broken contact. Averaging over mutations without removing them gives 4.03 Å; this value is a property of the self-contact convention rather than of how far the structural change propagates.

| Definition | Broken contacts | Unchanged contacts |
|---|---|---|
| All broken contacts, per-mutation mean | 4.03 Å (95% CI 3.11–5.42) | 15.70 Å (13.95–18.18) |
| **Self-contacts removed, per-mutation mean (n = 234)** | **7.93 Å (95% CI 6.63–10.11)** | **15.93 Å (13.93–18.86)** |
| Self-contacts removed, per-contact mean | 11.93 Å | — |
| Self-contacts removed, per-mutation median | 5.88 Å | — |

The main text reports the self-contacts-removed per-mutation mean. The qualitative conclusion is unchanged under every definition: broken contacts lie markedly closer to the mutated residue than contacts that persist, and 30.1% of non-self broken contacts lie within 5 Å. The distance is not a sufficient statistic for the effect of a contact on ΔΔG, which is why the association analysis in Section 4.5 uses counts rather than distances.
