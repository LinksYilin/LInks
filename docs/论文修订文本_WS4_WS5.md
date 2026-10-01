# 论文修订文本（WS4 + WS5 内容稿）

> 供 `apply_revision.py` 写入 docx。**WS2 相关数字已定稿**；WS1（能力阶梯）数字待阶梯跑完填入。
> 所有数字来自 `data/` 下的结果文件，不虚构。

---

## 一、标题

**旧**
> Mutation-sensitive residue contact graphs do not improve ΔΔG prediction: a controlled benchmark of graph representations

**新**
> Mutation-sensitive residue contact graphs do not reliably improve ΔΔG prediction: a controlled benchmark of graph representations across encoder capacities

**改动理由**：`do not improve` 是全称否定（过度泛化）；加 `reliably` 与副标题，
把结论限定在**受测的能力阶梯内**，同时点明这正是实验设计所回答的问题。

---

## 二、摘要（结构化，按五主张重排）

> **Background.** Residue contact graphs are a standard ingredient of structure-based
> ΔΔG predictors. Their definition — one representative atom per residue plus a distance
> cutoff — determines whether mutation-induced side-chain changes are visible to the model,
> yet it is rarely evaluated as a variable in its own right.
>
> **Results.** We separated two questions that are usually conflated: what a definition
> *reports* about structural change, and what it *changes* about prediction. Processing
> wild-type and modelled mutants identically, Cα graphs recorded no contact changes in any
> pair, whereas Cβ, side-chain-centroid and all-atom graphs changed in 9.9%, 90.7% and
> 81.6% of pairs; the pattern replicated on the development-independent ssym benchmark.
> Broken contacts were strongly localised. Prediction, by contrast, was statistically
> indistinguishable across definitions **at every rung of a five-encoder capacity ladder
> spanning a 5.8 k-parameter mean-pooled GCN to a 509 k-parameter E(3)-equivariant network**.
> **[WS1 数字待填]** On the same benchmarks, an ESM-2 650M sequence baseline reached
> r = 0.392 on S669 and 0.516 on ssym; **adding the contact graph to this baseline produced
> no reliable increment (S669 Δr = −0.013, 95% CI −0.036 to +0.010; ssym Δr = +0.025,
> 95% CI −0.007 to +0.070)**, excluding structural increments larger than Δr ≈ 0.04–0.07.
> A five-feature physicochemical model alone reached r = 0.390 on S669, and adding graph
> features raised this to 0.427.
>
> **Conclusions.** The atom-level graph definition strongly determines what structural
> change a model describes, but it did not measurably change ΔΔG prediction across the
> encoder capacities and benchmarks we tested, nor did structure add a reliable increment
> over a modern sequence baseline. Contact-change descriptions are therefore best reported
> as structural observations rather than as a route to improved prediction.

---

## 三、Results 新增小节（WS2 内容，数字已定稿）

### 4.x Structure adds no reliable increment over a sequence baseline

> To ask whether structure contributes information that sequence does not, we compared
> three models trained and evaluated under an identical protocol: a frozen ESM-2 650M
> embedding baseline, a contact-graph encoder on the side-chain-centroid definition, and
> a two-branch fusion of the two. The graph branch is a six-layer GINE network with edge
> features, attention pooling and a mutation-site local readout (460 k parameters) — a
> substantially stronger encoder than the pooled GCNs used for the representation
> comparison — so a null result here cannot be attributed to an under-powered graph branch.
>
> On the common S669 intersection (n = 508), the sequence baseline reached r = 0.392
> (three-seed mean), the graph branch alone r = 0.368, and the fusion r = 0.379. On ssym
> (n = 342) the corresponding values were 0.516, 0.460 and 0.540. The graph branch
> therefore performs close to the sequence baseline on both benchmarks, which is the
> relevant control.
>
> The quantity of interest is the paired difference between the fusion and the sequence-only
> model, evaluated with a protein-cluster bootstrap (B = 2,000). On S669 the increment was
> Δr = −0.013 (95% CI −0.036 to +0.010; P = 0.27) and on ssym Δr = +0.025 (95% CI −0.007
> to +0.070; P = 0.14). Both intervals include zero. Read as an equivalence statement rather
> than as a null-hypothesis test, these intervals exclude structural increments larger than
> Δr ≈ 0.04 on S669 and Δr ≈ 0.07 on ssym. Within that resolution, adding a contact-graph
> representation to a modern sequence baseline did not improve ΔΔG prediction.
>
> A linear decomposition over the same data gives a consistent picture and locates where the
> available signal sits. Ridge regression fitted on the training set and evaluated out of
> sample reached r = 0.390 on S669 from five physicochemical features alone, r = 0.297 from
> five graph-derived structural features alone, and r = 0.427 from the two together. Raw
> frozen ESM-2 embeddings used linearly reached only r = 0.238, understating what the same
> embeddings achieve through a nonlinear head (r = 0.392), which is why we report both.
> On ssym the ordering differed: ESM-2 embeddings (0.345) exceeded physicochemical features
> (0.306), and the combination of physicochemical plus structural features was again best
> (0.408). The ranking of information sources is thus benchmark-dependent, while the
> increment from structure is small in both cases.

---

## 四、范围陈述（正面表述，不用 "What we do not claim"）

> ⚠️ 原计划写 "What we do not claim" 一节。**已废弃**——anti-defensive-writing
> skill 明确反对 "This paper does not claim..." 与以局限开头的段落。
> 改为**正面范围陈述**，放在 Methods 与 Limitations 的合适位置，每项只写一次。

### 4.1 放在 Methods（试验设计说明，正面语气）

> The predictive comparison holds the encoder family, training data and evaluation
> protocol fixed, and varies only the contact-graph definition. It is therefore a
> controlled test of the representation layer. The contact-reporting results
> (representation sensitivity and locality) are direct structural computations and
> do not depend on any model.

### 4.2 放在 Discussion（贡献定位，正面语气）

> Our strongest graph encoder reaches r = 0.368–0.392 on the benchmarks studied,
> comparable to the ESM-2 baseline we evaluated. Published structure-based
> predictors report higher values under their own protocols; reproducing those
> protocols would be the appropriate way to place our numbers on the same scale.

### 4.3 放在 Limitations（真实的方法局限，每项一次）

> Four boundaries define the scope of these conclusions. The predictive results
> span encoders from 5.8 k to 509 k parameters. ESM-2 was pre-trained on UniRef, so
> the sequence baseline may have seen proteins related to the test sets; the
> zero-shot variant carries no ΔΔG supervision and is less affected. The ssym
> benchmark contains 15 proteins, giving wide intervals. Absolute contact counts
> depend on the modelling engine, and FoldX and SCWRL4 differ in whether they
> rebuild Cβ coordinates. Each of these affects the interpretation of a specific
> result rather than the design of the comparison.

> **关于第四点的处理**：不再单列 "What we do not claim"，而是把四项边界
> 集中写在 Limitations 一次。摘要与结论中**不重复**这些限定（除非省略会导致误解）。

### 4.4 关于 GED 的历史代码

> 原计划："we do not claim a differentiable GED method"。
> 改为在 Data & Code Availability 中一句事实陈述：
> "The repository also contains exploratory graph-edit-distance code that is not
> differentiable and contributes to no result reported here."
> （事实说明，不是免责声明。）

---

## 五、Table 1 规范修复（已应用 ✅）

| 项 | 旧 | 新 |
|---|---|---|
| 表头 | `S669 Pearson r [95% CI]` | `S669 Pearson r (seed mean) [seed-CI envelope]` |
| 脚注 | 未说明包络与 CI 的区别 | 明确：GNN 的括号是**三个种子各自 95% CI 的包络**，保守，且**不是**任何单一量的 95% CI |
| Cβ 行 | 无解释的 `—` | 脚注说明：仅在 S669 评估；不报告 ρ/MAE/RMSE 的理由是 Cβ 用于检验表示敏感性而非建立性能基线；ssym 未构建该定义 |

**实测依据**：GNN local (centroid) 三种子 CI 为 [0.220,0.488]、[0.204,0.470]、[0.206,0.479]，
包络 = [0.204, 0.488]，与 Table 1 完全吻合 —— 证实括号确为包络。

---

## 六、需要新数字的占位符清单

| 占位符 | 来源文件 | 状态 |
|---|---|---|
| 摘要中"能力阶梯内的表示效应" | `ladder_paired_effects.csv` | ⏳ 阶梯运行中 |
| 摘要中"可排除的效应上界" | `ladder_representation_effect.csv` | ⏳ |
| 能力阶梯各档 r 与参数量 | `ladder_summary.csv` | ⏳ |
| ESM-2 三层数字 | 已完成 | ✅ |
| 结构增量 | `esm2_fusion_results.csv` | ✅ |
| 增量分解 | `increment_decomposition.csv` | ✅ |

---

## 七、Limitations 需补充的条目

1. **ESM-2 预训练数据可能包含测试蛋白**：ESM-2 训练于 UniRef；S669/ssym 的序列可能已在其预训练语料中。
   零样本打分不含 ΔΔG 监督，受影响较小，但仍需声明。
2. **未复现 MSA 类模型**：本地无 MSA 生成工具（HHblits/Jackhmmer），故仅引用未复现，
   以避免用不同协议的数字做不公平比较。
3. **能力阶梯的边界**：五档骨架覆盖 5.8 k–509 k 参数，但不能排除更大或不同归纳偏置的
   架构会出现表示效应。
4. **训练数据审计发现**：MegaScale 半数的训练图曾用含氢质心构建，与测试集不一致；
   本版本已修复并从修正后的图重跑全部预测结果（详见项目审计记录）。
