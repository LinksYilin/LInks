# WS4：主张精确限定方案

> 目的：解决"第 1 点要加强 vs 第 4 点要收窄"的矛盾——**限定词要放在正确的主张上**，
> 而不是笼统收窄。
> 状态：措辞框架已定，**具体数字待新实验完成后填入**。

---

## 一、五主张表（论文的核心逻辑结构）

| # | 主张 | 证据类型 | 依赖模型能力？ | 依赖 FoldX？ | **正确的限定词** |
|---|---|---|---|---|---|
| **A1** | 不同定义**报告的**接触变化差异巨大 | 结构直接计算 | ❌ 不依赖 | ✅ 依赖 | "under **FoldX BuildModel** side-chain repacking" |
| **A2** | 断边强局部 | 结构直接计算 | ❌ 不依赖 | ✅ 依赖 | 同上 |
| **A3** | 跨引擎（SCWRL4）可复现 A1/A2 的定性结论 | 结构直接计算 | ❌ 不依赖 | ✅ 依赖双引擎 | "across two side-chain modelling engines" |
| **B** | 换图定义**不改变预测** | 受控实验 | ✅ **强依赖** | ❌ **不依赖** | "**across a capacity ladder** from a 5.8 K-parameter mean-pooled GCN to a 509 K-parameter E(3)-equivariant network" |
| **C** | 结构特征在**序列基线之上无增量** | 受控实验 | ✅ 依赖 | ❌ 不依赖 | "with the tested training data and the tested sequence baseline" |
| **D** | 接触计数特征与 ΔΔG 弱相关 | 相关性 | ❌ 不依赖 | ❌ 不依赖 | "on the quality-controlled mutation pairs" |

## 关键推论

### 推论 1：A 是结构事实，不受"弱模型"质疑
"你的 GNN 太弱"**只威胁 B 和 C**，不威胁 A。→ **论文应把 A 单独成节**，与模型能力解耦。

### 推论 2："FoldX repacking" **不能**加在 B/C 上
ΔΔG 预测**只用野生型结构**（`contact_graphs_*`），完全不碰 FoldX 突变体。
把它加在 B/C 上会让审稿人误以为预测依赖 FoldX —— **制造新错误**。

### 推论 3：B 的限定词是"能力阶梯"而非"某个家族"
用"across a capacity ladder from X to Y"限定，**既是诚实的限定，也是证据的陈述**——
比 "within the tested GNN family" 更强，因为阶梯本身就是回答"弱模型"质疑的实验。

---

## 二、标题修改

| | 措辞 | 问题 |
|---|---|---|
| ❌ 现标题 | "Mutation-sensitive residue contact graphs **do not improve** ΔΔG prediction: a controlled benchmark of graph representations" | `do not improve` 是**全称否定**，过度泛化 |
| ✅ 建议 | "Mutation-sensitive residue contact graphs **do not reliably improve** ΔΔG prediction: a controlled benchmark of graph representations across encoder capacities" | 加 `reliably` + 副标题点明能力阶梯 |

> 备选：把 "do not reliably improve" 换成 "provide no reliable increment over sequence baselines"
> ——但这取决于 C 的结果，待定。

---

## 三、需要新增的一节："What we do NOT claim"

**位置**：Discussion 末尾或 Limitations 之后。

**草案**（数字待填）：

> **What we do not claim.**
> First, we do not claim that contact-graph representations are uninformative in general.
> Our comparison varies the representation while holding the encoder, training data, and
> evaluation protocol fixed; a representation that is uninformative under one encoder
> could in principle be informative under another. We therefore scope claim (B) to the
> capacity ladder we evaluated rather than to all possible models.
>
> Second, we do not claim that our graph encoders are state of the art, nor that the
> correlations we report represent an upper bound on what structure-based methods can
> achieve. Our strongest encoder reaches r = 【填入】 on S669, comparable to a
> 【ESM-2 zero-shot/supervised】 sequence baseline at r = 【填入】; published
> structure-based predictors report higher values under their own protocols, which we
> did not reproduce and which would not be directly comparable.
>
> Third, we do not claim a differentiable graph-edit-distance method. The GED-related
> code in this repository implements soft matching and edit attribution for diagnostic
> purposes only; it is not differentiable and does not contribute to any result in this
> paper. We state this explicitly because an earlier version of this project explored
> that direction and we do not want the current results attributed to it.
>
> Fourth, we do not claim that absolute contact-change counts are engine-independent.
> FoldX and SCWRL4 differ in whether they rebuild Cβ coordinates, and we show that this
> difference changes the reported edit rate for the Cβ definition.

---

## 四、摘要结构调整

**当前**：Motivation / Results / Conclusion / Availability

**建议**：
1. `Motivation` → `Background`（BMC 结构式摘要惯例）
2. `Availability` 从摘要移入 Declarations（BMC 要求数据可得性写在声明区）
3. Results 中**明确区分**：
   - A 类结论（结构事实，一句话）
   - B 类结论（能力阶梯内的受控实验）
   - C 类结论（对强序列基线的增量）

---

## 五、Results 结构建议

| 节 | 内容 | 依赖 |
|---|---|---|
| 4.1 | **表示决定"报告什么"**（A1/A2/A3，结构事实） | 不依赖模型 |
| 4.2 | **表示不改变"预测什么"**（B，能力阶梯） | 依赖阶梯 |
| 4.3 | **结构在序列之上的增量**（C，融合实验） | 依赖 ESM |
| 4.4 | FoldX 建模可复现性与跨引擎检验 | — |
| 4.5 | 案例研究 | — |

> **关键**：把 4.1 放在最前，让审稿人先看到**不依赖模型强度的结构事实**，
> 再看依赖模型的 B。

---

## 六、待填入的数字（新实验完成后）

| 位置 | 需要的量 | 来源 |
|---|---|---|
| 标题/摘要 | 能力阶梯范围 | `ladder_summary.csv` |
| 摘要 | ESM-2 zero-shot r | `esm2_zeroshot_*.csv` |
| 摘要 | ESM-2 有监督 r | `esm2_sup_*_results.csv` |
| C 节 | 结构在 ESM 之上的 Δr + CI | `esm2_fusion_results.csv` |
| B 节 | 各骨架内的表示效应 Δr + CI | `ladder_paired_effects.csv` |
| 等价性 | 可排除的 |Δr| 上界 | `ladder_representation_effect.csv` |
| A1 | 编辑率 0% → X% | `edits_corrected_summary.csv` |
| A2 | 断边/保持边距离 | `locality_corrected.csv` |

---

## 七、诚实边界（必须保留的声明）

| 边界 | 位置 |
|---|---|
| ESM-2 预训练数据可能包含测试蛋白（UniRef） | Limitations |
| ssym 仅 15 个蛋白 | Limitations |
| FoldX 与 SCWRL4 均只做刚性骨架侧链重排 | Limitations |
| 我们未复现 ThermoMPNN / Stability Oracle（其图定义焊死在架构内，无法只换表示） | Discussion |
| 绝对编辑数依赖建模引擎 | Limitations |
