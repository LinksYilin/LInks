# Zenodo 手动上传材料（备选方案）

> 若 GitHub 集成不便，可在 https://zenodo.org/uploads/new 手动上传。
> 上传文件用 GitHub 自动打包的 zip（无需本地打包，约 60 MB）：
>
> **https://github.com/LinksYilin/LInks/archive/refs/tags/v1.0.0.zip**

---

## 直接复制粘贴的字段

### Resource type
```
Software
```

### Title
```
Contact-graph definitions shape what a model describes but not what it predicts: analysis code and result tables
```

### Creators
```
Huang, Yilin
```
（其余作者待你补充；单位：Xi'an Jiaotong-Liverpool University）

### Description
```html
<p>Analysis code, result tables, publication figures and the audit trail for a controlled comparison of residue contact-graph definitions in protein &Delta;&Delta;G prediction.</p>

<p><strong>What the study tests.</strong> A residue contact graph is defined by choosing one representative atom per residue and a distance cutoff. The repository evaluates that choice by separating three questions that are usually conflated: (A) what structural change a definition <em>reports</em>; (B) whether the definition changes <em>prediction</em> across a five-rung encoder capacity ladder (5.8 k to 509 k parameters); and (C) whether structure adds any increment over a sequence baseline.</p>

<p><strong>Principal results.</strong> Processing wild-type and modelled mutant structures identically, C&alpha; graphs recorded no contact change in any of the 505 quality-controlled pairs at the reference 8 &Aring; cutoff, whereas C&beta;, side-chain-centroid and all-atom graphs changed in 9.9%, 90.7% and 81.6% of pairs. On the capacity ladder, one of thirteen paired definition comparisons was significant before correction and none survived Holm correction. Adding a side-chain-centroid graph to an ESM-2 650M baseline produced no increment above +0.009 on S669 or +0.070 on ssym.</p>

<p><strong>Contents.</strong> 148 analysis scripts; 85 result tables including per-sample predictions and seed-aware aggregates with Holm-corrected multiplicity; 5 publication figures with source data; 25 documents including the audit trail for two corrected pipeline defects.</p>

<p><strong>Reproducibility.</strong> Inputs (PDB, MegaScale, ThermoMutDB, S669, ssym) are public and are not redistributed. The pipeline resolves its own root directory, so it runs after cloning to any location; set <code>GED_ROOT</code> to keep data elsewhere. Statistical reporting limits, including the definition dependence of the locality statistic and the directional (asymmetric) structure-increment bounds, are documented in the supplementary material.</p>
```

### Keywords
```
protein stability
mutation effect prediction
contact graph
graph neural network
FoldX
side-chain centroid
benchmark
reproducibility
```

### License
```
MIT License
```

### Version
```
1.0.0
```

### Related identifiers
| Relation | Identifier | Scheme |
|---|---|---|
| is supplement to | *(论文 DOI，接收后填)* | DOI |
| is identical to | https://github.com/LinksYilin/LInks | URL |
| is source of | https://github.com/LinksYilin/LInks/tree/v1.0.0 | URL |

---

## 上传后

1. 点 **Publish** → Zenodo 生成两个 DOI：
   - **Version DOI**：`10.5281/zenodo.XXXXXXX`（指向 v1.0.0）
   - **Concept DOI**：`10.5281/zenodo.YYYYYYY`（始终指向最新版）
2. **论文里应引用 Concept DOI**（这样以后更新版本不必改论文）
3. 把 DOI 告诉我，我填进论文的 4 个位置并重跑全部验证

---

## 为什么这一步必须由你完成

DOI 由 Zenodo 铸造，Zenodo 要求用**你本人的账号**登录并点击发布。这一步代表你以作者身份同意将代码归档并公开，属于身份行为，程序无法代为完成。

我能做的（已全部完成）：
- ✅ 代码推送到 GitHub（289 文件）
- ✅ 创建 Release v1.0.0
- ✅ 准备上传包与全部元数据
- ⏳ 拿到 DOI 后填入论文并验证
