# 索引约定 bug（第 7 个问题）

> 发现日期：本次审计。影响面：S669 的 13 个突变（其中 10 个被错误丢弃、3 个写错位置）。

---

## 一、bug 本质

`align_and_filter.py`（第 78-89 行）这样定义 `_pdb_res_idx`：

```python
residues = chain_residues(pdb_path, a['chain'])   # 整条链的全部残基
offset = int(a['offset'])
idx = offset + pos                                 # ← 整链索引
...
new_row['_pdb_res_idx'] = idx
```

即 `_pdb_res_idx = offset + pos_in_label_sequence`，是**整条链**上的索引。

但 `build_contact_graph.build_graph` 会**按 offset/length 切片**：

```python
if offset > 0 or length is not None:
    residues = residues[offset:end]      # 只保留标签序列那一段
```

所以接触图的节点索引是 **0 … len(wt_seq)-1**（切片后），
正确节点索引应为 `_pdb_res_idx − offset`。

**但所有图加载代码都直接用了 `_pdb_res_idx`**：
`train_gnn_baseline.load_s669_samples`、`unified_eval.load_test`、
`benchmark_eval_full.load_benchmark`、`gedmut_train`、`full_eval_pipeline` 等。

## 二、影响量化（实测）

### S669：543 个突变中，13 个 offset>0（全部 offset=120）

| 蛋白 | 突变 | 整链索引 | 正确节点索引 | 后果 |
|---|---|---|---|---|
| 1BA3 | H461D | 456 | 336 | 越界 → **被丢弃** |
| 1BA3 | H489D/H489K/H489M | 484 | 364 | 越界 → **被丢弃** |
| 1F8I | F345A | 344 | 224 | 越界 → **被丢弃** |
| 1NM1 | D187N/D187Y | 175 | 55 | 未越界但**写到残基 P（应为 D）** |
| 1PRE | C159S | 157 | 37 | 未越界但**写到残基 R（应为 C）** |
| 1R2Y | R244E | 242 | 122 | 越界 → **被丢弃** |
| 3D2A | M134E/M137P/S163P | 131/134/160 | 11/14/40 | 越界 → **被丢弃** |
| 5OAQ | Y429H | 198 | 78 | 越界 → **被丢弃** |

- 被错误丢弃：**10 个**（本应有效）
- 写在错误残基上：**3 个**

### ssym：342 个突变全部 offset=0 → **不受影响**

## 三、连带影响

这解释了数据流审计里一直存在的疑问：
> "有接触图 511 条，但 idx 有效只有 501 条"

那"消失的 10 条"就是被这个 bug 错误丢弃的。

## 四、修复

1. `fix_index_bug.py`：给两个基准 CSV 增加 `_node_idx = _pdb_res_idx − offset` 列，并校验残基与声称 wt 一致
   - S669：**512/543 通过**（31 个未通过的全部是 3DV0 一个蛋白）
   - ssym：**342/342 通过**
2. 修补加载器（优先用 `_node_idx`）：`train_gnn_baseline.py`、`unified_eval.py`、`benchmark_eval_full.py`

## 五、修复后样本数

| 基准 | 有图 | 旧索引可用 | **新索引可用** |
|---|---|---|---|
| S669 | 511 | 501 | **511** |
| ssym | 342 | 342 | 342 |

**S669 共同交集从 501 → 511**（找回 10 个突变）。

## 六、3DV0 的处理

3DV0 有 31 个突变索引校验失败，且此前已知：
- WT 42 残基 vs FoldX MT 41 残基（残基数不匹配）
- 处于"无接触图"名单中

**结论：3DV0 应从所有分析中排除**，论文需明确说明。

## 七、教训

这是本次审计发现的**第 3 个实现层 bug**（前两个：氢原子不一致、WT/MT 残基数不匹配）。

三个 bug 的共同模式：
> **"索引/集合约定"在流程的不同环节被隐式假设，且没有任何断言去校验它。**

具体到这里：`align_and_filter` 用整链索引，`build_graph` 用切片索引，
两者之间的转换从未被显式写出，也没有断言检查"落在该索引上的残基是否等于声称的野生型残基"。

**防御性修复建议**：任何以索引取残基的地方，都应立即断言 `seq[idx] == claimed_wt_aa`，
失败即报错而不是静默跳过。若当初加了这一行断言，这个 bug 在第一次运行就会暴露。
