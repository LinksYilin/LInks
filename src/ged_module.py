"""
ged_module.py — 真实的可微局部图编辑距离 + 编辑归因（GEDMut 2.0 重构）
======================================================================
修正旧版三个错误：
  1. 匹配从「恒等映射」改为「真实软匹配」（Sinkhorn-Knopp，可微）。
  2. edge_cost_fn 真正参与计算（进入边编辑代价与总目标）。
  3. 边差分从 0/1 硬值改为「软归因 + 置信度」。

设计（见 GEDMut2.0_设计定稿.md）：
  - 单点突变 → 只在突变位点 k-hop 局部子图内做编辑匹配（可解）。
  - 代价矩阵 C[i,j] = 节点替换代价（+ 可叠加结构失配项）。
  - Sinkhorn 得到软分配 S（可微，梯度可回传）。
  - 软编辑归因：sub_prob / break_prob / form_prob ∈ [0,1]。
  - 可微 GED 目标 = 对齐代价 + λ·边编辑代价，作为下游输入与正则。

本文件为纯 numpy 参考实现，训练时换 torch（Sinkhorn 与代价函数均可微）。
"""
import numpy as np


# ----------------------------------------------------------------------
# 可微软匹配
# ----------------------------------------------------------------------
def sinkhorn(cost, tau=0.1, n_iter=100, tol=1e-6):
    """
    Sinkhorn-Knopp 软分配：把「代价矩阵 cost」转成近似双随机矩阵 S。
    S[i,j] ≈ 节点 wt_i 匹配到 mt_j 的软概率。
    可微：对 cost 的梯度可经迭代传播（numpy 参考版输出 S，torch 版由 autograd 回传）。
    """
    n = cost.shape[0]
    K = np.exp(-cost / tau)
    K = K / (K.sum(axis=1, keepdims=True) + 1e-12)
    for _ in range(n_iter):
        K = K / (K.sum(axis=0, keepdims=True) + 1e-12)
        K = K / (K.sum(axis=1, keepdims=True) + 1e-12)
        if np.max(np.abs(K.sum(axis=1) - 1.0)) < tol:
            break
    return K


def hungarian(cost):
    """硬匹配（仅用于输出离散编辑序列，不参与梯度）。贪心近似，避免 scipy 依赖。

    对每一行（wt 节点），从未占用的列中选代价最小者。逐行处理，保证一一对应。
    """
    n = cost.shape[0]
    used_cols = set()
    perm = np.zeros(n, dtype=np.int64)
    for i in range(n):
        best = -1
        best_v = np.inf
        for j in range(n):
            if j in used_cols:
                continue
            if cost[i, j] < best_v:
                best_v = cost[i, j]
                best = j
        perm[i] = best
        used_cols.add(best)
    return perm


# ----------------------------------------------------------------------
# 局部子图提取
# ----------------------------------------------------------------------
def extract_local_subgraph(nodes, edges, mut_pos, k=1):
    """提取突变位点 mut_pos 的 k-hop 局部子图。返回 (sub_nodes, sub_edges, node_ids)。
    node_ids 把子图索引映射回全局索引。edges: (E,2) 无向边（含双向）。
    """
    n = nodes.shape[0]
    adj = {i: set() for i in range(n)}
    for (a, b) in edges:
        adj[a].add(b)
        adj[b].add(a)
    visited = {mut_pos}
    frontier = {mut_pos}
    for _ in range(k):
        nxt = set()
        for u in frontier:
            for v in adj.get(u, ()):
                if v not in visited:
                    visited.add(v)
                    nxt.add(v)
        frontier = nxt
    node_ids = sorted(visited)
    idx_map = {gid: i for i, gid in enumerate(node_ids)}
    sub_nodes = nodes[node_ids]
    sub_edges = []
    for (a, b) in edges:
        if a in idx_map and b in idx_map:
            sub_edges.append((idx_map[a], idx_map[b]))
    return sub_nodes, sub_edges, node_ids


# ----------------------------------------------------------------------
# 主函数：可微编辑距离 + 软归因
# ----------------------------------------------------------------------
def compute_edit_attribution(nodes_wt, nodes_mt, edges_wt, edges_mt,
                             aa_seq_wt, aa_seq_mt, mut_pos,
                             node_cost_fn, edge_attr_fn, edge_cost_fn,
                             k=1, tau=0.1, n_iter=100, lam_edge=1.0):
    """
    计算 WT→MT 的可微编辑距离与软编辑归因（局部子图内）。

    返回 dict：soft_align, hard_perm, sub_prob, sub_cost, break_prob, form_prob,
               ged_soft, n_local, n_edge_edits, edge_cost_used
    """
    # 1. 局部子图（以 WT 节点为准；单点突变不改节点集，仅改特征）
    sub_nodes_wt, sub_edges_wt, node_ids = extract_local_subgraph(nodes_wt, edges_wt, mut_pos, k)
    m = len(node_ids)
    if m == 0:
        return None

    # MT 的局部子图边：用相同的节点集（node_ids），从 MT 边集中提取
    idx_map = {gid: i for i, gid in enumerate(node_ids)}
    sub_edges_mt = []
    for (a, b) in edges_mt:
        if a in idx_map and b in idx_map:
            sub_edges_mt.append((idx_map[a], idx_map[b]))

    # 2. 代价矩阵（节点替换代价）
    C = np.zeros((m, m), dtype=np.float32)
    for i in range(m):
        gi = node_ids[i]
        for j in range(m):
            gj = node_ids[j]
            C[i, j] = node_cost_fn(aa_seq_wt[gi], aa_seq_mt[gj])

    # 3. 软匹配（真对齐，非恒等）
    S = sinkhorn(C, tau=tau, n_iter=n_iter)
    hard_perm = hungarian(C)

    # 4. 节点软归因
    sub_prob = 1.0 - np.clip(np.diag(S), 0.0, 1.0)
    sub_cost = np.array([C[i, hard_perm[i]] for i in range(m)], dtype=np.float32)

    # 5. 边编辑代价（真实参与）
    def build_edge_map(sub_edges, aa_seq):
        emap = {}
        for (a, b) in sub_edges:
            key = (min(a, b), max(a, b))
            ga, gb = node_ids[a], node_ids[b]
            emap[key] = edge_attr_fn(aa_seq[ga], aa_seq[gb])
        return emap

    wt_emap = build_edge_map(sub_edges_wt, aa_seq_wt)
    mt_emap = build_edge_map(sub_edges_mt, aa_seq_mt)  # 用 MT 自己的边集

    all_keys = set(wt_emap.keys()) | set(mt_emap.keys())
    break_prob = np.zeros((m, m), dtype=np.float32)
    form_prob = np.zeros((m, m), dtype=np.float32)
    edge_cost_sum = 0.0
    n_edge_edits = 0

    for key in all_keys:
        a, b = key
        in_wt = key in wt_emap
        in_mt = key in mt_emap
        p = float(S[a, a] * S[b, b])
        if in_wt and not in_mt:
            break_prob[a, b] = break_prob[b, a] = p
            c = edge_cost_fn(wt_emap[key], np.zeros_like(wt_emap[key]))
            edge_cost_sum += c
            n_edge_edits += 1
        elif (not in_wt) and in_mt:
            form_prob[a, b] = form_prob[b, a] = p
            c = edge_cost_fn(np.zeros_like(mt_emap[key]), mt_emap[key])
            edge_cost_sum += c
            n_edge_edits += 1
        else:
            c = edge_cost_fn(wt_emap[key], mt_emap[key])
            if c > 1e-6:
                edge_cost_sum += c
                n_edge_edits += 1

    ged_align = float(np.sum(C * S))
    ged_soft = ged_align + lam_edge * edge_cost_sum

    return {
        'soft_align': S,
        'hard_perm': hard_perm,
        'sub_prob': sub_prob,
        'sub_cost': sub_cost,
        'break_prob': break_prob,
        'form_prob': form_prob,
        'ged_soft': ged_soft,
        'n_local': m,
        'n_edge_edits': n_edge_edits,
        'edge_cost_used': n_edge_edits > 0,
    }


def summarize_attribution(attr, mut_pos):
    """汇总编辑归因（局部子图内）。"""
    n = attr['sub_prob'].shape[0]
    if mut_pos >= n:
        mut_pos = 0
    print(f"  局部子图节点数={attr['n_local']}, 边编辑数={attr['n_edge_edits']}, "
          f"GED_soft={attr['ged_soft']:.3f}")
    print(f"  突变位点替换概率={attr['sub_prob'][mut_pos]:.3f}, "
          f"替换代价={attr['sub_cost'][mut_pos]:.3f}")
    broken = np.argwhere(attr['break_prob'] > 0.5)
    formed = np.argwhere(attr['form_prob'] > 0.5)
    print(f"  高置信断边数={len(broken)//2}, 高置信成边数={len(formed)//2}")
    return


if __name__ == '__main__':
    print("ged_module.py（GEDMut 2.0）导入成功。请用 test_ged.py 做单元测试。")
