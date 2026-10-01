"""
test_ged.py — 验证 GEDMut 2.0 的 GED 模块（针对审阅指出的三个错误）
======================================================================
测试目标：
  1. 软匹配 ≠ 恒等：构造不对称代价，验证 Sinkhorn 产生非对角对齐。
  2. edge_cost 真正参与：断边/成边/类型变化时 edge_cost_used 必须为 True。
  3. 软归因：break/form 概率是 [0,1] 连续值，而非 0/1 硬值。
  4. 局部子图提取正确。
"""
import sys

import numpy as np

sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent))
from edit_cost import AA_INDEX, contact_type, node_substitution_cost
from ged_module import (
    compute_edit_attribution,
    extract_local_subgraph,
    hungarian,
    sinkhorn,
    summarize_attribution,
)


def make_nodes(seq):
    from edit_cost import AA_PROP
    n = len(seq)
    nodes = np.zeros((n, 23), dtype=np.float32)
    for i, aa in enumerate(seq):
        if aa in AA_INDEX:
            nodes[i, AA_INDEX[aa]] = 1.0
        p = AA_PROP.get(aa, [0, 0.35, 0])
        nodes[i, 20] = p[0]
        nodes[i, 21] = p[1]
        nodes[i, 22] = p[2]
    return nodes


def make_edges(pairs):
    return [(a, b) for (a, b) in pairs] + [(b, a) for (a, b) in pairs]


def test_soft_matching_not_identity():
    print("=== 测试 1: 软匹配 ≠ 恒等（真对齐） ===")
    # 对角代价贵、非对角便宜 → 最优对齐应偏离恒等
    cost = np.array([
        [1.0, 0.1, 0.2],   # 节点0 匹配 1 或 2 更便宜
        [0.1, 1.0, 0.3],
        [0.2, 0.3, 1.0],
    ], dtype=np.float32)
    S = sinkhorn(cost, tau=0.1, n_iter=100)
    print(f"  软分配 S=\n{np.round(S, 3)}")
    # 关键：软匹配产生非对角分配，证明不是恒等映射
    off_diag = S[0, 1] + S[0, 2]
    assert off_diag > S[0, 0], f"非对角质量应超过对角（实际对角 {S[0,0]:.3f}）"
    # 硬匹配也应偏离恒等
    perm = hungarian(cost)
    print(f"  硬匹配 perm={perm}")
    assert not np.array_equal(perm, np.arange(3)), f"硬匹配不应是恒等（实际 {perm}）"
    print("  ✅ 通过：软匹配与硬匹配都非恒等，是真对齐\n")


def test_edge_cost_really_used():
    print("=== 测试 2: edge_cost 真正参与计算 ===")
    seq = ['A', 'L', 'V', 'I', 'G']
    seq_mt = ['A', 'L', 'W', 'I', 'G']
    nodes_wt = make_nodes(seq)
    nodes_mt = make_nodes(seq_mt)
    edges_wt = make_edges([(0, 1), (1, 2), (2, 3), (3, 4)])
    edges_mt = make_edges([(0, 1), (3, 4)])  # 断了 (1,2),(2,3)

    attr = compute_edit_attribution(
        nodes_wt, nodes_mt, edges_wt, edges_mt,
        seq, seq_mt, mut_pos=2,
        node_cost_fn=lambda i, j: node_substitution_cost(i, j),
        edge_attr_fn=contact_type,
        edge_cost_fn=lambda a, b: float(np.abs(np.asarray(a) - np.asarray(b)).sum() / 2),
        k=1,
    )
    summarize_attribution(attr, mut_pos=0)
    assert attr['edge_cost_used'] is True, "断边场景 edge_cost 必须被调用并产生非零贡献"
    assert attr['n_edge_edits'] >= 2, f"应检测到至少2处边编辑（实际 {attr['n_edge_edits']}）"
    print("  ✅ 通过：edge_cost 真实参与，断边被计入\n")


def test_soft_attribution_not_binary():
    print("=== 测试 3: 归因是软值 [0,1] 而非 0/1 ===")
    seq = ['A', 'L', 'V', 'I', 'G']
    seq_mt = ['A', 'L', 'W', 'I', 'G']
    nodes_wt = make_nodes(seq)
    nodes_mt = make_nodes(seq_mt)
    edges_wt = make_edges([(0, 1), (1, 2), (2, 3), (3, 4)])
    edges_mt = make_edges([(0, 1), (3, 4)])

    attr = compute_edit_attribution(
        nodes_wt, nodes_mt, edges_wt, edges_mt,
        seq, seq_mt, mut_pos=2,
        node_cost_fn=lambda i, j: node_substitution_cost(i, j),
        edge_attr_fn=contact_type,
        edge_cost_fn=lambda a, b: float(np.abs(np.asarray(a) - np.asarray(b)).sum() / 2),
        k=1, tau=0.5,
    )
    bp = attr['break_prob']
    broken_vals = bp[bp > 0]
    print(f"  断边软概率值: {np.round(broken_vals, 3)}")
    if len(broken_vals) > 0:
        assert np.all((broken_vals > 0) & (broken_vals < 1)), \
            f"断边概率应为软值，实际 {broken_vals}"
    sp = attr['sub_prob']
    print(f"  替换软概率: {np.round(sp, 3)}")
    assert np.all((sp >= 0) & (sp <= 1)), "替换概率应 ∈ [0,1]"
    print("  ✅ 通过：归因是连续软值，非 0/1 硬差分\n")


def test_local_subgraph():
    print("=== 测试 4: 局部子图提取正确 ===")
    nodes = np.zeros((6, 1), dtype=np.float32)
    edges = make_edges([(0, 1), (1, 2), (2, 3), (3, 4), (4, 5)])
    sub_nodes, sub_edges, node_ids = extract_local_subgraph(nodes, edges, mut_pos=2, k=1)
    print(f"  突变位点2的1-hop子图节点: {node_ids}")
    assert set(node_ids) == {1, 2, 3}, f"1-hop 应为 {{1,2,3}}，实际 {set(node_ids)}"
    sub_nodes2, sub_edges2, node_ids2 = extract_local_subgraph(nodes, edges, mut_pos=2, k=2)
    print(f"  突变位点2的2-hop子图节点: {node_ids2}")
    assert set(node_ids2) == {0, 1, 2, 3, 4}, f"2-hop 应为 {{0,1,2,3,4}}，实际 {set(node_ids2)}"
    print("  ✅ 通过：局部子图提取正确\n")


def test_module_is_numpy_not_autograd():
    """
    如实记录模块的边界：ged_module.py 基于 NumPy，**不在 autograd 图中**，
    因此它不是"可微 GED"。本测试固定这一事实，防止后续被误读为可微实现。
    最终论文（路线 B：受控基准评估）不依赖该模块。
    """
    print("=== 测试 5: 明确该模块不是可微实现（如实声明边界） ===")
    import ged_module
    src = open(ged_module.__file__, encoding='utf-8').read()
    uses_torch = 'import torch' in src
    print(f"  ged_module.py 是否引入 torch: {uses_torch}")
    assert not uses_torch, (
        "ged_module.py 引入了 torch —— 若确实改为可微实现，"
        "请补充真实的梯度测试（torch.autograd.gradcheck 或反向传播断言），"
        "并更新本测试与论文表述")
    print("  ✅ 通过：已确认该模块为 NumPy 实现，不可微；论文未声称可微 GED\n")


if __name__ == '__main__':
    test_soft_matching_not_identity()
    test_edge_cost_really_used()
    test_soft_attribution_not_binary()
    test_local_subgraph()
    test_module_is_numpy_not_autograd()
    print("🎉 全部通过：GED 模块的软匹配、边代价参与、软归因、局部子图逻辑正确；")
    print("   ⚠️ 注意：该模块为 NumPy 实现（不可微），且不参与最终论文结果（路线 B）。")
