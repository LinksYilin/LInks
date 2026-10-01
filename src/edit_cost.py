"""
edit_cost.py — 编辑代价函数（生物先验显式化）
==============================================
节点替换代价：BLOSUM62 派生（保守替换代价低）。
边编辑代价：接触类型变化。

正式运行时若装了 biopython 可用其 BLOSUM62；这里提供内置简化版，
用 20 种氨基酸的理化属性距离作 fallback，保证无依赖也能跑。
"""
import numpy as np

AA_ORDER = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
AA_INDEX = {aa: i for i, aa in enumerate(AA_ORDER)}

# [疏水性(Kyte-Doolittle), 体积, 电荷] 归一化后的理化属性
AA_PROP = {
    'A': [1.8, 0.354, 0], 'R': [-4.5, 0.694, 1], 'N': [-3.5, 0.456, 0],
    'D': [-3.5, 0.444, -1], 'C': [2.5, 0.434, 0], 'Q': [-3.5, 0.575, 0],
    'E': [-3.5, 0.554, -1], 'G': [-0.4, 0.240, 0], 'H': [-3.2, 0.613, 0],
    'I': [4.5, 0.667, 0], 'L': [3.8, 0.667, 0], 'K': [-3.9, 0.674, 1],
    'M': [1.9, 0.652, 0], 'F': [2.8, 0.760, 0], 'P': [-1.6, 0.451, 0],
    'S': [-0.8, 0.356, 0], 'T': [-0.7, 0.464, 0], 'W': [-0.9, 0.911, 0],
    'Y': [-1.3, 0.774, 0], 'V': [4.2, 0.560, 0],
}


def _physicochem_distance(aa_i, aa_j):
    """理化属性欧氏距离（无 biopython 时的 fallback）。"""
    p_i = np.array(AA_PROP.get(aa_i, [0, 0.35, 0]), dtype=np.float32)
    p_j = np.array(AA_PROP.get(aa_j, [0, 0.35, 0]), dtype=np.float32)
    # 电荷差加倍加权（电荷改变对稳定性影响大）
    d = p_i - p_j
    d[2] *= 2.0
    return float(np.linalg.norm(d))


def node_substitution_cost(aa_i, aa_j):
    """
    节点替换代价（越小越相似），范围约 [0,1]。
    优先用 BLOSUM62；否则用理化属性距离（归一化到 0~1）。
    """
    try:
        from Bio.Align import substitution_matrices
        blosum = substitution_matrices.load("BLOSUM62")
        s = blosum.get((aa_i, aa_j))
        if s is not None:
            # BLOSUM 分数高 = 保守替换 = 代价低；归一化
            smax = 11.0  # BLOSUM62 对角上限约 11
            return max(0.0, min(1.0, 1.0 - s / smax))
    except Exception:
        pass
    # fallback：除以尺度常数（最大理化距离约 6.6），让差异有区分度
    d = _physicochem_distance(aa_i, aa_j)
    return min(1.0, d / 7.0)  # A->V≈0.34、A->W≈0.39、A->R≈0.95


def edge_change_cost(edge_type_i, edge_type_j):
    """边编辑代价：接触类型变化（0=相同类型，1=完全改变）。"""
    if np.array_equal(edge_type_i, edge_type_j):
        return 0.0
    return 1.0


def contact_type(aa_i, aa_j):
    """根据两个残基理化属性判定接触类型。返回 one-hot 特征。"""
    p_i, p_j = AA_PROP.get(aa_i, [0,0,0]), AA_PROP.get(aa_j, [0,0,0])
    hydrophobic = 1.0 if (p_i[0] > 1.5 and p_j[0] > 1.5) else 0.0
    electrostatic = 1.0 if (p_i[2] * p_j[2] < 0) else 0.0
    return np.array([hydrophobic, electrostatic], dtype=np.float32)


if __name__ == '__main__':
    # 简单自测
    print("保守替换 A->V 代价:", round(node_substitution_cost('A', 'V'), 3))
    print("激进替换 A->W 代价:", round(node_substitution_cost('A', 'W'), 3))
    print("相同残基 A->A 代价:", round(node_substitution_cost('A', 'A'), 3))
