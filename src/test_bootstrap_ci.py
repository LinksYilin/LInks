"""
test_bootstrap_ci.py — 验证分层 paired bootstrap 的正确性
============================================================
关键要验证：
  1. 蛋白分层抽样：bootstrap 以蛋白为单元（不是以突变为单元）。
  2. paired Δ：正确计算 A−B 的指标差与置信区间。
  3. 已知"方法 A 恒优于 B"时，CI 应显著（不含 0）。
  4. protein-level 汇总避免大蛋白主导。
"""
import sys

import numpy as np

sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent))
from bootstrap_ci import (
    bootstrap_paired,
    group_by_protein,
    pearson,
    protein_level_summary,
)


def test_group_by_protein():
    print("=== 测试 1: 蛋白分层分组 ===")
    pids = ['p1', 'p1', 'p1', 'p2', 'p2']
    g = group_by_protein(pids)
    print(f"  分组: {g}")
    assert g['p1'] == [0, 1, 2]
    assert g['p2'] == [3, 4]
    print("  ✅ 通过：按蛋白正确分组\n")


def test_bootstrap_significant_when_better():
    print("=== 测试 2: A 恒优于 B 时 Δ 应显著 ===")
    # 构造 4 个蛋白，每个蛋白 5 个突变
    rng = np.random.default_rng(0)
    n_proteins = 4
    n_per = 5
    y_true = []
    y_a = []  # 方法 A（好）：与 y_true 高度相关
    y_b = []  # 方法 B（差）：基本是噪声
    pids = []
    for p in range(n_proteins):
        base = rng.normal(0, 1)
        for i in range(n_per):
            pid = f'p{p}'
            t = base + rng.normal(0, 0.3)
            a = t + rng.normal(0, 0.1)   # A 紧贴真实
            b = rng.normal(0, 1)          # B 是纯噪声
            y_true.append(t)
            y_a.append(a)
            y_b.append(b)
            pids.append(pid)

    bp = bootstrap_paired(y_true, y_a, y_b, pids, metric='pearson', B=500, seed=1)
    print(f"  Δ(A−B) = {bp['delta_point']:.3f}, CI = ({bp['ci_low']:.3f}, {bp['ci_high']:.3f})")
    assert bp['delta_point'] > 0, "A 应优于 B"
    assert bp['ci_low'] > 0, f"A 的 CI 下限应 > 0（显著优于 B），实际 {bp['ci_low']:.3f}"
    print("  ✅ 通过：A 优于 B 时 Δ 显著为正\n")


def test_bootstrap_insignificant_when_equal():
    print("=== 测试 3: A 与 B 相当时 Δ 不显著 ===")
    rng = np.random.default_rng(2)
    y_true = []
    y_a = []
    y_b = []
    pids = []
    for p in range(4):
        base = rng.normal(0, 1)
        for i in range(5):
            t = base + rng.normal(0, 0.3)
            a = t + rng.normal(0, 0.5)
            b = t + rng.normal(0, 0.5)  # 两者同质量噪声
            y_true.append(t); y_a.append(a); y_b.append(b)
            pids.append(f'p{p}')
    bp = bootstrap_paired(y_true, y_a, y_b, pids, metric='pearson', B=500, seed=1)
    print(f"  Δ(A−B) = {bp['delta_point']:.3f}, CI = ({bp['ci_low']:.3f}, {bp['ci_high']:.3f})")
    # A、B 同质量，Δ 的 CI 应包含 0
    assert bp['ci_low'] <= 0 <= bp['ci_high'], "CI 应包含 0（不显著）"
    print("  ✅ 通过：A≈B 时 CI 包含 0（不显著）\n")


def test_protein_level_not_dominated():
    print("=== 测试 4: protein-level 汇总避免大蛋白主导 ===")
    # 蛋白 p_big 有 100 个突变（预测很差），p_small 有 3 个（预测很好）
    y_true = list(np.linspace(0, 1, 103))
    y_pred = list(np.linspace(0, 1, 103))
    # p_big 占前 100，预测差；p_small 占后 3，预测好
    y_pred = [0.0] * 100 + list(np.linspace(0.8, 1.0, 3))
    pids = ['p_big'] * 100 + ['p_small'] * 3

    pl = protein_level_summary(y_true, y_pred, pids, metric='pearson')
    # mutation-level（全部 103 个一起算）会被 p_big 主导，r 应偏低
    ml = pearson(y_true, y_pred)
    print(f"  mutation-level r = {ml:.3f}, protein-level r = {pl:.3f}")
    assert pl > ml, "protein-level 应避免大蛋白主导，r 应高于 mutation-level"
    print("  ✅ 通过：protein-level 汇总避免大蛋白主导\n")


if __name__ == '__main__':
    test_group_by_protein()
    test_bootstrap_significant_when_better()
    test_bootstrap_insignificant_when_equal()
    test_protein_level_not_dominated()
    print("🎉 全部通过：分层 paired bootstrap 逻辑正确")
