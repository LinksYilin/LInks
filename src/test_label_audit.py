"""
test_label_audit.py — 验证 label_audit 的纯逻辑函数
====================================================
覆盖：单位换算、去重、冲突判定、符号归一。
"""
import sys

sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent))
from label_audit import (
    classify_conflict,
    normalize_sign,
    resolve_duplicates,
    unit_to_kcalmol,
)


def test_unit_conversion():
    print("=== 测试 1: 单位换算 ===")
    assert abs(unit_to_kcalmol(1.0, 'kcal/mol') - 1.0) < 1e-9
    assert abs(unit_to_kcalmol(4.184, 'kj/mol') - 1.0) < 1e-2  # 4.184 kJ ≈ 1 kcal
    print(f"  4.184 kJ/mol -> {unit_to_kcalmol(4.184, 'kj/mol'):.3f} kcal/mol")
    print("  ✅ 通过：kcal 原样、kJ 换算正确\n")


def test_resolve_duplicates():
    print("=== 测试 2: 重复值取中位数 ===")
    med, n, rng = resolve_duplicates([1.0, 1.2, 3.0])
    assert med == 1.2 and n == 3
    assert abs(rng - 2.0) < 1e-9
    print(f"  中位数={med}, 来源数={n}, 极差={rng}")
    print("  ✅ 通过：中位数/来源数/极差正确\n")


def test_conflict():
    print("=== 测试 3: 矛盾值判定 ===")
    # 符号相反且幅度大 → 矛盾
    assert classify_conflict(2.5, -2.8) is True
    # 符号相反但幅度小 → 不矛盾（可能是 0 附近噪声）
    assert classify_conflict(0.3, -0.4) is False
    # 同号 → 不矛盾
    assert classify_conflict(2.0, 3.0) is False
    print("  ✅ 通过：矛盾判定正确（符号相反且幅度>2）\n")


def test_sign_normalize():
    print("=== 测试 4: 符号归一 ===")
    # mt_minus_wt 约定：正值=失稳，原样返回
    assert normalize_sign(2.0, 'mt_minus_wt') == 2.0
    # wt_minus_mt 约定：需取负
    assert normalize_sign(2.0, 'wt_minus_mt') == -2.0
    print("  ✅ 通过：符号约定归一正确\n")


if __name__ == '__main__':
    test_unit_conversion()
    test_resolve_duplicates()
    test_conflict()
    test_sign_normalize()
    print("🎉 全部通过：label_audit 纯逻辑函数正确")
