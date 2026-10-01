"""
label_audit.py — 标签规范化与审计（证据包第 1 层）
====================================================
对 ΔΔG 标签做行级规范化、去重、单位/符号统一、WT 残基核对、剔除原因记录。
按 标签规范化协议.md 实现。

核心纯逻辑函数（可单元测试，无外部依赖）：
  - unit_to_kcalmol(value, unit) -> kcal/mol
  - resolve_duplicates(values) -> (median, n, range)
  - classify_conflict(v1, v2) -> bool  是否矛盾（符号相反且幅度大）
  - normalize_sign(ddg, convention) -> 统一符号约定

用法（完整流程，需 biopython 做 PDB 核对）：
  python label_audit.py --input raw.csv --pdb_dir structures/ \
      --output labels_audit.csv --exclusion_summary exclusion_summary.txt
"""
import argparse
import statistics

# ----------------------------------------------------------------------
# 纯逻辑函数（可单测）
# ----------------------------------------------------------------------
KCAL_PER_KJ = 0.239


def unit_to_kcalmol(value, unit):
    """把 ΔΔG 换算到 kcal/mol。unit: 'kcal/mol' 或 'kj/mol'。"""
    v = float(value)
    u = (unit or '').lower().replace(' ', '')
    if u in ('kcal/mol', 'kcalmol', 'kcal'):
        return v
    if u in ('kj/mol', 'kjmol', 'kj'):
        return v * KCAL_PER_KJ
    raise ValueError(f"未知单位: {unit!r}")


def resolve_duplicates(values):
    """同一突变的多个值：返回 (中位数, 来源数, 极差)。"""
    vals = [float(v) for v in values]
    med = statistics.median(vals)
    rng = max(vals) - min(vals)
    return med, len(vals), rng


def classify_conflict(v1, v2, tol=2.0):
    """判断两个值是否矛盾：符号相反且幅度差 > tol。"""
    v1, v2 = float(v1), float(v2)
    return (v1 * v2 < 0) and (abs(v1 - v2) > tol)


def normalize_sign(ddg, convention='mt_minus_wt'):
    """统一符号约定。convention 表示当前值的含义：
       'mt_minus_wt': ΔG(mt) - ΔG(wt)（标准，正值=失稳）
       'wt_minus_mt': ΔG(wt) - ΔG(mt)（需取负）
    返回标准约定（正值=失稳）下的值。"""
    ddg = float(ddg)
    if convention == 'mt_minus_wt':
        return ddg
    if convention == 'wt_minus_mt':
        return -ddg
    raise ValueError(f"未知符号约定: {convention!r}")


# ----------------------------------------------------------------------
# 主流程（完整版依赖 biopython 做 PDB 核对；此处给骨架）
# ----------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', required=True)
    ap.add_argument('--pdb_dir', default=None)
    ap.add_argument('--output', required=True)
    ap.add_argument('--exclusion_summary', required=True)
    ap.parse_args()  # 仅校验必需参数是否提供（骨架版本暂不用其返回值）

    # 占位：完整实现见 标签规范化协议.md 第五节接口约定
    print("label_audit.py 骨架已就绪。核心纯逻辑函数见本文件，可单元测试。")


if __name__ == '__main__':
    main()
