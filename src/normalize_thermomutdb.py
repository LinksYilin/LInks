"""
normalize_thermomutdb.py — ThermoMutDB 标签规范化
==================================================
按 标签规范化协议.md 处理 ThermoMutDB：
  1. 符号翻号：ThermoMutDB 负值=失稳，标准是正值=失稳 → 标准 ddg = -thermomutdb_ddg
  2. 只保留单点突变（mutation_code 不含逗号）
  3. 剔除 ddg 缺失/非数值
  4. 记录 pH、temperature，标记异常（temperature > 150°C 视为可疑）
  5. 产出规范化 CSV + 剔除统计

用法：
  python normalize_thermomutdb.py \
    --input data/raw/thermomutdb/thermomutdb.json \
    --out_clean data/thermomutdb_clean.csv \
    --out_excluded data/thermomutdb_excluded.csv
"""
import argparse
import csv
import json


def is_single_point(mutation_code):
    """单点突变：mutation_code 不含逗号（如 'E49M'，非 'W138Y,W126Y'）。"""
    return mutation_code and ',' not in mutation_code


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', required=True)
    ap.add_argument('--out_clean', required=True)
    ap.add_argument('--out_excluded', required=True)
    args = ap.parse_args()

    with open(args.input, encoding='utf-8') as f:
        data = json.load(f)

    kept, excluded = [], []
    reasons = {}
    for r in data:
        reason = None
        ddg = r.get('ddg')
        mc = r.get('mutation_code', '')
        temp = r.get('temperature')

        if not is_single_point(mc):
            reason = '多点突变或缺失mutation_code'
        elif not isinstance(ddg, (int, float)):
            reason = 'ddg缺失或非数值'
        elif temp is not None and (temp > 400 or temp < 250):
            # 温度单位是开尔文(K)：273K≈0°C, 373K≈100°C。正常 250–400K，超出视为可疑。
            reason = f'温度异常({temp}K)'
        else:
            # 符号翻号 + 保留关键字段
            row = {
                'pdb_id': r.get('PDB_wild', ''),
                'mutation_code': mc,
                'ddg': -float(ddg),  # 翻号：标准正值=失稳
                'ph': r.get('ph'),
                'temperature': temp,
                'method': r.get('method', ''),
                'effect': r.get('effect', ''),
                'source_pmid': r.get('PMID', ''),
                'measure': r.get('measure', ''),
            }
            kept.append(row)
            continue
        if reason:
            reasons[reason] = reasons.get(reason, 0) + 1
            ex = {'mutation_code': mc, 'ddg': ddg, 'reason': reason}
            excluded.append(ex)

    # 写 CSV
    with open(args.out_clean, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=['pdb_id', 'mutation_code', 'ddg', 'ph', 'temperature',
                                          'method', 'effect', 'source_pmid', 'measure'])
        w.writeheader()
        w.writerows(kept)
    with open(args.out_excluded, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=['mutation_code', 'ddg', 'reason'])
        w.writeheader()
        w.writerows(excluded)

    print(f'保留: {len(kept)}, 剔除: {len(excluded)}')
    print('剔除原因分布:')
    for k, v in sorted(reasons.items(), key=lambda x: -x[1]):
        print(f'  {k}: {v}')
    # 验证翻号后符号
    pos = sum(1 for r in kept if r['ddg'] > 0)
    neg = sum(1 for r in kept if r['ddg'] < 0)
    print(f'翻号后: 正值(失稳) {pos} 条, 负值(稳定) {neg} 条')


if __name__ == '__main__':
    main()
