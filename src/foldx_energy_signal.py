"""
foldx_energy_signal.py — FoldX 能量项 vs 实验 ΔΔG 的信号诊断
=============================================================
关键问题：FoldX 的能量分解（而非几何接触变化）与实验 ΔΔG 相关性多高？
这决定"能量编辑"（方案 B）是否值得做。
"""
import os

import numpy as np
import pandas as pd

from paths import DATA, FIGURES, PUB_FIGURES, ROOT, SRC, TOOLS  # 路径集中管理

SRC_PATH = str(SRC)
DATA_PATH = str(DATA)
ROOT_PATH = str(ROOT)
FIG_PATH = str(FIGURES)
PUB_PATH = str(PUB_FIGURES)
TOOLS_PATH = str(TOOLS)




DATA = r'D:\GED_mutation\data'


def parse_dif(path):
    """解析 FoldX Dif 文件，返回 total energy 和各能量项。"""
    lines = open(path).read().split('\n')
    # 找表头和数据行
    header = None
    data_row = None
    for i, l in enumerate(lines):
        if l.startswith('Pdb\t'):
            header = l.split('\t')
            data_row = lines[i + 1].split('\t')
            break
    if header is None or data_row is None:
        return None
    return dict(zip(header, data_row))


def main():
    s669 = pd.read_csv(os.path.join(DATA, 'benchmarks_s669_clean.csv'))
    rows = []
    for _, r in s669.iterrows():
        pid = r['pdb_id']
        mut = r['mut_info']
        dif = os.path.join(DATA, 'mutant_structures_s669', pid, mut, 'work', 'out',
                           f'Dif_pdb{pid.lower()}.fxout')
        if not os.path.exists(dif):
            continue
        d = parse_dif(dif)
        if d is None:
            continue
        try:
            total = float(d['total energy'])
        except (KeyError, ValueError):
            continue
        rows.append({
            'pdb_id': pid, 'mut_info': mut,
            'foldx_ddg': total,
            'ddg_exp': r['ddg'],
        })
    df = pd.DataFrame(rows)
    print(f'有 FoldX 能量的突变: {len(df)}')

    # FoldX ΔΔG vs 实验 ΔΔG 相关性
    r = np.corrcoef(df['foldx_ddg'], df['ddg_exp'])[0, 1]
    print(f'FoldX ΔΔG 与实验 ΔΔG 相关系数: r = {r:.3f}')

    # MAE
    mae = np.mean(np.abs(df['foldx_ddg'] - df['ddg_exp']))
    print(f'MAE = {mae:.2f} kcal/mol')

    # 符号方向检查
    print(f'FoldX ΔΔG 均值 {df["foldx_ddg"].mean():.2f}, 实验 ΔΔG 均值 {df["ddg_exp"].mean():.2f}')


if __name__ == '__main__':
    main()
