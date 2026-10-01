"""
foldx_reproducibility_v2.py — FoldX 建模可复现性重测（一致定义 + 新案例）
=========================================================================
原论文用 1BFM M35W（26 断边）报告 86.7% 一致率，
但该案例在修正定义下只有 0 条真实断边，不能再用。

改为对 3 个"修正后编辑数较多"的案例各跑 3 次 FoldX BuildModel，
用一致定义（排除氢）计算断边/成边的一致率（intersection / union）。
"""
import os
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from contact_graph_defs import extract_residues, pairwise_within, repr_coords
from paths import DATA, FIGURES, PUB_FIGURES, ROOT, SRC, TOOLS  # 路径集中管理

SRC_PATH = str(SRC)
DATA_PATH = str(DATA)
ROOT_PATH = str(ROOT)
FIG_PATH = str(FIGURES)
PUB_PATH = str(PUB_FIGURES)
TOOLS_PATH = str(TOOLS)




DATA = r'D:\GED_mutation\data'
FOLDX = r'D:\GED_mutation\tools\foldx\foldx_1_20270131.exe'
WORK = r'D:\GED_mutation\tools\foldx\repro_v2'
TH, MAXC = 8.0, 10.0

# (pdb_id, mut_info, foldx_mut, chain, offset, length)
CASES = [
    ('1R2Y', 'R244E', 'RA244E;', 'A', 120, 153),
    ('3O39', 'L32P', 'LA32P;', 'A', 0, None),
    ('1XZO', 'W36A', 'WA36A;', 'A', 0, None),
]

ss = pd.read_csv(os.path.join(DATA, 'benchmarks_s669_clean.csv'))
align = pd.read_csv(os.path.join(DATA, 'alignment_map_s669.csv'))


def edges_from(residues):
    coords = repr_coords(residues, 'centroid')
    return {k for k, d in pairwise_within(coords, MAXC).items() if d < TH}


def main():
    os.makedirs(WORK, exist_ok=True)
    results = []
    for pid, mut, fx_mut, chain, off, fixed_len in CASES:
        row = ss[(ss['pdb_id'] == pid) & (ss['mut_info'] == mut)]
        if len(row) == 0:
            print(f'{pid} {mut}: 不在干净集，跳过')
            continue
        row = row.iloc[0]
        L = fixed_len if fixed_len else len(row['wt_seq'])

        d = os.path.join(WORK, f'{pid}_{mut}')
        os.makedirs(d, exist_ok=True)
        pdb_src = os.path.join(DATA, 'structures', f'pdb{pid.lower()}.ent')
        pdb_dst = os.path.join(d, f'pdb{pid.lower()}.pdb')
        if not os.path.exists(pdb_dst):
            shutil.copy(pdb_src, pdb_dst)
        with open(os.path.join(d, 'individual_list.txt'), 'w') as f:
            f.write(fx_mut)
        os.makedirs(os.path.join(d, 'out'), exist_ok=True)

        print(f'--- {pid} {mut} (FoldX: {fx_mut}) ---')
        env = dict(os.environ)
        cmd = [FOLDX, '--command=BuildModel', f'--pdb=pdb{pid.lower()}.pdb',
               '--mutant-file=individual_list.txt', '--output-dir=out', '--numberOfRuns=3']
        p = subprocess.run(cmd, cwd=d, capture_output=True, text=True, env=env)  # noqa: PLW1510 (返回值由调用方检查)
        if p.returncode != 0:
            print(f'  FoldX 失败: {p.stderr[-300:]}')
            continue

        # 读 3 个 run 的接触图
        sets = []
        for k in range(3):
            fp = os.path.join(d, 'out', f'pdb{pid.lower()}_1_{k}.pdb')
            if not os.path.exists(fp):
                print(f'  缺少 {os.path.basename(fp)}')
                continue
            res, seq = extract_residues(fp, chain, off, L)
            if res is None:
                continue
            sets.append(edges_from(res))
        if len(sets) < 2:
            print('  可用 run 不足，跳过')
            continue

        wt_res, wt_seq = extract_residues(pdb_src, chain, off, L)
        W = edges_from(wt_res)
        # 各 run 相对 WT 的断边集合
        broken = [W - s for s in sets]
        formed = [s - W for s in sets]
        inter_b = set.intersection(*broken)
        union_b = set.union(*broken)
        inter_f = set.intersection(*formed)
        union_f = set.union(*formed)
        rec = {
            'pdb_id': pid, 'mut_info': mut, 'n_runs': len(sets),
            'wt_edges': len(W),
            'broken_mean': float(np.mean([len(b) for b in broken])),
            'formed_mean': float(np.mean([len(f) for f in formed])),
            'broken_intersection': len(inter_b), 'broken_union': len(union_b),
            'formed_intersection': len(inter_f), 'formed_union': len(union_f),
            'broken_agreement': len(inter_b) / len(union_b) if union_b else float('nan'),
            'formed_agreement': len(inter_f) / len(union_f) if union_f else float('nan'),
        }
        results.append(rec)
        print(f'  断边 均值 {rec["broken_mean"]:.1f}, 交集 {len(inter_b)}, 并集 {len(union_b)}, '
              f'一致率 {rec["broken_agreement"]*100:.1f}%')
        print(f'  成边 均值 {rec["formed_mean"]:.1f}, 交集 {len(inter_f)}, 并集 {len(union_f)}, '
              f'一致率 {rec["formed_agreement"]*100:.1f}%')

    if results:
        df = pd.DataFrame(results)
        out = os.path.join(DATA, 'foldx_reproducibility_v2.csv')
        df.to_csv(out, index=False)
        print(f'\n已保存 {out}')
        print(f'平均断边一致率: {df["broken_agreement"].mean()*100:.1f}%')
        print(f'平均成边一致率: {df["formed_agreement"].mean()*100:.1f}%')


if __name__ == '__main__':
    main()
