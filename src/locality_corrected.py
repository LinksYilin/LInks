"""
locality_corrected.py — 修正后的"局部性"分析
=============================================
用一致氢排除 + 质控后的数据，比较：
  - 断裂接触 到突变位点的距离
  - 保持（未变）接触 到突变位点的距离
并做蛋白簇 bootstrap 检验差异。
这直接影响论文 §3.3 的"断边不局部"结论。
"""
import os
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




DATA = DATA_PATH  # 来自 paths.py，可用 GED_ROOT 环境变量覆盖
TH = 8.0
MAXC = 10.0


def pearson(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 3 or a.std() == 0 or b.std() == 0:
        return float('nan')
    return float(np.corrcoef(a, b)[0, 1])


def boot_ci_mean(vals, pid, B=1000, seed=0):
    rng = np.random.default_rng(seed)
    prot = np.array(pid)
    uniq, inv = np.unique(prot, return_inverse=True)
    pidx = [np.where(inv == i)[0] for i in range(len(uniq))]
    st = []
    for _ in range(B):
        idx = np.concatenate([pidx[i] for i in rng.integers(0, len(uniq), len(uniq))])
        if len(idx) >= 3:
            st.append(np.mean(vals[idx]))
    return float(np.percentile(st, 2.5)), float(np.percentile(st, 97.5))


def main():
    df = pd.read_csv(os.path.join(DATA, 'benchmarks_s669_clean.csv'))
    align = pd.read_csv(os.path.join(DATA, 'alignment_map_s669.csv'))
    amap = {r['pdb_id']: r for _, r in align.iterrows()}
    mt_dir = os.path.join(DATA, 'mutant_structures_s669')

    cache = {}
    rows = []
    for _, row in df.iterrows():
        pid, mut = row['pdb_id'], row['mut_info']
        a = amap.get(pid)
        if a is None:
            continue
        base = f'pdb{pid.lower()}'
        mt_path = os.path.join(mt_dir, pid, mut, 'work', 'out', f'{base}_1.pdb')
        wt_path = os.path.join(DATA, 'structures', f'{base}.ent')
        if not (os.path.exists(mt_path) and os.path.exists(wt_path)):
            continue
        L = len(row['wt_seq'])
        ch, off = a['chain'], int(a['offset'])

        if pid not in cache:
            wr, sw = extract_residues(wt_path, ch, off, L)
            cache[pid] = (wr, sw)
        wr, sw = cache[pid]
        if wr is None:
            continue
        mr, sm = extract_residues(mt_path, ch, off, L)
        if mr is None or len(sw) != len(sm):
            continue

        wcoord = repr_coords(wr, 'centroid')
        mcoord = repr_coords(mr, 'centroid')
        W = {k for k, d in pairwise_within(wcoord, MAXC).items() if d < TH}
        M = {k for k, d in pairwise_within(mcoord, MAXC).items() if d < TH}
        broken = W - M
        kept = W & M
        mi = int(row['_node_idx'])  # 修正：用切片后索引
        if mi >= len(wcoord):
            continue
        mc = wcoord[mi]

        def dists(edges, wcoord=wcoord, mc=mc):  # 绑定循环变量（B023）
            out = []
            for (i, j) in edges:
                out.append(min(np.linalg.norm(wcoord[i] - mc), np.linalg.norm(wcoord[j] - mc)))
            return out

        db, dk = dists(broken), dists(kept)
        if not db:
            continue
        rows.append({
            'pdb_id': pid, 'mut_info': mut, 'ddg': row['ddg'],
            'n_broken': len(broken), 'n_kept': len(kept),
            'mean_d_broken': float(np.mean(db)),
            'mean_d_kept': float(np.mean(dk)) if dk else np.nan,
            'frac_broken_lt5': float(np.mean(np.array(db) < 5)),
            'median_d_broken': float(np.median(db)),
        })

    m = pd.DataFrame(rows)
    m.to_csv(os.path.join(DATA, 'locality_corrected.csv'), index=False)
    print(f'有效突变: {len(m)}')
    print()
    print('=== 断裂接触 vs 保持接触 到突变位点距离 ===')
    print(f'  平均断边距离: {m["mean_d_broken"].mean():.2f} Å')
    print(f'  平均保持边距离: {m["mean_d_kept"].mean():.2f} Å')
    lo1, hi1 = boot_ci_mean(m['mean_d_broken'].values, m['pdb_id'].values)
    lo2, hi2 = boot_ci_mean(m['mean_d_kept'].values, m['pdb_id'].values)
    print(f'  断边 95% CI: [{lo1:.2f}, {hi1:.2f}]')
    print(f'  保持边 95% CI: [{lo2:.2f}, {hi2:.2f}]')
    print()
    print(f'  断裂接触中 < 5 Å 的比例: {m["frac_broken_lt5"].mean()*100:.1f}%')
    print(f'  断裂接触距离中位数: {m["median_d_broken"].median():.2f} Å')
    print()
    print('=== 距离与 ΔΔG 的关联 ===')
    print(f'  平均断边距离 vs ΔΔG: r = {pearson(m["mean_d_broken"].values, m["ddg"].values):.3f}')
    print(f'  断边数 vs ΔΔG: r = {pearson(m["n_broken"].values, m["ddg"].values):.3f}')
    print()
    print('=== 断边距离分布 ===')
    bins = [(0, 3), (3, 5), (5, 8), (8, 12), (12, 100)]
    for lo, hi in bins:
        msk = (m['mean_d_broken'] >= lo) & (m['mean_d_broken'] < hi)
        print(f'  [{lo:>2},{hi:>3}) Å: n={msk.sum():>3} ({msk.mean()*100:>5.1f}%)')


if __name__ == '__main__':
    main()
