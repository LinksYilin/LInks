"""
corrected_edit_analysis_ssym.py — ssym 的接触编辑分析（4 定义 × 5 阈值）
=======================================================================
与 S669 版本同协议：一致氢排除 + WT/MT 残基数质控。
"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from contact_graph_defs import compute_all_defs, edges_at
from paths import DATA, FIGURES, PUB_FIGURES, ROOT, SRC, TOOLS  # 路径集中管理

SRC_PATH = str(SRC)
DATA_PATH = str(DATA)
ROOT_PATH = str(ROOT)
FIG_PATH = str(FIGURES)
PUB_PATH = str(PUB_FIGURES)
TOOLS_PATH = str(TOOLS)

DATA = DATA_PATH
ATOM_DEFS = ['ca', 'cb', 'centroid', 'allatom']
THRESHOLDS = [6.0, 7.0, 8.0, 9.0, 10.0]


def main():
    df = pd.read_csv(os.path.join(DATA, 'benchmarks_ssym_clean.csv'))
    align = pd.read_csv(os.path.join(DATA, 'alignment_map_ssym.csv'))
    amap = {r['pdb_id']: r for _, r in align.iterrows()}
    mt_dir = os.path.join(DATA, 'mutant_structures_ssym')

    wt_cache, rows = {}, []
    n_bad = 0
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
        chain, off = a['chain'], int(a['offset'])

        if pid not in wt_cache:
            wt_cache[pid] = compute_all_defs(wt_path, chain, off, L, max_cut=max(THRESHOLDS))
        dw, sw = wt_cache[pid]
        if dw is None:
            continue
        dm, sm = compute_all_defs(mt_path, chain, off, L, max_cut=max(THRESHOLDS))
        if dm is None or sm is None:
            continue
        if len(sw) != len(sm):
            n_bad += 1
            continue

        wt_aa, mt_aa = mut[0], mut[-1]
        gly = int((wt_aa == 'G') != (mt_aa == 'G'))
        for ad in ATOM_DEFS:
            for th in THRESHOLDS:
                W = edges_at(dw[ad], th)
                M = edges_at(dm[ad], th)
                rows.append({'pdb_id': pid, 'mut_info': mut, 'ddg': row['ddg'], 'gly_switch': gly,
                             'atom_def': ad, 'threshold': th,
                             'wt_edges': len(W), 'mt_edges': len(M),
                             'n_broken': len(W - M), 'n_formed': len(M - W),
                             'n_edit': len(W - M) + len(M - W)})

    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(DATA, 'edits_ssym_corrected.csv'), index=False)
    n_mut = out[['pdb_id', 'mut_info']].drop_duplicates().shape[0]
    print(f'ssym 有效突变对: {n_mut}（剔除残基数不匹配 {n_bad} 个）')

    summ = []
    for ad in ATOM_DEFS:
        for th in THRESHOLDS:
            s = out[(out['atom_def'] == ad) & (out['threshold'] == th)]
            if len(s) == 0:
                continue
            summ.append({'atom_def': ad, 'threshold': th, 'n': len(s),
                         'mean_broken': s['n_broken'].mean(), 'mean_formed': s['n_formed'].mean(),
                         'mean_edit': s['n_edit'].mean(),
                         'pct_any_change': (s['n_edit'] > 0).mean() * 100,
                         'median_edit': s['n_edit'].median()})
    sm = pd.DataFrame(summ)
    sm.to_csv(os.path.join(DATA, 'edits_ssym_corrected_summary.csv'), index=False)

    print('\n=== ssym 8 Å 各定义 ===')
    s8 = sm[sm['threshold'] == 8.0]
    print(f'{"定义":<10} {"n":>5} {"WT边":>7} {"断边":>7} {"成边":>7} {"中位":>6} {"≥1变化":>8}')
    for _, r in s8.iterrows():
        print(f'{r["atom_def"]:<10} {int(r["n"]):>5} {r["mean_broken"]*0:>0}'
              f'{out[(out.atom_def == r["atom_def"]) & (out.threshold == 8.0)]["wt_edges"].mean():>8.1f} '
              f'{r["mean_broken"]:>7.2f} {r["mean_formed"]:>7.2f} {r["median_edit"]:>6.1f} '
              f'{r["pct_any_change"]:>7.1f}%')

    print('\n=== ssym 阈值敏感性（≥1 变化 %）===')
    print(sm.pivot(index='atom_def', columns='threshold', values='pct_any_change').round(1).to_string())
    print('\n已保存 edits_ssym_corrected.csv / edits_ssym_corrected_summary.csv')


if __name__ == '__main__':
    main()
