"""
corrected_edit_analysis.py — 修正后的接触编辑分析（论文用最终数据）
====================================================================
修正三件事：
  1. 一致定义：所有原子定义排除氢（修正 WT含氢/MT无氢 的不一致）
  2. 质量控制：排除 WT/MT 残基数不匹配的突变对
  3. 甘氨酸分层：单列"涉及 GLY 的原子切换"造成的伪编辑

输出：
  data/edits_corrected.csv       逐突变 × 定义 × 阈值
  data/edits_corrected_summary.csv  汇总（论文用）
"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from contact_graph_defs import compute_all_defs, edges_at, extract_residues
from paths import DATA, FIGURES, PUB_FIGURES, ROOT, SRC, TOOLS  # 路径集中管理

SRC_PATH = str(SRC)
DATA_PATH = str(DATA)
ROOT_PATH = str(ROOT)
FIG_PATH = str(FIGURES)
PUB_PATH = str(PUB_FIGURES)
TOOLS_PATH = str(TOOLS)




DATA = DATA_PATH  # 来自 paths.py，可用 GED_ROOT 环境变量覆盖
ATOM_DEFS = ['ca', 'cb', 'centroid', 'allatom']
THRESHOLDS = [6.0, 7.0, 8.0, 9.0, 10.0]


def main():
    df = pd.read_csv(os.path.join(DATA, 'benchmarks_s669_clean.csv'))
    align = pd.read_csv(os.path.join(DATA, 'alignment_map_s669.csv'))
    amap = {r['pdb_id']: r for _, r in align.iterrows()}
    mt_dir = os.path.join(DATA, 'mutant_structures_s669')

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

        # --- 质控：残基数必须一致 ---
        if len(sw) != len(sm):
            n_bad += 1
            continue

        res, _ = extract_residues(wt_path, chain, off, L)
        wt_has_h = 1 if res and any(x.element == 'H' for r in res for x in r.get_atoms()) else 0
        wt_aa, mt_aa = mut[0], mut[-1]
        gly_switch = int((wt_aa == 'G') != (mt_aa == 'G'))

        for ad in ATOM_DEFS:
            for th in THRESHOLDS:
                W = edges_at(dw[ad], th)
                M = edges_at(dm[ad], th)
                rows.append({
                    'pdb_id': pid, 'mut_info': mut, 'mut_idx': int(row['_node_idx']),  # 修正：节点索引（切片后）
                    'ddg': row['ddg'], 'wt_has_h': wt_has_h, 'gly_switch': gly_switch,
                    'atom_def': ad, 'threshold': th,
                    'wt_edges': len(W), 'mt_edges': len(M),
                    'n_broken': len(W - M), 'n_formed': len(M - W),
                    'n_edit': len(W - M) + len(M - W),
                })

    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(DATA, 'edits_corrected.csv'), index=False)
    n_mut = out[['pdb_id', 'mut_info']].drop_duplicates().shape[0]
    print(f'质控后有效突变对: {n_mut}（剔除残基数不匹配 {n_bad} 个）')

    # 汇总
    summ = []
    for ad in ATOM_DEFS:
        for th in THRESHOLDS:
            s = out[(out['atom_def'] == ad) & (out['threshold'] == th)]
            if len(s) == 0:
                continue
            tot = s['n_edit']
            summ.append({
                'atom_def': ad, 'threshold': th, 'n_mutations': len(s),
                'mean_wt_edges': s['wt_edges'].mean(),
                'mean_broken': s['n_broken'].mean(),
                'mean_formed': s['n_formed'].mean(),
                'mean_edit': tot.mean(),
                'pct_any_change': (tot > 0).mean() * 100,
                'pct_gt5': (tot > 5).mean() * 100,
                'median_edit': tot.median(),
            })
    sm = pd.DataFrame(summ)
    sm.to_csv(os.path.join(DATA, 'edits_corrected_summary.csv'), index=False)

    print('\n=== 8 Å 下各定义（一致氢排除 + 质控后）===')
    s8 = sm[sm['threshold'] == 8.0]
    print(f'{"定义":<10} {"n":>5} {"WT边":>7} {"断边":>7} {"成边":>7} {"中位编辑":>8} {"≥1变化%":>9} {"编辑>5%":>8}')
    for _, r in s8.iterrows():
        print(f'{r["atom_def"]:<10} {int(r["n_mutations"]):>5} {r["mean_wt_edges"]:>7.1f} '
              f'{r["mean_broken"]:>7.2f} {r["mean_formed"]:>7.2f} {r["median_edit"]:>8.1f} '
              f'{r["pct_any_change"]:>8.1f}% {r["pct_gt5"]:>7.1f}%')

    print('\n=== 阈值敏感性（≥1 变化比例 %）===')
    piv = sm.pivot(index='atom_def', columns='threshold', values='pct_any_change')
    print(piv.round(1).to_string())
    print('\n=== 阈值敏感性（平均编辑数）===')
    piv2 = sm.pivot(index='atom_def', columns='threshold', values='mean_edit')
    print(piv2.round(2).to_string())

    print('\n=== 甘氨酸分层（8 Å，质心）===')
    c8 = out[(out['atom_def'] == 'centroid') & (out['threshold'] == 8.0)]
    for lab, sub in [('不含 GLY 切换', c8[c8['gly_switch'] == 0]), ('含 GLY 切换', c8[c8['gly_switch'] == 1])]:
        print(f'  {lab} (n={len(sub)}): 断边 {sub["n_broken"].mean():.2f}, 成边 {sub["n_formed"].mean():.2f}, '
              f'≥1变化 {(sub["n_edit"] > 0).mean()*100:.1f}%')

    print('\n已保存 data/edits_corrected.csv 与 data/edits_corrected_summary.csv')


if __name__ == '__main__':
    main()
