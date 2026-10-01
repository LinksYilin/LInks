"""
diag_cb_both_tools.py — 对比 FoldX 与 SCWRL4 的 CA / CB 位移
==============================================================
公平比较的关键：两个工具各自移动了多少骨架原子和 Cβ 原子。
"""
import os

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser

parser = PDBParser(QUIET=True)
from pathlib import Path as _P
DATA = str(_P(__file__).resolve().parent.parent / 'data')
def load_ids(path, chain='A'):
    st = parser.get_structure('p', path)
    chs = [c for c in st[0].get_chains() if c.id == chain]
    if not chs:
        return {}
    out = {}
    for r in chs[0]:
        if not r.has_id('CA'):
            continue
        d = {'CA': r['CA'].get_coord()}
        if r.has_id('CB'):
            d['CB'] = r['CB'].get_coord()
        out[r.id[1]] = d
    return out


def stats(W, M):
    ids = sorted(set(W) & set(M))
    if len(ids) < 3:
        return None
    dca = [np.linalg.norm(W[i]['CA'] - M[i]['CA']) for i in ids]
    cbi = [i for i in ids if 'CB' in W[i] and 'CB' in M[i]]
    dcb = [np.linalg.norm(W[i]['CB'] - M[i]['CB']) for i in cbi] if cbi else []
    return max(dca), (max(dcb) if dcb else np.nan), (np.mean(dcb) if dcb else np.nan), len(dcb)


df = pd.read_csv(os.path.join(DATA, 'benchmarks_s669_clean.csv'))
align = pd.read_csv(os.path.join(DATA, 'alignment_map_s669.csv'))
amap = {r['pdb_id']: r for _, r in align.iterrows()}

res_foldx, res_scwrl = [], []
scwrl_root = str(__import__('pathlib').Path(__file__).resolve().parent.parent / 'tools' / 'scwrl4' / 'runs')

for _, row in df.head(40).iterrows():
    pid, mut = row['pdb_id'], row['mut_info']
    a = amap.get(pid)
    if a is None:
        continue
    chain = a['chain']
    # --- FoldX ---
    fx = os.path.join(DATA, 'mutant_structures_s669', pid, mut, 'work', 'out', f'pdb{pid.lower()}_1.pdb')
    wt = os.path.join(DATA, 'structures', f'pdb{pid.lower()}.ent')
    if os.path.exists(fx) and os.path.exists(wt):
        try:
            W = load_ids(wt, chain)
            M = load_ids(fx, chain)
            s = stats(W, M)
            if s:
                res_foldx.append((pid, mut) + s)
        except Exception:
            pass
    # --- SCWRL4 ---
    dd = os.path.join(scwrl_root, pid, mut)
    if os.path.isdir(dd):
        cands = [f for f in os.listdir(dd) if f.startswith('mt_scwrl4_') and f.endswith('.pdb')]
        wt_seg = os.path.join(dd, 'wt_seg.pdb')
        if cands and os.path.exists(wt_seg):
            try:
                W = load_ids(wt_seg, chain)
                M = load_ids(os.path.join(dd, min(cands)), chain)
                s = stats(W, M)
                if s:
                    res_scwrl.append((pid, mut) + s)
            except Exception:
                pass

print('=' * 74)
print('CA / CB 位移对比（同一批突变，两个工具）')
print('=' * 74)
for tag, res in [('FoldX BuildModel', res_foldx), ('SCWRL4（局部重排）', res_scwrl)]:
    if not res:
        print(f'\n{tag}: 无数据')
        continue
    d = pd.DataFrame(res, columns=['pid', 'mut', 'ca_max', 'cb_max', 'cb_mean', 'n_cb'])
    print(f'\n{tag}（n={len(d)}）:')
    print(f'  CA 最大位移: 均值 {d["ca_max"].mean():.4f} Å, 最大 {d["ca_max"].max():.4f} Å')
    print(f'  CB 最大位移: 均值 {d["cb_max"].mean():.4f} Å, 最大 {d["cb_max"].max():.4f} Å')
    print(f'  CB 平均位移: 均值 {d["cb_mean"].mean():.4f} Å')

if res_foldx and res_scwrl:
    d1 = pd.DataFrame(res_foldx, columns=['pid', 'mut', 'ca_max', 'cb_max', 'cb_mean', 'n_cb'])
    d2 = pd.DataFrame(res_scwrl, columns=['pid', 'mut', 'ca_max', 'cb_max', 'cb_mean', 'n_cb'])
    m = d1.merge(d2, on=['pid', 'mut'], suffixes=('_fx', '_sc'))
    print('\n' + '=' * 74)
    print(f'配对比较（n={len(m)} 个共同突变）')
    print('=' * 74)
    if len(m):
        print(f'  CA 位移:  FoldX {m["ca_max_fx"].mean():.4f} Å  vs  SCWRL4 {m["ca_max_sc"].mean():.4f} Å')
        print(f'  CB 位移:  FoldX {m["cb_max_fx"].mean():.4f} Å  vs  SCWRL4 {m["cb_max_sc"].mean():.4f} Å')
        print()
        print('  → 若 FoldX 的 CB 位移也显著大于 0，则"FoldX 保留 CB"的说法不成立')
        print('  → 若 FoldX 的 CB 位移接近 0，则"FoldX 保留 CB"成立，两工具差异有明确来源')
