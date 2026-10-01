"""
diag_scwrl4_cb.py — 测量 SCWRL4 是否重建了 Cβ
================================================
若 SCWRL4 重建 Cβ，则 Cβ 坐标会发生小位移，
这解释了"SCWRL4 下 Cβ 编辑率 97.8% 而 FoldX 仅 9.9%"。
"""
import os

import numpy as np
from Bio.PDB import PDBParser

parser = PDBParser(QUIET=True)
ROOT = str(__import__('pathlib').Path(__file__).resolve().parent.parent / 'tools' / 'scwrl4' / 'runs')


def load(path, chain='A'):
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


print('对比 WT 输入 与 SCWRL4 输出 的 CA / CB 位移')
print('（若 CB 位移显著大于 CA，说明 SCWRL4 重建了 CB）')
print()
print(f'{"案例":<16} {"CA 最大位移":>12} {"CB 最大位移":>12} {"CB 平均位移":>12}')
print('-' * 60)

n_checked = 0
for case in sorted(os.listdir(ROOT))[:12]:
    d = os.path.join(ROOT, case)
    if not os.path.isdir(d):
        continue
    muts = [m for m in os.listdir(d) if os.path.isdir(os.path.join(d, m))]
    if not muts:
        continue
    mut = muts[0]
    dd = os.path.join(d, mut)
    wt_p = os.path.join(dd, 'wt_seg.pdb')
    cands = [f for f in os.listdir(dd) if f.startswith('mt_scwrl4_') and f.endswith('.pdb')]
    if not wt_p or not cands:
        continue
    mt_p = os.path.join(dd, min(cands))
    try:
        W = load(wt_p)
        M = load(mt_p)
    except Exception as e:
        print(f'{case} {mut}: 读取失败 {e}')
        continue
    ids = sorted(set(W) & set(M))
    if len(ids) < 3:
        continue
    dca = [np.linalg.norm(W[i]['CA'] - M[i]['CA']) for i in ids]
    cb_ids = [i for i in ids if 'CB' in W[i] and 'CB' in M[i]]
    dcb = [np.linalg.norm(W[i]['CB'] - M[i]['CB']) for i in cb_ids]
    if not dcb:
        continue
    print(f'{case+" "+mut:<16} {max(dca):>12.4f} {max(dcb):>12.4f} {np.mean(dcb):>12.4f}')
    n_checked += 1

print()
print(f'检查了 {n_checked} 个案例')
print()
print('解读：')
print('  CA 位移 ≈ 0 → SCWRL4 保持骨架（符合文档"identical backbone"）')
print('  CB 位移 > 0 → SCWRL4 重建了 CB 坐标 → 接触图随之变化')
