"""
diagnose_hydrogen_bias.py — 系统性核查"WT 含氢 / MT 无氢"导致的假编辑
======================================================================
对 S669 全部有突变体结构的突变，比较两种定义下的编辑统计：
  (A) 既有定义（质心含氢；WT 含氢 → MT 无氢）→ 复现论文数字
  (B) 一致定义（质心排除氢；WT/MT 都一样）
并统计 WT/MT 结构中氢原子的出现率。
"""
import os

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from scipy.spatial import cKDTree

from pathlib import Path as _P
DATA = str(_P(__file__).resolve().parent.parent / 'data')
BACKBONE = {'N', 'CA', 'C', 'O', 'OXT'}
parser = PDBParser(QUIET=True)


def load_chain(path, chain_id=None, offset=0, length=None):
    st = parser.get_structure('p', path)
    chs = list(st[0].get_chains())
    if chain_id is not None:
        chs = [c for c in chs if c.id == chain_id]
    if not chs:
        return None
    res = [r for r in chs[0] if r.has_id('CA')]
    if offset > 0 or length is not None:
        end = (offset + length) if length is not None else None
        res = res[offset:end]
    return res


def count_atoms(residues):
    n_all = sum(1 for r in residues for _ in r.get_atoms())
    n_h = sum(1 for r in residues for a in r.get_atoms() if a.element == 'H')
    return n_all, n_h


def sc_centroid(r, include_h):
    a = [x for x in r.get_atoms()
         if x.get_name() not in BACKBONE and (include_h or x.element != 'H')]
    if not a:
        return r['CA'].get_coord()
    return np.mean([x.get_coord() for x in a], axis=0)


def edge_set(residues, include_h, cut=8.0):
    coords = np.array([sc_centroid(r, include_h) for r in residues])
    tree = cKDTree(coords)
    pairs = tree.query_pairs(r=cut, output_type='ndarray')
    return {(min(int(i), int(j)), max(int(i), int(j))) for i, j in pairs}


def main():
    df = pd.read_csv(os.path.join(DATA, 'benchmarks_s669_clean.csv'))
    align = pd.read_csv(os.path.join(DATA, 'alignment_map_s669.csv'))
    amap = {r['pdb_id']: r for _, r in align.iterrows()}
    mt_dir = os.path.join(DATA, 'mutant_structures_s669')

    wt_cache = {}
    wt_h_stats = []
    mt_h_stats = []
    rows = []

    for _, row in df.iterrows():
        pid = row['pdb_id']
        mut = row['mut_info']
        a = amap.get(pid)
        if a is None:
            continue
        base = f'pdb{pid.lower()}'
        mt_path = os.path.join(mt_dir, pid, mut, 'work', 'out', f'{base}_1.pdb')
        if not os.path.exists(mt_path):
            continue
        wt_path = os.path.join(DATA, 'structures', f'{base}.ent')
        if not os.path.exists(wt_path):
            continue
        L = len(row['wt_seq'])
        chain, off = a['chain'], int(a['offset'])

        if pid not in wt_cache:
            wr = load_chain(wt_path, chain, off, L)
            if wr is None:
                wt_cache[pid] = None
            else:
                na, nh = count_atoms(wr)
                wt_h_stats.append((pid, na, nh))
                wt_cache[pid] = wr
        wr = wt_cache[pid]
        if wr is None:
            continue
        mr = load_chain(mt_path, chain, off, L)
        if mr is None:
            continue
        na, nh = count_atoms(mr)
        mt_h_stats.append((pid, na, nh))

        # (A) 既有定义：WT 含氢，MT 实际无氢
        WA = edge_set(wr, include_h=True)
        MA = edge_set(mr, include_h=True)
        # (B) 一致定义：都排除氢
        WB = edge_set(wr, include_h=False)
        MB = edge_set(mr, include_h=False)

        rows.append({
            'pdb_id': pid, 'mut_info': mut,
            'wt_edges_A': len(WA), 'mt_edges_A': len(MA),
            'broken_A': len(WA - MA), 'formed_A': len(MA - WA),
            'wt_edges_B': len(WB), 'mt_edges_B': len(MB),
            'broken_B': len(WB - MB), 'formed_B': len(MB - WB),
        })

    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(DATA, 'hydrogen_bias_audit.csv'), index=False)

    print(f'处理 {len(out)} 个突变对')
    print()
    print('=== 氢原子出现情况 ===')
    wt_with_h = sum(1 for _, na, nh in wt_h_stats if nh > 0)
    mt_with_h = sum(1 for _, na, nh in mt_h_stats if nh > 0)
    print(f'  WT 结构（去重后 {len(wt_h_stats)} 个）含氢的: {wt_with_h}')
    print(f'  MT 结构（{len(mt_h_stats)} 个）含氢的: {mt_with_h}')
    print()
    print('=== 编辑统计对比 ===')
    for tag, bcol, fcol in [('(A) 既有定义（WT含氢/MT无氢）', 'broken_A', 'formed_A'),
                            ('(B) 一致定义（均排除氢）', 'broken_B', 'formed_B')]:
        b = out[bcol]
        f = out[fcol]
        tot = b + f
        print(f'  {tag}:')
        print(f'    断边 均值 {b.mean():.2f} 中位 {b.median():.0f}')
        print(f'    成边 均值 {f.mean():.2f} 中位 {f.median():.0f}')
        print(f'    至少 1 个接触变化的突变比例: {(tot > 0).mean()*100:.1f}%')
        print(f'    总编辑数 > 5 的比例: {(tot > 5).mean()*100:.1f}%')
    print()
    print('论文报告: 94.3% 至少一个接触变化; 断边均值 21.0; 成边均值 7.3')
    print(f'实测 (A): {(out["broken_A"]+out["formed_A"] > 0).mean()*100:.1f}%, '
          f'{out["broken_A"].mean():.1f}, {out["formed_A"].mean():.1f}')
    print(f'实测 (B): {(out["broken_B"]+out["formed_B"] > 0).mean()*100:.1f}%, '
          f'{out["broken_B"].mean():.1f}, {out["formed_B"].mean():.1f}')


if __name__ == '__main__':
    main()
