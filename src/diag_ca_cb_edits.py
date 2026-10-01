"""
diag_ca_cb_edits.py — 查清 Cα/Cβ 非零编辑的来源
=================================================
假设 1：甘氨酸（无 CB）在突变后引入 CB → 代表原子从 CA 变成 CB → 产生"伪编辑"
假设 2：阈值边界噪声（FoldX 输出坐标微小变化使恰好 8.0 Å 的边翻转）
"""
import os

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from scipy.spatial import cKDTree

DATA = r'D:\GED_mutation\data'
BACKBONE = {'N', 'CA', 'C', 'O', 'OXT'}
parser = PDBParser(QUIET=True)


def load_chain(path, chain_id, offset, length):
    st = parser.get_structure('p', path)
    chs = [c for c in st[0].get_chains() if c.id == chain_id]
    if not chs:
        return None
    res = [r for r in chs[0] if r.has_id('CA')]
    if offset > 0 or length is not None:
        end = (offset + length) if length is not None else None
        res = res[offset:end]
    return res


def coords_for(residues, atom_def):
    out = []
    for r in residues:
        if atom_def == 'ca':
            out.append(r['CA'].get_coord())
        elif atom_def == 'cb':
            out.append(r['CB'].get_coord() if r.has_id('CB') else r['CA'].get_coord())
        else:
            a = [x for x in r.get_atoms() if x.get_name() not in BACKBONE and x.element != 'H']
            out.append(np.mean([x.get_coord() for x in a], axis=0) if a else r['CA'].get_coord())
    return np.array(out)


def edges(coords, cut=8.0):
    tree = cKDTree(coords)
    pr = tree.query_pairs(r=cut, output_type='ndarray')
    return {(min(int(i), int(j)), max(int(i), int(j))) for i, j in pr}


def main():
    df = pd.read_csv(os.path.join(DATA, 'benchmarks_s669_clean.csv'))
    align = pd.read_csv(os.path.join(DATA, 'alignment_map_s669.csv'))
    amap = {r['pdb_id']: r for _, r in align.iterrows()}
    mt_dir = os.path.join(DATA, 'mutant_structures_s669')

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
        chain, off = a['chain'], int(a['offset'])
        wr = load_chain(wt_path, chain, off, L)
        mr = load_chain(mt_path, chain, off, L)
        if wr is None or mr is None:
            continue
            wt_aa = mut[0]
        mt_aa = mut[-1]
        # 突变位点是否涉及甘氨酸（原子切换）
        gly_switch = (wt_aa == 'G' and mt_aa != 'G') or (wt_aa != 'G' and mt_aa == 'G')

        rec = {'pdb_id': pid, 'mut_info': mut, 'gly_switch': gly_switch,
               'wt_aa': wt_aa, 'mt_aa': mt_aa}
        for ad in ['ca', 'cb', 'centroid']:
            W = edges(coords_for(wr, ad))
            M = edges(coords_for(mr, ad))
            rec[f'broken_{ad}'] = len(W - M)
            rec[f'formed_{ad}'] = len(M - W)
        rows.append(rec)

    out = pd.DataFrame(rows)
    print(f'总突变对: {len(out)}')
    print(f'涉及甘氨酸原子切换的: {out["gly_switch"].sum()}')
    print()

    print('=== 按是否涉及甘氨酸分组 ===')
    for label, sub in [('不含 GLY 切换', out[~out['gly_switch']]),
                       ('含 GLY 切换', out[out['gly_switch']])]:
        print(f'{label} (n={len(sub)}):')
        for ad in ['ca', 'cb', 'centroid']:
            tot = sub[f'broken_{ad}'] + sub[f'formed_{ad}']
            print(f'   {ad:<9} 断边均值 {sub[f"broken_{ad}"].mean():>6.2f}  '
                  f'成边均值 {sub[f"formed_{ad}"].mean():>6.2f}  ≥1变化 {(tot > 0).mean()*100:>5.1f}%')
        print()

    # Cα 在"不含 GLY 切换"组中仍非零 → 说明是阈值边界噪声
    sub = out[~out['gly_switch']]
    n_ca_edit = ((sub['broken_ca'] + sub['formed_ca']) > 0).sum()
    print(f'=== Cα 在无 GLY 切换组中的非零编辑：{n_ca_edit}/{len(sub)} '
          f'({n_ca_edit/len(sub)*100:.1f}%) ===')
    print('   Cα 在 WT/MT 都是 Cα 原子，理论上应完全一致；')
    print('   非零说明 FoldX 输出坐标有微小变化导致恰好位于阈值附近的边翻转。')

    out.to_csv(os.path.join(DATA, 'ca_cb_edit_diagnosis.csv'), index=False)
    print(f'\n已保存 {os.path.join(DATA, "ca_cb_edit_diagnosis.csv")}')


if __name__ == '__main__':
    main()
