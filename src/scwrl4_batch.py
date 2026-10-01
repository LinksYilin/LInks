"""
scwrl4_batch.py — 用 SCWRL4 批量建模突变体，并与 FoldX 对比
==============================================================
流程（每个突变）：
  1. 从 WT 结构切出本流程使用的链片段，写成单链 PDB
  2. 用该片段的序列构造突变序列（突变位点换成新氨基酸，全大写=全局重排）
  3. 运行 SCWRL4（-h 不输出氢）
  4. 读回突变体结构，构建侧链质心接触图（排除氢）
  5. 与 WT 图对比，得到断边/成边
  6. 保存，供与 FoldX 结果比较

关键防坑：
  - SCWRL4 输出残基顺序可能与输入不同 → 一律按 (链, 残基号) 映射，不按位置
  - SCWRL4 会丢弃骨架不完整的残基 → 记录残基集合变化
"""
import argparse
import os
import subprocess
import sys

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa, protein_letters_3to1

sys.path.insert(0, os.path.dirname(__file__))
from contact_graph_defs import pairwise_within
from paths import DATA, FIGURES, PUB_FIGURES, ROOT, SRC, TOOLS  # 路径集中管理

SRC_PATH = str(SRC)
DATA_PATH = str(DATA)
ROOT_PATH = str(ROOT)
FIG_PATH = str(FIGURES)
PUB_PATH = str(PUB_FIGURES)
TOOLS_PATH = str(TOOLS)




SCWRL = r'D:\GED_mutation\tools\scwrl4\Scwrl4.exe'

WORKROOT = r'D:\GED_mutation\tools\scwrl4\runs'
TH, MAXC = 8.0, 10.0
BACKBONE = {'N', 'CA', 'C', 'O', 'OXT'}
parser = PDBParser(QUIET=True)


def chain_residues(path, chain_id):
    st = parser.get_structure('p', path)
    chs = [c for c in st[0].get_chains() if c.id == chain_id]
    if not chs:
        return None
    return [r for r in chs[0] if r.has_id('CA')]


def res_aa(res):
    nm = res.get_resname().strip()
    if nm == 'MSE':
        return 'M'
    if not is_aa(res, standard=True):
        return None
    return protein_letters_3to1.get(nm)


def write_segment(residues, dst, chain_id):
    """把残基列表写成单链 PDB（排除氢）。"""
    with open(dst, 'w') as f:
        k = 0
        for res in residues:
            aa = res_aa(res)
            if aa is None:
                continue
            for a in res.get_atoms():
                if a.element == 'H':
                    continue
                c = a.get_coord()
                f.write(f"ATOM  {k+1:>5} {a.get_name():<4} {res.get_resname():>3} "
                        f"{chain_id}{res.id[1]:>4}{str(res.id[2]).strip():<1}   "
                        f"{c[0]:>8.3f}{c[1]:>8.3f}{c[2]:>8.3f}  1.00  0.00          {a.element:>2}\n")
                k += 1
            f.write('TER\n')
        f.write('END\n')


def centroid_edges(residues):
    coords = []
    for r in residues:
        atoms = [a for a in r.get_atoms() if a.get_name() not in BACKBONE and a.element != 'H']
        coords.append(np.mean([a.get_coord() for a in atoms], axis=0) if atoms else r['CA'].get_coord())
    coords = np.array(coords)
    return {k for k, d in pairwise_within(coords, MAXC).items() if d < TH}, coords


def ca_edges(residues):
    coords = np.array([r['CA'].get_coord() for r in residues])
    return {k for k, d in pairwise_within(coords, MAXC).items() if d < TH}


def cb_edges(residues):
    coords = np.array([r['CB'].get_coord() if r.has_id('CB') else r['CA'].get_coord()
                       for r in residues])
    return {k for k, d in pairwise_within(coords, MAXC).items() if d < TH}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--limit', type=int, default=0)
    ap.add_argument('--out', default=os.path.join(DATA, 'edits_scwrl4.csv'))
    ap.add_argument('--local-radius', type=float, default=0.0,
                    help='>0 时仅重排突变位点该半径内的残基（与 FoldX 范围对齐）；0=全局重排')
    ap.add_argument('--label_csv', default=os.path.join(DATA, 'benchmarks_s669_clean.csv'))
    ap.add_argument('--align_csv', default=os.path.join(DATA, 'alignment_map_s669.csv'))
    args = ap.parse_args()

    df = pd.read_csv(args.label_csv)
    align = pd.read_csv(args.align_csv)
    amap = {r['pdb_id']: r for _, r in align.iterrows()}
    os.makedirs(WORKROOT, exist_ok=True)

    rows = []
    n_ok = n_fail = 0
    todo = df if args.limit <= 0 else df.head(args.limit)

    for _, row in todo.iterrows():
        pid, mut = row['pdb_id'], row['mut_info']
        a = amap.get(pid)
        if a is None:
            continue
        wt_path = os.path.join(DATA, 'structures', f'pdb{pid.lower()}.ent')
        if not os.path.exists(wt_path):
            continue
        L = len(row['wt_seq'])
        chain, off = a['chain'], int(a['offset'])
        wr = chain_residues(wt_path, chain)
        if wr is None:
            continue
        seg = wr[off:off + L] if (off > 0 or L) else wr
        seg = [r for r in seg if res_aa(r) is not None]
        if len(seg) != L:
            n_fail += 1
            continue
        wseq = ''.join(res_aa(r) for r in seg)
        if '_node_idx' not in df.columns or pd.isna(row['_node_idx']):
            raise ValueError(
                f"缺少已验证的 _node_idx（{(row['pdb_id'], row['mut_info'])}）；"
                "拒绝回退到 _pdb_res_idx，请先运行 fix_index_bug.py")
        mi = int(row['_node_idx'])
        if not (0 <= mi < len(wseq)):
            n_fail += 1
            continue
        # 构造突变序列
        mseq = list(wseq)
        mseq[mi] = mut[-1]
        if args.local_radius > 0:
            # 局部模式：仅突变位点及其邻域大写（重排），其余小写（保留原构象）
            mc = seg[mi]['CA'].get_coord()
            keep_upper = set()
            for k, r_ in enumerate(seg):
                # 用侧链质心判断邻域（与接触图定义一致）
                atoms = [a for a in r_.get_atoms()
                         if a.get_name() not in BACKBONE and a.element != 'H']
                cen = np.mean([a.get_coord() for a in atoms], axis=0) if atoms else r_['CA'].get_coord()
                if np.linalg.norm(cen - mc) <= args.local_radius:
                    keep_upper.add(k)
            keep_upper.add(mi)
            mseq = ''.join(c.upper() if k in keep_upper else c.lower()
                           for k, c in enumerate(mseq))
            # 小写残基的侧链必须完整，否则 SCWRL4 会报错 → 不完整的改回大写
            fixed = []
            for k, c in enumerate(mseq):
                if c.islower():
                    atoms = [a for a in seg[k].get_atoms()
                             if a.get_name() not in BACKBONE and a.element != 'H']
                    if len(atoms) < 1:
                        fixed.append(c.upper())
                    else:
                        fixed.append(c)
                else:
                    fixed.append(c)
            mseq = ''.join(fixed)
        else:
            mseq = ''.join(mseq).upper()

        d = os.path.join(WORKROOT, pid, mut)
        os.makedirs(d, exist_ok=True)
        wt_seg = os.path.join(d, 'wt_seg.pdb')
        write_segment(seg, wt_seg, chain)
        with open(os.path.join(d, 'seq_mut.txt'), 'w') as f:
            f.write(mseq + '\n')

        tag = 'global' if args.local_radius <= 0 else f'local{args.local_radius:g}A'
        out_pdb = os.path.join(d, f'mt_scwrl4_{tag}.pdb')
        r = subprocess.run([SCWRL, '-i', 'wt_seg.pdb', '-o', os.path.basename(out_pdb),  # noqa: PLW1510 (返回值由调用方检查)
                            '-s', 'seq_mut.txt', '-h'],
                           cwd=d, capture_output=True, text=True)
        if r.returncode != 0 or not os.path.exists(out_pdb):
            n_fail += 1
            continue

        mr = chain_residues(out_pdb, chain)
        if mr is None:
            n_fail += 1
            continue
        mr = [x for x in mr if res_aa(x) is not None]
        # 按残基号映射，保证与 WT 顺序对齐
        wt_by_id = {r.id[1]: r for r in seg}
        mt_by_id = {r.id[1]: r for r in mr}
        common_ids = sorted(set(wt_by_id) & set(mt_by_id))
        if len(common_ids) != len(seg):
            n_fail += 1
            continue
        wt_r = [wt_by_id[i] for i in common_ids]
        mt_r = [mt_by_id[i] for i in common_ids]

        W_sc, _ = centroid_edges(wt_r)
        M_sc, _ = centroid_edges(mt_r)
        W_ca = ca_edges(wt_r)
        M_ca = ca_edges(mt_r)
        W_cb = cb_edges(wt_r)
        M_cb = cb_edges(mt_r)

        rows.append({
            'pdb_id': pid, 'mut_info': mut, 'n_res': len(common_ids),
            'mode': 'global' if args.local_radius <= 0 else f'local{args.local_radius:g}A',
            'wt_edges_c': len(W_sc), 'mt_edges_c': len(M_sc),
            'broken_scwrl4': len(W_sc - M_sc), 'formed_scwrl4': len(M_sc - W_sc),
            'broken_ca': len(W_ca - M_ca), 'formed_ca': len(M_ca - W_ca),
            'broken_cb': len(W_cb - M_cb), 'formed_cb': len(M_cb - W_cb),
            'n_edit_ca': len(W_ca - M_ca) + len(M_ca - W_ca),
            'n_edit_cb': len(W_cb - M_cb) + len(M_cb - W_cb),
            'n_edit_centroid': len(W_sc - M_sc) + len(M_sc - W_sc),
        })
        n_ok += 1
        if n_ok % 50 == 0:
            print(f'  已完成 {n_ok}')

    out = pd.DataFrame(rows)
    out.to_csv(args.out, index=False)
    print(f'\n成功 {n_ok}, 失败/跳过 {n_fail}')
    print(f'已保存 {args.out}')
    if len(out):
        for ad, col in [('Cα', 'n_edit_ca'), ('Cβ', 'n_edit_cb'), ('侧链质心', 'n_edit_centroid')]:
            print(f'  {ad}: 断边均值 '
                  f'{out["broken_" + ("scwrl4" if ad == "侧链质心" else ad.lower().replace("α", "a").replace("β", "b"))].mean():.2f}, '
                  f'≥1变化 {(out[col] > 0).mean()*100:.1f}%')


if __name__ == '__main__':
    main()
