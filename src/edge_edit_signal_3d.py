"""
edge_edit_signal_3d.py — 用 3D 空间距离重做边编辑预测诊断
============================================================
关键修正：之前用"图距离"（hop 数）发现断边不集中（AUROC 0.486），
但图距离 ≠ 3D 距离。改用 3D 空间距离重新诊断。
"""
import ast
import os
import sys

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, Polypeptide

sys.path.insert(0, os.path.dirname(__file__))
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

from pathlib import Path as _P
DATA = str(_P(__file__).resolve().parent.parent / 'data')
VOL = {'A':88.6,'R':173.4,'N':114.1,'D':111.1,'C':108.5,'Q':143.8,'E':138.4,'G':60.1,'H':153.2,
       'I':166.7,'L':166.7,'K':168.6,'M':162.9,'F':189.9,'P':112.7,'S':89.0,'T':116.1,'W':227.8,'Y':193.6,'V':140.0}
BACKBONE = {'N','CA','C','O','OXT'}


def parse_edges(s):
    return ast.literal_eval(s) if isinstance(s, str) and s else []


def sc_centroid(res):
    atoms = [a for a in res.get_atoms() if a.get_name() not in BACKBONE]
    if atoms:
        return np.mean([a.get_coord() for a in atoms], axis=0)
    if res.has_id('CA'):
        return res['CA'].get_coord()
    return None


def get_coords(pdb_path):
    parser = PDBParser(QUIET=True)
    s = parser.get_structure('x', pdb_path)
    m = s[0]
    coords = []
    for c in m.get_chains():
        for r in c.get_residues():
            if Polypeptide.is_aa(r, standard=True):
                sc = sc_centroid(r)
                if sc is not None:
                    coords.append(sc)
        break
    return np.array(coords)


def main():
    ed = pd.read_csv(os.path.join(DATA, 'true_edits_s669_sc.csv'))
    X, y = [], []
    for _, row in ed.iterrows():
        pid = row['pdb_id']
        mi = int(row['mut_idx'])
        wt_pdb = os.path.join(DATA, 'structures', f'pdb{pid.lower()}.ent')
        if not os.path.exists(wt_pdb):
            continue
        coords = get_coords(wt_pdb)
        if mi >= len(coords):
            continue
        mut_coord = coords[mi]

        broken = {tuple(sorted(e)) for e in parse_edges(row['broken_edges'])}
        gpath = os.path.join(DATA, 'contact_graphs_s669_sc', f'{pid}.npz')
        if not os.path.exists(gpath):
            continue
        wt_ei = np.load(gpath)['edge_index']
        wt_edges = {tuple(sorted((int(a), int(b)))) for a, b in wt_ei.T}

        vol_diff = abs(VOL.get(str(row['mut_info'])[0], 100) - VOL.get(str(row['mut_info'])[-1], 100))

        for (a, b) in wt_edges:
            if a >= len(coords) or b >= len(coords):
                continue
            # 3D 距离：边端点到突变位点的最小距离
            d3d = min(np.linalg.norm(coords[a] - mut_coord), np.linalg.norm(coords[b] - mut_coord))
            X.append([d3d, vol_diff])
            y.append(1 if (a, b) in broken else 0)

    X = np.array(X)
    y = np.array(y)
    print(f'边级别样本: {len(y)} 条（断 {y.sum()}，保持 {len(y)-y.sum()}）')

    # 单特征 AUROC
    for i, name in [(0, '3D距离'), (1, '体积变化')]:
        auc = roc_auc_score(y, X[:, i])
        print(f'  单特征 {name} AUROC: {auc:.3f}')

    # 逻辑回归
    from sklearn.model_selection import cross_val_score
    clf = LogisticRegression(max_iter=1000)
    aucs = cross_val_score(clf, X, y, cv=5, scoring='roc_auc')
    print(f'逻辑回归（3D距离+体积变化）AUROC: {aucs.mean():.3f} ± {aucs.std():.3f}')


if __name__ == '__main__':
    main()
