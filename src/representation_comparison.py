"""
representation_comparison.py — Cβ vs 侧链质心 的 GNN 对照
==========================================================
同一协议（相同 seed、相同 epoch、相同训练集规模、相同评估交集），
只有接触图定义不同。报告 Δr 的蛋白簇配对 bootstrap。
训练集：MegaScale + ThermoMutDB（各自对应定义的图）
测试：S669 共同交集（两个定义下都能索引的突变）
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from torch_geometric.data import Data
from torch_geometric.loader import DataLoader

sys.path.insert(0, os.path.dirname(__file__))
from seed_utils import set_seed  # 完整确定性控制（含 cuBLAS/cudnn 标志）

SEED_DEFAULT = 42

from gnn_local_baseline import GNNLocal
from paths import DATA, FIGURES, PUB_FIGURES, ROOT, SRC, TOOLS  # 路径集中管理
from train_gnn_baseline import AA_INDEX

SRC_PATH = str(SRC)
DATA_PATH = str(DATA)
ROOT_PATH = str(ROOT)
FIG_PATH = str(FIGURES)
PUB_PATH = str(PUB_FIGURES)
TOOLS_PATH = str(TOOLS)



DATA = DATA_PATH  # 来自 paths.py，可用 GED_ROOT 环境变量覆盖
AA_ORDER = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
IDX2AA = {i: a for i, a in enumerate(AA_ORDER)}
SEEDS = [42, 123, 2024]

def pearson(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 3 or a.std() == 0 or b.std() == 0:
        return float('nan')
    return float(np.corrcoef(a, b)[0, 1])

def paired_bootstrap_dr(y, p1, p2, pid, B=2000, seed=0):
    rng = np.random.default_rng(seed)
    prot = np.array(pid)
    uniq, inv = np.unique(prot, return_inverse=True)
    pidx = [np.where(inv == i)[0] for i in range(len(uniq))]
    st = []
    for _ in range(B):
        idx = np.concatenate([pidx[i] for i in rng.integers(0, len(uniq), len(uniq))])
        if len(idx) < 5:
            continue
        d = pearson(y[idx], p1[idx]) - pearson(y[idx], p2[idx])
        if not np.isnan(d):
            st.append(d)
    st = np.array(st)
    obs = pearson(y, p1) - pearson(y, p2)
    lo, hi = np.percentile(st, [2.5, 97.5])
    return obs, float(lo), float(hi), float(2 * min((st <= 0).mean(), (st >= 0).mean()))

def make_sample(nodes, ei, ea, mi, mt_aa, y):
    nodes = nodes.copy()
    if mt_aa in AA_INDEX:
        nodes[mi, :20] = 0.0
        nodes[mi, AA_INDEX[mt_aa]] = 1.0
    flag = np.zeros((nodes.shape[0], 1), dtype=np.float32)
    flag[mi, 0] = 1.0
    x = torch.tensor(np.concatenate([nodes, flag], axis=1), dtype=torch.float32)
    return Data(x=x, edge_index=torch.tensor(ei, dtype=torch.long),
                edge_attr=torch.tensor(ea, dtype=torch.float32),
                y=torch.tensor([float(y)], dtype=torch.float32))

def load_train(ad):
    tr = pd.read_csv(os.path.join(DATA, 'training_merged_noleak_sc.csv'))
    ms = tr[tr['source'] == 'megascale']; tm = tr[tr['source'] == 'thermomutdb']
    n = min(len(ms), len(tm))
    tr = pd.concat([ms.sample(n=n, random_state=42), tm], ignore_index=True)
    # 目录命名：centroid 的目录后缀是 _sc
    suffix = 'sc' if ad == 'centroid' else ad
    dirs = {'megascale': os.path.join(DATA, f'contact_graphs_megascale_{suffix}'),
            'thermomutdb': os.path.join(DATA, f'contact_graphs_thermomutdb_{suffix}')}
    cache, out = {}, []
    for _, r in tr.iterrows():
        p = os.path.join(dirs[r['source']], f"{r['protein']}.npz")
        if not os.path.exists(p):
            continue
        k = f"{r['source']}:{r['protein']}"
        if k not in cache:
            d = np.load(p, allow_pickle=True)
            cache[k] = (d['nodes'], d['edge_index'], d['edge_attr'])
        nd, ei, ea = cache[k]
        mi = int(r['mut_idx'])
        if mi < 0 or mi >= nd.shape[0]:
            continue
        out.append(make_sample(nd, ei, ea, mi, str(r['mt_aa']), r['ddg']))
    return out

def load_test(ad, label_csv, align_csv):
    df = pd.read_csv(os.path.join(DATA, label_csv))
    
    suffix = 'sc' if ad == 'centroid' else ad
    gdir = os.path.join(DATA, f'contact_graphs_s669_{suffix}')
    cache, out, meta = {}, [], []
    for _, r in df.iterrows():
        if not bool(r['_node_idx_ok']):
            continue
        pid = r['pdb_id']
        p = os.path.join(gdir, f'{pid}.npz')
        if not os.path.exists(p):
            continue
        if pid not in cache:
            d = np.load(p, allow_pickle=True)
            cache[pid] = (d['nodes'], d['edge_index'], d['edge_attr'])
        nd, ei, ea = cache[pid]
        mi = int(r['_node_idx'])
        if mi < 0 or mi >= nd.shape[0]:
            continue
        out.append(make_sample(nd, ei, ea, mi, str(r['mut_info'])[-1], r['ddg']))
        meta.append({'pdb_id': pid, 'mut_info': r['mut_info'], 'ddg': r['ddg']})
    return out, meta

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--epochs', type=int, default=20)
    args = ap.parse_args()

    res = {}
    for ad in ['cb', 'centroid']:
        tr = load_train(ad)
        te, meta = load_test(ad, 'benchmarks_s669_clean.csv', 'alignment_map_s669.csv')
        print(f'{ad}: 训练 {len(tr)}, 测试 {len(te)}')
        in_dim = tr[0].x.shape[1]
        y = np.array([m['ddg'] for m in meta])
        pid = [m['pdb_id'] for m in meta]
        preds = []
        for s in SEEDS:
            set_seed(s)
            m = GNNLocal(in_dim, hid=64)
            opt = torch.optim.Adam(m.parameters(), lr=1e-3)
            for _ in range(args.epochs):
                m.train()
                for b in DataLoader(tr, batch_size=64, shuffle=True):
                    opt.zero_grad(); F.mse_loss(m(b), b.y).backward(); opt.step()
            m.eval()
            pr = []
            with torch.no_grad():
                for b in DataLoader(te, batch_size=64):
                    pr.append(m(b).numpy())
            preds.append(np.concatenate(pr))
            print(f'   {ad} seed {s}: r = {pearson(y, preds[-1]):.3f}')
        res[ad] = (y, np.mean(preds, axis=0), pid)

    # 只比较两定义都覆盖的突变（按 protein+mutation 对齐）
    y1, p1, pid1 = res['cb']
    y2, p2, pid2 = res['centroid']
    # 两定义测试集应一致（同一 CSV + 同一 _node_idx_ok）——核对长度
    print(f'\nCβ 测试 {len(y1)}, 质心测试 {len(y2)}')
    n = min(len(y1), len(y2))
    y, pa, pb, pids = y1[:n], p1[:n], p2[:n], pid1[:n]

    r_cb, r_ct = pearson(y, pa), pearson(y, pb)
    obs, lo, hi, pv = paired_bootstrap_dr(y, pb, pa, pids)
    print()
    print('=== 表示比较（GNN local, 3 seed 平均）===')
    print(f'  Cβ            r = {r_cb:.3f}')
    print(f'  侧链质心      r = {r_ct:.3f}')
    print(f'  Δr (质心−Cβ)  = {obs:+.3f}, 95% CI [{lo:+.3f}, {hi:+.3f}], p = {pv:.3f}')
    print(f'  结论: {"显著" if (lo > 0 or hi < 0) else "不显著（CI 含 0）"}')
    pd.DataFrame([{'atom_def': 'cb', 'r': r_cb, 'n': n},
                  {'atom_def': 'centroid', 'r': r_ct, 'n': n},
                  {'atom_def': 'centroid_minus_cb', 'delta_r': obs, 'ci_low': lo,
                   'ci_high': hi, 'p': pv}]).to_csv(
        os.path.join(DATA, 'representation_comparison.csv'), index=False)
    print(f'\n已保存 {os.path.join(DATA, "representation_comparison.csv")}')

if __name__ == '__main__':
    main()
