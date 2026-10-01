# -*- coding: utf-8 -*-
"""
test_trep.py — 候选方法 TREP（按边类型分解的边缘池化）对照实验
================================================================
假设：边类型信息在常规节点级消息传递中被稀释（gnn_edge 不优于 gnn_local）。
      若在突变位点邻域**按边类型分别池化**，能否恢复这部分信息？

对照（同一协议、同一种子）：
  gnn_local   基线：节点级局部池化
  deep_gine   强基线：6 层 GINE + 注意力池化（无类型分解）
  TREP        候选：4 层 GINE + 按类型边缘池化

判定：TREP 相对 deep_gine / gnn_local 的配对 Δr（蛋白簇 bootstrap）
"""
import os
import sys
import time

import numpy as np
import torch
import torch.nn.functional as F
from torch_geometric.loader import DataLoader

sys.path.insert(0, os.path.dirname(__file__))
from ladder_common import load_test, load_train
from seed_utils import set_seed
from strong_backbones import DeepGINE, TypePoolGINE
from gnn_local_baseline import GNNLocal

DEV = 'cuda' if torch.cuda.is_available() else 'cpu'


def pearson(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return float(np.corrcoef(a, b)[0, 1]) if a.std() > 0 and b.std() > 0 else float('nan')


def cluster_ci(y, p, pid, B=1000, seed=0):
    rng = np.random.default_rng(seed)
    uniq, inv = np.unique(np.array(pid), return_inverse=True)
    pidx = [np.where(inv == i)[0] for i in range(len(uniq))]
    st = []
    for _ in range(B):
        idx = np.concatenate([pidx[i] for i in rng.integers(0, len(uniq), len(uniq))])
        if len(idx) >= 3:
            v = pearson(y[idx], p[idx])
            if not np.isnan(v):
                st.append(v)
    return (float(np.percentile(st, 2.5)), float(np.percentile(st, 97.5))) if st else (np.nan, np.nan)


def paired(y, p1, p2, pid, B=2000, seed=0):
    rng = np.random.default_rng(seed)
    uniq, inv = np.unique(np.array(pid), return_inverse=True)
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
    p = 2 * min((st <= 0).mean(), (st >= 0).mean())
    return obs, float(lo), float(hi), float(min(p, 1.0))


def run(tag, factory, tr, te, y, pid, seeds, epochs=20, lr=1e-3):
    ps = []
    for s in seeds:
        set_seed(s)
        m = factory().to(DEV)
        opt = torch.optim.Adam(m.parameters(), lr=lr)
        dl = DataLoader(tr, batch_size=64, shuffle=True)
        t0 = time.time()
        for _ in range(epochs):
            m.train()
            for b in dl:
                b = b.to(DEV)
                opt.zero_grad()
                F.mse_loss(m(b), b.y).backward()
                opt.step()
        m.eval()
        with torch.no_grad():
            p = np.concatenate([m(b.to(DEV)).cpu().numpy()
                                for b in DataLoader(te, batch_size=128)])
        ps.append(p)
        print(f'    seed {s}: r={pearson(y, p):+.4f} ({time.time()-t0:.0f}s)', flush=True)
    pm = np.mean(ps, axis=0)
    r = pearson(y, pm)
    lo, hi = cluster_ci(y, pm, pid)
    npar = sum(x.numel() for x in factory().parameters())
    print(f'  {tag:<12} 参数 {npar/1000:>7.1f}K  种子均值 r={r:.4f} [{lo:.3f},{hi:.3f}]')
    return pm, r


def main():
    ad = os.environ.get('AD', 'centroid')
    seeds = [42, 123, 2024]
    tr = load_train(ad, with_coords=False)
    te, meta = load_test(ad, 's669', with_coords=False)
    y = np.array([m['ddg'] for m in meta])
    pid = [m['pdb_id'] for m in meta]
    in_dim = tr[0].x.shape[1]
    print(f'定义 {ad}: 训练 {len(tr)}, 测试 {len(te)}, 种子 {seeds}\n', flush=True)

    print('=== 对照实验 ===')
    p_loc, _ = run('gnn_local', lambda: GNNLocal(in_dim, hid=64), tr, te, y, pid, seeds)
    p_dg, _ = run('deep_gine', lambda: DeepGINE(in_dim, hid=128, n_layers=6, k_hop=2),
                  tr, te, y, pid, seeds)
    p_trep, r_trep = run('TREP', lambda: TypePoolGINE(in_dim, hid=128, n_layers=4, k_hop=2),
                         tr, te, y, pid, seeds)

    print('\n=== 配对检验 ===')
    for a, b, pa, pb in [('TREP', 'deep_gine', p_trep, p_dg),
                         ('TREP', 'gnn_local', p_trep, p_loc),
                         ('deep_gine', 'gnn_local', p_dg, p_loc)]:
        obs, lo, hi, pv = paired(y, pa, pb, pid)
        sig = '显著' if (lo > 0 or hi < 0) else '不显著'
        print(f'  {a} − {b}: Δr={obs:+.4f} [{lo:+.4f},{hi:+.4f}] p={pv:.3f}  {sig}')

    print(f'\n结论：TREP 相对强基线的增量为 '
          f'{paired(y, p_trep, p_dg, pid)[0]:+.4f} '
          f'（{"成立" if paired(y, p_trep, p_dg, pid)[1] > 0 else "不成立"}）')


if __name__ == '__main__':
    main()
