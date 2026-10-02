# -*- coding: utf-8 -*-
"""within_arch_ladder.py — 同架构容量阶梯（回应 R2-M3）

审稿人 R2-M3：原阶梯的四个稳定档位跨 GCN → GINE → EGNN，**容量与架构混淆**，
因此"表示效应不随容量增长"可能只是架构差异，无法排除"编码器能力不足"。

本脚本在**单一架构（DeepGINE）内**只变宽度与层数，构造容量阶梯：

  hid=64,  layers=2    ~50k
  hid=96,  layers=4   ~185k
  hid=128, layers=6   ~460k
  hid=160, layers=6   ~715k

每个档位 × 3 种子 × {ca, centroid}，其余训练设置（数据、划分、轮数）与原阶梯一致。
输出：data/within_arch_ladder_results.csv / _predictions.csv


STATUS: not run to completion. Only 10 of the 24 planned runs finished
before this experiment was stopped, and no result from it appears in the
manuscript, the supplementary material or any released table. The script is
kept because it is the direct test of reviewer concern R2-M3 (the capacity
ladder confounds capacity with architecture); running it to completion would
settle that concern. Do not cite partial output of this script.
"""
import argparse
import os
import sys
import time

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from torch_geometric.loader import DataLoader

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ladder_common import build_model, load_test, load_train  # noqa: E402
from seed_utils import set_seed  # noqa: E402
from strong_backbones import DeepGINE  # noqa: E402

DEV = 'cuda' if torch.cuda.is_available() else 'cpu'
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'data')

# 同架构内的容量档位（只变 hid 与 n_layers）
RUNGS = [
    ('wine_50k', dict(hid=64, n_layers=2, k_hop=1)),
    ('wine_185k', dict(hid=96, n_layers=4, k_hop=2)),
    ('wine_460k', dict(hid=128, n_layers=6, k_hop=2)),
    ('wine_715k', dict(hid=160, n_layers=6, k_hop=2)),
]
DEFS = ['ca', 'centroid']


def pearson(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 3 or a.std() == 0 or b.std() == 0:
        return float('nan')
    return float(np.corrcoef(a, b)[0, 1])


def cluster_ci(y, p, pid, B=1000, seed=0):
    rng = np.random.default_rng(seed)
    prot = np.array(pid)
    uniq, inv = np.unique(prot, return_inverse=True)
    pidx = [np.where(inv == i)[0] for i in range(len(uniq))]
    st = []
    for _ in range(B):
        idx = np.concatenate([pidx[i] for i in rng.integers(0, len(uniq), len(uniq))])
        v = pearson(y[idx], p[idx])
        if not np.isnan(v):
            st.append(v)
    if not st:
        return float('nan'), float('nan')
    return float(np.percentile(st, 2.5)), float(np.percentile(st, 97.5))


def train_eval(rung, cfg, ad, seed, epochs, tr_cache, te_cache):
    in_dim = tr_cache[ad][0][0].x.shape[1]
    set_seed(seed)
    model = DeepGINE(in_dim, edge_dim=2, out=1, **cfg).to(DEV)
    opt = torch.optim.Adam(model.parameters(), lr=1e-3)
    tr = tr_cache[ad][0]
    loader = DataLoader(tr, batch_size=64, shuffle=True)
    for _ in range(epochs):
        model.train()
        for b in loader:
            b = b.to(DEV)
            opt.zero_grad()
            F.mse_loss(model(b).view(-1), b.y.view(-1)).backward()
            opt.step()
    model.eval()
    te, meta = te_cache[ad]
    te_loader = DataLoader(te, batch_size=128)
    preds = []
    with torch.no_grad():
        for b in te_loader:
            preds.append(model(b.to(DEV)).view(-1).cpu().numpy())
    preds = np.concatenate(preds)
    y = np.array([m['ddg'] for m in meta])
    pid = [m['pdb_id'] for m in meta]
    r = pearson(y, preds)
    lo, hi = cluster_ci(y, preds, pid)
    return r, lo, hi, y, preds, meta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--epochs', type=int, default=20)
    ap.add_argument('--seeds', default='42,123,2024')
    ap.add_argument('--defs', default=','.join(DEFS))
    ap.add_argument('--rungs', default=','.join(r for r, _ in RUNGS))
    args = ap.parse_args()
    seeds = [int(s) for s in args.seeds.split(',')]
    defs = args.defs.split(',')
    want = set(args.rungs.split(','))

    print(f'device={DEV} | epochs={args.epochs} | seeds={seeds} | defs={defs}', flush=True)
    t0 = time.time()
    tr_cache, te_cache = {}, {}
    for ad in defs:
        tr_cache[ad] = (load_train(ad), None)
        te_cache[ad] = load_test(ad, 's669')
        print(f'  {ad}: train {len(tr_cache[ad][0])}, test {len(te_cache[ad][0])} '
              f'({time.time()-t0:.0f}s)', flush=True)

    rows, prows = [], []
    for rung, cfg in RUNGS:
        if rung not in want:
            continue
        for ad in defs:
            for seed in seeds:
                ts = time.time()
                r, lo, hi, y, preds, meta = train_eval(
                    rung, cfg, ad, seed, args.epochs, tr_cache, te_cache)
                n_par = sum(p.numel() for p in
                            DeepGINE(21, edge_dim=2, out=1, **cfg).parameters())
                rows.append({'rung': rung, 'atom_def': ad, 'seed': seed,
                             'n_params': n_par, 'n_test': len(y), 'r': r,
                             'ci_lo': lo, 'ci_hi': hi})
                for m, p in zip(meta, preds):
                    prows.append({'rung': rung, 'atom_def': ad, 'seed': seed,
                                  'protein_id': m['pdb_id'], 'mutation_id': m['mut_info'],
                                  'y_true': m['ddg'], 'y_pred': float(p)})
                print(f'  {rung:<10} {ad:<9} s{seed:<5} params={n_par:>7,} '
                      f'r={r:+.4f}  [{time.time()-ts:.0f}s]', flush=True)
                pd.DataFrame(rows).to_csv(
                    os.path.join(OUT, 'within_arch_ladder_results.csv'), index=False)
                pd.DataFrame(prows).to_csv(
                    os.path.join(OUT, 'within_arch_ladder_predictions.csv'), index=False)

    print(f'\n完成，用时 {(time.time()-t0)/60:.1f} 分钟', flush=True)
    res = pd.DataFrame(rows)
    print(res.groupby(['rung', 'n_params'])['r']
          .agg(['mean', 'std', 'count']).to_string())


if __name__ == '__main__':
    main()
