# -*- coding: utf-8 -*-
"""
run_ladder.py — 能力阶梯主实验（WS1）
=======================================
对每个（骨架 × 图定义 × 种子）训练并评估，回答：
  "表示效应是否随模型能力增强而出现？"

骨架（弱 → 强）：
  gnn_global   GCN + 均值池化              （最弱）
  gnn_local    GCN + 突变位点局部池化
  deep_gine    6 层 GINE + 注意力池化 + 局部读出
  egnn         E(3) 等变网络（用 3D 坐标）

定义：ca / cb / centroid / allatom
输出：
  data/ladder_results.csv       每个组合的指标
  data/ladder_predictions.csv   逐样本预测（供配对检验）
"""
import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from torch_geometric.loader import DataLoader

sys.path.insert(0, os.path.dirname(__file__))
from ladder_common import DEFS, build_model, load_test, load_train
from seed_utils import describe, set_seed

DEV = 'cuda' if torch.cuda.is_available() else 'cpu'
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'data')

MODELS = {
    'gnn_global': dict(kind='gcn_global'),        # 最弱：GCN + 均值池化
    'gnn_local': dict(kind='gcn_local'),          # GCN + 突变位点局部池化
    'gnn_edge': dict(kind='gine_edge'),           # GINEConv 边感知（浅层）
    'deep_gine': dict(kind='deep_gine'),          # 6 层 GINE + 注意力池化
    'egnn': dict(kind='egnn'),                    # E(3) 等变网络
}
NEEDS_COORDS = {'egnn'}


def pearson(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 3 or a.std() == 0 or b.std() == 0:
        return float('nan')
    return float(np.corrcoef(a, b)[0, 1])


def spearman(a, b):
    from scipy.stats import spearmanr
    return float(spearmanr(a, b).statistic)


def cluster_ci(y, p, pid, B=1000, seed=0):
    rng = np.random.default_rng(seed)
    prot = np.array(pid)
    uniq, inv = np.unique(prot, return_inverse=True)
    pidx = [np.where(inv == i)[0] for i in range(len(uniq))]
    st = []
    for _ in range(B):
        idx = np.concatenate([pidx[i] for i in rng.integers(0, len(uniq), len(uniq))])
        if len(idx) >= 3:
            v = pearson(y[idx], p[idx])
            if not np.isnan(v):
                st.append(v)
    return (float(np.percentile(st, 2.5)), float(np.percentile(st, 97.5))) if st else (np.nan, np.nan)


def make_model(name, in_dim):
    kind = MODELS[name]['kind']
    if kind == 'gcn_global':
        # ★ 真正的全局均值池化（最弱档）
        from train_gnn_baseline import GNNRegressor
        return GNNRegressor(in_dim, hid=64)
    if kind == 'gcn_local':
        # 突变位点局部池化
        from gnn_local_baseline import GNNLocal
        return GNNLocal(in_dim, hid=64)
    if kind == 'gine_edge':
        # 边感知 GINEConv + 局部池化
        from strong_backbones import DeepGINE
        return DeepGINE(in_dim, edge_dim=2, hid=64, n_layers=2, k_hop=1)
    if kind == 'deep_gine':
        from strong_backbones import DeepGINE
        return DeepGINE(in_dim, edge_dim=2, hid=128, n_layers=6, k_hop=2)
    if kind == 'egnn':
        from strong_backbones import EGNN
        # Explicit scaled=False reproduces the 508,934-parameter configuration reported
        # in the manuscript. The strong_backbones default (scaled=True) adds four
        # coordinate-scale parameters (508,938) and is only the stability-control variant.
        return EGNN(in_dim, edge_dim=2, hid=128, n_layers=4, k_hop=2,
                    update_coords=True, scaled=False)
    raise ValueError(name)


def train_eval(model_name, samples, test_samples, meta, seed, epochs, lr=1e-3):
    set_seed(seed)
    in_dim = samples[0].x.shape[1]
    m = make_model(model_name, in_dim).to(DEV)
    n_par = sum(p.numel() for p in m.parameters())
    opt = torch.optim.Adam(m.parameters(), lr=lr)
    dl = DataLoader(samples, batch_size=64, shuffle=True)
    t0 = time.time()
    for _ in range(epochs):
        m.train()
        for b in dl:
            b = b.to(DEV)
            opt.zero_grad()
            F.mse_loss(m(b), b.y).backward()
            opt.step()
    train_t = time.time() - t0
    m.eval()
    with torch.no_grad():
        p = np.concatenate([m(b.to(DEV)).cpu().numpy()
                            for b in DataLoader(test_samples, batch_size=128)])
    return p, n_par, train_t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--models', default=','.join(MODELS))
    ap.add_argument('--defs', default=','.join(DEFS))
    ap.add_argument('--seeds', default='42,123,2024')
    ap.add_argument('--epochs', type=int, default=20)
    ap.add_argument('--bench', default='s669')
    ap.add_argument('--tag', default='ladder')
    args = ap.parse_args()

    models = [x for x in args.models.split(',') if x]
    defs = [x for x in args.defs.split(',') if x]
    seeds = [int(x) for x in args.seeds.split(',')]

    print(f'设备 {DEV} | {describe()}')
    print(f'骨架 {models}')
    print(f'定义 {defs}')
    print(f'种子 {seeds} | epochs {args.epochs} | 基准 {args.bench}')
    print()

    rows, preds = [], []
    # 缓存训练/测试数据（按定义）
    cache = {}
    for ad in defs:
        need_c = any(m in NEEDS_COORDS for m in models)
        t0 = time.time()
        tr = load_train(ad, with_coords=need_c)
        te, meta = load_test(ad, args.bench, with_coords=need_c)
        cache[ad] = (tr, te, meta)
        print(f'  [{ad}] 训练 {len(tr)}, 测试 {len(te)}  ({time.time()-t0:.0f}s)')
    print()

    y = None
    pid = None
    for ad in defs:
        tr, te, meta = cache[ad]
        y = np.array([m['ddg'] for m in meta])
        pid = [m['pdb_id'] for m in meta]
        for model in models:
            for seed in seeds:
                t0 = time.time()
                try:
                    p, n_par, train_t = train_eval(model, tr, te, meta, seed, args.epochs)
                except Exception as e:
                    print(f'  ❌ {model}/{ad}/s{seed}: {type(e).__name__}: {str(e)[:120]}')
                    continue
                r = pearson(y, p)
                rho = spearman(y, p)
                lo, hi = cluster_ci(y, p, pid)
                rows.append({
                    'benchmark': args.bench, 'atom_def': ad, 'model': model, 'seed': seed,
                    'n': len(y), 'n_proteins': len(set(pid)),
                    'r': r, 'ci_low': lo, 'ci_high': hi, 'spearman': rho,
                    'mae': float(np.mean(np.abs(y - p))),
                    'rmse': float(np.sqrt(np.mean((y - p) ** 2))),
                    'n_params': n_par, 'train_seconds': train_t,
                })
                for i in range(len(y)):
                    preds.append({'benchmark': args.bench, 'atom_def': ad, 'model': model,
                                  'seed': seed, 'protein_id': pid[i],
                                  'mutation_id': meta[i]['mut_info'],
                                  'y_true': y[i], 'y_pred': float(p[i])})
                print(f'  {model:<11} {ad:<9} s{seed:<5} r={r:.4f} [{lo:.3f},{hi:.3f}] '
                      f'({time.time()-t0:.0f}s)', flush=True)
                # ★ 增量保存（长任务需要断点）
                pd.DataFrame(rows).to_csv(os.path.join(OUT, f'{args.tag}_results.csv'), index=False)
                pd.DataFrame(preds).to_csv(os.path.join(OUT, f'{args.tag}_predictions.csv'),
                                           index=False)

    rf = os.path.join(OUT, f'{args.tag}_results.csv')
    pf = os.path.join(OUT, f'{args.tag}_predictions.csv')
    pd.DataFrame(rows).to_csv(rf, index=False)
    pd.DataFrame(preds).to_csv(pf, index=False)
    print(f'\n已保存 {rf}')
    print(f'已保存 {pf}')

    # 汇总
    if rows:
        df = pd.DataFrame(rows)
        print('\n=== 各骨架 × 定义 的 r（种子均值）===')
        piv = df.pivot_table(index='model', columns='atom_def', values='r', aggfunc='mean')
        print(piv.round(4).to_string())


if __name__ == '__main__':
    main()
