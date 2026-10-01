# -*- coding: utf-8 -*-
"""
test_egnn_stability.py — EGNN 稳定性验证（改进 vs 原版）
==========================================================
问题：原 EGNN 在多数定义下出现"好/坏双峰"（如 centroid: 0.393 vs −0.081），
      这使最高能力档无法作为"强编码器"的证据。
改进：坐标更新乘以可学习 scale，初始为 0 → 初始坐标完全固定，逐步学习位移。
本脚本用多个种子对比两种设置的稳定性。
"""
import argparse
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
from strong_backbones import DeepGINE, EGNN

DEV = 'cuda' if torch.cuda.is_available() else 'cpu'


def pearson(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return float(np.corrcoef(a, b)[0, 1]) if a.std() > 0 and b.std() > 0 else float('nan')


def run(tag, factory, tr, te, y, seeds, epochs=20, lr=1e-3):
    rs = []
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
        r = pearson(y, p)
        rs.append(r)
        print(f'    seed {s}: r={r:+.4f}  ({time.time()-t0:.0f}s)', flush=True)
    rs = np.array(rs)
    print(f'  {tag}: mean={rs.mean():+.4f}  sd={rs.std():.4f}  '
          f'min={rs.min():+.4f}  max={rs.max():+.4f}  '
          f'崩塌(<-0.05)={int((rs < -0.05).sum())}/{len(rs)}')
    return rs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--atom_def', default='centroid')
    ap.add_argument('--seeds', default='42,123,2024,7,99')
    ap.add_argument('--epochs', type=int, default=20)
    args = ap.parse_args()
    seeds = [int(x) for x in args.seeds.split(',')]

    tr = load_train(args.atom_def, with_coords=True)
    te, meta = load_test(args.atom_def, 's669', with_coords=True)
    y = np.array([m['ddg'] for m in meta])
    in_dim = tr[0].x.shape[1]
    print(f'定义 {args.atom_def}: 训练 {len(tr)}, 测试 {len(te)}, 种子 {seeds}\n', flush=True)

    print('=== 原版 EGNN（坐标更新无 scale）===')
    run('EGNN 原版', lambda: EGNN(in_dim, hid=128, n_layers=4, update_coords=True, scaled=False),
        tr, te, y, seeds, args.epochs)

    print('\n=== 改进版 EGNN（coord_scale 初始 0）===')
    run('EGNN 改进', lambda: EGNN(in_dim, hid=128, n_layers=4, update_coords=True, scaled=True),
        tr, te, y, seeds, args.epochs)

    print('\n=== 对照：深层 GINE（无坐标）===')
    if os.environ.get('SKIP_DEEPGINE') == '1':
        print('  （已跳过）')
    else:
        run('DeepGINE', lambda: DeepGINE(in_dim, hid=128, n_layers=6, k_hop=2),
            tr, te, y, seeds, args.epochs)


if __name__ == '__main__':
    main()
