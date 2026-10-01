"""
complexity_report.py — 计算复杂度与耗时实测（第 5 项）
========================================================
报告：
  1. 各模型的参数量
  2. Ridge vs GNN 的训练/推理时间（实测）
  3. FoldX 建模每个突变的平均耗时（从文件时间戳估算）
  4. 接触图构建耗时
"""
import os
import sys
import time

import numpy as np
import torch
import torch.nn.functional as F
from torch_geometric.loader import DataLoader

sys.path.insert(0, os.path.dirname(__file__))
from gnn_edge_baseline import GNNEdge
from gnn_local_baseline import GNNLocal
from train_gnn_baseline import GNNRegressor

DATA = r'D:\GED_mutation\data'
OUT = r'D:\GED_mutation\数据与代码说明.md'


def count_params(m):
    return sum(p.numel() for p in m.parameters())


def main():
    lines = []
    print('=== 1. 模型参数量 ===')
    # 用真实特征维度
    import glob
    f = sorted(glob.glob(os.path.join(DATA, 'contact_graphs_megascale_sc', '*.npz')))[0]
    d = np.load(f, allow_pickle=True)
    in_dim = d['nodes'].shape[1] + 1
    edge_dim = d['edge_attr'].shape[1]
    models = {
        'GNN global pooling': GNNRegressor(in_dim, hid=64),
        'GNN local pooling': GNNLocal(in_dim, hid=64),
        'GNN local + BLOSUM': GNNLocal(in_dim + 1, hid=64),
        'GNN edge-aware (GINEConv)': GNNEdge(in_dim, edge_dim, hid=64),
    }
    for name, m in models.items():
        n = count_params(m)
        print(f'  {name:<28} {n:>8,} 参数')
        lines.append(f'| {name} | {n:,} |')
    print(f'  {"Ridge (5 特征)":<28} {6:>8,} 参数')
    lines.append('| Ridge (5 特征) | 6 |')

    print()
    print('=== 2. 训练/推理时间实测（小样本）===')
    # 取 500 个训练样本做基准
    import glob as g2

    from torch_geometric.data import Data
    files = sorted(g2.glob(os.path.join(DATA, 'contact_graphs_thermomutdb_sc', '*.npz')))[:50]
    samples = []
    for fp in files:
        x = np.load(fp, allow_pickle=True)
        nodes = x['nodes']
        flag = np.zeros((nodes.shape[0], 1), dtype=np.float32); flag[0, 0] = 1.0
        samples.append(Data(x=torch.tensor(np.concatenate([nodes, flag], axis=1), dtype=torch.float32),
                            edge_index=torch.tensor(x['edge_index'], dtype=torch.long),
                            edge_attr=torch.tensor(x['edge_attr'], dtype=torch.float32),
                            y=torch.tensor([0.5], dtype=torch.float32)))
    print(f'  基准样本数: {len(samples)}（ThermoMutDB 前 50 个蛋白的图）')
    for name, m in [('GNN local', GNNLocal(in_dim, hid=64))]:
        loader = DataLoader(samples, batch_size=16, shuffle=True)
        opt = torch.optim.Adam(m.parameters(), lr=1e-3)
        t0 = time.time()
        for _ in range(3):
            m.train()
            for b in loader:
                opt.zero_grad(); F.mse_loss(m(b), b.y).backward(); opt.step()
        t1 = time.time()
        m.eval()
        with torch.no_grad():
            t2 = time.time()
            for b in DataLoader(samples, batch_size=16):
                m(b)
            t3 = time.time()
        per_epoch = (t1 - t0) / 3
        print(f'  {name}: 训练 {per_epoch*1000:.0f} ms/epoch（{len(samples)} 图）, '
              f'推理 {(t3-t2)*1000/len(samples):.2f} ms/图')

    print()
    print('=== 3. FoldX 建模耗时（从文件时间戳估算）===')
    mt_dir = os.path.join(DATA, 'mutant_structures_s669')
    stamps = []
    for root, dirs, fs in os.walk(mt_dir):
        for fn in fs:
            if fn.endswith('_1.pdb') and fn.startswith('pdb'):
                stamps.append(os.path.getmtime(os.path.join(root, fn)))
    stamps = sorted(stamps)
    if len(stamps) > 10:
        diffs = np.diff(stamps)
        diffs = diffs[(diffs > 0) & (diffs < 600)]
        print(f'  FoldX 建模突变体总数: {len(stamps)}')
        print(f'  相邻完成的平均间隔: {diffs.mean():.1f} s（中位 {np.median(diffs):.1f} s）')
        print(f'  估算单突变 FoldX 耗时: ≈{np.median(diffs):.1f} s')
        lines.append(f'| FoldX BuildModel 单突变 | ≈{np.median(diffs):.0f} s |')

    print()
    print('=== 4. 接触图构建耗时 ===')
    from contact_graph_defs import compute_all_defs
    fp = os.path.join(DATA, 'structures', 'pdb1bfm.ent')
    t0 = time.time()
    for _ in range(20):
        compute_all_defs(fp, 'A', 0, 69, max_cut=10.0)
    t1 = time.time()
    print(f'  4 种定义一次解析（69 残基）: {(t1-t0)/20*1000:.1f} ms')
    lines.append(f'| 4 种定义接触图（69 残基） | {(t1-t0)/20*1000:.0f} ms |')

    print()
    print('（参数量与耗时已汇总，可写入论文 Methods 的资源说明）')


if __name__ == '__main__':
    main()
