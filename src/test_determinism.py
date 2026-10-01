"""
test_determinism.py — 验证训练在同一 seed 下是否真正可复现
============================================================
审查指出：原脚本未控制 CUDA/cuDNN/DataLoader generator。
本测试用同一 seed 训练两次，比较逐样本预测是否完全一致。
"""
import os
import sys

import numpy as np
import torch
import torch.nn.functional as F
from torch_geometric.data import Data
from torch_geometric.loader import DataLoader

sys.path.insert(0, os.path.dirname(__file__))
from gnn_local_baseline import GNNLocal
from seed_utils import describe, make_generator, set_seed, worker_init_fn

from pathlib import Path as _P
DATA = str(_P(__file__).resolve().parent.parent / 'data')
SEED = 42
EPOCHS = 3


def load(d):
    df = np.load(d, allow_pickle=True)
    return df['nodes'], df['edge_index'], df['edge_attr']


def build():
    graphs = []
    gdir = os.path.join(DATA, 'contact_graphs_megascale_sc')
    names = sorted(f for f in os.listdir(gdir) if f.endswith('.npz'))[:120]
    for f in names:
        n, e, a = load(os.path.join(gdir, f))
        graphs.append(Data(x=torch.tensor(np.concatenate(
            [n, np.zeros((n.shape[0], 1), dtype=np.float32)], axis=1), dtype=torch.float32),
            edge_index=torch.tensor(e, dtype=torch.long),
            edge_attr=torch.tensor(a, dtype=torch.float32),
            y=torch.tensor([0.5], dtype=torch.float32)))
    return graphs


def train_once(graphs, deterministic=True, device='cpu'):
    set_seed(SEED, deterministic=deterministic)
    m = GNNLocal(graphs[0].x.shape[1], hid=32).to(device)
    opt = torch.optim.Adam(m.parameters(), lr=1e-3)
    dl = DataLoader(graphs, batch_size=16, shuffle=True,
                    generator=make_generator(SEED),
                    worker_init_fn=worker_init_fn(SEED)) if deterministic else \
        DataLoader(graphs, batch_size=16, shuffle=True)
    for _ in range(EPOCHS):
        m.train()
        for b in dl:
            b = b.to(device)
            opt.zero_grad(); F.mse_loss(m(b), b.y).backward(); opt.step()
    m.eval()
    with torch.no_grad():
        return np.concatenate([m(b.to(device)).cpu().numpy()
                               for b in DataLoader(graphs, batch_size=16)])


def main():
    print('环境:', describe())
    graphs = build()
    print(f'测试图数量: {len(graphs)}\n')

    print('--- 测试 A: 同 seed、开启完整确定性、两次运行 ---')
    a1 = train_once(graphs, deterministic=True, device='cpu')
    a2 = train_once(graphs, deterministic=True, device='cpu')
    same = np.array_equal(a1, a2)
    print(f'  逐元素完全相同: {same}')
    if not same:
        print(f'  最大差: {np.abs(a1 - a2).max():.3e}')

    print('\n--- 测试 B: 不开启确定性（模拟原脚本）、两次运行 ---')
    b1 = train_once(graphs, deterministic=False, device='cpu')
    b2 = train_once(graphs, deterministic=False, device='cpu')
    same_b = np.array_equal(b1, b2)
    print(f'  逐元素完全相同: {same_b}')

    print('\n--- 测试 C: CPU vs GPU（若可用） ---')
    if torch.cuda.is_available():
        c1 = train_once(graphs, deterministic=True, device='cuda')
        diff_cpu_gpu = np.abs(a1 - c1).max()
        print(f'  CPU vs GPU 最大差: {diff_cpu_gpu:.3e}')
        c2 = train_once(graphs, deterministic=True, device='cuda')
        print(f'  GPU 两次完全相同: {np.array_equal(c1, c2)}')
    else:
        print('  CUDA 不可用，跳过')

    print('\n' + '=' * 60)
    verdict = []
    if same:
        verdict.append('✅ 开启确定性后，同 seed 结果可完全复现')
    else:
        verdict.append('❌ 开启确定性后仍不可复现 —— 必须排查')
    if not same_b:
        verdict.append('⚠️ 不开启确定性时不可复现 → 说明原脚本存在该风险')
    else:
        verdict.append('ℹ️ 本机不开启确定性也可复现（但不保证其他机器）')
    for v in verdict:
        print(v)


if __name__ == '__main__':
    main()
