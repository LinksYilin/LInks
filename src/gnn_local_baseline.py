"""
gnn_local_baseline.py — 更强的 GNN 基线（突变位点局部池化）
============================================================
相比 gnn_baseline 的全局 mean pooling（会稀释突变位点信号），
本基线显式提取突变位点的节点表示，concat 到全局表示上，
让模型能真正利用突变位点附近的结构上下文。

这是"结构信息真实上限"的确认实验。
"""
import argparse
import os
import sys

import numpy as np
import torch
import torch.nn.functional as F
from torch_geometric.loader import DataLoader
from torch_geometric.nn import GCNConv, global_mean_pool

sys.path.insert(0, os.path.dirname(__file__))
from paths import DATA, FIGURES, PUB_FIGURES, ROOT, SRC, TOOLS  # 路径集中管理
from train_gnn_baseline import load_s669_samples, load_training_samples, pearson

SRC_PATH = str(SRC)
DATA_PATH = str(DATA)
ROOT_PATH = str(ROOT)
FIG_PATH = str(FIGURES)
PUB_PATH = str(PUB_FIGURES)
TOOLS_PATH = str(TOOLS)


class GNNLocal(torch.nn.Module):
    """2 层 GCN + 全局池化 + 突变位点局部表示 + MLP。"""
    def __init__(self, in_dim, hid=64, out=1):
        super().__init__()
        self.conv1 = GCNConv(in_dim, hid)
        self.conv2 = GCNConv(hid, hid)
        # 输入 = 全局表示(hid) + 突变位点表示(hid) + 突变位点原始理化差(hid 用不到，先不加)
        self.fc = torch.nn.Linear(hid * 2, out)

    def forward(self, data):
        x = self.conv1(data.x, data.edge_index).relu()
        x = self.conv2(x, data.edge_index).relu()
        # 全局表示
        glob = global_mean_pool(x, data.batch)
        # 突变位点表示：最后一列是 flag（标记突变位点），用 mask 取出
        flag = data.x[:, -1]  # 突变位点 flag（0/1）
        # 每个图取突变位点节点的表示
        # 用 scatter 方式：对每个图，找 flag==1 的节点
        local_repr = []
        for i in range(data.num_graphs):
            node_mask = (data.batch == i) & (flag > 0.5)
            if node_mask.sum() > 0:
                local_repr.append(x[node_mask].mean(0))  # 突变位点（若有多个取平均）
            else:
                local_repr.append(glob[i])
        local = torch.stack(local_repr)
        feat = torch.cat([glob, local], dim=1)
        return self.fc(feat).squeeze(-1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--training_csv', default=os.path.join(DATA_PATH, 'training_merged.csv'))
    ap.add_argument('--ms_graph_dir', default=os.path.join(DATA_PATH, 'contact_graphs_megascale'))
    ap.add_argument('--tm_graph_dir', default=os.path.join(DATA_PATH, 'contact_graphs_thermomutdb'))
    ap.add_argument('--s669_graph_dir', default=os.path.join(DATA_PATH, 'contact_graphs_s669'))
    ap.add_argument('--s669_label', default=os.path.join(DATA_PATH, 'benchmarks_s669_clean.csv'))
    ap.add_argument('--epochs', type=int, default=20)
    ap.add_argument('--hid', type=int, default=64)
    ap.add_argument('--max_samples', type=int, default=0)
    ap.add_argument('--seed', type=int, default=42)
    args = ap.parse_args()

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)

    graph_dirs = {'megascale': args.ms_graph_dir, 'thermomutdb': args.tm_graph_dir}
    train_samples, _ = load_training_samples(
        args.training_csv, graph_dirs,
        max_samples=args.max_samples if args.max_samples > 0 else None)
    print(f'训练样本: {len(train_samples)}')

    s669_samples = load_s669_samples(args.s669_graph_dir, args.s669_label)
    print(f'S669 样本: {len(s669_samples)}')

    train_loader = DataLoader(train_samples, batch_size=64, shuffle=True)
    s669_loader = DataLoader(s669_samples, batch_size=64)

    in_dim = train_samples[0].x.shape[1]
    model = GNNLocal(in_dim, hid=args.hid)
    opt = torch.optim.Adam(model.parameters(), lr=1e-3)

    for epoch in range(args.epochs):
        model.train()
        for batch in train_loader:
            opt.zero_grad()
            pred = model(batch)
            loss = F.mse_loss(pred, batch.y)
            loss.backward()
            opt.step()
        if (epoch + 1) % 4 == 0 or epoch == args.epochs - 1:
            model.eval()
            preds, trues = [], []
            with torch.no_grad():
                for batch in s669_loader:
                    preds.append(model(batch))
                    trues.append(batch.y)
            r = pearson(torch.cat(trues), torch.cat(preds)).item()
            print(f'Epoch {epoch+1}: S669 r={r:.3f}')

    model.eval()
    preds, trues = [], []
    with torch.no_grad():
        for batch in s669_loader:
            preds.append(model(batch))
            trues.append(batch.y)
    r = pearson(torch.cat(trues), torch.cat(preds)).item()
    print(f'\n最终 S669 r = {r:.3f}（突变位点局部池化基线）')


if __name__ == '__main__':
    main()
