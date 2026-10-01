"""
gnn_edge_baseline.py — 支持边特征的 GNN 基线
==============================================
发现：GCNConv 完全忽略 edge_attr（接触类型），结构信息的主要载体被浪费。
本脚本用 GINEConv（支持边特征的消息传递）替换 GCNConv，
让模型能区分"疏水接触 / 氢键 / 静电"等不同接触类型。

这是"结构信息真实上限"的进一步探索。
"""
import argparse
import os
import sys

import numpy as np
import torch
import torch.nn.functional as F
from torch_geometric.loader import DataLoader
from torch_geometric.nn import GINEConv, global_mean_pool

sys.path.insert(0, os.path.dirname(__file__))
from paths import DATA, FIGURES, PUB_FIGURES, ROOT, SRC, TOOLS  # 路径集中管理
from train_gnn_baseline import load_s669_samples, load_training_samples, pearson

SRC_PATH = str(SRC)
DATA_PATH = str(DATA)
ROOT_PATH = str(ROOT)
FIG_PATH = str(FIGURES)
PUB_PATH = str(PUB_FIGURES)
TOOLS_PATH = str(TOOLS)


class GNNEdge(torch.nn.Module):
    """2 层 GINEConv（用边特征）+ 全局池化 + 突变位点局部表示 + MLP。"""
    def __init__(self, in_dim, edge_dim, hid=64, out=1):
        super().__init__()
        # GINEConv 需要一个 MLP 把边特征嵌入后与节点特征拼接
        self.edge_mlp1 = torch.nn.Sequential(
            torch.nn.Linear(edge_dim, hid), torch.nn.ReLU(), torch.nn.Linear(hid, hid))
        self.conv1 = GINEConv(torch.nn.Sequential(
            torch.nn.Linear(in_dim, hid), torch.nn.ReLU(), torch.nn.Linear(hid, hid)), edge_dim=hid)
        self.edge_mlp2 = torch.nn.Sequential(
            torch.nn.Linear(edge_dim, hid), torch.nn.ReLU(), torch.nn.Linear(hid, hid))
        self.conv2 = GINEConv(torch.nn.Sequential(
            torch.nn.Linear(hid, hid), torch.nn.ReLU(), torch.nn.Linear(hid, hid)), edge_dim=hid)
        self.fc = torch.nn.Linear(hid * 2, out)

    def forward(self, data):
        edge_attr1 = self.edge_mlp1(data.edge_attr)
        x = self.conv1(data.x, data.edge_index, edge_attr1).relu()
        edge_attr2 = self.edge_mlp2(data.edge_attr)
        x = self.conv2(x, data.edge_index, edge_attr2).relu()
        glob = global_mean_pool(x, data.batch)
        flag = data.x[:, -1]
        local_repr = []
        for i in range(data.num_graphs):
            mask = (data.batch == i) & (flag > 0.5)
            if mask.sum() > 0:
                local_repr.append(x[mask].mean(0))
            else:
                local_repr.append(glob[i])
        local = torch.stack(local_repr)
        feat = torch.cat([glob, local], dim=1)
        return self.fc(feat).squeeze(-1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--training_csv', default=os.path.join(DATA_PATH, 'training_merged_noleak.csv'))
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
    edge_dim = train_samples[0].edge_attr.shape[1]
    model = GNNEdge(in_dim, edge_dim, hid=args.hid)
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
    print(f'\n最终 S669 r = {r:.3f}（边特征 GINEConv 基线）')


if __name__ == '__main__':
    main()
