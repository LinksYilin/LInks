"""
gnn_baseline.py — 容量匹配的 GNN 对照基线（消融变体 B）
========================================================
纯接触图 GNN 回归：只吃 WT 接触图 + 突变位点标记，不做编辑建模。
这是 GEDMut 的对照（证明 GED 编辑建模的增量）。

架构：2-3 层 GCN/GAT + 全局池化 + MLP 回归头。
输入：接触图（节点特征 + 边）+ 突变位点（把突变残基的节点特征增强一个 flag）。

用法：
  python gnn_baseline.py \
    --graph_dir data/contact_graphs_s669 \
    --label_csv data/benchmarks_s669_clean.csv \
    --epochs 100 --seed 42
"""
import argparse
import json
import os

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from torch_geometric.data import Data
from torch_geometric.loader import DataLoader
from torch_geometric.nn import GCNConv, global_mean_pool


def load_graphs(graph_dir, label_csv):
    """加载接触图 + 标签，返回 (samples, protein_ids)。
    protein_ids 用于按蛋白分层划分（防止同蛋白突变跨 train/val）。
    """
    df = pd.read_csv(label_csv)
    graphs = {}
    for pid in df['pdb_id'].unique():
        p = os.path.join(graph_dir, f'{pid}.npz')
        if not os.path.exists(p):
            continue
        d = np.load(p)
        graphs[pid] = {
            'nodes': d['nodes'],
            'edge_index': d['edge_index'],
            'edge_attr': d['edge_attr'],
        }

    samples = []
    protein_ids = []
    for _, row in df.iterrows():
        pid = row['pdb_id']
        if pid not in graphs:
            continue
        g = graphs[pid]
        mut_idx = int(row['_pdb_res_idx']) if '_pdb_res_idx' in row and not pd.isna(row['_pdb_res_idx']) else None
        nodes = g['nodes'].copy()
        flag = np.zeros((nodes.shape[0], 1), dtype=np.float32)
        if mut_idx is not None and 0 <= mut_idx < nodes.shape[0]:
            flag[mut_idx, 0] = 1.0
        nodes = np.concatenate([nodes, flag], axis=1)

        x = torch.tensor(nodes, dtype=torch.float32)
        ei = torch.tensor(g['edge_index'], dtype=torch.long)
        ea = torch.tensor(g['edge_attr'], dtype=torch.float32)
        y = torch.tensor([float(row['ddg'])], dtype=torch.float32)
        data = Data(x=x, edge_index=ei, edge_attr=ea, y=y)
        samples.append(data)
        protein_ids.append(pid)
    return samples, protein_ids


class GNNRegressor(torch.nn.Module):
    """2 层 GCN + 全局均值池化 + MLP 回归头。"""
    def __init__(self, in_dim, hid=64, out=1):
        super().__init__()
        self.conv1 = GCNConv(in_dim, hid)
        self.conv2 = GCNConv(hid, hid)
        self.fc = torch.nn.Linear(hid, out)

    def forward(self, data):
        x = self.conv1(data.x, data.edge_index).relu()
        x = self.conv2(x, data.edge_index).relu()
        x = global_mean_pool(x, data.batch)  # 用 PyG 自动生成的 batch 索引
        return self.fc(x).squeeze(-1)


def pearson(y_true, y_pred):
    t = y_true - y_true.mean()
    p = y_pred - y_pred.mean()
    return (t * p).sum() / (t.norm() * p.norm() + 1e-12)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--graph_dir', required=True)
    ap.add_argument('--label_csv', required=True)
    ap.add_argument('--epochs', type=int, default=100)
    ap.add_argument('--seed', type=int, default=42)
    ap.add_argument('--hid', type=int, default=64)
    ap.add_argument('--n_folds', type=int, default=5)
    args = ap.parse_args()

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)

    samples, protein_ids = load_graphs(args.graph_dir, args.label_csv)
    print(f'加载 {len(samples)} 个样本, {len(set(protein_ids))} 个蛋白')

    # 按蛋白分层的 GroupKFold
    from sklearn.model_selection import GroupKFold
    gkf = GroupKFold(n_splits=args.n_folds)
    idx = np.arange(len(samples))
    groups = np.array(protein_ids)

    in_dim = samples[0].x.shape[1]
    fold_rs = []
    for fold, (train_idx, val_idx) in enumerate(gkf.split(idx, groups=groups)):
        train = [samples[i] for i in train_idx]
        val = [samples[i] for i in val_idx]
        train_loader = DataLoader(train, batch_size=32, shuffle=True)
        val_loader = DataLoader(val, batch_size=32)

        model = GNNRegressor(in_dim, hid=args.hid)
        opt = torch.optim.Adam(model.parameters(), lr=1e-3)
        for epoch in range(args.epochs):
            model.train()
            for batch in train_loader:
                opt.zero_grad()
                pred = model(batch)
                loss = F.mse_loss(pred, batch.y)
                loss.backward()
                opt.step()
        model.eval()
        preds, trues = [], []
        with torch.no_grad():
            for batch in val_loader:
                preds.append(model(batch))
                trues.append(batch.y)
        pred = torch.cat(preds)
        true = torch.cat(trues)
        r = pearson(true, pred).item()
        fold_rs.append(r)
        print(f'  Fold {fold}: val r={r:.3f} (val蛋白数={len(set(groups[val_idx]))})')

    print(f'\n各折 r: {[f"{r:.3f}" for r in fold_rs]}')
    print(f'平均 r: {np.mean(fold_rs):.3f} ± {np.std(fold_rs):.3f}')


if __name__ == '__main__':
    main()
