"""
train_gnn_baseline.py — GNN 基线（统一训练集训练 + S669 评估）
================================================================
正确的数据流：
  训练：统一训练集（MegaScale + ThermoMutDB，见 merge_training.py 产出）
  最终评估：S669（一次性，不调参）

容量匹配的纯 GNN 回归（消融变体 B），作为 GEDMut 的对照。
"""
import argparse
import os

import numpy as np

from seed_utils import set_seed  # 完整确定性控制（含 cuBLAS/cudnn 标志）

SEED_DEFAULT = 42

import pandas as pd
import torch
import torch.nn.functional as F
from torch_geometric.data import Data
from torch_geometric.loader import DataLoader
from torch_geometric.nn import GCNConv, global_mean_pool

from paths import DATA, FIGURES, PUB_FIGURES, ROOT, SRC, TOOLS  # 路径集中管理

SRC_PATH = str(SRC)
DATA_PATH = str(DATA)
ROOT_PATH = str(ROOT)
FIG_PATH = str(FIGURES)
PUB_PATH = str(PUB_FIGURES)
TOOLS_PATH = str(TOOLS)




# 氨基酸 → one-hot 索引（与 build_contact_graph 的 AA_ORDER 一致）
AA_ORDER = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
AA_INDEX = {aa: i for i, aa in enumerate(AA_ORDER)}


def load_training_samples(training_csv, graph_dirs, max_samples=None, balance=True, tm_frac=0.5):
    """加载统一训练集样本。

    training_csv: merge_training.py 产出，列 = [protein, mut_idx, ddg, source]
    graph_dirs: dict {source: graph_dir}
    balance: True 时按 source 平衡采样（tm_frac 控制 ThermoMutDB 比例）
    """
    df = pd.read_csv(training_csv)
    if balance:
        ms = df[df['source'] == 'megascale']
        tm = df[df['source'] == 'thermomutdb']
        n_tm = len(tm)
        n_ms = int(n_tm / tm_frac * (1 - tm_frac)) if tm_frac < 1 else 0
        n_ms = min(n_ms, len(ms))
        ms = ms.sample(n=n_ms, random_state=42)
        df = pd.concat([ms, tm], ignore_index=True)
        print(f'  平衡采样: MegaScale {len(ms)} + ThermoMutDB {len(tm)} = {len(df)}')

    graph_cache = {}
    samples = []
    protein_ids = []
    for _, row in df.iterrows():
        protein = row['protein']
        source = row['source']
        gdir = graph_dirs.get(source)
        if gdir is None:
            continue
        p = os.path.join(gdir, f'{protein}.npz')
        if not os.path.exists(p):
            continue
        cache_key = f'{source}:{protein}'
        if cache_key not in graph_cache:
            d = np.load(p, allow_pickle=True)
            graph_cache[cache_key] = (d['nodes'], d['edge_index'], d['edge_attr'])
        nodes, ei, ea = graph_cache[cache_key]
        mut_idx = int(row['mut_idx'])
        if mut_idx < 0 or mut_idx >= nodes.shape[0]:
            continue
        # 关键修复：把突变位点的 one-hot 改为突变后氨基酸（mt_aa），
        # 让模型感知"这个位点变成了什么"（而非只加 0/1 flag）
        nodes = nodes.copy()
        mt_aa = str(row.get('mt_aa', ''))
        if mt_aa in AA_INDEX:
            nodes[mut_idx, :20] = 0.0
            nodes[mut_idx, AA_INDEX[mt_aa]] = 1.0
        # 追加突变位点标记（保留 flag，作为"这是突变位点"的定位信息）
        flag = np.zeros((nodes.shape[0], 1), dtype=np.float32)
        flag[mut_idx, 0] = 1.0
        x = torch.tensor(np.concatenate([nodes, flag], axis=1), dtype=torch.float32)
        edge_index = torch.tensor(ei, dtype=torch.long)
        edge_attr = torch.tensor(ea, dtype=torch.float32)
        y = torch.tensor([float(row['ddg'])], dtype=torch.float32)
        samples.append(Data(x=x, edge_index=edge_index, edge_attr=edge_attr, y=y))
        protein_ids.append(cache_key)
        if max_samples and len(samples) >= max_samples:
            break
    return samples, protein_ids


def load_s669_samples(graph_dir, label_csv):
    """加载 S669 评估样本（同样的构造方式：突变位点 one-hot 改为 mt_aa）。"""
    df = pd.read_csv(label_csv)
    graph_cache = {}
    samples = []
    for _, row in df.iterrows():
        pid = row['pdb_id']
        p = os.path.join(graph_dir, f'{pid}.npz')
        if not os.path.exists(p):
            continue
        if pid not in graph_cache:
            d = np.load(p, allow_pickle=True)
            graph_cache[pid] = (d['nodes'], d['edge_index'], d['edge_attr'])
        nodes, ei, ea = graph_cache[pid]
        # ★ 索引修正：_pdb_res_idx 是整链索引，图节点是切片后索引 → 优先用 _node_idx
        if '_node_idx' not in row.index or pd.isna(row['_node_idx']):
            raise ValueError(
                f"缺少已验证的 _node_idx（{(row['pdb_id'], row['mut_info'])}）；"
                "拒绝回退到有问题的 _pdb_res_idx，请先运行 fix_index_bug.py")
        mut_idx = int(row['_node_idx'])
        if mut_idx < 0 or mut_idx >= nodes.shape[0]:
            continue
        # 突变位点 one-hot 改为 mt_aa（mut_info 最后一个字符是 mt_aa）
        nodes = nodes.copy()
        mt_aa = str(row['mut_info'])[-1]
        if mt_aa in AA_INDEX:
            nodes[mut_idx, :20] = 0.0
            nodes[mut_idx, AA_INDEX[mt_aa]] = 1.0
        flag = np.zeros((nodes.shape[0], 1), dtype=np.float32)
        flag[mut_idx, 0] = 1.0
        x = torch.tensor(np.concatenate([nodes, flag], axis=1), dtype=torch.float32)
        edge_index = torch.tensor(ei, dtype=torch.long)
        edge_attr = torch.tensor(ea, dtype=torch.float32)
        y = torch.tensor([float(row['ddg'])], dtype=torch.float32)
        samples.append(Data(x=x, edge_index=edge_index, edge_attr=edge_attr, y=y))
    return samples


class GNNRegressor(torch.nn.Module):
    def __init__(self, in_dim, hid=64, out=1):
        super().__init__()
        self.conv1 = GCNConv(in_dim, hid)
        self.conv2 = GCNConv(hid, hid)
        self.fc = torch.nn.Linear(hid, out)

    def forward(self, data):
        x = self.conv1(data.x, data.edge_index).relu()
        x = self.conv2(x, data.edge_index).relu()
        x = global_mean_pool(x, data.batch)
        return self.fc(x).squeeze(-1)


def pearson(y_true, y_pred):
    t = y_true - y_true.mean()
    p = y_pred - y_pred.mean()
    return (t * p).sum() / (t.norm() * p.norm() + 1e-12)


def evaluate(model, loader):
    model.eval()
    preds, trues = [], []
    with torch.no_grad():
        for batch in loader:
            preds.append(model(batch))
            trues.append(batch.y)
    if not preds:
        return float('nan')
    pred = torch.cat(preds)
    true = torch.cat(trues)
    return pearson(true, pred).item()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--training_csv', default=r'D:\GED_mutation\data\training_merged.csv')
    ap.add_argument('--ms_graph_dir', default=r'D:\GED_mutation\data\contact_graphs_megascale')
    ap.add_argument('--tm_graph_dir', default=r'D:\GED_mutation\data\contact_graphs_thermomutdb')
    ap.add_argument('--s669_graph_dir', default=r'D:\GED_mutation\data\contact_graphs_s669')
    ap.add_argument('--s669_label', default=r'D:\GED_mutation\data\benchmarks_s669_clean.csv')
    ap.add_argument('--epochs', type=int, default=5)
    ap.add_argument('--hid', type=int, default=64)
    ap.add_argument('--max_samples', type=int, default=0, help='0=不限制')
    ap.add_argument('--seed', type=int, default=42)
    args = ap.parse_args()

    set_seed(args.seed)

    graph_dirs = {'megascale': args.ms_graph_dir, 'thermomutdb': args.tm_graph_dir}
    print('加载统一训练集...')
    train_samples, _ = load_training_samples(
        args.training_csv, graph_dirs,
        max_samples=args.max_samples if args.max_samples > 0 else None)
    print(f'  训练样本: {len(train_samples)}')

    print('加载 S669 评估样本...')
    s669_samples = load_s669_samples(args.s669_graph_dir, args.s669_label)
    print(f'  S669 样本: {len(s669_samples)}')

    train_loader = DataLoader(train_samples, batch_size=64, shuffle=True)
    s669_loader = DataLoader(s669_samples, batch_size=64)

    in_dim = train_samples[0].x.shape[1]
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
        if (epoch + 1) % 2 == 0 or epoch == args.epochs - 1:
            r = evaluate(model, s669_loader)
            print(f'Epoch {epoch+1}: S669 r={r:.3f}')

    r = evaluate(model, s669_loader)
    print(f'\n最终 S669 Pearson r = {r:.3f}（统一训练集）')


if __name__ == '__main__':
    main()
