"""
add_cb_to_predictions.py — 训练 Cβ 的 GNN 并把逐样本预测追加到 predictions_unified.csv
======================================================================================
Figure 1 之前缺少 Cβ 家族（表中有数据，图中未画）。
本脚本用与其它定义完全相同的协议训练 Cβ GNN（3 seed），
把预测按 figure 脚本期望的格式追加进 predictions_unified.csv。
"""
import os
import sys

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from torch_geometric.data import Data
from torch_geometric.loader import DataLoader

sys.path.insert(0, os.path.dirname(__file__))
from seed_utils import set_seed  # 完整确定性控制（含 cuBLAS/cudnn 标志）

SEED_DEFAULT = 42

from gnn_local_baseline import GNNLocal
from paths import DATA, FIGURES, PUB_FIGURES, ROOT, SRC, TOOLS  # 路径集中管理
from train_gnn_baseline import AA_INDEX

SRC_PATH = str(SRC)
DATA_PATH = str(DATA)
ROOT_PATH = str(ROOT)
FIG_PATH = str(FIGURES)
PUB_PATH = str(PUB_FIGURES)
TOOLS_PATH = str(TOOLS)




DATA = DATA_PATH  # 来自 paths.py，可用 GED_ROOT 环境变量覆盖
AA_ORDER = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
IDX2AA = {i: a for i, a in enumerate(AA_ORDER)}
SEEDS = [42, 123, 2024]


def make_sample(nodes, ei, ea, mi, mt, y):
    nodes = nodes.copy()
    if mt in AA_INDEX:
        nodes[mi, :20] = 0.0
        nodes[mi, AA_INDEX[mt]] = 1.0
    flag = np.zeros((nodes.shape[0], 1), dtype=np.float32)
    flag[mi, 0] = 1.0
    return Data(x=torch.tensor(np.concatenate([nodes, flag], axis=1), dtype=torch.float32),
                edge_index=torch.tensor(ei, dtype=torch.long),
                edge_attr=torch.tensor(ea, dtype=torch.float32),
                y=torch.tensor([float(y)], dtype=torch.float32))


def load_train():
    tr = pd.read_csv(os.path.join(DATA, 'training_merged_noleak_sc.csv'))
    ms = tr[tr['source'] == 'megascale']; tm = tr[tr['source'] == 'thermomutdb']
    n = min(len(ms), len(tm))
    tr = pd.concat([ms.sample(n=n, random_state=42), tm], ignore_index=True)
    dirs = {'megascale': os.path.join(DATA, 'contact_graphs_megascale_cb'),
            'thermomutdb': os.path.join(DATA, 'contact_graphs_thermomutdb_cb')}
    cache, out = {}, []
    for _, r in tr.iterrows():
        p = os.path.join(dirs[r['source']], f"{r['protein']}.npz")
        if not os.path.exists(p):
            continue
        k = f"{r['source']}:{r['protein']}"
        if k not in cache:
            d = np.load(p, allow_pickle=True)
            cache[k] = (d['nodes'], d['edge_index'], d['edge_attr'])
        nd, ei, ea = cache[k]
        mi = int(r['mut_idx'])
        if mi < 0 or mi >= nd.shape[0]:
            continue
        out.append(make_sample(nd, ei, ea, mi, str(r['mt_aa']), r['ddg']))
    return out


def load_test():
    df = pd.read_csv(os.path.join(DATA, 'benchmarks_s669_clean.csv'))
    gdir = os.path.join(DATA, 'contact_graphs_s669_cb')
    cache, out, meta = {}, [], []
    for _, r in df.iterrows():
        if '_node_idx_ok' in df.columns and not bool(r['_node_idx_ok']):
            continue
        pid = r['pdb_id']
        p = os.path.join(gdir, f'{pid}.npz')
        if not os.path.exists(p):
            continue
        if pid not in cache:
            d = np.load(p, allow_pickle=True)
            cache[pid] = (d['nodes'], d['edge_index'], d['edge_attr'])
        nd, ei, ea = cache[pid]
        if '_node_idx' not in df.columns or pd.isna(r['_node_idx']):
            raise ValueError(
                f"缺少已验证的 _node_idx（{(r['pdb_id'], r['mut_info'])}）；"
                "拒绝回退到 _pdb_res_idx，请先运行 fix_index_bug.py")
        mi = int(r['_node_idx'])
        if mi < 0 or mi >= nd.shape[0]:
            continue
        out.append(make_sample(nd, ei, ea, mi, str(r['mut_info'])[-1], r['ddg']))
        meta.append({'pdb_id': pid, 'mut_info': r['mut_info'], 'ddg': r['ddg']})
    return out, meta


def pearson(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return float(np.corrcoef(a, b)[0, 1]) if a.std() > 0 and b.std() > 0 else float('nan')


def main():
    tr = load_train()
    te, meta = load_test()
    print(f'Cβ 训练 {len(tr)}，测试 {len(te)}')
    in_dim = tr[0].x.shape[1]
    y = np.array([m['ddg'] for m in meta])

    rows = []
    preds = []
    for s in SEEDS:
        set_seed(s)
        m = GNNLocal(in_dim, hid=64)
        opt = torch.optim.Adam(m.parameters(), lr=1e-3)
        for _ in range(20):
            m.train()
            for b in DataLoader(tr, batch_size=64, shuffle=True):
                opt.zero_grad(); F.mse_loss(m(b), b.y).backward(); opt.step()
        m.eval()
        pr = []
        with torch.no_grad():
            for b in DataLoader(te, batch_size=64):
                pr.append(m(b).numpy())
        p = np.concatenate(pr)
        preds.append(p)
        print(f'  seed {s}: r = {pearson(y, p):.3f}')
        for i in range(len(y)):
            rows.append({'y_true': y[i], 'y_pred': float(p[i]),
                         'protein_id': meta[i]['pdb_id'], 'mutation_id': meta[i]['mut_info'],
                         'model': f'gnn_local_cb_s{s}', 'seed': s})
    print(f'  3 seed 平均 r = {pearson(y, np.mean(preds, axis=0)):.3f}')

    out_p = os.path.join(DATA, 'predictions_unified.csv')
    base = pd.read_csv(out_p)
    n_before = base['model'].nunique()
    base = base[~base['model'].str.startswith('gnn_local_cb')]   # 幂等：避免重复追加
    combined = pd.concat([base, pd.DataFrame(rows)], ignore_index=True)
    combined.to_csv(out_p, index=False)
    print(f'\npredictions_unified.csv: 模型 {n_before} -> {combined["model"].nunique()} 个')
    print(f'  新增: {sorted(set(combined["model"]) - set(base["model"]))}')


if __name__ == '__main__':
    main()
