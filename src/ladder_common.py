# -*- coding: utf-8 -*-
"""
ladder_common.py — 能力阶梯实验的公共数据加载与模型工厂
==========================================================
统一 4 种图定义的数据加载（含 coords），并提供骨架工厂。
"""
import os
import sys

import numpy as np
import pandas as pd
import torch
from torch_geometric.data import Data

sys.path.insert(0, os.path.dirname(__file__))
from paths import DATA
from seed_utils import set_seed  # noqa: F401  （完整确定性控制）
from train_gnn_baseline import AA_INDEX
from gnn_local_baseline import GNNLocal
from strong_backbones import DeepGINE, EGNN

DATA_PATH = str(DATA)

# 定义 → 图目录后缀（与既有约定一致）
SUFFIX = {'ca': 'ca', 'cb': 'cb', 'centroid': 'sc', 'allatom': 'allatom'}
DEFS = ['ca', 'cb', 'centroid', 'allatom']


def graph_dir(split, ad, cutoff=None):
    """split: megascale / thermomutdb / s669 / ssym"""
    sfx = SUFFIX[ad]
    extra = '' if cutoff in (None, 8.0) else f'{cutoff:g}A'
    if ad == 'centroid':
        # centroid 的阈值图目录名为 _centroid{th}A（历史原因）
        return os.path.join(DATA_PATH, f'contact_graphs_{split}_centroid{extra}' if extra
                            else f'contact_graphs_{split}_sc')
    return os.path.join(DATA_PATH, f'contact_graphs_{split}_{sfx}{extra}')


def make_sample(nodes, ei, ea, mi, mt_aa, y, coords=None):
    nodes = nodes.copy()
    if mt_aa in AA_INDEX:
        nodes[mi, :20] = 0.0
        nodes[mi, AA_INDEX[mt_aa]] = 1.0
    flag = np.zeros((nodes.shape[0], 1), dtype=np.float32)
    flag[mi, 0] = 1.0
    x = torch.tensor(np.concatenate([nodes, flag], axis=1), dtype=torch.float32)
    kw = {}
    if coords is not None:
        kw['coords'] = torch.tensor(np.asarray(coords, dtype=np.float32))
    return Data(x=x, edge_index=torch.tensor(ei, dtype=torch.long),
                edge_attr=torch.tensor(ea, dtype=torch.float32),
                y=torch.tensor([float(y)], dtype=torch.float32), **kw)


def _load_npz(path, with_coords):
    z = np.load(path, allow_pickle=True)
    if with_coords:
        if 'coords' not in z.files:
            raise RuntimeError(f'{path} 缺 coords；EGNN 需要坐标')
        return z['nodes'], z['edge_index'], z['edge_attr'], z['coords']
    return z['nodes'], z['edge_index'], z['edge_attr'], None


def load_train(ad, cutoff=None, with_coords=False, seed=42):
    """MegaScale + ThermoMutDB 按 1:1 采样，返回 PyG Data 列表。"""
    tr = pd.read_csv(os.path.join(DATA_PATH, 'training_merged_noleak_sc.csv'))
    ms = tr[tr['source'] == 'megascale']
    tm = tr[tr['source'] == 'thermomutdb']
    n = min(len(ms), len(tm))
    tr = pd.concat([ms.sample(n=n, random_state=seed), tm], ignore_index=True)
    dirs = {'megascale': graph_dir('megascale', ad, cutoff),
            'thermomutdb': graph_dir('thermomutdb', ad, cutoff)}
    cache, out = {}, []
    for _, r in tr.iterrows():
        p = os.path.join(dirs[r['source']], f"{r['protein']}.npz")
        if not os.path.exists(p):
            continue
        k = f"{r['source']}:{r['protein']}"
        if k not in cache:
            cache[k] = _load_npz(p, with_coords)
        nd, ei, ea, xy = cache[k]
        mi = int(r['mut_idx'])
        if mi < 0 or mi >= nd.shape[0]:
            continue
        out.append(make_sample(nd, ei, ea, mi, str(r['mt_aa']), r['ddg'], xy))
    return out


def load_test(ad, bench='s669', cutoff=None, with_coords=False):
    """返回 (samples, meta)；meta 含 y_true / pdb_id / mut_info。"""
    label = 'benchmarks_s669_clean.csv' if bench == 's669' else 'benchmarks_ssym_clean.csv'
    df = pd.read_csv(os.path.join(DATA_PATH, label))
    gdir = graph_dir(bench, ad, cutoff)
    cache, out, meta = {}, [], []
    for _, r in df.iterrows():
        if '_node_idx_ok' in df.columns and not bool(r['_node_idx_ok']):
            continue
        pid = r['pdb_id']
        p = os.path.join(gdir, f'{pid}.npz')
        if not os.path.exists(p):
            continue
        if pid not in cache:
            cache[pid] = _load_npz(p, with_coords)
        nd, ei, ea, xy = cache[pid]
        if '_node_idx' not in r.index or pd.isna(r['_node_idx']):
            raise ValueError(f'缺少 _node_idx: {pid} {r["mut_info"]}')
        mi = int(r['_node_idx'])
        if mi < 0 or mi >= nd.shape[0]:
            continue
        out.append(make_sample(nd, ei, ea, mi, str(r['mut_info'])[-1], r['ddg'], xy))
        meta.append({'pdb_id': pid, 'mut_info': r['mut_info'], 'ddg': float(r['ddg'])})
    return out, meta


# --------------------------------------------------------------------- 模型
def build_model(name, in_dim, **kw):
    if name == 'gnn_global':
        from gnn_baseline import GNNBaseline  # noqa: F401
        return GNNLocal(in_dim, hid=64)          # 兼容旧名，实际由 runner 区分
    if name == 'gnn_local':
        return GNNLocal(in_dim, hid=64)
    if name == 'deep_gine':
        return DeepGINE(in_dim, edge_dim=2, hid=128, n_layers=6, k_hop=2, **kw)
    if name == 'egnn':
        return EGNN(in_dim, edge_dim=2, hid=128, n_layers=4, k_hop=2, **kw)
    raise ValueError(name)
