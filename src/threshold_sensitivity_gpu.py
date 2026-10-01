# -*- coding: utf-8 -*-
"""
threshold_sensitivity_gpu.py — 阈值敏感性（GPU 版）
====================================================
原 threshold_sensitivity_prediction.py 在 CPU 上训练（30 次），需数小时。
本版：GPU + 增量保存 + 与主阶梯一致的 GNNLocal 配置。
"""
import os
import sys
import time

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from torch_geometric.loader import DataLoader

sys.path.insert(0, os.path.dirname(__file__))
from paths import DATA
from seed_utils import set_seed
from gnn_local_baseline import GNNLocal

DEV = 'cuda' if torch.cuda.is_available() else 'cpu'
OUT = os.path.join(str(DATA), 'threshold_sensitivity_prediction.csv')
THRESHOLDS = [6.0, 7.0, 8.0, 9.0, 10.0]
SEEDS = [42, 123]


def suffix(th):
    return 'sc' if abs(th - 8.0) < 1e-6 else f'centroid{th:g}A'


def pearson(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return float(np.corrcoef(a, b)[0, 1]) if a.std() > 0 and b.std() > 0 else float('nan')


def cluster_ci(y, p, pid, B=1000, seed=0):
    rng = np.random.default_rng(seed)
    uniq, inv = np.unique(np.array(pid), return_inverse=True)
    pidx = [np.where(inv == i)[0] for i in range(len(uniq))]
    st = []
    for _ in range(B):
        idx = np.concatenate([pidx[i] for i in rng.integers(0, len(uniq), len(uniq))])
        if len(idx) >= 3:
            v = pearson(y[idx], p[idx])
            if not np.isnan(v):
                st.append(v)
    return (float(np.percentile(st, 2.5)), float(np.percentile(st, 97.5))) if st else (np.nan, np.nan)


def main():
    from ladder_common import load_train, load_test

    rows = []
    for bench, prefix in [('S669', 's669'), ('ssym', 'ssym')]:
        for th in THRESHOLDS:
            t0 = time.time()
            tr = load_train(prefix if False else None, None) if False else None
            # 训练集按阈值取对应目录
            tr = _load_train_th(th)
            te, meta = _load_test_th(th, prefix)
            if not tr or not te:
                print(f'{bench}/{th}: 数据缺失', flush=True)
                continue
            y = np.array([m['ddg'] for m in meta])
            pid = [m['pdb_id'] for m in meta]
            in_dim = tr[0].x.shape[1]
            preds = []
            for s in SEEDS:
                set_seed(s)
                m_ = GNNLocal(in_dim, hid=64).to(DEV)
                opt = torch.optim.Adam(m_.parameters(), lr=1e-3)
                for _ in range(20):
                    m_.train()
                    for b in DataLoader(tr, batch_size=64, shuffle=True):
                        b = b.to(DEV)
                        opt.zero_grad()
                        F.mse_loss(m_(b), b.y).backward()
                        opt.step()
                m_.eval()
                with torch.no_grad():
                    pr = [m_(b.to(DEV)).cpu().numpy()
                          for b in DataLoader(te, batch_size=128)]
                preds.append(np.concatenate(pr))
            pm = np.mean(preds, axis=0)
            r = pearson(y, pm)
            lo, hi = cluster_ci(y, pm, pid)
            rows.append({'benchmark': bench, 'threshold': th, 'n': len(y),
                         'mean_edges': float(np.mean([s.edge_index.shape[1] // 2 for s in te])),
                         'r': r, 'ci_low': lo, 'ci_high': hi})
            pd.DataFrame(rows).to_csv(OUT, index=False)   # 增量保存
            print(f'  {bench} {th:g} Å: n={len(y)}, r={r:.3f} [{lo:.3f},{hi:.3f}] '
                  f'({time.time()-t0:.0f}s)', flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(OUT, index=False)
    print(f'\n已保存 {OUT}')


def _load_train_th(th):
    """按阈值加载训练集。样本构造复用 ladder_common.make_sample，
    确保突变氨基酸 one-hot 与突变标志列与主流程完全一致。"""
    from ladder_common import make_sample

    main_csv = os.path.join(str(DATA), 'training_merged_noleak_sc.csv')
    tr = pd.read_csv(main_csv)
    ms, tm = tr[tr['source'] == 'megascale'], tr[tr['source'] == 'thermomutdb']
    n = min(len(ms), len(tm))
    tr = pd.concat([ms.sample(n=n, random_state=42), tm], ignore_index=True)
    sfx = suffix(th)
    dirs = {'megascale': os.path.join(str(DATA), f'contact_graphs_megascale_{sfx}'),
            'thermomutdb': os.path.join(str(DATA), f'contact_graphs_thermomutdb_{sfx}')}
    out = []
    for _, r in tr.iterrows():
        src, prot = r['source'], str(r['protein'])
        p = os.path.join(dirs[src], f'{prot}.npz')
        if not os.path.exists(p):
            continue
        z = np.load(p, allow_pickle=True)
        mi = int(r['mut_idx'])
        if not (0 <= mi < z['nodes'].shape[0]):
            continue
        out.append(make_sample(z['nodes'], z['edge_index'], z['edge_attr'],
                               mi, str(r['mt_aa']), float(r['ddg'])))
    return out


def _load_test_th(th, prefix):
    """按阈值加载测试集，同样复用 make_sample（突变氨基酸 + 标志列）。"""
    from ladder_common import make_sample

    sfx = suffix(th)
    gd = os.path.join(str(DATA), f'contact_graphs_{prefix}_{sfx}')
    label = f'benchmarks_{prefix}_clean.csv'
    df = pd.read_csv(os.path.join(str(DATA), label))
    if '_node_idx_ok' in df.columns:
        df = df[df['_node_idx_ok']]
    out, meta = [], []
    for _, r in df.iterrows():
        pid = r['pdb_id']
        p = os.path.join(gd, f'{pid}.npz')
        if not os.path.exists(p):
            continue
        z = np.load(p, allow_pickle=True)
        mi = int(r['_node_idx'])
        if not (0 <= mi < z['nodes'].shape[0]):
            continue
        out.append(make_sample(z['nodes'], z['edge_index'], z['edge_attr'],
                               mi, str(r['mut_info'])[-1], float(r['ddg'])))
        meta.append({'pdb_id': pid, 'mut_info': r['mut_info'], 'ddg': float(r['ddg'])})
    return out, meta


if __name__ == '__main__':
    main()
