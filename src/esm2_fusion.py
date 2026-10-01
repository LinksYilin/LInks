# -*- coding: utf-8 -*-
"""
esm2_fusion.py — ESM-2 + 接触图融合（WS2c，核心增量检验）
===========================================================
回答："在强序列基线之上，结构接触图能否加出增量？"

三个模型在同一划分上训练与评估：
  A) ESM-only      : 冻结的 ESM-2 嵌入 → MLP
  B) Graph-only    : 接触图 GNN（本文骨架）
  C) ESM + Graph   : 双分支拼接 → MLP   ← 关键

比较 C 与 A 的配对 Δr（蛋白簇 bootstrap），给出结构增量的证据。

用法：
  python esm2_fusion.py --model esm2_650m --defs centroid --seeds 42,123,2024
"""
import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from torch_geometric.loader import DataLoader

sys.path.insert(0, os.path.dirname(__file__))
from ladder_common import graph_dir, load_train
from paths import DATA
from seed_utils import set_seed
from strong_backbones import DeepGINE

CACHE = os.path.join(str(DATA), 'esm_emb_cache')
OUT = str(DATA)
DEV = 'cuda' if torch.cuda.is_available() else 'cpu'


def pearson(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 3 or a.std() == 0 or b.std() == 0:
        return float('nan')
    return float(np.corrcoef(a, b)[0, 1])


def clusters(pid):
    prot = np.array(pid)
    uniq, inv = np.unique(prot, return_inverse=True)
    return [np.where(inv == i)[0] for i in range(len(uniq))], len(uniq)


def paired_dr(y, p1, p2, pidx, n_uniq, B=2000, seed=0):
    rng = np.random.default_rng(seed)
    st = []
    for _ in range(B):
        idx = np.concatenate([pidx[i] for i in rng.integers(0, n_uniq, n_uniq)])
        if len(idx) < 5:
            continue
        d = pearson(y[idx], p1[idx]) - pearson(y[idx], p2[idx])
        if not np.isnan(d):
            st.append(d)
    st = np.array(st)
    obs = pearson(y, p1) - pearson(y, p2)
    lo, hi = np.percentile(st, [2.5, 97.5])
    p = 2 * min((st <= 0).mean(), (st >= 0).mean())
    return obs, float(lo), float(hi), float(min(p, 1.0))


# ------------------------------------------------------------------ 模型
class GraphBranch(torch.nn.Module):
    def __init__(self, in_dim, hid=128, out=128):
        super().__init__()
        self.enc = DeepGINE(in_dim, edge_dim=2, hid=hid, n_layers=4, k_hop=2)
        # ★ DeepGINE.embed 返回 (glob | loc)，维度为 2*hid
        self.head = torch.nn.Linear(hid * 2, out)

    def forward(self, g):
        return self.head(self.enc.embed(g))


class Fusion(torch.nn.Module):
    def __init__(self, esm_dim, in_dim, use_esm=True, use_graph=True, hid=256):
        super().__init__()
        self.use_esm, self.use_graph = use_esm, use_graph
        d = 0
        if use_esm:
            self.esm_head = torch.nn.Sequential(
                torch.nn.Linear(esm_dim, hid), torch.nn.ReLU(), torch.nn.Dropout(0.1))
            d += hid
        if use_graph:
            self.gb = GraphBranch(in_dim, hid=128, out=hid)
            d += hid
        self.fc = torch.nn.Sequential(
            torch.nn.Linear(d, hid), torch.nn.ReLU(), torch.nn.Dropout(0.1),
            torch.nn.Linear(hid, 1))

    def forward(self, esm, g):
        parts = []
        if self.use_esm:
            parts.append(self.esm_head(esm))
        if self.use_graph:
            parts.append(self.gb(g))
        return self.fc(torch.cat(parts, dim=1)).squeeze(-1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--model', default='esm2_650m')
    ap.add_argument('--defs', default='centroid')
    ap.add_argument('--seeds', default='42,123,2024')
    ap.add_argument('--epochs', type=int, default=20)
    args = ap.parse_args()

    emb_file = os.path.join(CACHE, f'emb_{args.model}.npz')
    if not os.path.exists(emb_file):
        raise SystemExit(f'缺少嵌入缓存 {emb_file}；请先运行 esm2_supervised.py --stage cache')
    z = np.load(emb_file, allow_pickle=True)
    E, keys = z['emb'], list(z['key'])
    idx = {k: i for i, k in enumerate(keys)}
    print(f'ESM 嵌入 {E.shape}')

    # 训练集图 + 键映射（与嵌入缓存一致）
    tr = pd.read_csv(os.path.join(OUT, 'training_merged_noleak_sc.csv'))
    ms, tm = tr[tr['source'] == 'megascale'], tr[tr['source'] == 'thermomutdb']
    n = min(len(ms), len(tm))
    tr = pd.concat([ms.sample(n=n, random_state=42), tm], ignore_index=True)

    results, preds_out = [], []
    for ad in args.defs.split(','):
        samples = load_train(ad, with_coords=False)
        # load_train 内部顺序与 tr 一致（同样采样），重建键
        tr_keys = []
        for _, r in tr.iterrows():
            tr_keys.append(f'{r["source"]}:{r["protein"]}:{int(r["mut_idx"])}:{r["mt_aa"]}')
        keep = [i for i, k in enumerate(tr_keys) if k in idx and i < len(samples)]
        g_sub = [samples[i] for i in keep]
        e_sub = torch.tensor(np.stack([E[idx[tr_keys[i]]] for i in keep]))
        y_sub = torch.tensor(np.array([float(tr.iloc[i]['ddg']) for i in keep], dtype=np.float32))
        print(f'[{ad}] 图 {len(g_sub)}, ESM {tuple(e_sub.shape)}')

        # 测试集
        test_g, test_e, test_y, test_pid, test_mut = {}, {}, {}, {}, {}
        for bench in ['s669', 'ssym']:
            label = 'benchmarks_s669_clean.csv' if bench == 's669' else 'benchmarks_ssym_clean.csv'
            df = pd.read_csv(os.path.join(OUT, label))
            if '_node_idx_ok' in df.columns:
                df = df[df['_node_idx_ok']]
            gd = graph_dir(bench, ad)
            gs, es, ys, ps, ms_ = [], [], [], [], []
            for _, r in df.iterrows():
                pid = r['pdb_id']
                p = os.path.join(gd, f'{pid}.npz')
                k = f'{bench}:{pid}:{r["mut_info"]}'
                if not os.path.exists(p) or k not in idx:
                    continue
                zz = np.load(p, allow_pickle=True)
                mi = int(r['_node_idx'])
                from ladder_common import make_sample
                gs.append(make_sample(zz['nodes'], zz['edge_index'], zz['edge_attr'],
                                      mi, str(r['mut_info'])[-1], r['ddg']))
                es.append(E[idx[k]])
                ys.append(float(r['ddg']))
                ps.append(pid)
                ms_.append(str(r['mut_info']))
            test_g[bench] = gs
            test_e[bench] = torch.tensor(np.stack(es)) if es else None
            test_y[bench] = np.array(ys)
            test_pid[bench] = ps
            test_mut[bench] = ms_

        in_dim = g_sub[0].x.shape[1]
        esm_dim = e_sub.shape[1]

        for tag, use_esm, use_graph in [('ESM-only', True, False),
                                        ('Graph-only', False, True),
                                        ('ESM+Graph', True, True)]:
            P = {b: [] for b in test_g}
            for seed in [int(x) for x in args.seeds.split(',')]:
                set_seed(seed)
                m = Fusion(esm_dim, in_dim, use_esm, use_graph).to(DEV)
                opt = torch.optim.Adam(m.parameters(), lr=1e-3)
                dl = DataLoader(list(range(len(g_sub))), batch_size=64, shuffle=True)
                for _ in range(args.epochs):
                    m.train()
                    for bidx in dl:
                        bidx = bidx.tolist()
                        gg = [g_sub[i].to(DEV) for i in bidx]
                        from torch_geometric.data import Batch
                        gb = Batch.from_data_list(gg)
                        eb = e_sub[bidx].to(DEV)
                        yb = y_sub[bidx].to(DEV)
                        opt.zero_grad()
                        F.mse_loss(m(eb, gb), yb).backward()
                        opt.step()
                m.eval()
                from torch_geometric.data import Batch
                for bench in test_g:
                    if not test_g[bench]:
                        continue
                    with torch.no_grad():
                        pv = m(test_e[bench].to(DEV),
                               Batch.from_data_list([g.to(DEV) for g in test_g[bench]])).cpu().numpy()
                    P[bench].append(pv)
            for bench in test_g:
                if not P[bench]:
                    continue
                pm = np.mean(P[bench], axis=0)
                yv = test_y[bench]
                r = pearson(yv, pm)
                results.append({'benchmark': bench, 'atom_def': ad, 'model': tag,
                                'n': len(yv), 'r': r})
                for i in range(len(yv)):
                    preds_out.append({'benchmark': bench, 'atom_def': ad, 'model': tag,
                                      'protein_id': test_pid[bench][i],
                                      'mutation_id': test_mut[bench][i],
                                      'y_true': yv[i], 'y_pred': float(pm[i])})
                print(f'  [{ad}] {tag:<11} {bench}: n={len(yv)} r={r:.4f}', flush=True)

        # 配对增量：ESM+Graph vs ESM-only
        for bench in test_g:
            a = [p for p in preds_out if p['atom_def'] == ad and p['benchmark'] == bench
                 and p['model'] == 'ESM+Graph']
            b = [p for p in preds_out if p['atom_def'] == ad and p['benchmark'] == bench
                 and p['model'] == 'ESM-only']
            if not a or not b:
                continue
            ya = np.array([x['y_true'] for x in a])
            pa = np.array([x['y_pred'] for x in a])
            pb = np.array([x['y_pred'] for x in b])
            pidv = [x['protein_id'] for x in a]
            pidx, nu = clusters(pidv)
            obs, lo, hi, pv = paired_dr(ya, pa, pb, pidx, nu)
            print(f'  ★ [{ad}] {bench}: 结构增量 Δr={obs:+.4f} [{lo:+.4f},{hi:+.4f}] p={pv:.3f}')
            results.append({'benchmark': bench, 'atom_def': ad,
                            'model': 'increment(ESM+Graph vs ESM-only)',
                            'n': len(ya), 'r': obs, 'ci_low': lo, 'ci_high': hi, 'p': pv})

    pd.DataFrame(results).to_csv(os.path.join(OUT, 'esm2_fusion_results.csv'), index=False)
    pd.DataFrame(preds_out).to_csv(os.path.join(OUT, 'esm2_fusion_predictions.csv'), index=False)
    print('\n已保存 esm2_fusion_results.csv / esm2_fusion_predictions.csv')


if __name__ == '__main__':
    main()
