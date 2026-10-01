# -*- coding: utf-8 -*-
"""
esm2_fusion_fast.py — ESM-2 + 接触图融合（优化版，WS2c）
=========================================================
相比 esm2_fusion.py 的关键改动：
  1. 把 ESM 嵌入直接存入 PyG Data，用 PyG 原生 DataLoader 批处理
     （原版每个 batch 手工 Batch.from_data_list，CPU 成为瓶颈）
  2. 每完成一个 (模型, 种子) 立即写盘 → 可随时查看部分结果
  3. 评估集批量前向一次算完

三个模型：ESM-only / Graph-only / ESM+Graph
核心产出：ESM+Graph 相对 ESM-only 的配对 Δr（蛋白簇 bootstrap）
"""
import argparse
import os
import sys
import time

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from torch_geometric.data import Data, DataLoader
from torch_geometric.loader import DataLoader as PyGLoader

sys.path.insert(0, os.path.dirname(__file__))
from ladder_common import graph_dir, make_sample
from paths import DATA
from seed_utils import set_seed
from strong_backbones import DeepGINE

CACHE = os.path.join(str(DATA), 'esm_emb_cache')
OUT = str(DATA)
DEV = 'cuda' if torch.cuda.is_available() else 'cpu'
AA20 = 'ACDEFGHIKLMNPQRSTVWY'


def pearson(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 3 or a.std() == 0 or b.std() == 0:
        return float('nan')
    return float(np.corrcoef(a, b)[0, 1])


def paired_dr(y, p1, p2, pid, B=2000, seed=0):
    rng = np.random.default_rng(seed)
    prot = np.array(pid)
    uniq, inv = np.unique(prot, return_inverse=True)
    pidx = [np.where(inv == i)[0] for i in range(len(uniq))]
    st = []
    for _ in range(B):
        idx = np.concatenate([pidx[i] for i in rng.integers(0, len(uniq), len(uniq))])
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


class GraphBranch(torch.nn.Module):
    def __init__(self, in_dim, hid=128, out=128):
        super().__init__()
        self.enc = DeepGINE(in_dim, edge_dim=2, hid=hid, n_layers=4, k_hop=2)
        self.head = torch.nn.Linear(hid * 2, out)

    def forward(self, g, esm=None):
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

    def forward(self, data):
        parts = []
        if self.use_esm:
            parts.append(self.esm_head(data.esm))
        if self.use_graph:
            parts.append(self.gb(data))
        return self.fc(torch.cat(parts, dim=1)).squeeze(-1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--model', default='esm2_650m')
    ap.add_argument('--defs', default='centroid')
    ap.add_argument('--seeds', default='42,123,2024')
    ap.add_argument('--epochs', type=int, default=20)
    args = ap.parse_args()

    z = np.load(os.path.join(CACHE, f'emb_{args.model}.npz'), allow_pickle=True)
    E_all = z['emb']
    key2i = {k: i for i, k in enumerate(list(z['key']))}
    print(f'ESM 嵌入 {E_all.shape}', flush=True)

    tr = pd.read_csv(os.path.join(OUT, 'training_merged_noleak_sc.csv'))
    ms, tm = tr[tr['source'] == 'megascale'], tr[tr['source'] == 'thermomutdb']
    n = min(len(ms), len(tm))
    tr = pd.concat([ms.sample(n=n, random_state=42), tm], ignore_index=True)

    res_f, pred_f = os.path.join(OUT, 'esm2_fusion_results.csv'), \
        os.path.join(OUT, 'esm2_fusion_predictions.csv')
    results, preds_out = [], []

    for ad in args.defs.split(','):
        # ---------- 训练集 ----------
        trs = []
        for _, r in tr.iterrows():
            src, prot = r['source'], str(r['protein'])
            p = os.path.join(graph_dir(src, ad), f'{prot}.npz')
            k = f'{src}:{prot}:{int(r["mut_idx"])}:{r["mt_aa"]}'
            if not os.path.exists(p) or k not in key2i:
                continue
            zz = np.load(p, allow_pickle=True)
            mi = int(r['mut_idx'])
            if not (0 <= mi < zz['nodes'].shape[0]):
                continue
            d = make_sample(zz['nodes'], zz['edge_index'], zz['edge_attr'], mi,
                            str(r['mt_aa']), r['ddg'])
            d.esm = torch.tensor(E_all[key2i[k]]).unsqueeze(0)
            trs.append(d)
        in_dim, esm_dim = trs[0].x.shape[1], trs[0].esm.shape[1]
        print(f'[{ad}] 训练 {len(trs)}, in_dim={in_dim}, esm_dim={esm_dim}', flush=True)

        # ---------- 测试集 ----------
        tests = {}
        for bench in ['s669', 'ssym']:
            label = f'benchmarks_{bench}_clean.csv'
            df = pd.read_csv(os.path.join(OUT, label))
            if '_node_idx_ok' in df.columns:
                df = df[df['_node_idx_ok']]
            gd = graph_dir(bench, ad)
            xs, ys, ps, ms_ = [], [], [], []
            for _, r in df.iterrows():
                pid = r['pdb_id']
                p = os.path.join(gd, f'{pid}.npz')
                k = f'{bench}:{pid}:{r["mut_info"]}'
                if not os.path.exists(p) or k not in key2i:
                    continue
                zz = np.load(p, allow_pickle=True)
                mi = int(r['_node_idx'])
                if not (0 <= mi < zz['nodes'].shape[0]):
                    continue
                d = make_sample(zz['nodes'], zz['edge_index'], zz['edge_attr'], mi,
                                str(r['mut_info'])[-1], r['ddg'])
                d.esm = torch.tensor(E_all[key2i[k]]).unsqueeze(0)
                xs.append(d)
                ys.append(float(r['ddg']))
                ps.append(pid)
                ms_.append(str(r['mut_info']))
            tests[bench] = (xs, np.array(ys), ps, ms_)
            print(f'  [{ad}] {bench}: {len(xs)}', flush=True)

        for tag, use_esm, use_graph in [('ESM-only', True, False),
                                        ('Graph-only', False, True),
                                        ('ESM+Graph', True, True)]:
            P = {b: [] for b in tests}
            for seed in [int(x) for x in args.seeds.split(',')]:
                set_seed(seed)
                m = Fusion(esm_dim, in_dim, use_esm, use_graph).to(DEV)
                opt = torch.optim.Adam(m.parameters(), lr=1e-3)
                dl = PyGLoader(trs, batch_size=64, shuffle=True)
                t0 = time.time()
                for _ in range(args.epochs):
                    m.train()
                    for b in dl:
                        b = b.to(DEV)
                        opt.zero_grad()
                        F.mse_loss(m(b), b.y).backward()
                        opt.step()
                m.eval()
                for bench, (xs, ys, ps, ms_) in tests.items():
                    if not xs:
                        continue
                    with torch.no_grad():
                        p = np.concatenate([m(b.to(DEV)).cpu().numpy()
                                            for b in PyGLoader(xs, batch_size=128)])
                    P[bench].append(p)
                print(f'  [{ad}] {tag:<11} seed{seed} 完成 ({time.time()-t0:.0f}s)', flush=True)

            for bench, (xs, ys, ps, ms_) in tests.items():
                if not P[bench]:
                    continue
                pm = np.mean(P[bench], axis=0)
                r = pearson(ys, pm)
                results.append({'benchmark': bench, 'atom_def': ad, 'model': tag,
                                'n': len(ys), 'r': r})
                for i in range(len(ys)):
                    preds_out.append({'benchmark': bench, 'atom_def': ad, 'model': tag,
                                      'protein_id': ps[i], 'mutation_id': ms_[i],
                                      'y_true': ys[i], 'y_pred': float(pm[i])})
                print(f'  [{ad}] {tag:<11} {bench}: n={len(ys)} r={r:.4f}', flush=True)
            # ★ 增量保存
            pd.DataFrame(results).to_csv(res_f, index=False)
            pd.DataFrame(preds_out).to_csv(pred_f, index=False)

        # 增量：ESM+Graph vs ESM-only
        for bench, (xs, ys, ps, ms_) in tests.items():
            a = [x for x in preds_out if x['atom_def'] == ad and x['benchmark'] == bench
                 and x['model'] == 'ESM+Graph']
            b = [x for x in preds_out if x['atom_def'] == ad and x['benchmark'] == bench
                 and x['model'] == 'ESM-only']
            if not a or not b:
                continue
            ya = np.array([x['y_true'] for x in a])
            pa = np.array([x['y_pred'] for x in a])
            pb = np.array([x['y_pred'] for x in b])
            obs, lo, hi, pv = paired_dr(ya, pa, pb, [x['protein_id'] for x in a])
            print(f'  ★ [{ad}] {bench}: 结构增量 Δr={obs:+.4f} [{lo:+.4f},{hi:+.4f}] p={pv:.3f}',
                  flush=True)
            results.append({'benchmark': bench, 'atom_def': ad,
                            'model': 'increment(ESM+Graph vs ESM-only)',
                            'n': len(ya), 'r': obs, 'ci_low': lo, 'ci_high': hi, 'p': pv})
        pd.DataFrame(results).to_csv(res_f, index=False)
        pd.DataFrame(preds_out).to_csv(pred_f, index=False)

    print('\n完成', flush=True)


if __name__ == '__main__':
    main()
