# -*- coding: utf-8 -*-
"""
esm2_supervised.py — ESM-2 有监督基线（WS2b）
===============================================
流程：
  1. 从图元数据取每个蛋白的序列
  2. 对每个突变构造突变序列，用 ESM-2 取
       - 全局表示：全序列嵌入的均值
       - 局部表示：突变位点的上下文嵌入
     （ESM-2 权重冻结，不反传；嵌入缓存到磁盘）
  3. 在缓存嵌入上训练 MLP 回归头（同一训练集划分）
  4. 在 S669 / ssym 上评估，保存逐样本预测

用法：
  python esm2_supervised.py --stage cache --model esm2_650m
  python esm2_supervised.py --stage train --model esm2_650m
"""
import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd
import torch

sys.path.insert(0, os.path.dirname(__file__))
from paths import DATA
from seed_utils import set_seed

os.environ.setdefault('HF_HOME', str(__import__('pathlib').Path(__file__).resolve().parent.parent / '.hf_cache'))
CACHE = os.path.join(str(DATA), 'esm_emb_cache')
os.makedirs(CACHE, exist_ok=True)

MODELS = {
    'esm2_150m': 'facebook/esm2_t30_150M_UR50D',
    'esm2_650m': 'facebook/esm2_t33_650M_UR50D',
}
SUFFIX = {'centroid': 'sc'}


def seq_map(split):
    """从图目录的 metadata 取 protein -> 序列。"""
    from ladder_common import graph_dir
    d = graph_dir(split, 'centroid')
    out = {}
    if not os.path.isdir(d):
        return out
    for f in os.listdir(d):
        if not f.endswith('.npz'):
            continue
        md = json.loads(np.load(os.path.join(d, f), allow_pickle=True)['metadata'][0])
        if 'seq' in md:
            out[f[:-4]] = md['seq']
    return out


def acid_index(mut_info):
    """'A131G' -> (0-based 位置由调用方给, wt, mt)"""
    return str(mut_info)[0], str(mut_info)[-1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--stage', default='all', choices=['cache', 'train', 'all'])
    ap.add_argument('--model', default='esm2_650m', choices=list(MODELS))
    ap.add_argument('--epochs', type=int, default=40)
    ap.add_argument('--seeds', default='42,123,2024')
    ap.add_argument('--max_train', type=int, default=0)
    args = ap.parse_args()

    dev = 'cuda' if torch.cuda.is_available() else 'cpu'
    emb_file = os.path.join(CACHE, f'emb_{args.model}.npz')
    meta_file = os.path.join(CACHE, f'emb_{args.model}_meta.csv')

    # ---------------- 1) 构建任务清单 ----------------
    tr = pd.read_csv(os.path.join(str(DATA), 'training_merged_noleak_sc.csv'))
    ms = tr[tr['source'] == 'megascale']
    tm = tr[tr['source'] == 'thermomutdb']
    n = min(len(ms), len(tm))
    tr = pd.concat([ms.sample(n=n, random_state=42), tm], ignore_index=True)
    if args.max_train:
        tr = tr.head(args.max_train)

    seqs = {}
    seqs.update({('megascale', k): v for k, v in seq_map('megascale').items()})
    seqs.update({('thermomutdb', k): v for k, v in seq_map('thermomutdb').items()})
    print(f'序列库: megascale {len(seq_map("megascale"))}, thermomutdb {len(seq_map("thermomutdb"))}')

    jobs = []  # (kind, source, protein, pos, wt, mt, y, key)
    for _, r in tr.iterrows():
        key = (r['source'], str(r['protein']))
        s = seqs.get(key)
        if not s:
            continue
        pos = int(r['mut_idx'])
        if not (0 <= pos < len(s)):
            continue
        jobs.append(('train', r['source'], str(r['protein']), pos, s[pos], str(r['mt_aa']),
                     float(r['ddg']), f'{r["source"]}:{r["protein"]}:{pos}:{r["mt_aa"]}'))
    print(f'训练任务 {len(jobs)}')

    from ladder_common import load_test
    for bench in ['s669', 'ssym']:
        te, meta = load_test('centroid', bench)
        for m in meta:
            # 测试序列取自 s669 图 metadata
            pass
    test_jobs = []
    for bench, label in [('s669', 'benchmarks_s669_clean.csv'), ('ssym', 'benchmarks_ssym_clean.csv')]:
        df = pd.read_csv(os.path.join(str(DATA), label))
        if '_node_idx_ok' in df.columns:
            df = df[df['_node_idx_ok']]
        s_map = seq_map(bench)
        for _, r in df.iterrows():
            pid = r['pdb_id']
            s = s_map.get(pid)
            if not s:
                continue
            pos = int(r['_node_idx'])
            if not (0 <= pos < len(s)):
                continue
            wt, mt = acid_index(r['mut_info'])
            if s[pos] != wt:
                # 序列与标注不一致则跳过（并计数）
                continue
            test_jobs.append((bench, pid, pos, s[pos], mt, float(r['ddg']),
                              f'{bench}:{pid}:{r["mut_info"]}'))
    print(f'测试任务 {len(test_jobs)}')

    # ---------------- 2) 缓存嵌入 ----------------
    if args.stage in ('cache', 'all'):
        from transformers import AutoModel, AutoTokenizer
        name = MODELS[args.model]
        print(f'加载 {name} ...')
        tok = AutoTokenizer.from_pretrained(name)
        model = AutoModel.from_pretrained(name).eval()
        if dev == 'cuda':
            model = model.half().cuda()

        all_jobs = jobs + test_jobs
        embs, keep = [], []
        t0 = time.time()
        for i, j in enumerate(all_jobs):
            if j[0] == 'train':
                _, src, prot, pos, wt, mt, y, key = j
            else:
                bench, prot, pos, wt, mt, y, key = j
            # 构造突变序列
            base = seqs.get(('megascale', prot)) or seqs.get(('thermomutdb', prot)) \
                or seq_map('s669').get(prot) or seq_map('ssym').get(prot)
            if base is None:
                continue
            if pos >= len(base):
                continue
            mut_seq = base[:pos] + mt + base[pos + 1:]
            enc = tok(mut_seq, return_tensors='pt')
            enc = {k: v.to(dev) for k, v in enc.items()}
            with torch.no_grad():
                out = model(**enc).last_hidden_state[0].float()
            glob = out[1:-1].mean(0)          # 去掉 <cls>/<eos>
            loc = out[min(pos + 1, out.shape[0] - 2)]
            embs.append(torch.cat([glob, loc]).cpu().numpy().astype(np.float32))
            keep.append(key)
            if (i + 1) % 1000 == 0:
                print(f'  嵌入 {i+1}/{len(all_jobs)}  {time.time()-t0:.0f}s')
        E = np.stack(embs)
        np.savez(emb_file, emb=E, key=np.array(keep))
        print(f'嵌入完成: {E.shape} -> {emb_file}  ({time.time()-t0:.0f}s)')

    # ---------------- 3) 训练回归头 ----------------
    if args.stage in ('train', 'all'):
        z = np.load(emb_file, allow_pickle=True)
        E, keys = z['emb'], list(z['key'])
        idx = {k: i for i, k in enumerate(keys)}

        Xtr, ytr = [], []
        for j in jobs:
            k = j[7]
            if k in idx:
                Xtr.append(E[idx[k]])
                ytr.append(j[6])
        Xtr = torch.tensor(np.stack(Xtr))
        ytr = torch.tensor(np.array(ytr, dtype=np.float32))
        print(f'训练张量 {Xtr.shape}')

        results, preds_out = [], []
        for seed in [int(x) for x in args.seeds.split(',')]:
            set_seed(seed)
            head = torch.nn.Sequential(
                torch.nn.Linear(Xtr.shape[1], 256), torch.nn.ReLU(), torch.nn.Dropout(0.1),
                torch.nn.Linear(256, 64), torch.nn.ReLU(), torch.nn.Linear(64, 1)).to(dev)
            opt = torch.optim.Adam(head.parameters(), lr=1e-3, weight_decay=1e-5)
            Xd, yd = Xtr.to(dev), ytr.to(dev)
            for ep in range(args.epochs):
                head.train()
                perm = torch.randperm(len(Xd), device=dev)
                for b in range(0, len(Xd), 256):
                    sel = perm[b:b + 256]
                    opt.zero_grad()
                    loss = torch.nn.functional.mse_loss(
                        head(Xd[sel]).squeeze(-1), yd[sel])
                    loss.backward()
                    opt.step()
            head.eval()
            for bench in ['s669', 'ssym']:
                Xte, yte, pidv, kk = [], [], [], []
                for j in test_jobs:
                    if j[0] != bench:
                        continue
                    k = j[6]
                    if k in idx:
                        Xte.append(E[idx[k]])
                        yte.append(j[5])
                        pidv.append(j[1])
                        kk.append(k)
                if not Xte:
                    continue
                with torch.no_grad():
                    p = head(torch.tensor(np.stack(Xte)).to(dev)).squeeze(-1).cpu().numpy()
                yv = np.array(yte)
                r = float(np.corrcoef(yv, p)[0, 1])
                results.append({'benchmark': bench, 'model': f'esm2_sup_{args.model}',
                                'seed': seed, 'n': len(yv), 'r': r})
                for i in range(len(yv)):
                    preds_out.append({'benchmark': bench, 'model': f'esm2_sup_{args.model}',
                                      'seed': seed, 'protein_id': pidv[i],
                                      'mutation_id': kk[i].split(':')[-1],
                                      'y_true': yv[i], 'y_pred': float(p[i])})
                print(f'  seed {seed} {bench}: n={len(yv)} r={r:.4f}')
        pd.DataFrame(results).to_csv(os.path.join(str(DATA), f'esm2_sup_{args.model}_results.csv'),
                                     index=False)
        pd.DataFrame(preds_out).to_csv(
            os.path.join(str(DATA), f'esm2_sup_{args.model}_predictions.csv'), index=False)
        print('已保存结果与预测')


if __name__ == '__main__':
    main()
