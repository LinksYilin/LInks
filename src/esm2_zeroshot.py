# -*- coding: utf-8 -*-
"""
esm2_zeroshot.py — ESM-2 零样本 ΔΔG 打分（WS2a）
==================================================
方法：掩码边缘（masked-marginal）
  对每个突变，把 WT 序列中该位点替换为 <mask>，前向一次，
  取该位点的 log-softmax，计算
      score = log P(mt | context) − log P(wt | context)
  该 score 越小表示突变越不利 → 与"正值=失稳"的 ΔΔG 应负相关。

输出 data/esm2_zeroshot.csv：pdb_id, mut_info, ddg, score, model
"""
import argparse
import os
import sys
import time

import numpy as np
import pandas as pd
import torch

sys.path.insert(0, os.path.dirname(__file__))
from paths import DATA

os.environ.setdefault('HF_HOME', str(__import__('pathlib').Path(__file__).resolve().parent.parent / '.hf_cache'))

AA20 = 'ACDEFGHIKLMNPQRSTVWY'
MODELS = {
    'esm2_150m': 'facebook/esm2_t30_150M_UR50D',
    'esm2_650m': 'facebook/esm2_t33_650M_UR50D',
}


def pearson(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if a.std() == 0 or b.std() == 0:
        return float('nan')
    return float(np.corrcoef(a, b)[0, 1])


def spearman(a, b):
    from scipy.stats import spearmanr
    return float(spearmanr(a, b).statistic)


def cluster_ci(y, p, pid, B=1000, seed=0):
    rng = np.random.default_rng(seed)
    prot = np.array(pid)
    uniq, inv = np.unique(prot, return_inverse=True)
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
    ap = argparse.ArgumentParser()
    ap.add_argument('--model', default='esm2_650m', choices=list(MODELS))
    ap.add_argument('--bench', default='s669')
    ap.add_argument('--batch_proteins', type=int, default=1)
    args = ap.parse_args()

    from transformers import AutoModelForMaskedLM, AutoTokenizer

    name = MODELS[args.model]
    dev = 'cuda' if torch.cuda.is_available() else 'cpu'
    print(f'模型 {name} | 设备 {dev}')
    tok = AutoTokenizer.from_pretrained(name)
    model = AutoModelForMaskedLM.from_pretrained(name).eval()
    if dev == 'cuda':
        model = model.half().cuda()
    print('  加载完成')

    mask_id = tok.mask_token_id
    aa_ids = {a: tok.convert_tokens_to_ids(a) for a in AA20}

    label = 'benchmarks_s669_clean.csv' if args.bench == 's669' else 'benchmarks_ssym_clean.csv'
    df = pd.read_csv(os.path.join(str(DATA), label))
    if '_node_idx_ok' in df.columns:
        df = df[df['_node_idx_ok']].reset_index(drop=True)
    print(f'  突变 {len(df)} 条, 蛋白 {df["pdb_id"].nunique()} 个')

    # 按蛋白分组，逐蛋白逐突变打分
    rows = []
    t0 = time.time()
    n_done = 0
    for pid, grp in df.groupby('pdb_id'):
        seq = grp['wt_seq'].iloc[0]
        enc = tok(seq, return_tensors='pt')
        ids = enc['input_ids'].to(dev)
        attn = enc['attention_mask'].to(dev)
        # 第 i 个残基对应 token 位置 i+1（前面有 <cls>）
        for _, r in grp.iterrows():
            pos = int(r['pos'])
            wt_aa, mt_aa = str(r['mut_info'])[0], str(r['mut_info'])[-1]
            if wt_aa not in aa_ids or mt_aa not in aa_ids:
                rows.append({'pdb_id': pid, 'mut_info': r['mut_info'], 'ddg': r['ddg'],
                             'score': np.nan, 'model': args.model})
                continue
            tok_pos = pos + 1
            if tok_pos >= ids.shape[1] - 1:
                rows.append({'pdb_id': pid, 'mut_info': r['mut_info'], 'ddg': r['ddg'],
                             'score': np.nan, 'model': args.model})
                continue
            ids_m = ids.clone()
            ids_m[0, tok_pos] = mask_id
            with torch.no_grad():
                out = model(input_ids=ids_m, attention_mask=attn)
            logp = torch.log_softmax(out.logits[0, tok_pos].float(), dim=-1)
            score = float(logp[aa_ids[mt_aa]] - logp[aa_ids[wt_aa]])
            rows.append({'pdb_id': pid, 'mut_info': r['mut_info'], 'ddg': float(r['ddg']),
                         'score': score, 'model': args.model})
        n_done += 1
        if n_done % 10 == 0:
            print(f'  已处理 {n_done} 个蛋白  ({time.time()-t0:.0f}s, {len(rows)} 条)', flush=True)

    out = pd.DataFrame(rows)
    p = os.path.join(str(DATA), f'esm2_zeroshot_{args.bench}_{args.model}.csv')
    out.to_csv(p, index=False)

    sub = out.dropna(subset=['score'])
    y, s, pidv = sub['ddg'].values, sub['score'].values, sub['pdb_id'].values
    # score 越小越不利 → 与失稳正相关应为 -score
    r_pos = pearson(y, s)
    r_neg = pearson(y, -s)
    print(f'\n=== {args.bench} / {args.model} （n={len(sub)}）===')
    print(f'  r(ΔΔG,  score) = {r_pos:+.4f}')
    print(f'  r(ΔΔG, −score) = {r_neg:+.4f}   ← 若为正，说明 −score 是合理的失稳方向')
    print(f'  Spearman(ΔΔG, −score) = {spearman(y, -s):+.4f}')
    lo, hi = cluster_ci(y, -s, pidv)
    print(f'  95% CI（蛋白簇 bootstrap） = [{lo:+.3f}, {hi:+.3f}]')
    print(f'\n已保存 {p}')


if __name__ == '__main__':
    main()
