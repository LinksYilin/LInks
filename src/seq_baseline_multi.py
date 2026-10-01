# -*- coding: utf-8 -*-
"""
seq_baseline_multi.py — 多序列模型零样本基线（稳健性检验）
=============================================================
目的：检验"结构无增量"的结论是否只在 ESM-2 上成立。
     换用不同预训练语料的序列模型，若结论一致则更稳。

支持模型（按预训练语料区分）：
  esm2_150m  facebook/esm2_t30_150M_UR50D   UniRef50 (2022)
  esm2_650m  facebook/esm2_t33_650M_UR50D   UniRef50 (2022)
  esm1b      facebook/esm1b_t33_650M_UR50S  UniRef50 (2019)
  protbert   Rostlab/prot_bert               UniRef100 (BERT)

输出 data/seq_zeroshot_{bench}_{model}.csv（含 score 与 ddg）
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
    'esm1b': 'facebook/esm1b_t33_650M_UR50S',
    'protbert': 'Rostlab/prot_bert',
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
    ap = argparse.ArgumentParser()
    ap.add_argument('--model', default='esm1b', choices=list(MODELS))
    ap.add_argument('--bench', default='s669')
    args = ap.parse_args()

    from transformers import AutoModelForMaskedLM, AutoTokenizer

    name = MODELS[args.model]
    dev = 'cuda' if torch.cuda.is_available() else 'cpu'
    print(f'模型 {name} | 设备 {dev}', flush=True)

    # ProtBert 用空格分隔的序列
    spaced = args.model == 'protbert'
    tok = AutoTokenizer.from_pretrained(name)
    model = AutoModelForMaskedLM.from_pretrained(name).eval()
    if dev == 'cuda':
        model = model.half().cuda()
    print('  加载完成', flush=True)

    mask_id = tok.mask_token_id
    aa_ids = {a: tok.convert_tokens_to_ids(a) for a in AA20}
    if any(v is None or v < 0 for v in aa_ids.values()):
        # ProtBert 的氨基酸 token 是单字母
        aa_ids = {a: tok.convert_tokens_to_ids(a) for a in AA20}
    print(f'  mask_id={mask_id}, 示例 token id: {list(aa_ids.items())[:4]}', flush=True)

    label = 'benchmarks_s669_clean.csv' if args.bench == 's669' else 'benchmarks_ssym_clean.csv'
    df = pd.read_csv(os.path.join(str(DATA), label))
    if '_node_idx_ok' in df.columns:
        df = df[df['_node_idx_ok']].reset_index(drop=True)
    print(f'  突变 {len(df)}', flush=True)

    rows, t0, nd = [], time.time(), 0
    for pid, grp in df.groupby('pdb_id'):
        seq = grp['wt_seq'].iloc[0]
        text = ' '.join(seq) if spaced else seq
        enc = tok(text, return_tensors='pt')
        ids = enc['input_ids'].to(dev)
        attn = enc.get('attention_mask')
        attn = attn.to(dev) if attn is not None else None
        # ProtBert 有空格时 token 位置 = 1 + 2*pos
        step = 2 if spaced else 1
        for _, r in grp.iterrows():
            pos = int(r['pos'])
            wt, mt = str(r['mut_info'])[0], str(r['mut_info'])[-1]
            tok_pos = 1 + step * pos
            if (wt not in aa_ids or mt not in aa_ids
                    or tok_pos >= ids.shape[1] - 1):
                rows.append({'pdb_id': pid, 'mut_info': r['mut_info'], 'ddg': float(r['ddg']),
                             'score': np.nan, 'model': args.model})
                continue
            ids_m = ids.clone()
            ids_m[0, tok_pos] = mask_id
            with torch.no_grad():
                out = model(input_ids=ids_m, attention_mask=attn)
            logp = torch.log_softmax(out.logits[0, tok_pos].float(), dim=-1)
            rows.append({'pdb_id': pid, 'mut_info': r['mut_info'], 'ddg': float(r['ddg']),
                         'score': float(logp[aa_ids[mt]] - logp[aa_ids[wt]]),
                         'model': args.model})
        nd += 1
        if nd % 20 == 0:
            print(f'    {nd} 蛋白 ({time.time()-t0:.0f}s)', flush=True)

    out = pd.DataFrame(rows)
    p = os.path.join(str(DATA), f'seq_zeroshot_{args.bench}_{args.model}.csv')
    out.to_csv(p, index=False)

    sub = out.dropna(subset=['score'])
    y, s, pidv = sub['ddg'].values, sub['score'].values, sub['pdb_id'].values
    print(f'\n=== {args.bench} / {args.model} (n={len(sub)}) ===')
    print(f'  r(ΔΔG, −score) = {pearson(y, -s):+.4f}')
    print(f'  Spearman       = {spearman(y, -s):+.4f}')
    lo, hi = cluster_ci(y, -s, pidv)
    print(f'  95% CI         = [{lo:+.3f}, {hi:+.3f}]')
    print(f'已保存 {p}')


if __name__ == '__main__':
    main()
