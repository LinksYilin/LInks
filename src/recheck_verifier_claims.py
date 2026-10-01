# -*- coding: utf-8 -*-
"""recheck_verifier_claims.py — 复核独立验证者的关键指控"""
import os

import numpy as np
import pandas as pd

D = r'D:\GED_mutation\data'
print('=' * 78)
print('复核独立验证者的指控')
print('=' * 78)

# --- 1) 甘氨酸率：74.6%/6.8% 还是 72.4%/1.8%？ ---
e = pd.read_csv(os.path.join(D, 'edits_corrected.csv'))
s = e[(e.threshold == 8.0) & (e.atom_def == 'cb')]
if 'gly_switch' in s.columns:
    g = s[s.gly_switch.astype(bool)]
    ng = s[~s.gly_switch.astype(bool)]
else:
    g = s[(s.wt_aa == 'G') | (s.mt_aa == 'G')]
    ng = s[(s.wt_aa != 'G') & (s.mt_aa != 'G')]
print(f'\n[1] 甘氨酸（505 对, 8Å, Cβ）:')
print(f'    含 Gly {len(g)} 对, 有编辑率 {100*(g.n_edit>0).mean():.2f}%')
print(f'    不含 Gly {len(ng)} 对, 有编辑率 {100*(ng.n_edit>0).mean():.2f}%')
for t in [6.0, 7.0, 9.0, 10.0]:
    st = e[(e.threshold == t) & (e.atom_def == 'cb')]
    gg = st[(st.wt_aa == 'G') | (st.mt_aa == 'G')] if 'gly_switch' not in st.columns else st[st.gly_switch.astype(bool)]
    nn = st[~st.index.isin(gg.index)]
    print(f'    {t:g}Å: Gly {100*(gg.n_edit>0).mean():.2f}%  非Gly {100*(nn.n_edit>0).mean():.2f}%')
c = pd.read_csv(os.path.join(D, 'ca_cb_edit_diagnosis.csv'))
print(f'    旧诊断文件 ca_cb_edit_diagnosis.csv: {len(c)} 行')

# --- 2) 案例研究距离 ---
loc = pd.read_csv(os.path.join(D, 'locality_corrected.csv'))
print(f'\n[2] 案例研究（locality_corrected.csv）:')
for pid, mut in [('1R2Y', 'R244E'), ('3O39', 'L32P'), ('1XZO', 'W36A')]:
    r = loc[(loc.pdb_id == pid) & (loc.mut_info.astype(str).str.contains(mut[:1]))]
    if len(r):
        print(f'    {pid} {mut}: mean_d_broken = {r.mean_d_broken.iloc[0]:.2f} Å')

# --- 3) 0.362 是 Cβ 还是 centroid？ ---
sm = pd.read_csv(os.path.join(D, 'ladder_seed_summary_audited.csv'))
dg = sm[sm.model == 'deep_gine']
print(f'\n[3] deep_gine 各定义 seed_mean_r:')
for _, r in dg.iterrows():
    print(f'    {r.atom_def:<10} {r.seed_mean_r:.3f}')

# --- 4) 集合嵌套 ---
from itertools import combinations
print(f'\n[4] 集合嵌套关系:')
print(f'    edits_corrected (505) 是否 ⊂ 508/511? 见下方集合差')

# --- 5) 平均 |Δr|：共有 vs 全部 ---
eff = pd.read_csv(os.path.join(D, 'ladder_paired_effects_audited.csv'))
COMMON = ['ca vs centroid', 'cb vs centroid']
print(f'\n[5] 各档平均 |Δr|:')
for m, g in eff.groupby('model'):
    allm = g.delta_r.abs().mean()
    com = g[g.comparison.isin(COMMON)].delta_r.abs().mean()
    print(f'    {m:<12} 全部比较 {allm:.4f}  共有对 {com:.4f}')

# --- 6) 训练集蛋白数 ---
tr = pd.read_csv(os.path.join(D, 'training_merged_noleak_sc.csv'))
ms, tm = tr[tr.source == 'megascale'], tr[tr.source == 'thermomutdb']
n = min(len(ms), len(tm))
v = pd.concat([ms.sample(n=n, random_state=42), tm])
print(f'\n[6] 训练集: 总蛋白 {v.protein.nunique()}, '
      f'MegaScale {v[v.source=="megascale"].protein.nunique()}, '
      f'ThermoMutDB {v[v.source=="thermomutdb"].protein.nunique()}')

# --- 7) ssym ridge ---
inc = pd.read_csv(os.path.join(D, 'increment_decomposition.csv'))
print(f'\n[7] ridge P 值: {inc.columns.tolist()[:6]}')
print(inc.head(3).to_string(index=False))
