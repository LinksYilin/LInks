# -*- coding: utf-8 -*-
"""
fix_table1_structure.py — 修复 Table 1 的统计标注规范（WS5）
==============================================================
问题（已实测确认）：
  1. 表头写 "Pearson r [95% CI]"，但 GNN 行的括号实际是
     "三个种子各自 95% CI 的**包络**"（已验证：包络 = [min(ci_low), max(ci_high)]
     = [0.204, 0.488]，与 Table 1 完全一致）。
     包络**不是任何量的 95% CI**，且把优化变异与抽样不确定性混在一起。
  2. Cβ 行格式与其它行断裂（无括号、无 ρ/MAE/RMSE、无 ssym、脚注无解释）。

修复：
  - 表头明确写 "(seed mean) [95% CI]"，并注明括号为包络
  - 脚注完整说明三种区间的区别
  - Cβ 行补 "—" 的说明
"""
import argparse
import os
import shutil
import sys

sys.path.insert(0, os.path.dirname(__file__))
from paths import ROOT

DOCX = os.path.join(str(ROOT), 'manuscript_polished_final_20260930.docx')
BACKUP = os.path.join(str(ROOT), 'manuscript_polished_final_20260930_preWS5.docx')

HDR = {
    'S669 Pearson r [95% CI]': 'S669 Pearson r (seed mean) [seed-CI envelope]',
    'ssym Pearson r [95% CI]': 'ssym Pearson r (seed mean) [seed-CI envelope]',
}

FOOTNOTE = (
    'Table 1 | Predictive performance on the common 511-mutation S669 intersection and on ssym. '
    'MAE and RMSE are in kcal mol⁻¹. All methods were evaluated on the same mutation intersection. '
    'For Ridge, brackets give the 95% protein-cluster bootstrap confidence interval (1,000 resamples '
    'over proteins). For GNN families, the three seeds are optimisation replicates rather than '
    'independent samples: the reported interval is the **envelope** of the three per-seed 95% '
    'confidence intervals (that is, the minimum lower bound and the maximum upper bound across '
    'seeds), which is deliberately conservative and is NOT itself a 95% confidence interval for any '
    'single quantity. Point estimates and their intervals for individual seeds are plotted in '
    'Figure 1a. The Cβ family (row 3) was evaluated on S669 only; its Spearman, MAE and RMSE are not '
    'reported because the Cβ comparison exists to test representational sensitivity to side-chain '
    'repacking rather than to establish a performance baseline, and the corresponding ssym contact '
    'graphs were not constructed for that definition.'
)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--apply', action='store_true')
    args = ap.parse_args()

    if not args.apply:
        print('[dry-run] 将要做的修改：')
        for a, b in HDR.items():
            print(f'  表头: "{a}"\n     -> "{b}"')
        print(f'  脚注改为:\n    {FOOTNOTE[:200]}...')
        return

    if not os.path.exists(BACKUP):
        shutil.copy2(DOCX, BACKUP)
        print(f'已备份 -> {os.path.basename(BACKUP)}')

    from docx import Document
    doc = Document(DOCX)
    n_hdr = 0
    n_note = 0
    for t in doc.tables:
        for row in t.rows:
            for c in row.cells:
                for p in c.paragraphs:
                    old = p.text.strip()
                    if old in HDR:
                        if p.runs:
                            p.runs[0].text = HDR[old]
                            for r in p.runs[1:]:
                                r.text = ''
                        n_hdr += 1
    for p in doc.paragraphs:
        if p.text.strip().startswith('Table 1 |'):
            if p.runs:
                p.runs[0].text = FOOTNOTE
                for r in p.runs[1:]:
                    r.text = ''
            n_note += 1
    doc.save(DOCX)
    print(f'表头替换 {n_hdr} 处, 脚注替换 {n_note} 处')


if __name__ == '__main__':
    main()
