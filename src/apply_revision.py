# -*- coding: utf-8 -*-
"""
apply_revision.py — 一次性应用 WS4（主张限定）+ WS5（图表规范）+ 新数字
=========================================================================
从新结果文件读取数字，改写论文：
  1. 标题加 reliably
  2. 摘要按五主张结构重写（数字自动填入）
  3. Table 1 表头/脚注修正（消除 envelope vs 95% CI 矛盾）
  4. Table 1 的 Cβ 行补全或显式说明
  5. Figure 1 图注与表对应
  6. 新增 "What we do not claim" 段

用法：python apply_revision.py [--dry-run]
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from paths import DATA, ROOT

DOCX = os.path.join(str(ROOT), 'manuscript_polished_final_20260930.docx')
OUT_DOCX = os.path.join(str(ROOT), 'manuscript_revised_v3.docx')


def rd(name):
    p = os.path.join(str(DATA), name)
    if not os.path.exists(p):
        return None
    try:
        df = pd.read_csv(p)
        return df if len(df) else None
    except Exception:
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dry-run', action='store_true')
    args = ap.parse_args()

    print('=== 读取新结果 ===')
    lad = rd('ladder_results.csv')
    eff = rd('ladder_paired_effects.csv')
    fus = rd('esm2_fusion_results.csv')
    zs_s = rd('esm2_zeroshot_s669_esm2_650m.csv')
    zs_y = rd('esm2_zeroshot_ssym_esm2_650m.csv')
    sup = rd('esm2_sup_esm2_650m_results.csv')
    inc = rd('increment_decomposition.csv')

    def show(tag, df):
        print(f'  {tag}: {"缺失" if df is None else f"{len(df)} 行"}')
    for t, d in [('ladder_results', lad), ('ladder_paired_effects', eff),
                 ('esm2_fusion_results', fus), ('zeroshot_s669', zs_s),
                 ('zeroshot_ssym', zs_y), ('esm2_sup_results', sup),
                 ('increment_decomposition', inc)]:
        show(t, d)

    # ---- 汇总关键数字 ----
    print('\n=== 关键数字 ===')
    if lad is not None:
        piv = lad.pivot_table(index='model', columns='atom_def', values='r', aggfunc='mean')
        print('  阶梯 r（骨架 × 定义）:')
        print(piv.round(4).to_string())
        print(f'  参数量: {lad.groupby("model")["n_params"].first().to_dict()}')

    if fus is not None:
        for _, r in fus.iterrows():
            if str(r['model']).startswith('increment'):
                print(f'  ★ {r["benchmark"]}: {r["model"]} Δr={r["r"]:+.4f} '
                      f'[{r.get("ci_low", float("nan")):+.4f},{r.get("ci_high", float("nan")):+.4f}] '
                      f'p={r.get("p", float("nan")):.3f}')

    if inc is not None:
        print('  增量分解:')
        print(inc.round(4).to_string(index=False))

    if args.dry_run:
        print('\n[dry-run] 不修改文档')
        return

    # ---- 实际改写 ----
    from docx import Document
    doc = Document(DOCX)
    paras = doc.paragraphs

    def find(pred):
        for p in paras:
            if pred(p.text.strip()):
                return p
        return None

    def set_all(p, text):
        r = p.runs
        if not r:
            p.add_run(text)
            return
        r[0].text = text
        for x in r[1:]:
            x.text = ''

    # 1) 标题
    p = find(lambda t: t.startswith('Mutation-sensitive residue contact graphs'))
    if p:
        set_all(p, 'Mutation-sensitive residue contact graphs do not reliably improve '
                   'ΔΔG prediction: a controlled benchmark of graph representations '
                   'across encoder capacities')
        print('\n✅ 标题已改（加 reliably + 副标题）')

    print('\n下一步：编号填写需要在所有实验完成后进行')
    doc.save(OUT_DOCX)
    print(f'已保存 {OUT_DOCX}')


if __name__ == '__main__':
    main()
