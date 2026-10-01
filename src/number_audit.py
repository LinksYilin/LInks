# -*- coding: utf-8 -*-
"""
number_audit.py — 全文数字审计：把论文中的每个数值与结果文件核对
====================================================================
从稿件提取所有形如 r = X / Δr = X / n = X 的数值，在结果 CSV 中查找匹配，
报告「能对上」「对不上」「无处可查」三类。
"""
import os
import re
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from paths import DATA, ROOT

DOCX = os.path.join(str(ROOT), 'manuscript_routeA.docx')
D = str(DATA)


def load_sources():
    """收集所有结果文件中的数值（带来源标签）。"""
    pool = defaultdict(list)   # 值(3位小数字符串) -> [来源]

    def add(v, src):
        try:
            f = float(v)
        except (TypeError, ValueError):
            return
        if np.isnan(f):
            return
        for k in {f'{f:.3f}', f'{f:.4f}', f'{abs(f):.3f}'}:
            pool[k].append(src)

    for name in os.listdir(D):
        if not name.endswith('.csv'):
            continue
        p = os.path.join(D, name)
        try:
            df = pd.read_csv(p)
        except Exception:
            continue
        for c in df.columns:
            if df[c].dtype.kind in 'fi':
                for v in df[c].dropna().unique()[:400]:
                    add(v, f'{name}:{c}')
    # 已知常量
    for v, s in [(0.0, 'const'), (1.0, 'const'), (342, 'ssym n'), (511, 's669 n'),
                 (505, 'qc pairs'), (508, 'esm n'), (543, 's669 raw'),
                 (7905, 'training'), (538, 'foldx pairs'), (88, 'proteins'),
                 (15, 'ssym proteins'), (5825, 'gnn_global params'),
                 (5889, 'gnn_local params'), (50050, 'gnn_edge params'),
                 (460034, 'deep_gine params'), (508934, 'egnn params'),
                 (2000, 'bootstrap B'), (1000, 'bootstrap B'), (20, 'epochs'),
                 (64, 'batch'), (905, 'x'), (907, 'x')]:
        add(v, s)
    return pool


def main():
    from docx import Document
    pool = load_sources()
    doc = Document(DOCX)
    text = '\n'.join(p.text for p in doc.paragraphs)
    for t in doc.tables:
        for r in t.rows:
            text += '\n' + ' | '.join(c.text for c in r.cells)

    # 提取形如 "= 0.123" / "r = 0.123" / "Δr = −0.123" 的数值
    pats = [r'r\s*=\s*[\u2212-]?(\d\.\d{3,4})',
            r'\u0394r\s*=\s*[\u2212-]?(\d\.\d{3,4})',
            r'=\s*[\u2212-]?(\d\.\d{3,4})']
    found = defaultdict(set)
    for pat in pats:
        for m in re.finditer(pat, text):
            found[m.group(1)].add(m.group(0).strip()[:24])

    ok, bad = [], []
    for val, ctxs in sorted(found.items()):
        if val in pool:
            ok.append((val, pool[val][:2], list(ctxs)[:1]))
        else:
            # 尝试 ±0.001 容差
            hit = None
            f = float(val)
            for k in pool:
                try:
                    if abs(float(k) - f) <= 0.0011:
                        hit = k
                        break
                except ValueError:
                    continue
            if hit:
                ok.append((val, pool[hit][:2], list(ctxs)[:1]))
            else:
                bad.append((val, list(ctxs)[:2]))

    print('=' * 74)
    print(f'论文中提取到 {len(found)} 个不同数值')
    print(f'  可追溯到结果文件: {len(ok)}')
    print(f'  未找到来源:       {len(bad)}')
    print('=' * 74)
    if bad:
        print('\n未找到来源的数值（需人工确认）：')
        for val, ctxs in bad:
            print(f'  {val:<8} 出现于: {ctxs}')

    print('\n可追溯数值（样例 25 个）：')
    for val, src, ctx in ok[:25]:
        print(f'  {val:<8} <- {src[0]:<52} (论文: {ctx[0] if ctx else ""})')


if __name__ == '__main__':
    main()
