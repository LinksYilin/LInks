# -*- coding: utf-8 -*-
"""
consistency_sweep.py — 全文术语与数字一致性扫描（nature-polishing 要求）
==========================================================================
按 nature-polishing 的 Terminology Ledger 原则，检查同一概念是否用了多种写法，
以及数字/单位/符号是否漂移。
"""
import os
import re
from collections import Counter

from docx import Document

DOC = r'D:\GED_mutation\manuscript_routeA.docx'

# (规范形式, [可接受变体], 变体正则)
TERMS = [
    ('side-chain-centroid', [r'side-chain-centroid', r'side-chain centroid',
                             r'sidechain centroid', r'centroid graph']),
    ('\u0394\u0394G', [r'\u0394\u0394G', r'ddG', r'\u0394G\u2020', r'delta-delta G']),
    ('S669', [r'\bS669\b', r'\bs669\b', r'S-669']),
    ('ssym', [r'\bssym\b', r'\bSsym\b', r'\bSSYM\b', r's-ssym']),
    ('contact graph', [r'contact graph', r'contact-graph', r'contactgraphs']),
    ('wild-type', [r'wild-type', r'wild type', r'wildtype']),
    ('FoldX', [r'\bFoldX\b', r'\bFoldx\b', r'\bfoldx\b']),
    ('SCWRL4', [r'\bSCWRL4\b', r'\bScwrl4\b']),
    ('protein-cluster bootstrap', [r'protein-cluster bootstrap', r'protein cluster bootstrap',
                                   r'cluster bootstrap']),
    ('Pearson r', [r'Pearson r', r'Pearson correlation', r'pearson r']),
    ('ESM-2', [r'\bESM-2\b', r'\bESM2\b', r'\besm-2\b']),
    ('capacity ladder', [r'capacity ladder', r'capability ladder', r'encoder ladder']),
]

UNITS = [
    ('kcal mol', [r'kcal mol\u207b\u00b9', r'kcal/mol', r'kcal mol-1', r'kcal\u00b7mol']),
    ('\u00c5', [r'\u00c5', r'\bAngstrom', r'\bangstrom']),
]

NUMS = [
    ('0.907 / 90.7', [r'90\.7\s*%', r'0\.907']),
    ('0.099 / 9.9', [r'9\.9\s*%', r'0\.099']),
    ('0.816 / 81.6', [r'81\.6\s*%', r'0\.816']),
    ('511', [r'\b511\b', r'\b512\b']),
    ('342', [r'\b342\b']),
    ('7,905', [r'7,905', r'7905']),
    ('505', [r'\b505\b']),
    ('538', [r'\b538\b']),
]


def main():
    d = Document(DOC)
    paras = [p.text for p in d.paragraphs]
    full = '\n'.join(paras)

    print('=' * 74)
    print('术语一致性（同一概念的不同写法）')
    print('=' * 74)
    issues = 0
    for canon, variants in TERMS:
        counts = []
        for v in variants:
            n = len(re.findall(v, full))
            if n:
                counts.append((v, n))
        if len(counts) > 1:
            print(f'\n⚠ {canon}: {len(counts)} 种写法')
            for v, n in counts:
                print(f'    {v:<34} {n:>4} 次')
            issues += 1
        else:
            print(f'  ✓ {canon:<28} {counts[0][1] if counts else 0:>4} 次（单一写法）')

    print()
    print('=' * 74)
    print('单位与数字')
    print('=' * 74)
    for canon, variants in UNITS + NUMS:
        counts = [(v, len(re.findall(v, full))) for v in variants]
        counts = [(v, n) for v, n in counts if n]
        if len(counts) > 1:
            print(f'  ⚠ {canon}: ' + ', '.join(f'{v}={n}' for v, n in counts))
            issues += 1
        elif counts:
            print(f'  ✓ {canon:<20} {counts[0][0]} ({counts[0][1]} 次)')

    # 图/表引用完整性
    print()
    print('=' * 74)
    print('图/表引用')
    print('=' * 74)
    caps = re.findall(r'(Figure \d+) \|', full)
    refs = sorted(set(re.findall(r'Figure \d+[a-c]?', full)))
    cap_nums = sorted(set(caps), key=lambda x: int(x.split()[1]))
    print(f'  图注: {cap_nums}')
    print(f'  正文引用: {refs}')
    missing = [r for r in refs if r.split()[0] + ' ' + r.split()[1][0] not in
               [c.split()[0] + ' ' + c.split()[1][0] for c in cap_nums]]
    if missing:
        print(f'  ⚠ 引用但无图注: {missing}')
        issues += 1
    else:
        print('  ✓ 所有 Figure 引用都有对应图注')

    tabs = sorted(set(re.findall(r'Table \d+', full)))
    print(f'  Table 引用: {tabs}')

    print()
    print(f'合计 {issues} 类一致性问题')


if __name__ == '__main__':
    main()
