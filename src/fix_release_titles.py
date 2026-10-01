# -*- coding: utf-8 -*-
"""fix_release_titles.py — 更新发布包公开材料中的论文标题

CITATION.cff 与 README.md 仍用旧标题；src/ 与 docs/ 中的旧标题是历史修订记录，
属于正常留痕，不修改。
"""
import os
import re

REL = r'D:\GED_mutation\release'
OLD_A = ('Mutation-sensitive residue contact graphs: a controlled benchmark of graph '
         'representations for \u0394\u0394G prediction')
OLD_B = ('Mutation-sensitive residue contact graphs do not improve \u0394\u0394G prediction: '
         'a controlled benchmark of graph representations')
OLD_C = ('Mutation-sensitive residue contact graphs do not reliably improve \u0394\u0394G prediction: '
         'a controlled benchmark of graph representations')
NEW = ('Contact-graph definitions shape what a model describes but not what it predicts: '
       'a controlled comparison across five encoder capacities')


def main():
    targets = [os.path.join(REL, 'CITATION.cff'), os.path.join(REL, 'README.md')]
    n = 0
    for p in targets:
        if not os.path.exists(p):
            print(f'  缺失: {p}')
            continue
        t = open(p, encoding='utf-8').read()
        orig = t
        for old in (OLD_A, OLD_B, OLD_C):
            if old in t:
                t = t.replace(old, NEW)
                n += 1
        if t != orig:
            open(p, 'w', encoding='utf-8').write(t)
            print(f'  已更新 {os.path.basename(p)}')

    # 校验
    print()
    for p in targets:
        t = open(p, encoding='utf-8').read()
        bad = sum(t.count(x) for x in (OLD_A, OLD_B, OLD_C))
        good = t.count(NEW)
        print(f'  {os.path.basename(p)}: 新标题 {good} 处, 旧标题 {bad} 处')

    # CITATION.cff 的 doi/url 占位提醒
    c = open(os.path.join(REL, 'CITATION.cff'), encoding='utf-8').read()
    for marker in ['【USER】', '【REPO】']:
        if marker in c:
            print(f'  ⚠ CITATION.cff 仍含占位符 {marker}（需填入 GitHub 用户名/仓库名）')
    return n


if __name__ == '__main__':
    main()
