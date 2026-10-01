# -*- coding: utf-8 -*-
"""fill_author_info.py — 把已知的作者信息填入标题页（待用户确认/补充）

已知（来自 CITATION.cff）：
  Yilin Huang, Xi'an Jiaotong-Liverpool University, Yilin.Huang24@student.xjtlu.edu.cn
"""
import os
import re

ROOT = r'D:\GED_mutation'
TP = os.path.join(ROOT, '标题页与投稿清单_TitlePage_Checklist.md')

OLD_TABLE = """| Order | Full name | Affiliation | ORCID | Email |
|---|---|---|---|---|
| 1 | *to be supplied* | *to be supplied* | *to be supplied* | *to be supplied* |
| 2 | *to be supplied* | *to be supplied* | *to be supplied* | *to be supplied* |"""

NEW_TABLE = """> **Status: single author on record. Confirm before submission.**
> The details below are taken from `CITATION.cff`. If the manuscript has
> additional authors, add their rows and update `CITATION.cff` to match.

| Order | Full name | Affiliation | ORCID | Email |
|---|---|---|---|---|
| 1 | Yilin Huang | Xi'an Jiaotong-Liverpool University | *to be supplied* (optional but recommended) | Yilin.Huang24@student.xjtlu.edu.cn |"""

OLD_CORR = """**Corresponding author**
Yilin Huang — Yilin.Huang24@student.xjtlu.edu.cn

*Replace the placeholder rows above with the final author list, affiliations and ORCIDs before submission.*"""

NEW_CORR = """**Corresponding author**
Yilin Huang — Yilin.Huang24@student.xjtlu.edu.cn
Xi'an Jiaotong-Liverpool University

*If this is a single-author manuscript, no further author information is needed.
An ORCID iD is optional but recommended; add it above and in `CITATION.cff`.*"""


def main():
    t = open(TP, encoding='utf-8').read()
    n = 0
    for a, b in [(OLD_TABLE, NEW_TABLE), (OLD_CORR, NEW_CORR)]:
        if a in t:
            t = t.replace(a, b)
            n += 1
    open(TP, 'w', encoding='utf-8').write(t)
    print(f'标题页: {n} 处')

    # 校验
    t2 = open(TP, encoding='utf-8').read()
    print(f'  残留 "to be supplied" 作者行: {t2.count("| *to be supplied* |")}')
    print(f'  含作者名: {"Yilin Huang" in t2}')
    print(f'  含单位: {"Xi\'an Jiaotong-Liverpool" in t2}')

    # 同步到发布包
    import shutil
    dst = os.path.join(ROOT, 'release', 'docs', '标题页与投稿清单_TitlePage_Checklist.md')
    if os.path.isdir(os.path.dirname(dst)):
        shutil.copy2(TP, dst)
        print('  发布包副本已同步')


if __name__ == '__main__':
    main()
