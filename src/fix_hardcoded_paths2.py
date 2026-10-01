# -*- coding: utf-8 -*-
"""fix_hardcoded_paths2.py — 第二批：ROOT/OUT 常量、CLI 默认值、双引号 sys.path"""
import os
import re
import shutil

SRC = r'D:\GED_mutation\src'
BAK = r'D:\GED_mutation\_pathfix_backup2'

# 顺序有意义：先长后短
RULES = [
    # 双引号版 sys.path
    (re.compile(r'sys\.path\.insert\(0,\s*r"D:[\\/][^"]*[\\/]src"\s*\)'),
     "sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))"),
    # ROOT / OUT 常量
    (re.compile(r"^ROOT\s*=\s*r'D:\\GED_mutation'\s*$", re.M),
     "ROOT = str(__import__('pathlib').Path(__file__).resolve().parent.parent)"),
    (re.compile(r"^OUT\s*=\s*r'D:\\GED_mutation\\figures\\publication'\s*$", re.M),
     "OUT = str(__import__('pathlib').Path(__file__).resolve().parent.parent"
     " / 'figures' / 'publication')"),
    (re.compile(r"^OUT\s*=\s*r'D:\\GED_mutation\\([^']+)'\s*$", re.M),
     r"OUT = str(__import__('pathlib').Path(__file__).resolve().parent.parent / '\1')"),
    (re.compile(r"os\.path\.join\(r'D:\\GED_mutation\\src',\s*f\)"),
     "os.path.join(os.path.dirname(os.path.abspath(__file__)), f)"),
    # train_gnn_baseline CLI 默认值
    (re.compile(r"default=r'D:\\GED_mutation\\data\\training_merged\.csv'"),
     "default=str(__import__('pathlib').Path(__file__).resolve().parent.parent"
     " / 'data' / 'training_merged.csv')"),
    (re.compile(r"default=r'D:\\GED_mutation\\data\\contact_graphs_megascale'"),
     "default=str(__import__('pathlib').Path(__file__).resolve().parent.parent"
     " / 'data' / 'contact_graphs_megascale')"),
    (re.compile(r"default=r'D:\\GED_mutation\\data\\contact_graphs_thermomutdb'"),
     "default=str(__import__('pathlib').Path(__file__).resolve().parent.parent"
     " / 'data' / 'contact_graphs_thermomutdb')"),
    (re.compile(r"default=r'D:\\GED_mutation\\data\\contact_graphs_s669'"),
     "default=str(__import__('pathlib').Path(__file__).resolve().parent.parent"
     " / 'data' / 'contact_graphs_s669')"),
    (re.compile(r"default=r'D:\\GED_mutation\\data\\benchmarks_s669_clean\.csv'"),
     "default=str(__import__('pathlib').Path(__file__).resolve().parent.parent"
     " / 'data' / 'benchmarks_s669_clean.csv')"),
]


def main():
    os.makedirs(BAK, exist_ok=True)
    changed, rolled = [], []
    for f in sorted(os.listdir(SRC)):
        if not f.endswith('.py'):
            continue
        p = os.path.join(SRC, f)
        orig = open(p, encoding='utf-8').read()
        t = orig
        hits = []
        for pat, rep in RULES:
            if pat.search(t):
                t = pat.sub(rep, t)
                hits.append(pat.pattern[:38])
        if t == orig:
            continue
        if not os.path.exists(os.path.join(BAK, f)):
            shutil.copy2(p, os.path.join(BAK, f))
        open(p, 'w', encoding='utf-8').write(t)
        try:
            compile(t, p, 'exec')
            changed.append((f, len(hits)))
        except SyntaxError as e:
            shutil.copy2(os.path.join(BAK, f), p)
            rolled.append((f, str(e)))

    print(f'修改 {len(changed)} 个文件')
    for f, n in changed:
        print(f'  {f:<44} {n} 条规则')
    if rolled:
        print(f'\n回滚 {len(rolled)}:')
        for f, e in rolled:
            print(f'  {f}: {e}')


if __name__ == '__main__':
    main()
