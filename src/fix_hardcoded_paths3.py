# -*- coding: utf-8 -*-
"""fix_hardcoded_paths3.py — 第三批：data/ figures/ tools/ 前缀的项目内路径"""
import os
import re
import shutil

SRC = r'D:\GED_mutation\src'
BAK = r'D:\GED_mutation\_pathfix_backup3'

P = "__import__('pathlib').Path(__file__).resolve().parent.parent"

RULES = [
    # 项目内目录常量
    (re.compile(r"r'D:\\GED_mutation\\tools\\scwrl4\\runs'"), f"str({P} / 'tools' / 'scwrl4' / 'runs')"),
    (re.compile(r"r'D:\\GED_mutation\\tools\\scwrl4\\Scwrl4\.exe'"),
     f"str({P} / 'tools' / 'scwrl4' / 'Scwrl4.exe')"),
    (re.compile(r"r'D:\\GED_mutation\\tools\\foldx\\repro_v2'"),
     f"str({P} / 'tools' / 'foldx' / 'repro_v2')"),
    (re.compile(r"r'D:\\GED_mutation\\tools\\foldx\\foldx_1_20270131\.exe'"),
     "os.environ.get('FOLDX_BIN', str(" + P + " / 'tools' / 'foldx' / 'foldx_1_20270131.exe'))"),
    (re.compile(r"r'D:\\GED_mutation\\tools\\blast\\ncbi-blast-2\.17\.0\+\\bin'"),
     "os.environ.get('BLAST_BIN', str(" + P + " / 'tools' / 'blast' / 'ncbi-blast-2.17.0+' / 'bin'))"),
    # 具体数据文件
    (re.compile(r"r'D:\\GED_mutation\\data\\structures\\pdb1bfm\.ent'"),
     f"str({P} / 'data' / 'structures' / 'pdb1bfm.ent')"),
    (re.compile(r"r'D:\\GED_mutation\\data\\mutant_structures_s669\\1BFM\\M35W\\work\\out\\pdb1bfm_1\.pdb'"),
     f"str({P} / 'data' / 'mutant_structures_s669' / '1BFM' / 'M35W' / 'work' / 'out' / 'pdb1bfm_1.pdb')"),
    (re.compile(r"r'D:\\GED_mutation\\data\\contact_graphs_s669_sc\\1BFM\.npz'"),
     f"str({P} / 'data' / 'contact_graphs_s669_sc' / '1BFM.npz')"),
    # 目录前缀
    (re.compile(r"os\.path\.join\(r'D:\\GED_mutation\\figures',"), f"os.path.join(str({P} / 'figures'),"),
    (re.compile(r"os\.path\.join\(r'D:\\GED_mutation',"), f"os.path.join(str({P}),"),
    (re.compile(r"open\(r'D:\\GED_mutation\\manuscript_draft_v0\.md'"),
     f"open(str({P} / 'manuscript_draft_v0.md')"),
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
        n = 0
        for pat, rep in RULES:
            c = len(pat.findall(t))
            if c:
                t = pat.sub(rep, t)
                n += c
        if t == orig:
            continue
        if not os.path.exists(os.path.join(BAK, f)):
            shutil.copy2(p, os.path.join(BAK, f))
        open(p, 'w', encoding='utf-8').write(t)
        try:
            compile(t, p, 'exec')
            changed.append((f, n))
        except SyntaxError as e:
            shutil.copy2(os.path.join(BAK, f), p)
            rolled.append((f, str(e)))
    print(f'修改 {len(changed)} 个文件, 共 {sum(n for _, n in changed)} 处')
    for f, n in changed:
        print(f'  {f:<42} {n} 处')
    if rolled:
        print(f'\n回滚 {len(rolled)}:')
        for f, e in rolled:
            print(f'  {f}: {e}')


if __name__ == '__main__':
    main()
