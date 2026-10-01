# -*- coding: utf-8 -*-
"""
fix_hardcoded_paths.py — 让发布包脚本可在他人机器上运行

只处理三种明确、安全的模式：
  1. sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent))   →  os.path.dirname(__file__)
  2. DATA = r'D:\\...\\data'                →  若上方已 from paths import DATA 则删除该行
                                              否则替换为 Path(__file__).resolve().parent.parent/'data'
  3. os.environ.setdefault('HF_HOME', str(__import__('pathlib').Path(__file__).resolve().parent.parent / '.hf_cache'))  →  相对本文件的路径
其余硬编码（工具目录、外部程序）保留并记录，因为它们指向不随仓库分发的二进制。
每个文件修改后立即做语法校验；失败则回滚该文件。
"""
import os
import re
import shutil

SRC = r'D:\GED_mutation\src'
BAKDIR = r'D:\GED_mutation\_pathfix_backup'

PAT_SYSPATH = re.compile(r"sys\.path\.insert\(0,\s*r'[A-Za-z]:\\[^']*\\src'\s*\)")
PAT_DATA = re.compile(r"^DATA\s*=\s*r'([A-Za-z]:\\[^']*)\\data'\s*$", re.M)
PAT_HF = re.compile(r"os\.environ\.setdefault\('HF_HOME',\s*r'([A-Za-z]:\\[^']*)'\)")

REPL_SYSPATH = "sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent))"
PORTABLE_DATA = ("from pathlib import Path as _P\n"
                 "DATA = str(_P(__file__).resolve().parent.parent / 'data')")


def main():
    os.makedirs(BAKDIR, exist_ok=True)
    changed, deleted_lines, skipped = [], [], []
    for f in sorted(os.listdir(SRC)):
        if not f.endswith('.py'):
            continue
        p = os.path.join(SRC, f)
        orig = open(p, encoding='utf-8').read()
        t = orig
        hits = []

        # 1) sys.path
        if PAT_SYSPATH.search(t):
            t = PAT_SYSPATH.sub(REPL_SYSPATH, t)
            hits.append('sys.path')

        # 2) HF_HOME / 其它缓存目录
        m = PAT_HF.search(t)
        if m:
            t = PAT_HF.sub(
                "os.environ.setdefault('HF_HOME', str(__import__('pathlib').Path(__file__)"
                ".resolve().parent.parent / '.hf_cache'))", t)
            hits.append('HF_HOME')

        # 3) DATA
        m = PAT_DATA.search(t)
        if m:
            has_import = re.search(r'^from paths import .*\bDATA\b', t, re.M)
            if has_import:
                # 冗余覆盖：删除该行
                t = PAT_DATA.sub('', t)
                hits.append('DATA(删除冗余)')
                deleted_lines.append(f)
            else:
                t = PAT_DATA.sub(PORTABLE_DATA, t)
                hits.append('DATA(改为相对)')

        if t == orig:
            continue

        if not os.path.exists(os.path.join(BAKDIR, f)):
            shutil.copy2(p, os.path.join(BAKDIR, f))
        open(p, 'w', encoding='utf-8').write(t)
        try:
            compile(t, p, 'exec')
            changed.append((f, hits))
        except SyntaxError as e:
            shutil.copy2(os.path.join(BAKDIR, f), p)
            skipped.append((f, str(e)))

    print(f'修改成功 {len(changed)} 个文件')
    for f, h in changed:
        print(f'  {f:<44} {h}')
    if skipped:
        print()
        print(f'回滚 {len(skipped)} 个（语法错误）:')
        for f, e in skipped:
            print(f'  {f}: {e}')
    print()
    return changed


if __name__ == '__main__':
    main()
