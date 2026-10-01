# -*- coding: utf-8 -*-
"""fix_hardcoded_paths4.py — 最后一批 + 渲染工具路径可配置"""
import os
import re
import shutil

SRC = r'D:\GED_mutation\src'
BAK = r'D:\GED_mutation\_pathfix_backup4'
P = "__import__('pathlib').Path(__file__).resolve().parent.parent"

RULES = [
    # check_submission.py: os.listdir(r'D:\...\figures')
    (re.compile(r"os\.listdir\(r'D:\\GED_mutation\\figures'\)"),
     f"os.listdir(str({P} / 'figures'))"),
    # test_scwrl4.py
    (re.compile(r"r'D:\\GED_mutation\\tools\\scwrl4\\test_mut'"),
     f"str({P} / 'tools' / 'scwrl4' / 'test_mut')"),
    # submission_gate.py: 渲染工具改为环境变量可配置
    (re.compile(
        r"^(\s*)cli = r'D:\\新建文件夹\\resources\\app\.asar\.unpacked\\dsh\\node_modules"
        r"\\@deepseek-ai\\libreoffice-kit\\lib\\cli\.js'\s*$", re.M),
     r"\1cli = os.environ.get('DSH_LO_CLI', r'D:\新建文件夹\resources\app.asar.unpacked"
     r"\dsh\node_modules\@deepseek-ai\libreoffice-kit\lib\cli.js')"),
    (re.compile(
        r"^(\s*)node = r'D:\\新建文件夹\\resources\\runtime\\primary-runtime\\dependencies"
        r"\\node\\bin\\node\.exe'\s*$", re.M),
     r"\1node = os.environ.get('DSH_NODE', r'D:\新建文件夹\resources\runtime"
     r"\primary-runtime\dependencies\node\bin\node.exe')"),
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
                # lambda 避免替换串中的反斜杠被当作转义序列
                t = pat.sub(lambda m, r=rep: r, t)
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
        print(f'  {f:<40} {n} 处')
    if rolled:
        for f, e in rolled:
            print(f'  回滚 {f}: {e}')


if __name__ == '__main__':
    main()
