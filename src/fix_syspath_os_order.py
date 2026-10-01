# -*- coding: utf-8 -*-
"""fix_syspath_os_order.py — 修正 sys.path 替换引入的 os 依赖问题

若脚本在 import os 之前就调用 sys.path.insert(...os.path.dirname...)，会 NameError。
统一改成不依赖 os 的写法。
"""
import os
import re
import shutil

SRC = r'D:\GED_mutation\src'
BAK = r'D:\GED_mutation\_pathfix_backup3'

OLD = re.compile(r"sys\.path\.insert\(0,\s*os\.path\.dirname\(os\.path\.abspath\(__file__\)\)\)")
NEW = "sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent))"


def main():
    os.makedirs(BAK, exist_ok=True)
    fixed = []
    for f in sorted(os.listdir(SRC)):
        if not f.endswith('.py'):
            continue
        p = os.path.join(SRC, f)
        t = open(p, encoding='utf-8').read()
        if not OLD.search(t):
            continue
        # 检查 import os 是否出现在第一次使用之前
        first_use = OLD.search(t).start()
        m = re.search(r'^\s*(?:import os\b|import .*\bos\b)', t[:first_use], re.M)
        if m:
            continue  # os 已可用
        if not os.path.exists(os.path.join(BAK, f)):
            shutil.copy2(p, os.path.join(BAK, f))
        new = OLD.sub(NEW, t)
        open(p, 'w', encoding='utf-8').write(new)
        try:
            compile(new, p, 'exec')
            fixed.append(f)
        except SyntaxError as e:
            shutil.copy2(os.path.join(BAK, f), p)
            print(f'  回滚 {f}: {e}')
    print(f'修正 {len(fixed)} 个文件（os 未提前导入）:')
    for f in fixed:
        print(f'  {f}')


if __name__ == '__main__':
    main()
