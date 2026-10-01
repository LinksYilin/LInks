# -*- coding: utf-8 -*-
"""validate_release.py — 发布包可移植性与完整性验证

重点检查硬编码绝对路径，因为发布包要在他人机器上运行。
"""
import os
import re
import sys

_HERE = __import__('pathlib').Path(__file__).resolve().parent
_ROOT = _HERE.parent
# 两种布局都要支持：
#   工作副本  <root>/src/validate_release.py，发布包在 <root>/release/
#   发布包内  <root>/src/validate_release.py，本身就是发布包（无 <root>/release/）
REL = str(_ROOT / 'release') if (_ROOT / 'release').is_dir() else str(_ROOT)
ROOT = str(_ROOT)
# 允许出现的路径引用（相对或环境变量）
OK_PATTERNS = [
    r'os\.path\.dirname',
    r'__file__',
    r'GED_ROOT',
    r'os\.environ',
    r'Path\(',
]


def main():
    issues = []
    warnings = []

    # 1) Python 语法
    bad_syntax = []
    for dirpath, _, files in os.walk(os.path.join(REL, 'src')):
        for f in files:
            if not f.endswith('.py'):
                continue
            p = os.path.join(dirpath, f)
            try:
                compile(open(p, encoding='utf-8').read(), p, 'exec')
            except SyntaxError as e:
                bad_syntax.append(f'{f}: {e}')
    if bad_syntax:
        issues.append(f'语法错误 {len(bad_syntax)} 个: {bad_syntax[:3]}')

    # 2) 硬编码绝对路径（会导致他人无法运行）
    hard = {}
    for dirpath, _, files in os.walk(os.path.join(REL, 'src')):
        for f in files:
            if not f.endswith('.py'):
                continue
            # 路径修复工具自身把这些模式作为正则字符串，属正常
            if f.startswith('fix_hardcoded_paths') or f.startswith('check_hardcoded'):
                continue
            p = os.path.join(dirpath, f)
            t = open(p, encoding='utf-8').read()
            for m in re.finditer(r'[\'"]([A-Za-z]:\\\\[^\'"]{3,})[\'"]', t):
                hard.setdefault(f, []).append(m.group(1))
    if hard:
        warnings.append(f'{len(hard)} 个脚本含硬编码盘符路径')
        for f, paths in list(hard.items())[:8]:
            warnings.append(f'    {f}: {paths[0]}')

    # 3) CSV 可读性
    import pandas as pd
    bad_csv = []
    n_csv = 0
    for dirpath, _, files in os.walk(os.path.join(REL, 'results')):
        for f in files:
            if not f.endswith('.csv'):
                continue
            n_csv += 1
            p = os.path.join(dirpath, f)
            try:
                pd.read_csv(p, nrows=5)
            except Exception as e:
                bad_csv.append(f'{f}: {type(e).__name__}')
    if bad_csv:
        issues.append(f'不可读 CSV {len(bad_csv)} 个: {bad_csv[:5]}')

    # 4) 图表文件
    figdir = os.path.join(REL, 'figures', 'publication')
    figs = sorted(os.listdir(figdir)) if os.path.isdir(figdir) else []
    for want in ['figure1_study_design.png', 'figure2_contact_rewiring.png',
                 'figure3_capacity_ladder.png', 'figure4_information_budget.png',
                 'figure5_contact_diagnostics.png']:
        if want not in figs:
            issues.append(f'缺少图: {want}')

    # 5) 顶层必备文件
    for f in ['README.md', 'LICENSE', 'CITATION.cff', 'requirements.lock', '.gitignore']:
        if not os.path.exists(os.path.join(REL, f)):
            issues.append(f'缺少顶层文件: {f}')

    # 6) README 中的路径引用是否存在
    readme = open(os.path.join(REL, 'README.md'), encoding='utf-8').read()
    for m in re.finditer(r'`([a-zA-Z0-9_./\-]+\.(?:py|csv|md|png))`', readme):
        ref = m.group(1)
        if ref.startswith(('src/', 'results/', 'figures/', 'docs/')):
            if not os.path.exists(os.path.join(REL, ref)):
                warnings.append(f'README 引用了不存在的文件: {ref}')

    # 7) 路径解析：paths.py 必须能解析到一个存在的 DATA 目录
    #    （发布包用 results/，完整工作副本用 data/）
    try:
        import importlib.util
        src_paths = os.path.join(REL, 'src', 'paths.py')
        if os.path.exists(src_paths):
            spec = importlib.util.spec_from_file_location('_rel_paths', src_paths)
            mod = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mod)
            if not mod.DATA.exists():
                issues.append(f'paths.DATA 解析到不存在的目录: {mod.DATA}')
            else:
                n_in_data = len([f for f in os.listdir(mod.DATA) if f.endswith('.csv')])
                if n_in_data == 0:
                    issues.append(f'paths.DATA 目录中没有 CSV: {mod.DATA}')
    except Exception as e:
        issues.append(f'paths.py 无法导入: {type(e).__name__}: {e}')

    # 8) 本地模块导入完整性（缺失的本地模块会使脚本无法运行）
    import ast as _ast
    local = {f[:-3] for f in os.listdir(os.path.join(REL, 'src')) if f.endswith('.py')}
    KNOWN = set(getattr(sys, 'stdlib_module_names', ())) | {
        'numpy', 'pandas', 'scipy', 'torch', 'torch_geometric', 'sklearn', 'docx',
        'matplotlib', 'Bio', 'PIL', 'requests', 'statsmodels', 'networkx', 'lxml',
        'openpyxl', 'joblib', 'tqdm', 'yaml', 'h5py', 'seaborn', 'transformers',
        'edit_cost', 'ged_module', 'gnn_baseline', 'gnn_local_baseline',
        'gnn_edge_baseline', 'train_gnn_baseline', 'contact_graph_defs', 'paths',
        'seed_utils', 'strong_backbones', 'ladder_common', 'filter_leakage',
    }
    broken = []
    for f in sorted(os.listdir(os.path.join(REL, 'src'))):
        if not f.endswith('.py') or f.startswith(('fix_hardcoded', 'check_hardcoded')):
            continue
        p = os.path.join(REL, 'src', f)
        try:
            tree = _ast.parse(open(p, encoding='utf-8').read())
        except SyntaxError:
            continue
        imported = set()
        for n in _ast.walk(tree):
            if isinstance(n, _ast.Import):
                imported |= {a.name.split('.')[0] for a in n.names}
            elif isinstance(n, _ast.ImportFrom) and n.module and n.level == 0:
                imported.add(n.module.split('.')[0])
        miss = sorted(m for m in imported
                      if m not in local and m not in KNOWN and not m.startswith('_'))
        if miss:
            broken.append(f'{f} -> {miss}')
    if broken:
        issues.append(f'{len(broken)} 个脚本引用了缺失的本地模块')
        for b in broken[:6]:
            warnings.append(f'    {b}')

    print('=' * 74)
    print('发布包验证')
    print('=' * 74)
    print(f'脚本 {len([f for _,_,fs in os.walk(os.path.join(REL,"src")) for f in fs if f.endswith(".py")])} 个')
    print(f'CSV {n_csv} 个（不可读 {len(bad_csv)}）')
    print(f'图表 {len(figs)} 个')
    print()
    if issues:
        print('❌ 必须修复:')
        for x in issues:
            print(f'   {x}')
    else:
        print('✅ 无阻断性问题')
    if warnings:
        print()
        print('⚠ 需注意:')
        for x in warnings[:14]:
            print(f'   {x}')
    return 1 if issues else 0


if __name__ == '__main__':
    raise SystemExit(main())
