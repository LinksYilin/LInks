# -*- coding: utf-8 -*-
"""check_release_runnability.py — 检查发布包脚本能否导入（修正版）

只标记「既不是标准库、也不是已安装的第三方库、也不在发布包内」的模块。
"""
import ast
import os
import sys

REL = r'D:\GED_mutation\release'
SRC = os.path.join(REL, 'src')


def stdlib_names():
    return set(getattr(sys, 'stdlib_module_names', ())) | {
        'numpy', 'pandas', 'scipy', 'torch', 'torch_geometric', 'sklearn', 'docx',
        'matplotlib', 'Bio', 'transformers', 'PIL', 'requests', 'statsmodels', 'networkx',
        'lxml', 'openpyxl', 'joblib', 'tqdm', 'yaml', 'h5py', 'seaborn',
    }


def main():
    local = {f[:-3] for f in os.listdir(SRC) if f.endswith('.py')}
    known = stdlib_names()
    print(f'本地模块 {len(local)} 个 | 已知外部模块 {len(known)} 个')
    print()

    problems = []
    for f in sorted(os.listdir(SRC)):
        if not f.endswith('.py'):
            continue
        p = os.path.join(SRC, f)
        try:
            tree = ast.parse(open(p, encoding='utf-8').read())
        except SyntaxError as e:
            problems.append((f, 'SYNTAX', str(e)))
            continue
        imported = set()
        for n in ast.walk(tree):
            if isinstance(n, ast.Import):
                for a in n.names:
                    imported.add(a.name.split('.')[0])
            elif isinstance(n, ast.ImportFrom):
                if n.module and n.level == 0:
                    imported.add(n.module.split('.')[0])
        # 只关心：不是本地、也不是已知外部 → 可能是缺失的本地模块
        suspicious = sorted(m for m in imported
                            if m not in local and m not in known and not m.startswith('_'))
        if suspicious:
            problems.append((f, 'IMPORT', suspicious))

    print(f'有问题的脚本: {len(problems)} 个')
    for f, kind, detail in problems:
        print(f'  ❌ [{kind}] {f}: {detail}')

    # 检查 3.7 关键流程脚本是否可导入
    print()
    print('=== 论文关键流程脚本的存在性 ===')
    key = ['run_ladder.py', 'ladder_common.py', 'strong_backbones.py',
           'analyze_ladder_3seed.py', 'locality_no_self.py', 'esm2_fusion_fast.py',
           'seed_utils.py', 'paths.py', 'contact_graph_defs.py']
    for k in key:
        print(f'  {"✅" if os.path.exists(os.path.join(SRC, k)) else "❌"} {k}')

    return len(problems)


if __name__ == '__main__':
    sys.exit(0 if main() == 0 else 1)
