# -*- coding: utf-8 -*-
"""audit_methods_code.py — 检查 Methods 描述与发布代码是否一致"""
import os
import re

import pandas as pd
from docx import Document

DOC = r'D:\GED_mutation\manuscript_routeA.docx'
SRC = r'D:\GED_mutation\release\src'

d = Document(DOC)
ps = [p.text.strip() for p in d.paragraphs]
iref = next((i for i, t in enumerate(ps) if t == 'References'), len(ps))
methods = '\n'.join(t for t in ps[:iref])


def code(name):
    p = os.path.join(SRC, name)
    return open(p, encoding='utf-8').read() if os.path.exists(p) else ''


print('=' * 76)
print('Methods 描述 ↔ 代码 一致性检查')
print('=' * 76)

checks = []

# 1) 节点特征维度
lc = code('ladder_common.py')
m = re.search(r'nodes\[mi, :20\]', lc)
checks.append(('节点特征：20 维 one-hot', bool(m), 'ladder_common.make_sample'))
m2 = re.search(r'np\.concatenate\(\[nodes,\s*flag\]', lc)
checks.append(('节点特征：拼接 1 维突变标记', bool(m2), 'ladder_common.make_sample'))
print('\n[3.5 节点特征]')
print(f'  论文称 24 维 = 20 one-hot + 3 理化 + 1 标记')
for label, ok, where in checks:
    print(f'  {"OK " if ok else "FAIL"} {label:<36} ({where})')

# 2) 边特征维度
sb = code('strong_backbones.py')
rng = code('run_ladder.py')
m3 = re.search(r'edge_dim\s*=\s*2', rng + sb)
print(f'\n[3.5 边特征]')
print(f'  {"OK " if m3 else "FAIL"} 边特征 2 维（疏水/静电）          (run_ladder/strong_backbones)')

# 3) 五个编码器参数
print('\n[3.5 五个编码器]')
import subprocess
out = subprocess.run(
    [r'D:\GED_mutation\.venv\Scripts\python.exe', '-c',
     'import sys; sys.path.insert(0,r"' + SRC + '");'
     'from run_ladder import make_model;'
     'print("|".join(f"{n}={sum(p.numel() for p in make_model(n,24).parameters())}"'
     ' for n in ["gnn_global","gnn_local","gnn_edge","deep_gine","egnn"]))'],
    capture_output=True, text=True, cwd=SRC)
print('  ' + (out.stdout.strip() or out.stderr.strip()[:200]))

# 4) 训练设置
print('\n[3.6 训练设置]')
for pat, label in [(r'epochs', 'epochs 可配置'), (r'lr\s*=\s*1e-3', 'lr=1e-3'),
                   (r'batch_size\s*=\s*64', 'batch=64'), (r'Adam', 'Adam 优化器')]:
    found = bool(re.search(pat, rng))
    print(f'  {"OK " if found else "FAIL"} {label:<28} (run_ladder.py)')

# 5) bootstrap 重采样单位
print('\n[3.6 bootstrap]')
la = code('analyze_ladder_3seed.py')
for pat, label in [(r'protein', '按蛋白重采样'), (r'seed', '同时重采样种子'),
                   (r'percentile', '百分位区间')]:
    found = bool(re.search(pat, la, re.I))
    print(f'  {"OK " if found else "FAIL"} {label:<28} (analyze_ladder_3seed.py)')

# 6) 多重比较
print('\n[3.6 多重比较]')
print(f'  {"OK " if "holm" in la.lower() or "Holm" in la else "FAIL"} Holm step-down 校正          (analyze_ladder_3seed.py)')

# 7) 四个图定义
print('\n[3.1 四个图定义]')
cg = code('contact_graph_defs.py')
for ad in ['ca', 'cb', 'centroid', 'allatom']:
    found = ad in cg or ad in lc
    print(f'  {"OK " if found else "FAIL"} {ad}')
th = re.search(r'THRESHOLDS\s*=\s*\[([^\]]+)\]', code('compute_edits_by_definition.py'))
print(f'  阈值: {th.group(1) if th else "未找到"}')

# 8) 确定性
su = code('seed_utils.py')
print('\n[3.6 确定性]')
for pat, label in [(r'CUBLAS_WORKSPACE_CONFIG', 'CUBLAS 工作区设置'),
                   (r'use_deterministic_algorithms', '确定性算法启用')]:
    found = bool(re.search(pat, su))
    print(f'  {"OK " if found else "FAIL"} {label:<28} (seed_utils.py)')

print()
print('=' * 76)
