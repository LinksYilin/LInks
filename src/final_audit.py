"""
final_audit.py — 投稿前最终全面检查
=====================================
检查：文件完整性、数据一致性、论文数字、图表、代码、遗留问题。
"""
import os

import pandas as pd

ROOT = str(__import__('pathlib').Path(__file__).resolve().parent.parent)
DATA = os.path.join(ROOT, 'data')

print('=' * 70)
print('1. 论文文件')
print('=' * 70)
for f in ['manuscript_revised_20260930.docx', 'manuscript_revised_20260930_v2.docx']:
    p = os.path.join(ROOT, f)
    if os.path.exists(p):
        print(f'  {"OK " if os.path.exists(p) else ""} {f}: {os.path.getsize(p)/1024:.0f} KB, '
              f'修改 {pd.Timestamp(os.path.getmtime(p), unit="s").strftime("%Y-%m-%d %H:%M")}')

print()
print('=' * 70)
print('2. 图表文件（publication）')
print('=' * 70)
pub = os.path.join(ROOT, 'figures', 'publication')
for f in sorted(os.listdir(pub)):
    sz = os.path.getsize(os.path.join(pub, f))
    print(f'  {f:<48} {sz/1024:>8.1f} KB')

print()
print('=' * 70)
print('3. 接触图数据（修正后）')
print('=' * 70)
for d in ['contact_graphs_s669_sc', 'contact_graphs_ssym_sc', 'contact_graphs_thermomutdb_sc',
          'contact_graphs_megascale_sc', 'contact_graphs_s669_cb', 'contact_graphs_s669_ca']:
    p = os.path.join(DATA, d)
    if os.path.exists(p):
        n = len([f for f in os.listdir(p) if f.endswith('.npz')])
        print(f'  {d:<36} {n:>4} 图')
    else:
        print(f'  {d:<36} (不存在)')

print()
print('=' * 70)
print('4. 关键结果文件')
print('=' * 70)
for f in ['benchmark_summary.csv', 'benchmark_predictions_final.csv', 'paired_comparison.csv',
          'representation_comparison.csv', 'edits_corrected.csv', 'edits_corrected_summary.csv',
          'mechanism_corrected.csv', 'locality_corrected.csv', 'foldx_reproducibility_v2.csv',
          'index_bug_audit_s669.csv', 'hydrogen_bias_audit.csv']:
    p = os.path.join(DATA, f)
    if os.path.exists(p):
        try:
            d = pd.read_csv(p)
            print(f'  {f:<40} {len(d):>6} 行, {len(d.columns)} 列')
        except Exception as e:
            print(f'  {f:<40} 读取失败 {e}')
    else:
        print(f'  {f:<40} (缺失)')

print()
print('=' * 70)
print('5. 数据一致性交叉核对')
print('=' * 70)
# 编辑统计
ec = pd.read_csv(os.path.join(DATA, 'edits_corrected.csv'))
c8 = ec[(ec['atom_def'] == 'centroid') & (ec['threshold'] == 8.0)]
print(f'  edits_corrected (centroid, 8Å): {len(c8)} 突变对')
print(f'    断边均值 {c8["n_broken"].mean():.2f}, 成边均值 {c8["n_formed"].mean():.2f}, '
      f'≥1变化 {(c8["n_edit"] > 0).mean()*100:.1f}%')

# 预测结果
bs = pd.read_csv(os.path.join(DATA, 'benchmark_summary.csv'))
for bench in bs['benchmark'].unique():
    sub = bs[bs['benchmark'] == bench]
    n = sub['n'].iloc[0]
    print(f'  {bench}: n={n}, 模型 {sub["model"].nunique()} 个')

# 论文中的关键数字是否与结果文件一致
man = open(os.path.join(ROOT, 'v2_全文导出.md'), encoding='utf-8').read() if os.path.exists(
    os.path.join(ROOT, 'v2_全文导出.md')) else ''
print()
print('  论文 vs 结果文件一致性:')
pairs = [
    ('511', 'benchmark_summary S669 n=511', int(bs[bs['benchmark'] == 'S669']['n'].iloc[0]) == 511),
    ('90.7%', 'centroid ≥1变化', abs((c8['n_edit'] > 0).mean() * 100 - 90.7) < 0.1),
    ('2.11', 'centroid 断边均值', abs(c8['n_broken'].mean() - 2.11) < 0.01),
]
for label, desc, ok in pairs:
    print(f'    {"OK " if ok else "MISMATCH"} {label} <=> {desc}')

print()
print('=' * 70)
print('6. 测试套件')
print('=' * 70)
