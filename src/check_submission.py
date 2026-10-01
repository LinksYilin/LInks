"""投稿就绪性检查"""
import os

from PIL import Image

print('=== 图表状态 ===')
for f in sorted(os.listdir(str(__import__('pathlib').Path(__file__).resolve().parent.parent / 'figures'))):
    if f.endswith('.png'):
        img = Image.open(os.path.join(str(__import__('pathlib').Path(__file__).resolve().parent.parent / 'figures'), f))
        dpi = img.info.get('dpi', (72, 72))
        print(f'  {f}: {img.size[0]}x{img.size[1]}, DPI={dpi[0]:.0f}')

print()
print('=== 论文关键要素 ===')
txt = open(str(__import__('pathlib').Path(__file__).resolve().parent.parent / 'manuscript_draft_v0.md'), encoding='utf-8').read()
checks = {
    '标题（诚实版）': 'Side-Chain-Centroid Contact Graphs Reveal' in txt,
    'Abstract': '## Abstract' in txt,
    'Introduction': '## 1. Introduction' in txt,
    'Methods 3.1-3.7': all(f'### 3.{i}' in txt for i in range(1, 8)),
    'Results 4.1-4.5': all(f'### 4.{i}' in txt for i in range(1, 6)),
    'Discussion 5.1-5.4': all(f'### 5.{i}' in txt for i in range(1, 5)),
    'References': '## References' in txt,
    '邮箱': 'Yilin.Huang24' in txt,
    '代码声明': 'available from the corresponding author' in txt,
    'bootstrap CI 数字': '0.254, 0.529' in txt,
    'protein-level': 'Fisher-z' in txt,
}
for k, v in checks.items():
    print(f'  {k}: {"OK" if v else "MISSING"}')

print()
print('=== 交付物清单 ===')
for f in ['manuscript_draft_v0.docx', 'cover_letter.md', '投稿材料清单.md',
          'requirements.lock', 'data_hashes.md', '审计报告.md', '修改报告.md',
          '统一交集统计结果.md', '数据流审计表.md']:
    p = os.path.join(str(__import__('pathlib').Path(__file__).resolve().parent.parent), f)
    print(f'  {"OK " if os.path.exists(p) else "MISS "}{f}')
