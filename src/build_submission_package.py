# -*- coding: utf-8 -*-
"""build_submission_package.py — 生成可直接上传的投稿包

BMC Bioinformatics 投稿通常需要：
  1. 主稿件（DOCX 或 PDF）
  2. 补充材料（单独文件，DOCX/PDF 均可）
  3. Cover letter（在线表单粘贴或单独文件）
本脚本把 Markdown 材料转为 DOCX，并整理到 submission/ 目录。
"""
import os
import re
import shutil

from docx import Document
from docx.shared import Pt, Inches
from docx.enum.text import WD_ALIGN_PARAGRAPH

ROOT = r'D:\GED_mutation'
OUT = os.path.join(ROOT, 'submission')


def md_to_docx(md_path, docx_path, title=None):
    """把简单的 Markdown（标题/段落/表格/列表）转为 DOCX。"""
    lines = open(md_path, encoding='utf-8').read().split('\n')
    d = Document()
    st = d.styles['Normal']
    st.font.name = 'Times New Roman'
    st.font.size = Pt(11)

    i = 0
    while i < len(lines):
        line = lines[i].rstrip()
        # 表格
        if line.startswith('|') and i + 1 < len(lines) and set(lines[i + 1].replace('|', '').strip()) <= set('-: '):
            header = [c.strip() for c in line.strip('|').split('|')]
            rows = []
            i += 2
            while i < len(lines) and lines[i].startswith('|'):
                rows.append([c.strip() for c in lines[i].strip('|').split('|')])
                i += 1
            t = d.add_table(rows=1, cols=len(header))
            t.style = 'Table Grid'
            for j, h in enumerate(header):
                cell = t.rows[0].cells[j]
                cell.text = ''
                run = cell.paragraphs[0].add_run(re.sub(r'\*\*(.+?)\*\*', r'\1', h))
                run.bold = True
            for r in rows:
                cells = t.add_row().cells
                for j, c in enumerate(r[:len(header)]):
                    cells[j].text = re.sub(r'\*\*(.+?)\*\*', r'\1', c)
            d.add_paragraph()
            continue
        # 标题
        m = re.match(r'^(#{1,4})\s+(.*)$', line)
        if m:
            lvl = len(m.group(1))
            txt = re.sub(r'\*\*(.+?)\*\*', r'\1', m.group(2))
            d.add_heading(txt, level=min(lvl, 3))
            i += 1
            continue
        # 列表
        m = re.match(r'^(\s*)[-*]\s+(.*)$', line)
        if m:
            txt = re.sub(r'\*\*(.+?)\*\*', r'\1', m.group(2))
            d.add_paragraph(txt, style='List Bullet')
            i += 1
            continue
        m = re.match(r'^(\s*)\d+\.\s+(.*)$', line)
        if m:
            txt = re.sub(r'\*\*(.+?)\*\*', r'\1', m.group(2))
            d.add_paragraph(txt, style='List Number')
            i += 1
            continue
        # 普通段落
        if line.strip() and not line.strip().startswith('>'):
            txt = re.sub(r'\*\*(.+?)\*\*', r'\1', line)
            d.add_paragraph(txt)
        i += 1

    if title:
        # 在开头插入标题
        first = d.paragraphs[0]
        h = d.add_heading(title, level=0)
        first._p.addprevious(h._p)
    d.save(docx_path)
    return docx_path


def main():
    os.makedirs(OUT, exist_ok=True)

    # 1) 主稿件
    shutil.copy2(os.path.join(ROOT, 'manuscript_routeA.docx'),
                 os.path.join(OUT, '01_manuscript.docx'))
    print('  ✅ 01_manuscript.docx')

    # 2) 补充材料 -> DOCX
    md_to_docx(os.path.join(ROOT, '补充材料_Supplementary.md'),
               os.path.join(OUT, '02_supplementary.docx'),
               title='Supplementary Material')
    print('  ✅ 02_supplementary.docx')

    # 3) 投稿信
    md_to_docx(os.path.join(ROOT, '投稿信_CoverLetter.md'),
               os.path.join(OUT, '03_cover_letter.docx'),
               title='Cover Letter')
    print('  ✅ 03_cover_letter.docx')

    # 4) 标题页
    md_to_docx(os.path.join(ROOT, '标题页与投稿清单_TitlePage_Checklist.md'),
               os.path.join(OUT, '04_title_page.docx'),
               title='Title Page')
    print('  ✅ 04_title_page.docx')

    # 5) 图（单独文件，部分期刊要求）
    figdir = os.path.join(OUT, 'figures')
    os.makedirs(figdir, exist_ok=True)
    pd = os.path.join(ROOT, 'figures', 'publication')
    for n, f in [(1, 'figure1_study_design.png'), (2, 'figure2_contact_rewiring.png'),
                 (3, 'figure3_capacity_ladder.png'), (4, 'figure4_information_budget.png'),
                 (5, 'figure5_contact_diagnostics.png')]:
        src = os.path.join(pd, f)
        if os.path.exists(src):
            shutil.copy2(src, os.path.join(figdir, f'Figure{n}.png'))
    print(f'  ✅ figures/ (5 个 PNG)')

    # 6) 说明文件
    with open(os.path.join(OUT, 'SUBMISSION_README.txt'), 'w', encoding='utf-8') as f:
        f.write(
            'SUBMISSION PACKAGE\n'
            '==================\n\n'
            'Manuscript title:\n'
            '  Contact-graph definitions shape what a model describes but not what it\n'
            '  predicts: a controlled comparison across five encoder capacities\n\n'
            'Target journal: BMC Bioinformatics\n\n'
            'Files\n'
            '-----\n'
            '  01_manuscript.docx      main manuscript (17 pages, 5 figures, 1 table)\n'
            '  02_supplementary.docx   supplementary material (S1-S10)\n'
            '  03_cover_letter.docx    cover letter\n'
            '  04_title_page.docx      title page and submission checklist\n'
            '  figures/                the five figures as separate PNG files (336-450 dpi)\n\n'
            'Code and data\n'
            '-------------\n'
            '  GitHub: https://github.com/LinksYilin/LInks\n'
            '  Zenodo DOI: 10.5281/zenodo.23086956 (concept, always latest)\n'
            '             10.5281/zenodo.23086957 (v1.0.1)\n\n'
            'Author\n'
            '------\n'
            '  Yilin Huang, Xi\'an Jiaotong-Liverpool University\n'
            '  Yilin.Huang24@student.xjtlu.edu.cn (corresponding)\n\n'
            'Upload order suggested\n'
            '----------------------\n'
            '  1. 01_manuscript.docx as the main document\n'
            '  2. 02_supplementary.docx as supplementary material\n'
            '  3. figures/ if the system asks for figures separately\n'
            '  4. paste 03_cover_letter.docx content into the cover-letter field\n'
            '  5. copy author details from 04_title_page.docx\n')
    print('  ✅ SUBMISSION_README.txt')

    print()
    print(f'投稿包: {OUT}')
    total = sum(os.path.getsize(os.path.join(dp, f))
                for dp, _, fs in os.walk(OUT) for f in fs)
    print(f'总体积 {total/1024/1024:.1f} MB')
    for f in sorted(os.listdir(OUT)):
        p = os.path.join(OUT, f)
        if os.path.isfile(p):
            print(f'  {f:<28} {os.path.getsize(p)/1024:>8.0f} KB')


if __name__ == '__main__':
    main()
