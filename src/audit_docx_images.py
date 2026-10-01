# -*- coding: utf-8 -*-
"""audit_docx_images.py — 定位文档中每个图片段落及其对应媒体文件"""
import hashlib
import os
import zipfile

from docx import Document
from docx.text.paragraph import Paragraph
from docx.oxml.ns import qn

DOC = r'D:\GED_mutation\manuscript_routeA.docx'
FIG = r'D:\GED_mutation\figures\publication'
R = '{http://schemas.openxmlformats.org/officeDocument/2006/relationships}'

want = {}
for f in ['figure1_study_design.png', 'figure2_contact_rewiring.png',
          'figure3_capacity_ladder.png', 'figure4_information_budget.png',
          'figure5_contact_diagnostics.png']:
    fp = os.path.join(FIG, f)
    if os.path.exists(fp):
        want[hashlib.sha256(open(fp, 'rb').read()).hexdigest()] = f

with zipfile.ZipFile(DOC) as z:
    part = z.read('word/document.xml').decode('utf-8')
    # 用 python-docx 建立 rId -> partname 映射
    doc = Document(DOC)
    rels = doc.part.rels

    def rid_of(blip):
        return blip.get(R + 'embed')

    print(f'{"段":>4}  {"rId":<8} {"媒体":<20} {"匹配":<30} 上下文')
    for i, p in enumerate(doc.paragraphs):
        blips = p._p.findall('.//' + qn('a:blip'))
        for b in blips:
            rid = rid_of(b)
            target = rels[rid].target_ref if rid in rels else '?'
            media = target.split('/')[-1]
            try:
                data = rels[rid].target_part.blob
                sha = hashlib.sha256(data).hexdigest()
            except Exception:
                sha = ''
            match = want.get(sha, 'STALE/DUPLICATE')
            ctx = p.text.strip()[:40]
            print(f'{i:>4}  {rid:<8} {media:<20} {match:<30} {ctx}')
