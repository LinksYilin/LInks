# -*- coding: utf-8 -*-
"""clean_orphan_media.py — 移除 DOCX 中未被任何关系引用的媒体文件

安全做法：重建 zip，跳过孤立媒体条目，并保留其余所有条目原样。
重建后校验：段落数、表格数、内嵌图数、可打开性均不变。
"""
import hashlib
import os
import re
import shutil
import zipfile

from docx import Document
from docx.oxml.ns import qn

DOC = r'D:\GED_mutation\manuscript_routeA.docx'
BAK = r'D:\GED_mutation\manuscript_routeA_predocxclean.docx'


def snapshot(path):
    d = Document(path)
    with zipfile.ZipFile(path) as z:
        names = z.namelist()
        docxml = z.read('word/document.xml').decode('utf-8')
    return {
        'paras': len(d.paragraphs),
        'tables': len(d.tables),
        'shapes': len(d.inline_shapes),
        'media': sorted(n for n in names if n.startswith('word/media/')),
        'sha_doc': hashlib.sha256(docxml.encode()).hexdigest(),
    }


def main():
    before = snapshot(DOC)
    print('清理前:', {k: v for k, v in before.items() if k != 'media'})
    print(f'  媒体 {len(before["media"])} 个')

    with zipfile.ZipFile(DOC) as z:
        names = z.namelist()
        docxml = z.read('word/document.xml').decode('utf-8')
        rels = z.read('word/_rels/document.xml.rels').decode('utf-8')
        data = {n: z.read(n) for n in names}

    used = set(re.findall(r'r:embed="(rId\d+)"', docxml))
    rid2target = dict(re.findall(r'Id="(rId\d+)"[^>]*Target="([^"]+)"', rels))
    used_media = {rid2target[r].replace('media/', 'word/media/') for r in used if r in rid2target}
    orphan = sorted(n for n in names if n.startswith('word/media/') and n not in used_media)
    print(f'  孤立媒体: {[os.path.basename(o) for o in orphan]}')

    if not orphan:
        print('无需清理')
        return 0

    # 同时删除指向孤立媒体的 Relationship 条目
    orphan_rids = [rid for rid, tgt in rid2target.items()
                   if tgt.replace('media/', 'word/media/') in orphan
                   and (tgt.replace('media/', 'word/media/') not in used_media)]
    new_rels = rels
    for rid in orphan_rids:
        new_rels = re.sub(r'<Relationship[^>]*Id="' + rid + r'"[^>]*/>', '', new_rels)
    print(f'  删除关系条目 {len(orphan_rids)} 个')

    shutil.copy2(DOC, BAK)
    tmp = DOC + '.tmp'
    with zipfile.ZipFile(tmp, 'w', zipfile.ZIP_DEFLATED) as out:
        for n in names:
            if n in orphan:
                continue
            out.writestr(n, new_rels.encode() if n.endswith('document.xml.rels') else data[n])
    os.replace(tmp, DOC)

    after = snapshot(DOC)
    print('清理后:', {k: v for k, v in after.items() if k != 'media'})
    print(f'  媒体 {len(after["media"])} 个')

    ok = (before['paras'] == after['paras'] and before['tables'] == after['tables']
          and before['shapes'] == after['shapes']
          and before['sha_doc'] == after['sha_doc'])
    print()
    print('✅ 校验通过：段落/表格/图片/document.xml 均未变' if ok else '❌ 校验失败，已保留备份')
    if not ok:
        shutil.copy2(BAK, DOC)
        return 1
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
