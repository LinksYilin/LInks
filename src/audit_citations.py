# -*- coding: utf-8 -*-
"""audit_citations.py — 逐条核实参考文献是否真实存在

对每条文献提取 DOI/arXiv 编号（若有），并通过 Crossref / arXiv API 核验。
无法自动核验的标记为待人工确认。
"""
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request

from docx import Document

DOC = r'D:\GED_mutation\manuscript_routeA.docx'
UA = {'User-Agent': 'manuscript-audit/1.0 (mailto:Yilin.Huang24@student.xjtlu.edu.cn)'}

d = Document(DOC)
refs = {}
for p in d.paragraphs:
    m = re.match(r'^\[(\d+)\]\s+(.+)$', p.text.strip())
    if m:
        refs[int(m.group(1))] = m.group(2)

print(f'参考文献 {len(refs)} 条')


def crossref_title(title, author=None):
    """用标题查 Crossref，返回最匹配记录的标题/年/DOI。"""
    q = urllib.parse.urlencode({'query.bibliographic': title[:300], 'rows': 3})
    url = f'https://api.crossref.org/works?{q}'
    try:
        req = urllib.request.Request(url, headers=UA)
        with urllib.request.urlopen(req, timeout=30) as r:
            j = json.loads(r.read().decode())
    except Exception as e:
        return None, str(e)
    items = j.get('message', {}).get('items', [])
    for it in items:
        t = (it.get('title') or [''])[0]
        yr = None
        for k in ('published-print', 'published-online', 'issued'):
            if it.get(k, {}).get('date-parts'):
                yr = it[k]['date-parts'][0][0]
                break
        return {'title': t, 'year': yr, 'doi': it.get('DOI')}, None
    return None, 'no items'


def arxiv_lookup(arxid):
    url = f'http://export.arxiv.org/api/query?id_list={arxid}'
    try:
        req = urllib.request.Request(url, headers=UA)
        with urllib.request.urlopen(req, timeout=30) as r:
            x = r.read().decode()
    except Exception as e:
        return None, str(e)
    m = re.search(r'<entry>.*?<title>(.*?)</title>', x, re.S)
    if not m:
        return None, 'not found'
    return {'title': re.sub(r'\s+', ' ', m.group(1)).strip()}, None


def norm(s):
    return re.sub(r'[^a-z0-9]', '', (s or '').lower())


results = []
for n in sorted(refs):
    t = refs[n]
    # 抽取标题（第一个句点前的部分通常含标题）
    m = re.match(r'^(.*?)\((\d{4})\)\s*(.+?)\.\s', t)
    title = m.group(3) if m else t[:120]
    year = m.group(2) if m else '?'
    ax = re.search(r'arXiv[,\s]*(\d{4}\.\d{4,5})', t)
    print(f'\n[{n}] {title[:78]}')
    print(f'     声称年份 {year}')
    if ax:
        info, err = arxiv_lookup(ax.group(1))
        if info:
            ok = norm(title)[:40] in norm(info["title"]) or norm(info["title"])[:40] in norm(title)
            print(f'     arXiv {ax.group(1)} -> {info["title"][:70]}')
            print(f'     {"✅ 匹配" if ok else "⚠ 标题需人工核对"}')
            results.append((n, ok))
        else:
            print(f'     arXiv 查询失败: {err}')
            results.append((n, None))
        time.sleep(1)
        continue
    info, err = crossref_title(title)
    if info:
        ok = norm(title)[:38] in norm(info['title']) or norm(info['title'])[:38] in norm(title)
        yr_ok = (str(info['year']) == year) if info['year'] else None
        print(f'     Crossref -> {info["title"][:70]}')
        print(f'     DOI {info["doi"]}  年份 {info["year"]}')
        print(f'     {"✅ 标题匹配" if ok else "⚠ 标题需人工核对"}'
              f'{"" if yr_ok in (True, None) else "  ⚠ 年份不符"}')
        results.append((n, ok))
    else:
        print(f'     Crossref 查询失败: {err}')
        results.append((n, None))
    time.sleep(1)

print()
print('=' * 76)
n_ok = sum(1 for _, v in results if v)
n_bad = sum(1 for _, v in results if v is False)
n_un = sum(1 for _, v in results if v is None)
print(f'引用核实：匹配 {n_ok} / 需人工 {n_bad} / 无法查询 {n_un}  共 {len(results)}')
print('=' * 76)
