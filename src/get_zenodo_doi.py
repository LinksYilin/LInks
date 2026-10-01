# -*- coding: utf-8 -*-
"""get_zenodo_doi.py — 查询 Zenodo 记录的版本 DOI 与概念 DOI"""
import json
import urllib.request

for q in ['LinksYilin', 'contact-graph definitions encoder capacities']:
    url = 'https://zenodo.org/api/records?size=10&q=' + urllib.parse.quote(q) if False else \
          'https://zenodo.org/api/records?size=10&q=' + __import__('urllib.parse', fromlist=['x']).quote(q)
    try:
        with urllib.request.urlopen(url, timeout=30) as r:
            j = json.loads(r.read().decode())
    except Exception as e:
        print(f'[{q}] 失败: {e}')
        continue
    for h in j.get('hits', {}).get('hits', []):
        md = h.get('metadata', {})
        title = md.get('title', '')
        if 'contact-graph' not in title.lower():
            continue
        print('=' * 72)
        print(f"版本 DOI   : {h.get('doi')}")
        print(f"概念 DOI   : {h.get('conceptdoi')}")
        print(f"概念记录ID : {h.get('conceptrecid')}")
        print(f"版本       : {md.get('version')}")
        print(f"发布时间   : {md.get('publication_date')}")
        print(f"许可证     : {md.get('license', {}).get('id')}")
        print(f"标题       : {title}")
        print(f"作者       : {[c.get('name') for c in md.get('creators', [])]}")
        print(f"记录 URL   : {h.get('links', {}).get('self_html')}")
        files = h.get('files', [])
        print(f"文件       : {len(files)} 个")
        for f in files[:5]:
            print(f"   {f.get('key')}  {f.get('size', 0)/1e6:.1f} MB")
        print()
