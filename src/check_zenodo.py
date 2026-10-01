# -*- coding: utf-8 -*-
"""check_zenodo.py — 查询 Zenodo 是否已归档本仓库"""
import json
import urllib.parse
import urllib.request

QUERIES = ['LinksYilin', 'LInks contact-graph', 'contact-graph definitions encoder capacities']
for q in QUERIES:
    url = ('https://zenodo.org/api/records?size=5&q=' + urllib.parse.quote(q))
    try:
        with urllib.request.urlopen(url, timeout=25) as r:
            j = json.loads(r.read().decode())
    except Exception as e:
        print(f'[{q}] 查询失败: {e}')
        continue
    hits = j.get('hits', {}).get('hits', [])
    total = j.get('hits', {}).get('total', 0)
    print(f'[{q}] 命中 {total}')
    for h in hits[:3]:
        md = h.get('metadata', {})
        print(f'    {h.get("doi")}  {md.get("title","")[:70]}')
    print()
