# -*- coding: utf-8 -*-
"""trigger_zenodo_archive.py — 创建新 Release 以触发 Zenodo 归档

Zenodo 的 GitHub 集成只归档「新」Release。若集成在 v1.0.0 之后才开启，
需再发一个 Release。默认发 v1.0.1（内容与 v1.0.0 相同）。
"""
import argparse
import json
import subprocess
import sys
import time
import urllib.error
import urllib.request

REPO = 'LinksYilin/LInks'


def token():
    p = subprocess.run(['git', 'credential', 'fill'],
                       input='protocol=https\nhost=github.com\n\n',
                       capture_output=True, text=True)
    return next((l.split('=', 1)[1] for l in p.stdout.splitlines()
                 if l.startswith('password=')), None)


def api(path, method='GET', payload=None, tok=None):
    hdr = {'Authorization': f'Bearer {tok}', 'Accept': 'application/vnd.github+json',
           'User-Agent': 'dsh', 'Content-Type': 'application/json'}
    req = urllib.request.Request(
        'https://api.github.com' + path, method=method,
        data=json.dumps(payload).encode() if payload is not None else None, headers=hdr)
    try:
        with urllib.request.urlopen(req, timeout=30) as r:
            return r.status, json.loads(r.read().decode())
    except urllib.error.HTTPError as e:
        return e.code, e.read().decode()[:300]


def zenodo_doi_query():
    import urllib.parse
    url = ('https://zenodo.org/api/records?size=10&q='
           + urllib.parse.quote('LinksYilin OR "contact-graph definitions"'))
    try:
        with urllib.request.urlopen(url, timeout=30) as r:
            j = json.loads(r.read().decode())
        return j.get('hits', {}).get('hits', [])
    except Exception:
        return []


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--tag', default='v1.0.1')
    ap.add_argument('--wait', type=int, default=180,
                    help='触发后轮询 Zenodo 的秒数')
    args = ap.parse_args()

    tok = token()
    if not tok:
        print('❌ 未取到 GitHub token')
        return 1

    s, j = api(f'/repos/{REPO}/releases', tok=tok)
    if s != 200:
        print(f'❌ 查询 release 失败: {s} {j}')
        return 1
    existing = {r['tag_name'] for r in j}
    print(f'现有 release: {sorted(existing)}')

    if args.tag in existing:
        print(f'{args.tag} 已存在，直接进入轮询')
    else:
        body = (
            'Same content as v1.0.0, released to trigger archival after the Zenodo '
            'GitHub integration was enabled.\n\n'
            'Contact-graph definitions shape what a model describes but not what it '
            'predicts: a controlled comparison across five encoder capacities.')
        s, j = api(f'/repos/{REPO}/releases', 'POST',
                   {'tag_name': args.tag, 'name': args.tag, 'body': body,
                    'draft': False, 'prerelease': False, 'target_commitish': 'main'},
                   tok=tok)
        print(f'创建 {args.tag} -> HTTP {s}')
        if s == 201:
            print(f'  {j["html_url"]}')
        else:
            print(f'  响应: {j}')
            return 1

    print(f'\n轮询 Zenodo（最多 {args.wait}s）…')
    found = []
    for i in range(args.wait // 15):
        time.sleep(15)
        hits = zenodo_doi_query()
        for h in hits:
            doi = h.get('doi', '')
            title = h.get('metadata', {}).get('title', '')[:60]
            if doi and (doi, title) not in found:
                found.append((doi, title))
        if found:
            print('\n🎉 找到归档：')
            for doi, title in found:
                print(f'  DOI: {doi}')
                print(f'  {title}')
            return 0
        print(f'  [{15 * (i + 1)}s] 尚未归档…')

    print('\n⚠ 轮询超时。Zenodo 归档通常需 1–5 分钟；稍后重跑本脚本即可查看。')
    return 2


if __name__ == '__main__':
    sys.exit(main())
