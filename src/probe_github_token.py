# -*- coding: utf-8 -*-
"""probe_github_token.py — 只读探测 GitHub token 的权限范围"""
import json
import subprocess
import urllib.error
import urllib.request

# 先取出 token
inp = "protocol=https\nhost=github.com\n\n"
p = subprocess.run(['git', 'credential', 'fill'], input=inp, capture_output=True, text=True)
token = None
for line in p.stdout.splitlines():
    if line.startswith('password='):
        token = line.split('=', 1)[1]
if not token:
    print('未取到 token')
    raise SystemExit(1)
print(f'token 长度 {len(token)}，前缀 {token[:4]}…')

hdr = {'Authorization': f'Bearer {token}', 'Accept': 'application/vnd.github+json',
       'User-Agent': 'dsh'}


def api(path):
    req = urllib.request.Request('https://api.github.com' + path, headers=hdr)
    try:
        with urllib.request.urlopen(req, timeout=20) as r:
            return r.status, json.loads(r.read().decode()), dict(r.headers)
    except urllib.error.HTTPError as e:
        return e.code, None, dict(e.headers)
    except Exception as e:
        return None, str(e), {}


s, j, h = api('/user')
print(f'\n/user -> {s}')
if j:
    print(f'  登录名: {j.get("login")}')

print(f'  token scopes: {h.get("X-OAuth-Scopes", "(未提供)")}')
print(f'  accepted scopes: {h.get("X-Accepted-OAuth-Scopes", "(未提供)")}')

s2, j2, _ = api('/repos/LinksYilin/LInks')
print(f'\n/repos/LinksYilin/LInks -> {s2}')
if j2:
    print(f'  permissions: {j2.get("permissions")}')
    print(f'  visibility: {j2.get("visibility")}')
    print(f'  默认分支: {j2.get("default_branch")}')

s3, j3, _ = api('/repos/LinksYilin/LInks/releases')
print(f'\n/releases -> {s3}  已有 release 数: {len(j3) if isinstance(j3, list) else "?"}')
if isinstance(j3, list):
    for r in j3[:3]:
        print(f'    tag={r.get("tag_name")}  name={r.get("name")}')
