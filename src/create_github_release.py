# -*- coding: utf-8 -*-
"""create_github_release.py — 用缓存凭据创建 v1.0.0 Release

Release 本身是可引用版本，也是 Zenodo GitHub 集成的触发点。
"""
import json
import subprocess
import urllib.error
import urllib.request

REPO = 'LinksYilin/LInks'
TAG = 'v1.0.0'
TITLE = 'v1.0.0 — Contact-graph definitions and ΔΔG prediction'
BODY = """First archived release accompanying the manuscript:

**Contact-graph definitions shape what a model describes but not what it predicts: a controlled comparison across five encoder capacities**

Contents
- `src/` — 148 analysis scripts: graph construction, the five-rung capacity ladder, ESM-2 baselines and fusion, statistics and verification
- `results/` — 85 result tables, including per-sample predictions and the audited seed-aware aggregates
- `figures/` — 5 publication figures (PNG/PDF) with source data
- `docs/` — 25 documents, including the audit trail for two corrected pipeline defects (a hydrogen inconsistency in training graphs, and an in-sample fit in an exploratory script) and the submission materials

Notes
- Inputs (PDB structures, MegaScale, ThermoMutDB, S669, ssym) are public and are not redistributed; the pipeline regenerates every graph and result from them.
- The pipeline resolves its root from its own location, so it runs after cloning to any directory. Set `GED_ROOT` to keep data elsewhere.
- Reproducibility caveats, including the locality statistic's definition dependence, are documented in the supplementary material.
"""

inp = 'protocol=https\nhost=github.com\n\n'
p = subprocess.run(['git', 'credential', 'fill'], input=inp, capture_output=True, text=True)
token = next((l.split('=', 1)[1] for l in p.stdout.splitlines()
              if l.startswith('password=')), None)
if not token:
    raise SystemExit('未取到 token')
HDR = {'Authorization': f'Bearer {token}', 'Accept': 'application/vnd.github+json',
       'User-Agent': 'dsh', 'Content-Type': 'application/json'}


def api(path, method='GET', payload=None):
    req = urllib.request.Request(
        'https://api.github.com' + path, method=method,
        data=json.dumps(payload).encode() if payload is not None else None, headers=HDR)
    try:
        with urllib.request.urlopen(req, timeout=30) as r:
            return r.status, json.loads(r.read().decode())
    except urllib.error.HTTPError as e:
        return e.code, e.read().decode()[:400]


# 1) 已有 release？
s, j = api(f'/repos/{REPO}/releases')
if isinstance(j, list) and any(r.get('tag_name') == TAG for r in j):
    print(f'{TAG} 已存在，跳过')
    raise SystemExit(0)
print(f'现有 release: {len(j) if isinstance(j, list) else j}')

# 2) 创建 tag + release
s, j = api(f'/repos/{REPO}/releases', 'POST', {
    'tag_name': TAG, 'name': TITLE, 'body': BODY,
    'draft': False, 'prerelease': False, 'target_commitish': 'main'})
print(f'\n创建 release -> HTTP {s}')
if s == 201:
    print(f'  tag: {j["tag_name"]}')
    print(f'  url: {j["html_url"]}')
    print(f'  作者: {j["author"]["login"]}')
else:
    print(f'  响应: {j}')
