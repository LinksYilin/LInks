"""
paths.py — 项目路径集中管理（回应审查：脚本硬编码 D:\\GED_mutation）
=====================================================================
原则：
  - 默认从本文件位置自动推导项目根目录 → 克隆到任何路径都能直接运行
  - 允许用环境变量 GED_ROOT 覆盖（便于在别处存放数据）
  - 不假设盘符、不假设用户名

用法:
    from paths import DATA, FIGURES, TOOLS, ROOT
    df = pd.read_csv(DATA / 'benchmarks_s669_clean.csv')
"""
import os
from pathlib import Path

# 本文件位于 <ROOT>/src/paths.py → 根目录是上一级
_HERE = Path(__file__).resolve().parent
ROOT = Path(os.environ.get('GED_ROOT', _HERE.parent))

SRC = ROOT / 'src'
DATA = ROOT / 'data'
FIGURES = ROOT / 'figures'
PUB_FIGURES = FIGURES / 'publication'
TOOLS = ROOT / 'tools'
DOCS = ROOT / 'docs'

# 外部程序（可用环境变量覆盖，便于其他机器）
FOLDX_BIN = os.environ.get('FOLDX_BIN', str(TOOLS / 'foldx' / 'foldx_1_20270131.exe'))
SCWRL4_BIN = os.environ.get('SCWRL4_BIN', str(TOOLS / 'scwrl4' / 'Scwrl4.exe'))

# 论文与关键文件
MANUSCRIPT_V2 = ROOT / 'manuscript_revised_20260930_v2.docx'
MANUSCRIPT_ORIG = ROOT / 'manuscript_revised_20260930.docx'


def ensure_dirs():
    for d in (DATA, FIGURES, PUB_FIGURES):
        d.mkdir(parents=True, exist_ok=True)


def describe() -> dict:
    return {
        'ROOT': str(ROOT),
        'DATA': str(DATA),
        'FIGURES': str(PUB_FIGURES),
        'GED_ROOT_env': os.environ.get('GED_ROOT'),
        'FOLDX_exists': Path(FOLDX_BIN).exists(),
        'SCWRL4_exists': Path(SCWRL4_BIN).exists(),
    }


if __name__ == '__main__':
    for k, v in describe().items():
        print(f'{k:16} = {v}')
