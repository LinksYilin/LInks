"""
download_pdbs.py — 批量下载 PDB 结构并校验
============================================
从 RCSB PDB 按 ID 列表下载结构文件，逐文件校验完整性（%PDF→应为 PDB 文本头、能否被 BioPython 解析）。

用法：
  python download_pdbs.py --ids 1BFM,1D5G,... --out_dir data/structures/ --format pdb

依赖：biopython（已装）。
"""
import argparse
import os
import time

try:
    from Bio.PDB import PDBList, PDBParser
except ImportError:
    raise SystemExit("请先安装 biopython: pip install biopython")


def download_one(pdbl, pdb_id, out_dir, file_format='pdb', max_retries=3):
    """下载单个结构，返回本地文件路径或 None。"""
    for attempt in range(max_retries):
        try:
            path = pdbl.retrieve_pdb_file(
                pdb_id, pdir=out_dir, file_format=file_format, overwrite=False)
            return path
        except Exception as e:
            print(f"  [重试 {attempt+1}] {pdb_id}: {e}")
            time.sleep(2)
    return None


def validate_pdb(path, pdb_id):
    """校验 PDB 文件完整性：文件存在 + 非空 + 能被 BioPython 解析出至少一个残基。"""
    if path is None or not os.path.exists(path):
        return False, "文件不存在"
    size = os.path.getsize(path)
    if size == 0:
        return False, "文件为空"
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure(pdb_id, path)
        model = structure[0]
        n_res = sum(1 for _ in model.get_residues())
        if n_res == 0:
            return False, "解析出 0 残基"
        return True, f"OK ({n_res} 残基, {size} 字节)"
    except Exception as e:
        return False, f"解析失败: {e}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--ids', required=True, help='逗号分隔的 PDB ID（4字符）')
    ap.add_argument('--out_dir', required=True)
    ap.add_argument('--format', default='pdb', choices=['pdb', 'mmCif', 'cif'])
    args = ap.parse_args()

    ids = [i.strip() for i in args.ids.split(',') if i.strip()]
    os.makedirs(args.out_dir, exist_ok=True)
    pdbl = PDBList()

    ok, fail = [], []
    for pdb_id in ids:
        path = download_one(pdbl, pdb_id, args.out_dir, args.format)
        valid, msg = validate_pdb(path, pdb_id)
        if valid:
            ok.append(pdb_id)
            print(f"  ✅ {pdb_id}: {msg}")
        else:
            fail.append((pdb_id, msg))
            print(f"  ❌ {pdb_id}: {msg}")

    print(f"\n=== 汇总：成功 {len(ok)}/{len(ids)}，失败 {len(fail)} ===")
    if fail:
        print("失败清单:")
        for pid, msg in fail:
            print(f"  {pid}: {msg}")


if __name__ == '__main__':
    main()
