"""
data_flow_audit.py — 数据流审计：各阶段样本数 + 共同测试交集
============================================================
统一各方法的测试样本，输出逐样本的跳过记录（pdb_id, mutation, 原因, 阶段, 是否进入最终分析）。
"""
import os

import pandas as pd

from paths import DATA, FIGURES, PUB_FIGURES, ROOT, SRC, TOOLS  # 路径集中管理

SRC_PATH = str(SRC)
DATA_PATH = str(DATA)
ROOT_PATH = str(ROOT)
FIG_PATH = str(FIGURES)
PUB_PATH = str(PUB_FIGURES)
TOOLS_PATH = str(TOOLS)

DATA = r'D:\GED_mutation\data'


def main():
    # 1. 原始 S669
    raw = pd.read_csv(os.path.join(DATA, 'raw', 'benchmarks', 's669.csv'))
    raw['_key'] = raw['pdb_id'] + ':' + raw['mut_info']
    print(f'1. 原始 S669: {len(raw)} 条')

    # 2. 标签审计后（干净集）
    clean = pd.read_csv(os.path.join(DATA, 'benchmarks_s669_clean.csv'))
    clean['_key'] = clean['pdb_id'] + ':' + clean['mut_info']
    print(f'2. 标签审计后（clean）: {len(clean)} 条')

    # 3. 对齐后（有 alignment map 的蛋白）
    align = pd.read_csv(os.path.join(DATA, 'alignment_map_s669.csv'))
    aligned_pids = set(align[align['status'] == 'OK']['pdb_id'])
    clean['_aligned'] = clean['pdb_id'].isin(aligned_pids)
    print(f'3. 序列-结构对齐后: {clean["_aligned"].sum()} 条')

    # 4. 有接触图
    cg_dir = os.path.join(DATA, 'contact_graphs_s669_sc')
    cg_pids = {f.replace('.npz', '') for f in os.listdir(cg_dir) if f.endswith('.npz')}
    clean['_has_cg'] = clean['pdb_id'].isin(cg_pids)
    print(f'4. 有接触图: {clean["_has_cg"].sum()} 条')

    # 5. mutation index 有效（_pdb_res_idx 在范围内）
    def idx_valid(row):
        pid = row['pdb_id']
        if pid not in cg_pids:
            return False
        import numpy as np



        d = np.load(os.path.join(cg_dir, f'{pid}.npz'), allow_pickle=True)
        n = d['nodes'].shape[0]
        idx = row['_pdb_res_idx']
        return 0 <= idx < n
    clean['_idx_valid'] = clean.apply(idx_valid, axis=1)
    print(f'5. mutation index 有效: {clean["_idx_valid"].sum()} 条')

    # 6. 有 FoldX 突变体
    mt_dir = os.path.join(DATA, 'mutant_structures_s669')
    def has_mutant(row):
        p = os.path.join(mt_dir, row['pdb_id'], row['mut_info'], 'work', 'out',
                         f'pdb{row["pdb_id"].lower()}_1.pdb')
        return os.path.exists(p)
    clean['_has_mt'] = clean.apply(has_mutant, axis=1)
    print(f'6. 有 FoldX 突变体: {clean["_has_mt"].sum()} 条')

    # 7. 有真实编辑
    true_edits = pd.read_csv(os.path.join(DATA, 'true_edits_s669_sc.csv'))
    te_keys = set(true_edits['pdb_id'] + ':' + true_edits['mut_info'])
    clean['_has_te'] = clean['_key'].isin(te_keys)
    print(f'7. 有真实编辑: {clean["_has_te"].sum()} 条')

    # 8. 共同测试交集（GNN 用的：有接触图 + idx 有效）
    clean['_common_gnn'] = clean['_has_cg'] & clean['_idx_valid']
    print(f'8. GNN 共同测试交集（有接触图 + idx 有效）: {clean["_common_gnn"].sum()} 条')

    # 输出逐样本跳过记录
    rows = []
    for _, r in clean.iterrows():
        stages = []
        if not r['_aligned']:
            stages.append('对齐失败')
        if not r['_has_cg']:
            stages.append('无接触图')
        if r['_has_cg'] and not r['_idx_valid']:
            stages.append('idx越界')
        if r['_common_gnn'] and not r['_has_mt']:
            stages.append('无FoldX突变体')
        if r['_has_mt'] and not r['_has_te']:
            stages.append('无真实编辑')
        if stages:
            rows.append({
                'pdb_id': r['pdb_id'],
                'mutation': r['mut_info'],
                'reason': '; '.join(stages),
                'in_final_gnn_analysis': r['_common_gnn'],
            })
    skip_df = pd.DataFrame(rows)
    skip_out = os.path.join(DATA, 'data_flow_skip_log.csv')
    skip_df.to_csv(skip_out, index=False)
    print(f'\n跳过记录（非空）: {len(skip_df)} 条，已写入 {skip_out}')


if __name__ == '__main__':
    main()
