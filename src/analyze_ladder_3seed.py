# -*- coding: utf-8 -*-
"""Audit the capacity ladder using complete, architecture-matched seed sets.

For each definition, first align mutations, then average predictions over common
seeds. Two intervals are reported: a protein-cluster bootstrap on the fixed seed
ensemble (sampling uncertainty), and a joint seed/protein bootstrap (also reflects
optimisation variability). The latter is descriptive with only three seeds.
EGNN seeds from different parameter counts are NEVER combined.
"""
import os
from itertools import permutations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from paths import DATA

D = str(DATA)
KEY = ['benchmark', 'atom_def', 'model', 'seed', 'protein_id', 'mutation_id']
STABLE = ['gnn_global', 'gnn_local', 'gnn_edge', 'deep_gine']


def corr(y, p):
    if len(y) < 3 or np.std(y) == 0 or np.std(p) == 0:
        return np.nan
    return float(np.corrcoef(y, p)[0, 1])


def read_inputs():
    r0 = pd.read_csv(os.path.join(D, 'ladder_results.csv'))
    r1 = pd.read_csv(os.path.join(D, 'ladder_s2024_results.csv'))
    r2 = pd.read_csv(os.path.join(D, 'ladder_egnn_legacy_s2024_results.csv'))
    p0 = pd.read_csv(os.path.join(D, 'ladder_predictions.csv'))
    p1 = pd.read_csv(os.path.join(D, 'ladder_s2024_predictions.csv'))
    p2 = pd.read_csv(os.path.join(D, 'ladder_egnn_legacy_s2024_predictions.csv'))
    # Drop the mismatched 508938-parameter EGNN before adding the corrected seed.
    r1 = r1[~r1.model.eq('egnn')]
    p1 = p1[~p1.model.eq('egnn')]
    res = pd.concat([r0, r1, r2], ignore_index=True)
    pred = pd.concat([p0, p1, p2], ignore_index=True)
    if res.duplicated(['benchmark', 'atom_def', 'model', 'seed']).any() or pred.duplicated(KEY).any():
        raise ValueError('Duplicate result or prediction keys; resolve before aggregation')
    for m, g in res.groupby('model'):
        versions = sorted(g.n_params.dropna().unique())
        if len(versions) > 1:
            print(f'ARCHITECTURE MISMATCH: {m}: {versions}; excluding incompatible seeds')
            keep = int(r0[r0.model == m].n_params.iloc[0])
            bad = g[g.n_params != keep]
            pred = pred.merge(bad[['benchmark', 'atom_def', 'model', 'seed']].assign(exclude=1),
                              how='left', on=['benchmark', 'atom_def', 'model', 'seed'])
            pred = pred[pred.exclude.isna()].drop(columns='exclude')
            res = res[~res.index.isin(bad.index)]
    return res, pred


def aligned_arrays(sub, alt, ref):
    alt_g = sub[sub.atom_def == alt]
    ref_g = sub[sub.atom_def == ref]
    a_seeds = set(alt_g.seed.unique()); b_seeds = set(ref_g.seed.unique())
    seeds = sorted(a_seeds & b_seeds)
    if not seeds:
        raise ValueError(f'No common seeds {alt}/{ref}')
    index = ['protein_id', 'mutation_id']
    a = alt_g[alt_g.seed.isin(seeds)].pivot(index=index, columns='seed', values='y_pred')
    b = ref_g[ref_g.seed.isin(seeds)].pivot(index=index, columns='seed', values='y_pred')
    common = a.index.intersection(b.index)
    a = a.loc[common, seeds]; b = b.loc[common, seeds]
    if a.isna().any().any() or b.isna().any().any():
        raise ValueError(f'Incomplete seed x mutation matrix {alt}/{ref}')
    yy = alt_g.groupby(index).y_true.first().loc[common]
    yr = ref_g.groupby(index).y_true.first().loc[common]
    if not np.allclose(yy.to_numpy(), yr.to_numpy()):
        raise ValueError('Labels mismatch across definitions')
    return (yy.to_numpy(), a.to_numpy(), b.to_numpy(),
            np.asarray(common.get_level_values('protein_id')), seeds)


def paired_bootstrap(y, a, b, protein, n_boot=2000, seed=101):
    rng = np.random.default_rng(seed)
    labels, inv = np.unique(protein, return_inverse=True)
    groups = [np.flatnonzero(inv == j) for j in range(len(labels))]
    point = corr(y, a.mean(axis=1)) - corr(y, b.mean(axis=1))
    fixed, joint = [], []
    for _ in range(n_boot):
        ix = np.concatenate([groups[j] for j in rng.integers(len(groups), size=len(groups))])
        fixed.append(corr(y[ix], a[ix].mean(axis=1)) - corr(y[ix], b[ix].mean(axis=1)))
        ss = rng.integers(a.shape[1], size=a.shape[1])
        joint.append(corr(y[ix], a[ix][:, ss].mean(axis=1)) -
                     corr(y[ix], b[ix][:, ss].mean(axis=1)))
    fixed = np.asarray(fixed); joint = np.asarray(joint)
    if not np.isfinite(fixed).all() or not np.isfinite(joint).all():
        raise ValueError('Non-finite bootstrap differences')
    p_fixed = float(min(1, 2 * min(np.mean(fixed <= 0), np.mean(fixed >= 0))))
    p_joint = float(min(1, 2 * min(np.mean(joint <= 0), np.mean(joint >= 0))))
    return (point, np.percentile(fixed, [2.5, 97.5]),
            np.percentile(joint, [2.5, 97.5]), p_fixed, p_joint)


def main():
    res, pred = read_inputs()
    rows = []
    for m, g in res.groupby('model'):
        for ad, cell in g.groupby('atom_def'):
            if cell.seed.nunique() < 2:
                continue
            rows.append(dict(model=m, atom_def=ad, n_params=int(cell.n_params.iloc[0]),
                             n_seeds=cell.seed.nunique(), seed_mean_r=float(cell.r.mean()),
                             seed_sd_r=float(cell.r.std(ddof=1))))
    summary = pd.DataFrame(rows)
    summary.to_csv(os.path.join(D, 'ladder_seed_summary_audited.csv'), index=False)
    effects = []
    for m, sub in pred[pred.benchmark == 's669'].groupby('model'):
        if 'centroid' not in set(sub.atom_def):
            continue
        for ad in sorted(set(sub.atom_def) - {'centroid'}):
            y, a, b, prot, seeds = aligned_arrays(sub, ad, 'centroid')
            point, ci_fixed, ci_joint, p_fixed, p_joint = paired_bootstrap(y, a, b, prot)
            effects.append(dict(model=m, n_params=int(res[res.model == m].n_params.iloc[0]),
                                comparison=f'{ad} vs centroid', n_mutations=len(y),
                                n_proteins=len(set(prot)), n_seeds=len(seeds),
                                seeds=','.join(map(str, seeds)), delta_r=point,
                                fixed_lo=ci_fixed[0], fixed_hi=ci_fixed[1], fixed_p=p_fixed,
                                seed_aware_lo=ci_joint[0], seed_aware_hi=ci_joint[1],
                                seed_aware_p=p_joint))
    out = pd.DataFrame(effects)
    out.to_csv(os.path.join(D, 'ladder_paired_effects_audited.csv'), index=False)
    print('SEED SUMMARY\n', summary.round(4).to_string(index=False))
    print('\nPAIRED EFFECTS\n', out.round(4).to_string(index=False))
    if len(out):
        ordered = out.sort_values('seed_aware_p').copy()
        ordered['holm_p'] = np.minimum(1, ordered.seed_aware_p.to_numpy() *
                                      (len(ordered) - np.arange(len(ordered))))
        ordered['holm_p'] = np.maximum.accumulate(ordered.holm_p)
        ordered.to_csv(os.path.join(D, 'ladder_paired_effects_audited_holm.csv'), index=False)
        print('\nMultiplicity (Holm across all 13):\n',
              ordered[['model', 'comparison', 'seed_aware_p', 'holm_p']].round(4).to_string(index=False))
    stable = out[out.model.isin(STABLE)]
    # A capacity test needs a common definition-pair set; CA and CB exist at every rung.
    common = stable[stable.comparison.isin(['ca vs centroid', 'cb vs centroid'])]
    trend = common.groupby(['model', 'n_params']).delta_r.apply(lambda x: np.abs(x).mean()).reset_index(name='mean_abs_dr')
    trend = trend.sort_values('n_params')
    x = np.log10(trend.n_params.to_numpy()); vals = trend.mean_abs_dr.to_numpy()
    obs = float(spearmanr(x, vals).statistic)
    null = np.array([spearmanr(x, p).statistic for p in permutations(vals)])
    p_two = float(np.mean(np.abs(null) >= abs(obs) - 1e-12))
    print('\nCOMMON-DEFINITION TREND\n', trend.round(4).to_string(index=False))
    print(f'Spearman rho={obs:+.3f}, exact two-sided permutation p={p_two:.3f}, {len(null)} permutations')
    print('Caveat: four capacities and no pre-registered trend test; do not infer equivalence.')
    trend.assign(rho=obs, exact_p=p_two).to_csv(os.path.join(D, 'ladder_capacity_trend_audited.csv'), index=False)


if __name__ == '__main__':
    main()
