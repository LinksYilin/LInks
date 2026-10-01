# -*- coding: utf-8 -*-
"""Reproduce the original 508934-parameter EGNN at seed 2024.

The current EGNN default gained four coordinate-scale parameters (508938).
This runner explicitly uses scaled=False to match seeds 42 and 123. Never
combine the 508938-parameter candidate with this controlled three-seed ladder.
"""
import os
import time

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from torch_geometric.loader import DataLoader

from ladder_common import load_test, load_train
from seed_utils import set_seed
from strong_backbones import EGNN
from run_ladder import cluster_ci, pearson, spearman
from paths import DATA

DEV = 'cuda' if torch.cuda.is_available() else 'cpu'
SEED = 2024
OUTPUT = os.path.join(str(DATA), 'ladder_egnn_legacy_s2024')


def main():
    rows, predictions = [], []
    for atom_def in ('ca', 'cb', 'centroid'):
        train = load_train(atom_def, with_coords=True)
        test, metadata = load_test(atom_def, 's669', with_coords=True)
        set_seed(SEED)
        model = EGNN(train[0].x.shape[1], edge_dim=2, hid=128,
                     n_layers=4, k_hop=2, update_coords=True, scaled=False).to(DEV)
        n_params = sum(p.numel() for p in model.parameters())
        if n_params != 508934:
            raise AssertionError(f'Expected original EGNN 508934 parameters, got {n_params}')
        optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
        loader = DataLoader(train, batch_size=64, shuffle=True)
        t0 = time.time()
        for _ in range(20):
            model.train()
            for batch in loader:
                batch = batch.to(DEV)
                optimizer.zero_grad()
                F.mse_loss(model(batch), batch.y).backward()
                optimizer.step()
        model.eval()
        with torch.no_grad():
            pred = np.concatenate([model(b.to(DEV)).cpu().numpy()
                                   for b in DataLoader(test, batch_size=128)])
        y = np.array([m['ddg'] for m in metadata], dtype=float)
        pid = [m['pdb_id'] for m in metadata]
        lo, hi = cluster_ci(y, pred, pid)
        rows.append(dict(benchmark='s669', atom_def=atom_def, model='egnn', seed=SEED,
                         n=len(y), n_proteins=len(set(pid)), r=pearson(y, pred),
                         ci_low=lo, ci_high=hi, spearman=spearman(y, pred),
                         mae=float(np.mean(np.abs(y - pred))),
                         rmse=float(np.sqrt(np.mean((y - pred) ** 2))),
                         n_params=n_params, train_seconds=time.time()-t0))
        predictions.extend(dict(benchmark='s669', atom_def=atom_def, model='egnn', seed=SEED,
                                protein_id=pid[i], mutation_id=metadata[i]['mut_info'],
                                y_true=y[i], y_pred=float(v)) for i, v in enumerate(pred))
        pd.DataFrame(rows).to_csv(OUTPUT+'_results.csv', index=False)
        pd.DataFrame(predictions).to_csv(OUTPUT+'_predictions.csv', index=False)
        print(f'{atom_def}: n={len(y)}, r={rows[-1]["r"]:.4f}, params={n_params}, '
              f'time={rows[-1]["train_seconds"]:.0f}s', flush=True)
    print(f'Wrote {len(rows)} architecture-matched seed results', flush=True)


if __name__ == '__main__':
    main()
