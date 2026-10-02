# Per (cluster, dissection): mean log1p(counts per 10k) and fraction expressing, for groups of >= 20 nuclei
import sys, time, numpy as np, pandas as pd
sys.path.insert(0,'.')
from cxopen import open_c
import tiledbsoma as soma
c = open_c()
exp = c["census_data"]["homo_sapiens"]
var = pd.read_parquet('var.parquet'); G = int(var.soma_joinid.max()) + 1
cl = pd.read_parquet('cl_map.parquet').set_index('observation_joinid').cluster
dsets = pd.read_csv('hipamy_datasets.csv')
NC = int(cl.max()) + 1
import os
os.makedirs('agg3', exist_ok=True)
for k, row in dsets.iterrows():
    t0 = time.time()
    if os.path.exists(f'agg3/{k}.npz'): continue
    keys, means, fracs = [], [], []
    with exp.axis_query("RNA", obs_query=soma.AxisQuery(value_filter=f"dataset_id == '{row.dataset_id}'")) as q:
        obs = q.obs(column_names=["soma_joinid","observation_joinid","raw_sum"]).concat().to_pandas()
        obs['cl'] = obs.observation_joinid.map(cl).fillna(-1).astype(int)
        n = np.bincount(obs.cl[obs.cl >= 0], minlength=NC)
        cell_group = pd.Series(obs.cl.values, index=obs.soma_joinid.values)
        cell_size = pd.Series(obs.raw_sum.values, index=obs.soma_joinid.values)
        S = np.zeros(NC * G); F = np.zeros(NC * G)
        for tbl in q.X("raw").tables():
            ci = tbl["soma_dim_0"].to_numpy(); gj = tbl["soma_dim_1"].to_numpy(); v = tbl["soma_data"].to_numpy()
            grp = cell_group.reindex(ci).to_numpy(); ok = grp >= 0
            val = np.log1p(v[ok] / cell_size.reindex(ci[ok]).to_numpy() * 1e4)
            idx = grp[ok].astype(np.int64) * G + gj[ok]
            S += np.bincount(idx, weights=val, minlength=NC * G)
            F += np.bincount(idx, minlength=NC * G)
        S = S.reshape(NC, G); F = F.reshape(NC, G)
        for cc in np.where(n >= 20)[0]:
            keys.append((int(cc), k, int(n[cc]))); means.append((S[cc] / n[cc]).astype(np.float16)); fracs.append((F[cc] / n[cc]).astype(np.float16))
        del S, F
    np.savez(f'agg3/{k}.npz', keys=np.array(keys), means=np.stack(means), fracs=np.stack(fracs))
    print(k + 1, len(dsets), len(keys), f'{time.time()-t0:.0f}s', flush=True)
parts = [np.load(f'agg3/{k}.npz') for k in range(len(dsets))]
np.savez('agg3.npz', keys=np.concatenate([p['keys'] for p in parts]), means=np.concatenate([p['means'] for p in parts]), fracs=np.concatenate([p['fracs'] for p in parts]))
print('done', flush=True)
