import sys, time, numpy as np, pandas as pd
sys.path.insert(0,'.')
from cxopen import open_c
import tiledbsoma as soma
c = open_c()
exp = c["census_data"]["homo_sapiens"]
var = exp.ms["RNA"].var.read(column_names=["soma_joinid","feature_id","feature_name"]).concat().to_pandas()
G = int(var.soma_joinid.max()) + 1
cl = pd.read_parquet('cl_map.parquet').set_index('observation_joinid').cluster
dsets = pd.read_csv('hipamy_datasets.csv')
NC = int(cl.max()) + 1
S = np.zeros(NC * G); F = np.zeros(NC * G)
counts = np.zeros((NC, len(dsets)), dtype=np.int64)
for k, row in dsets.iterrows():
    t0 = time.time()
    with exp.axis_query("RNA", obs_query=soma.AxisQuery(value_filter=f"dataset_id == '{row.dataset_id}'")) as q:
        obs = q.obs(column_names=["soma_joinid","observation_joinid","raw_sum"]).concat().to_pandas()
        obs['cl'] = obs.observation_joinid.map(cl).fillna(-1).astype(int)
        np.add.at(counts[:, k], obs.cl[obs.cl >= 0].values, 1)
        cell_group = pd.Series(obs.cl.values, index=obs.soma_joinid.values)
        cell_size = pd.Series(obs.raw_sum.values, index=obs.soma_joinid.values)
        for tbl in q.X("raw").tables():
            ci = tbl["soma_dim_0"].to_numpy(); gj = tbl["soma_dim_1"].to_numpy(); v = tbl["soma_data"].to_numpy()
            grp = cell_group.reindex(ci).to_numpy(); ok = grp >= 0
            val = np.log1p(v[ok] / cell_size.reindex(ci[ok]).to_numpy() * 1e4)
            idx = grp[ok].astype(np.int64) * G + gj[ok]
            S += np.bincount(idx, weights=val, minlength=NC * G)
            F += np.bincount(idx, minlength=NC * G)
    print(k + 1, len(dsets), f'{time.time()-t0:.0f}s', flush=True)
n = counts.sum(1)
np.savez_compressed('agg2.npz', S=S.reshape(NC, G).astype(np.float32), F=F.reshape(NC, G).astype(np.float32), counts=counts, n=n)
var.to_parquet('var.parquet')
print('done', flush=True)
