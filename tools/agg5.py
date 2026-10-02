# Linear means: per dataset (dissection), sum of counts per 10k per gene for each supercluster
# (brief file) and each cluster (detailed file); same nuclei, same cell sizes as agg.py / agg2.py
import sys, os, time, numpy as np, pandas as pd
sys.path.insert(0,'.')
from cxopen import open_c
import tiledbsoma as soma
c = open_c()
exp = c["census_data"]["homo_sapiens"]
var = pd.read_parquet('var.parquet'); G = int(var.soma_joinid.max()) + 1
sc = pd.read_parquet('sc_map.parquet').set_index('observation_joinid').sc
cl = pd.read_parquet('cl_map.parquet').set_index('observation_joinid').cluster
NC = int(cl.max()) + 1
dsets = pd.read_csv('hipamy_datasets.csv')
os.makedirs('agg5', exist_ok=True)
LC = np.zeros(NC * G)
for k, row in dsets.iterrows():
    t0 = time.time()
    if os.path.exists(f'agg5/{k}.npz'):
        LC += np.load(f'agg5/{k}.npz')['LCk'].ravel(); continue
    with exp.axis_query("RNA", obs_query=soma.AxisQuery(value_filter=f"dataset_id == '{row.dataset_id}'")) as q:
        obs = q.obs(column_names=["soma_joinid","observation_joinid","raw_sum"]).concat().to_pandas()
        obs['sc'] = obs.observation_joinid.map(sc)
        groups = sorted(obs.sc.dropna().unique()); gi = {g:i for i,g in enumerate(groups)}
        sgrp = pd.Series(obs.sc.map(gi).fillna(-1).astype(int).values, index=obs.soma_joinid.values)
        cgrp = pd.Series(obs.observation_joinid.map(cl).fillna(-1).astype(int).values, index=obs.soma_joinid.values)
        size = pd.Series(obs.raw_sum.values, index=obs.soma_joinid.values)
        LS = np.zeros(len(groups) * G); LCk = np.zeros(NC * G)
        for tbl in q.X("raw").tables():
            ci = tbl["soma_dim_0"].to_numpy(); gj = tbl["soma_dim_1"].to_numpy(); v = tbl["soma_data"].to_numpy()
            val = v / size.reindex(ci).to_numpy() * 1e4
            s = sgrp.reindex(ci).to_numpy(); ok = s >= 0
            LS += np.bincount(s[ok].astype(np.int64) * G + gj[ok], weights=val[ok], minlength=len(groups) * G)
            cc = cgrp.reindex(ci).to_numpy(); ok = cc >= 0
            LCk += np.bincount(cc[ok].astype(np.int64) * G + gj[ok], weights=val[ok], minlength=NC * G)
    np.savez(f'agg5/{k}.npz', groups=np.array(groups), LS=LS.reshape(len(groups), G).astype(np.float32), LCk=LCk.reshape(NC, G).astype(np.float32))
    LC += LCk
    print(k + 1, len(dsets), f'{time.time()-t0:.0f}s', flush=True)
np.save('agg5_LC.npy', LC.reshape(NC, G).astype(np.float32))
print('done', flush=True)
