import sys, time, numpy as np, pandas as pd
sys.path.insert(0,'.')
from cxopen import open_c
import tiledbsoma as soma
c = open_c()
exp = c["census_data"]["homo_sapiens"]
var = exp.ms["RNA"].var.read(column_names=["soma_joinid","feature_id","feature_name"]).concat().to_pandas()
G = int(var.soma_joinid.max()) + 1
sc = pd.read_parquet('sc_map.parquet').set_index('observation_joinid').sc
dsets = pd.read_csv('hipamy_datasets.csv')
out = []
for k, row in dsets.iterrows():
    t0 = time.time()
    with exp.axis_query("RNA", obs_query=soma.AxisQuery(value_filter=f"dataset_id == '{row.dataset_id}'")) as q:
        obs = q.obs(column_names=["soma_joinid","observation_joinid","raw_sum"]).concat().to_pandas()
        obs['sc'] = obs.observation_joinid.map(sc)
        groups = sorted(obs.sc.dropna().unique())
        gi = {g:i for i,g in enumerate(groups)}
        cell_group = pd.Series(obs.sc.map(gi).fillna(-1).astype(int).values, index=obs.soma_joinid.values)
        cell_size = pd.Series(obs.raw_sum.values, index=obs.soma_joinid.values)
        S = np.zeros(len(groups)*G); F = np.zeros(len(groups)*G)
        nnz = 0
        for tbl in q.X("raw").tables():
            ci = tbl["soma_dim_0"].to_numpy(); gj = tbl["soma_dim_1"].to_numpy(); v = tbl["soma_data"].to_numpy()
            grp = cell_group.reindex(ci).to_numpy(); ok = grp >= 0
            val = np.log1p(v[ok] / cell_size.reindex(ci[ok]).to_numpy() * 1e4)
            idx = grp[ok].astype(np.int64) * G + gj[ok]
            S += np.bincount(idx, weights=val, minlength=len(groups)*G)
            F += np.bincount(idx, minlength=len(groups)*G)
            nnz += len(v)
        ncell = obs.sc.value_counts()
        for g, i in gi.items():
            out.append(dict(dissection=row.dataset_title.replace('Dissection: ',''), supercluster=g, n=int(ncell[g]),
                            mean=(S[i*G:(i+1)*G] / ncell[g]).astype(np.float32), frac=(F[i*G:(i+1)*G] / ncell[g]).astype(np.float32)))
    print(k+1, len(dsets), row.dataset_title[:60], len(obs), 'cells', nnz, 'nnz', f'{time.time()-t0:.0f}s', flush=True)
    pd.to_pickle((var, out), 'agg.pkl')
print('done', flush=True)
