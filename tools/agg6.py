# Per nucleus: assign to the closest Protein Atlas brain cluster (neurons c-0..c-94, non-neurons c-0..c-31)
# by Pearson r over marker genes (log1p counts per 10k vs the cluster's log1p(nCPM/100)); then per
# (cluster, dissection): sums of log1p(cp10k), cp10k and detections over all genes
import sys, os, time, numpy as np, pandas as pd, scipy.sparse as sp
sys.path.insert(0,'.')
from cxopen import open_c
import tiledbsoma as soma
c = open_c(); exp = c["census_data"]["homo_sapiens"]
var = pd.read_parquet('var.parquet'); G = int(var.soma_joinid.max()) + 1
ref = pd.read_parquet('hpa/ref_markers.parquet'); meta = pd.read_csv('hpa/meta.tsv', sep='\t')
idx = var.set_index('feature_id').soma_joinid
ref = ref[ref.index.isin(idx.index)]
mj = idx.loc[ref.index].values                      # soma_joinids of the markers
col = pd.Series(np.arange(len(mj)), index=mj)
R = ref.values; Rz = (R - R.mean(0)) / R.std(0)        # per cluster, z over the markers
m = len(mj); K = R.shape[1]
dsets = pd.read_csv('hipamy_datasets.csv')
os.makedirs('agg6', exist_ok=True)
for k, row in dsets.iterrows():
    if os.path.exists(f'agg6/{k}.npz'): continue
    t0 = time.time()
    with exp.axis_query("RNA", obs_query=soma.AxisQuery(value_filter=f"dataset_id == '{row.dataset_id}'")) as q:
        obs = q.obs(column_names=["soma_joinid","raw_sum"]).concat().to_pandas()
        pos = pd.Series(np.arange(len(obs)), index=obs.soma_joinid.values); size = obs.raw_sum.values
        rows_, cols_, vals_ = [], [], []
        S = np.zeros((K + 1) * G); Lc = np.zeros((K + 1) * G); F = np.zeros((K + 1) * G)
        tables = []
        for tbl in q.X("raw").tables():
            ci = pos.reindex(tbl["soma_dim_0"].to_numpy()).to_numpy(); gj = tbl["soma_dim_1"].to_numpy(); v = tbl["soma_data"].to_numpy()
            cp = v / size[ci] * 1e4
            ok = np.isin(gj, mj)
            rows_.append(ci[ok]); cols_.append(col.reindex(gj[ok]).to_numpy()); vals_.append(np.log1p(cp[ok]))
            tables.append((ci.astype(np.int32), gj.astype(np.int32), cp.astype(np.float32)))
        X = sp.csr_matrix((np.concatenate(vals_), (np.concatenate(rows_), np.concatenate(cols_))), shape=(len(obs), m))
        mu = np.asarray(X.sum(1)).ravel() / m; sq = np.asarray(X.multiply(X).sum(1)).ravel() / m
        sd = np.sqrt(np.maximum(sq - mu ** 2, 1e-12))
        r = (X @ Rz) / (m * sd[:, None])
        o = np.argsort(-r, 1); best = o[:, 0]; rb = r[np.arange(len(r)), best]; r2 = r[np.arange(len(r)), o[:, 1]]
        assign = np.where((rb >= 0.15), best, K)            # K = unassigned
        for ci, gj, cp in tables:
            g = assign[ci].astype(np.int64) * G + gj
            S += np.bincount(g, weights=np.log1p(cp), minlength=(K + 1) * G)
            Lc += np.bincount(g, weights=cp, minlength=(K + 1) * G)
            F += np.bincount(g, minlength=(K + 1) * G)
        n = np.bincount(assign, minlength=K + 1)
    np.savez(f'agg6/{k}.npz', S=S.reshape(K + 1, G).astype(np.float32), L=Lc.reshape(K + 1, G).astype(np.float32), F=F.reshape(K + 1, G).astype(np.float32), n=n, rb=rb.astype(np.float32), r2=r2.astype(np.float32), best=best)
    print(k + 1, len(dsets), len(obs), 'cells', f'unassigned {n[K]}', f'median r {np.median(rb):.2f}', f'{time.time()-t0:.0f}s', flush=True)
print('done', flush=True)
