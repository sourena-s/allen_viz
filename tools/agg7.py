# Whole-genome sums per (reclustered subfield, dissection) for the hippocampal principal neurons
# labelled by rc2.py: log1p(counts per 10k), counts per 10k, detections; per hippocampal dissection
import sys, os, time, numpy as np, pandas as pd
sys.path.insert(0,'.')
from cxopen import open_c
import tiledbsoma as soma
from scdefs import ORDER
OUT = os.environ.get('OUT', 'agg7')
LABS = ["CA1", "CA2", "CA3", "SUB", "DG", "none"]
lab = pd.read_parquet(os.environ.get('LABELS', 'rc_labels_1.0.parquet'))
code = pd.Series(lab.subfield.map({l: i for i, l in enumerate(LABS)}).values, index=lab.soma_joinid.values)
var = pd.read_parquet('var.parquet'); G = int(var.soma_joinid.max()) + 1; L = len(LABS)
dsets = pd.read_csv('hipamy_datasets.csv'); hip = {o[0] for o in ORDER[:11]}
c = open_c(); exp = c["census_data"]["homo_sapiens"]
os.makedirs(OUT, exist_ok=True)
for k, row in dsets.iterrows():
    title = row.dataset_title.replace('Dissection: ', '')
    if title not in hip or os.path.exists(f'{OUT}/{k}.npz'): continue
    t0 = time.time()
    with exp.axis_query("RNA", obs_query=soma.AxisQuery(value_filter=f"dataset_id == '{row.dataset_id}'")) as q:
        obs = q.obs(column_names=["soma_joinid", "raw_sum"]).concat().to_pandas()
        grp = pd.Series(code.reindex(obs.soma_joinid.values).fillna(-1).astype(int).values, index=obs.soma_joinid.values)
        size = pd.Series(obs.raw_sum.values, index=obs.soma_joinid.values)
        S = np.zeros(L * G); Lc = np.zeros(L * G); F = np.zeros(L * G)
        for tbl in q.X("raw").tables():
            ci = tbl["soma_dim_0"].to_numpy(); gj = tbl["soma_dim_1"].to_numpy(); v = tbl["soma_data"].to_numpy()
            g = grp.reindex(ci).to_numpy(); ok = g >= 0
            cp = v[ok] / size.reindex(ci[ok]).to_numpy() * 1e4; idx = g[ok].astype(np.int64) * G + gj[ok]
            S += np.bincount(idx, weights=np.log1p(cp), minlength=L * G); Lc += np.bincount(idx, weights=cp, minlength=L * G); F += np.bincount(idx, minlength=L * G)
        n = np.bincount(grp.values[grp.values >= 0], minlength=L)
    np.savez(f'{OUT}/{k}.npz', S=S.reshape(L, G).astype(np.float32), L=Lc.reshape(L, G).astype(np.float32), F=F.reshape(L, G).astype(np.float32), n=n, title=title)
    print(title, n.tolist(), f'{time.time()-t0:.0f}s', flush=True)
print('done', flush=True)
