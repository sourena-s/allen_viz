# Neocortex of Siletti et al. 2023: per neocortical dissection (area) and supercluster, whole-genome sums
# of log1p(cp10k), cp10k, detections, raw counts and nuclei with >= 2/4/8/16/32 counts (wbcx/<code>.npz).
# Neocortex: Cerebral cortex (Cx) dissections minus entorhinal, perirhinal, parahippocampal (TH-TL),
# retrosplenial (A29-A30) and prostriata, which are not isocortex.
import sys, os, time, numpy as np, pandas as pd
sys.path.insert(0, '.')
from cxopen import open_c_small as open_c
import tiledbsoma as soma
from wb_defs import SCS
CODES = ["M1C", "A46", "A44-A45", "A14", "A13", "A32", "A24", "A25", "A23", "Ig", "Idg", "FI", "S1C", "A5-A7", "A40", "A43",
         "A1C", "STG", "MTG", "ITG", "A38", "Temporal area TF", "V1C", "V2", "A19"]
def code_of(t): return t.split(' - ')[-1]
sc = pd.read_parquet('sc_map.parquet').set_index('observation_joinid').sc; scix = {s: i for i, s in enumerate(SCS)}; NS = len(SCS)
var = pd.read_parquet('var.parquet'); G = int(var.soma_joinid.max()) + 1; TH = [2, 4, 8, 16, 32]
ds = pd.read_csv('wb_datasets.csv'); os.makedirs('wbcx', exist_ok=True)
c = open_c(); exp = c["census_data"]["homo_sapiens"]
for k, row in ds.iterrows():
    t = row.dataset_title.replace('Dissection: ', '')
    if not t.startswith('Cerebral cortex (Cx)') or code_of(t) not in CODES: continue
    out = f'wbcx/{CODES.index(code_of(t)):02d}.npz'
    if os.path.exists(out): continue
    t0 = time.time()
    with exp.axis_query("RNA", obs_query=soma.AxisQuery(value_filter=f"dataset_id == '{row.dataset_id}'")) as q:
        obs = q.obs(column_names=["soma_joinid", "observation_joinid", "raw_sum"]).concat().to_pandas()
        s = obs.observation_joinid.map(sc).map(scix); u = np.where(s.notna(), s.fillna(0).astype(int), -1)
        base = int(obs.soma_joinid.min()); span = int(obs.soma_joinid.max()) - base + 1
        grpA = np.full(span, -1, np.int64); grpA[obs.soma_joinid.values - base] = u; sizeA = np.ones(span); sizeA[obs.soma_joinid.values - base] = obs.raw_sum.values
        N = np.bincount(u[u >= 0], minlength=NS)
        acc = {key: np.zeros(NS * G) for key in ['S', 'L', 'F', 'R']}; hacc = [np.zeros(NS * G) for _ in TH]
        for tbl in q.X("raw").tables():
            ci = tbl["soma_dim_0"].to_numpy(); gj = tbl["soma_dim_1"].to_numpy(); v = tbl["soma_data"].to_numpy()
            g = grpA[ci - base]; m = g >= 0
            g, gj, v, ci = g[m].astype(np.int64), gj[m], v[m], ci[m]
            idx = g * G + gj; cp = v / sizeA[ci - base] * 1e4
            acc['S'] += np.bincount(idx, weights=np.log1p(cp), minlength=NS * G); acc['L'] += np.bincount(idx, weights=cp, minlength=NS * G)
            acc['F'] += np.bincount(idx, minlength=NS * G); acc['R'] += np.bincount(idx, weights=v, minlength=NS * G)
            for i, th in enumerate(TH): sel = v >= th; hacc[i] += np.bincount(idx[sel], minlength=NS * G)
    np.savez(out, N=N, title=t, **{key: a.reshape(NS, G).astype(np.float32) for key, a in acc.items()}, H=np.stack([h.reshape(NS, G) for h in hacc]).astype(np.float32))
    print(t[:80], len(obs), f'{time.time()-t0:.0f}s', flush=True)
print('done', flush=True)
