# Neocortex at Siletti cluster level: the 25 neocortical dissections of wbcx.py pooled, per Siletti
# cluster (461 brain-wide), whole-genome sums of log1p(cp10k), cp10k, detections, raw counts and
# nuclei with >= 2/4/8/16/32 counts. Checkpoints in wbcl_sums.npz / wbcl_done.txt.
import sys, os, time, numpy as np, pandas as pd
sys.path.insert(0, '.')
from cxopen import open_c_small as open_c
import tiledbsoma as soma
CODES = ["M1C", "A46", "A44-A45", "A14", "A13", "A32", "A24", "A25", "A23", "Ig", "Idg", "FI", "S1C", "A5-A7", "A40", "A43",
         "A1C", "STG", "MTG", "ITG", "A38", "Temporal area TF", "V1C", "V2", "A19"]
def code_of(t): return t.split(' - ')[-1]
cl = pd.read_parquet('cl_map.parquet').set_index('observation_joinid').cluster
U = int(cl.max()) + 1; TH = [2, 4, 8, 16, 32]
var = pd.read_parquet('var.parquet'); G = int(var.soma_joinid.max()) + 1
ds = pd.read_csv('wb_datasets.csv')
done = set(open('wbcl_done.txt').read().split()) if os.path.exists('wbcl_done.txt') else set()
if os.path.exists('wbcl_sums.npz'):
    z = np.load('wbcl_sums.npz'); A = {k: z[k] for k in ['S', 'L', 'F', 'R', 'H', 'N']}
else:
    A = dict(S=np.zeros((U, G), np.float32), L=np.zeros((U, G), np.float32), F=np.zeros((U, G), np.float32), R=np.zeros((U, G), np.float32),
             H=np.zeros((len(TH), U, G), np.float32), N=np.zeros(U))
def save():
    np.savez('wbcl_tmp.npz', **A); os.replace('wbcl_tmp.npz', 'wbcl_sums.npz'); open('wbcl_done.txt', 'w').write('\n'.join(done))
c = open_c(); exp = c["census_data"]["homo_sapiens"]
for k, row in ds.iterrows():
    t = row.dataset_title.replace('Dissection: ', '')
    if not t.startswith('Cerebral cortex (Cx)') or code_of(t) not in CODES or row.dataset_id in done: continue
    t0 = time.time()
    with exp.axis_query("RNA", obs_query=soma.AxisQuery(value_filter=f"dataset_id == '{row.dataset_id}'")) as q:
        obs = q.obs(column_names=["soma_joinid", "observation_joinid", "raw_sum"]).concat().to_pandas()
        u = obs.observation_joinid.map(cl).fillna(-1).astype(np.int64).to_numpy()
        base = int(obs.soma_joinid.min()); span = int(obs.soma_joinid.max()) - base + 1
        grpA = np.full(span, -1, np.int64); grpA[obs.soma_joinid.values - base] = u; sizeA = np.ones(span); sizeA[obs.soma_joinid.values - base] = obs.raw_sum.values
        A['N'] += np.bincount(u[u >= 0], minlength=U)
        buf = []; nb = 0
        def flush():
            global buf, nb
            if not buf: return
            idx = np.concatenate([x[0] for x in buf]); cp = np.concatenate([x[1] for x in buf]); v = np.concatenate([x[2] for x in buf]); buf = []; nb = 0
            for key, w in (('S', np.log1p(cp)), ('L', cp), ('F', None), ('R', v)):
                A[key] += np.bincount(idx, weights=w, minlength=U * G).reshape(U, G).astype(np.float32)
            for i, th in enumerate(TH):
                sel = v >= th; A['H'][i] += np.bincount(idx[sel], minlength=U * G).reshape(U, G).astype(np.float32)
        for tbl in q.X("raw").tables():
            ci = tbl["soma_dim_0"].to_numpy(); gj = tbl["soma_dim_1"].to_numpy(); v = tbl["soma_data"].to_numpy()
            g = grpA[ci - base]; m = g >= 0
            g, gj, v, ci = g[m], gj[m], v[m], ci[m]
            buf.append((g * G + gj, (v / sizeA[ci - base] * 1e4).astype(np.float32), v.astype(np.float32))); nb += len(v)
            if nb > 60_000_000: flush()
        flush()
    done.add(row.dataset_id)
    if len(done) % 5 == 0: save()
    print(len(done), t[:70], len(obs), f'{time.time()-t0:.0f}s', flush=True)
save(); print('done', flush=True)
