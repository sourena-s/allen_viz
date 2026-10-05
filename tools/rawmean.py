# Mean raw counts per nucleus, per gene and row unit (groups of the brief files), for the cells
# window's raw ("flow-like") differential: <name>_raw.bin.gz = float32 per-gene maxima, then
# uint8 [gene][unit] = mean / max * 255. MODE=region (siletti_hipamy) or custom (siletti_custom).
import sys, os, time, gzip, json, struct, numpy as np, pandas as pd
sys.path.insert(0,'.')
from cxopen import open_c
import tiledbsoma as soma
from scdefs import ORDER
MODE = os.environ.get('MODE', 'region'); NAME = 'siletti_hipamy' if MODE == 'region' else 'siletti_custom'
b = gzip.open(f'/home/user/allen_viz/assets/{NAME}.bin.gz').read(); n = struct.unpack('<I', b[8:12])[0]; h = json.loads(b[12:12 + n])
dis = [o[0] for o in ORDER]; U = len(h['groups']); gid = {(d, c): u for u, (d, c, _) in enumerate(h['groups'])}
var = pd.read_parquet('var.parquet'); G = int(var.soma_joinid.max()) + 1
gcol = var.drop_duplicates('feature_name').set_index('feature_name').soma_joinid.reindex(h['genes']).to_numpy()
if MODE == 'region':
    sc = pd.read_parquet('sc_map.parquet').set_index('observation_joinid').sc; scix = {s: i for i, s in enumerate(h['superclusters'])}
else:
    lab = pd.read_parquet('rc_labels_v2.parquet'); keys = h['keys']
    code = pd.Series(lab.subfield.map({l: i for i, l in enumerate(keys)}).values, index=lab.soma_joinid.values)
S = np.zeros(U * G, np.float64); N = np.array([g[2] for g in h['groups']], np.float64)
dsets = pd.read_csv('hipamy_datasets.csv'); c = open_c(); exp = c["census_data"]["homo_sapiens"]
for k, row in dsets.iterrows():
    title = row.dataset_title.replace('Dissection: ', '')
    if title not in dis or (MODE == 'custom' and dis.index(title) >= 11): continue
    d = dis.index(title); t0 = time.time()
    with exp.axis_query("RNA", obs_query=soma.AxisQuery(value_filter=f"dataset_id == '{row.dataset_id}'")) as q:
        obs = q.obs(column_names=["soma_joinid", "observation_joinid"]).concat().to_pandas()
        cl = obs.observation_joinid.map(sc).map(scix) if MODE == 'region' else pd.Series(code.reindex(obs.soma_joinid.values).values)
        u = np.array([gid.get((d, int(x)), -1) if x == x else -1 for x in cl.values]); grp = pd.Series(u, index=obs.soma_joinid.values)
        for tbl in q.X("raw").tables():
            ci = tbl["soma_dim_0"].to_numpy(); gj = tbl["soma_dim_1"].to_numpy(); v = tbl["soma_data"].to_numpy()
            g = grp.reindex(ci).to_numpy(); m = g >= 0
            S += np.bincount(g[m].astype(np.int64) * G + gj[m], weights=v[m], minlength=U * G)
    print(title, f'{time.time()-t0:.0f}s', flush=True)
M = S.reshape(U, G) / N[:, None]
ok = ~np.isnan(gcol.astype(float)); R = np.zeros((len(h['genes']), U), np.float32); R[ok] = M[:, gcol[ok].astype(np.int64)].T
mx = R.max(1); mx[mx == 0] = 1
Q = np.round(R / mx[:, None] * 255).astype(np.uint8)
open(f'{NAME}_raw.bin.gz', 'wb').write(gzip.compress(mx.astype('<f4').tobytes() + Q.tobytes(), 6))
print('done', Q.shape, os.path.getsize(f'{NAME}_raw.bin.gz') / 1e6, 'MB', flush=True)
