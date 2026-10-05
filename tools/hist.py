# Count histograms for the cells window's count bars: per gene and row unit (group of the brief
# file), the fraction of nuclei with >= 2, 4, 8, 16 and 32 raw counts (>= 1 is in the file already).
# MODE=region: groups (dissection, Siletti supercluster) of siletti_hipamy; MODE=custom: groups
# (dissection, reclustered label) of siletti_custom. Output: <name>_hist.bin.gz, uint8 per
# [bin][gene][unit] in the brief file's gene order, value = fraction * 255.
import sys, os, time, gzip, json, struct, numpy as np, pandas as pd
sys.path.insert(0,'.')
from cxopen import open_c
import tiledbsoma as soma
from scdefs import ORDER
MODE = os.environ.get('MODE', 'region'); NAME = 'siletti_hipamy' if MODE == 'region' else 'siletti_custom'
b = gzip.open(f'/home/user/allen_viz/assets/{NAME}.bin.gz').read(); n = struct.unpack('<I', b[8:12])[0]; h = json.loads(b[12:12 + n])
dis = [o[0] for o in ORDER]; U = len(h['groups'])
gid = {(d, c): u for u, (d, c, _) in enumerate(h['groups'])}
var = pd.read_parquet('var.parquet'); G = int(var.soma_joinid.max()) + 1
jid = var.drop_duplicates('feature_name').set_index('feature_name').soma_joinid
gcol = jid.reindex(h['genes']).to_numpy()   # census joinid per file gene
if MODE == 'region':
    sc = pd.read_parquet('sc_map.parquet').set_index('observation_joinid').sc
    scix = {s: i for i, s in enumerate(h['superclusters'])}
else:
    lab = pd.read_parquet('rc_labels_v2.parquet'); keys = h['keys']
    code = pd.Series(lab.subfield.map({l: i for i, l in enumerate(keys)}).values, index=lab.soma_joinid.values)
TH = [2, 4, 8, 16, 32]; H = np.zeros((len(TH), U, G), np.float32); N = np.array([g[2] for g in h["groups"]], np.float64)
dsets = pd.read_csv('hipamy_datasets.csv'); c = open_c(); exp = c["census_data"]["homo_sapiens"]
for k, row in dsets.iterrows():
    title = row.dataset_title.replace('Dissection: ', '')
    if title not in dis or (MODE == 'custom' and dis.index(title) >= 11): continue
    d = dis.index(title); t0 = time.time()
    with exp.axis_query("RNA", obs_query=soma.AxisQuery(value_filter=f"dataset_id == '{row.dataset_id}'")) as q:
        obs = q.obs(column_names=["soma_joinid", "observation_joinid"]).concat().to_pandas()
        if MODE == 'region': cl = obs.observation_joinid.map(sc).map(scix)
        else: cl = pd.Series(code.reindex(obs.soma_joinid.values).values)
        u = np.array([gid.get((d, int(x)), -1) if x == x else -1 for x in cl.values])
        grp = pd.Series(u, index=obs.soma_joinid.values)
        for tbl in q.X("raw").tables():
            ci = tbl["soma_dim_0"].to_numpy(); gj = tbl["soma_dim_1"].to_numpy(); v = tbl["soma_data"].to_numpy()
            ok = v >= 2; g = grp.reindex(ci[ok]).to_numpy(); m = g >= 0
            g, gj2, v2 = g[m].astype(np.int64), gj[ok][m], v[ok][m]
            for i, t in enumerate(TH):
                s = v2 >= t; H[i] += np.bincount(g[s] * G + gj2[s], minlength=U * G).reshape(U, G)
    print(title, f'{time.time()-t0:.0f}s', flush=True)
Q = np.zeros((len(TH), len(h['genes']), U), np.uint8)
okg = ~np.isnan(gcol.astype(float))
for i in range(len(TH)):
    Q[i][okg] = np.round(H[i][:, gcol[okg].astype(np.int64)].T / N * 255).clip(0, 255).astype(np.uint8)
open(f'{NAME}_hist.bin.gz', 'wb').write(gzip.compress(Q.tobytes(), 6))
print('done', Q.shape, os.path.getsize(f'{NAME}_hist.bin.gz') / 1e6, 'MB', flush=True)
