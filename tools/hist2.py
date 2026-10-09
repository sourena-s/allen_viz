# Count histograms for every cells-window table in one pass over all 105 Siletti dissections: per gene and
# row unit, the share of the unit's nuclei with >= 2/4/8/16/32/64/128 raw counts (<file>_hist.bin.gz,
# uint8 [bin][gene][unit]); also mean raw counts per nucleus (<file>_raw.bin.gz) for the Siletti cluster
# and Protein Atlas tables, which had none.
import sys, os, time, gzip, json, struct, numpy as np, pandas as pd
sys.path.insert(0, '.')
from cxopen import open_c_small as open_c
import tiledbsoma as soma
from scdefs import ORDER
from wb_defs import REGIONS, SCS
A = '/home/user/allen_viz/assets/'
TH = [2, 4, 8, 16, 32, 64, 128]
def head(name):
    b = gzip.open(A + name + '.bin.gz').read(); n = struct.unpack('<I', b[8:12])[0]; return json.loads(b[12:12 + n])
var = pd.read_parquet('var.parquet'); G = int(var.soma_joinid.max()) + 1
jid = var.drop_duplicates('feature_name').set_index('feature_name').soma_joinid
dis = [o[0] for o in ORDER]
CODES = ["M1C", "A46", "A44-A45", "A14", "A13", "A32", "A24", "A25", "A23", "Ig", "Idg", "FI", "S1C", "A5-A7", "A40", "A43",
         "A1C", "STG", "MTG", "ITG", "A38", "Temporal area TF", "V1C", "V2", "A19"]
reg_of = {p: i for i, (_, ps) in enumerate(REGIONS) for p in ps}
F = {}
for name in ['siletti_hipamy', 'siletti_custom', 'siletti_clusters', 'siletti_hpa', 'siletti_brain', 'siletti_ncx', 'siletti_ncxcl']:
    h = head(name); F[name] = h
    h['U'] = len(h['clusters']) if 'clusters' in h and name in ('siletti_clusters', 'siletti_hpa') else len(h['groups'])
    h['H'] = np.zeros((len(TH), h['U'], G), np.float32); h['N'] = np.zeros(h['U']); h['R'] = np.zeros((h['U'], G), np.float32) if name in ('siletti_clusters', 'siletti_hpa') else None
    print(name, h['U'], flush=True)
key_region = {(d, c): u for u, (d, c, _) in enumerate(F['siletti_hipamy']['groups'])}; sc_region = {s: i for i, s in enumerate(F['siletti_hipamy']['superclusters'])}
key_custom = {(d, c): u for u, (d, c, _) in enumerate(F['siletti_custom']['groups'])}
key_brain = {(d, c): u for u, (d, c, _) in enumerate(F['siletti_brain']['groups'])}; sc_brain = {s: i for i, s in enumerate(F['siletti_brain']['superclusters'])}
key_ncx = {(d, c): u for u, (d, c, _) in enumerate(F['siletti_ncx']['groups'])}; sc_ncx = {s: i for i, s in enumerate(F['siletti_ncx']['superclusters'])}
cl_clusters = {c['id']: u for u, c in enumerate(F['siletti_clusters']['clusters'])}
cl_ncxcl = {c: u for u, c in enumerate(F['siletti_ncxcl']['clId'])}
hpa_unit = {c['id']: u for u, c in enumerate(F['siletti_hpa']['clusters'])}
sc = pd.read_parquet('sc_map.parquet').set_index('observation_joinid').sc
cl = pd.read_parquet('cl_map.parquet').set_index('observation_joinid').cluster
lab = pd.read_parquet('rc_labels_v2.parquet'); ckeys = F['siletti_custom']['keys']
code = pd.Series(lab.subfield.map({l: i for i, l in enumerate(ckeys)}).values, index=lab.soma_joinid.values)
hip_ds = pd.read_csv('hipamy_datasets.csv'); hip_k = {r.dataset_id: k for k, r in hip_ds.iterrows()}
Kh = len(pd.read_csv('hpa/meta.tsv', sep='\t'))
ds = pd.read_csv('wb_datasets.csv'); c = open_c(); exp = c["census_data"]["homo_sapiens"]
done = set(open('hist2_done.txt').read().split()) if os.path.exists('hist2_done.txt') else set()
if os.path.exists('hist2_sums.npz'):
    z = np.load('hist2_sums.npz', allow_pickle=True)
    for name, h in F.items():
        h['H'] = z[name + '_H']; h['N'] = z[name + '_N']
        if h['R'] is not None: h['R'] = z[name + '_R']
def save():
    np.savez('hist2_tmp.npz', **{f'{n}_{k}': h[k] for n, h in F.items() for k in ('H', 'N', 'R') if h[k] is not None})
    os.replace('hist2_tmp.npz', 'hist2_sums.npz'); open('hist2_done.txt', 'w').write('\n'.join(done))
for k, row in ds.iterrows():
    if row.dataset_id in done: continue
    t = row.dataset_title.replace('Dissection: ', ''); t0 = time.time()
    with exp.axis_query("RNA", obs_query=soma.AxisQuery(value_filter=f"dataset_id == '{row.dataset_id}'")) as q:
        obs = q.obs(column_names=["soma_joinid", "observation_joinid"]).concat().to_pandas()
        s = obs.observation_joinid.map(sc).values; cid = obs.observation_joinid.map(cl).values
        maps = {}
        n = len(obs); none = np.full(n, -1)
        if t in dis:   # hippocampus / PHG / amygdala dissections
            d = dis.index(t)
            maps['siletti_hipamy'] = np.array([key_region.get((d, sc_region.get(x, -9)), -1) for x in s])
            if d < 11:
                cc = code.reindex(obs.soma_joinid.values).values
                maps['siletti_custom'] = np.array([key_custom.get((d, int(x)), -1) if x == x else -1 for x in cc])
            maps['siletti_clusters'] = np.array([cl_clusters.get(int(x), -1) if x == x else -1 for x in cid])
            z6 = np.load(f'agg6/{hip_k[row.dataset_id]}.npz'); best = z6['best']; rb = z6['rb']
            if len(best) == n: maps['siletti_hpa'] = np.array([hpa_unit.get(int(b), -1) if r >= 0.15 else -1 for b, r in zip(best, rb)])
            else: print('hpa length mismatch', t, flush=True)
        r = reg_of.get(t.split(' - ')[0])
        if r is not None: maps['siletti_brain'] = np.array([key_brain.get((r, sc_brain.get(x, -9)), -1) for x in s])
        if t.startswith('Cerebral cortex (Cx)') and t.split(' - ')[-1] in CODES:
            a = CODES.index(t.split(' - ')[-1])
            maps['siletti_ncx'] = np.array([key_ncx.get((a, sc_ncx.get(x, -9)), -1) for x in s])
            maps['siletti_ncxcl'] = np.array([cl_ncxcl.get(int(x), -1) if x == x else -1 for x in cid])
        base = int(obs.soma_joinid.min()); span = int(obs.soma_joinid.max()) - base + 1
        UA = {}
        for name, u in maps.items():
            arr = np.full(span, -1, np.int64); arr[obs.soma_joinid.values - base] = u; UA[name] = arr
            F[name]['N'] += np.bincount(u[u >= 0], minlength=F[name]['U'])
        for tbl in q.X("raw").tables():
            ci = tbl["soma_dim_0"].to_numpy() - base; gj = tbl["soma_dim_1"].to_numpy(); v = tbl["soma_data"].to_numpy()
            big = v >= 2
            for name, arr in UA.items():
                h = F[name]; U = h['U']; g = arr[ci]
                if h['R'] is not None:
                    m = g >= 0; h['R'] += np.bincount(g[m] * G + gj[m], weights=v[m], minlength=U * G).reshape(U, G).astype(np.float32)
                m = big & (g >= 0); gg = g[m] * G + gj[m]; vv = v[m]
                for i, th in enumerate(TH):
                    sel = vv >= th
                    if not sel.any(): break
                    h['H'][i] += np.bincount(gg[sel], minlength=U * G).reshape(U, G).astype(np.float32)
    done.add(row.dataset_id)
    if len(done) % 5 == 0: save()
    print(len(done), '/', len(ds), t[:70], n, list(maps), f'{time.time()-t0:.0f}s', flush=True)
save()
# write the tables in each file's gene and unit order
for name, h in F.items():
    gcol = jid.reindex(h['genes']).to_numpy(); ok = ~np.isnan(gcol.astype(float)); gc = gcol[ok].astype(np.int64)
    N = np.maximum(h['N'], 1)
    Q = np.zeros((len(TH), len(h['genes']), h['U']), np.uint8)
    for i in range(len(TH)): Q[i][ok] = np.round(h['H'][i][:, gc].T / N * 255).clip(0, 255).astype(np.uint8)
    open(f'{name}_hist.bin.gz', 'wb').write(gzip.compress(Q.tobytes(), 6))
    if h['R'] is not None:
        Rm = np.zeros((len(h['genes']), h['U']), np.float32); Rm[ok] = (h['R'][:, gc] / N[:, None]).T
        mx = Rm.max(1); mx[mx == 0] = 1
        open(f'{name}_raw.bin.gz', 'wb').write(gzip.compress(mx.astype('<f4').tobytes() + np.round(Rm / mx[:, None] * 255).astype(np.uint8).tobytes(), 6))
    print(name, 'written', flush=True)
print('done', flush=True)
