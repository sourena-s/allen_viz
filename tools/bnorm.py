# Brain norm per donor: in one pass over all 105 Siletti dissections, each donor's whole-brain CP10k sum per gene
# and each cells-window unit's CP10k sums split by donor. Then per unit and gene: the mean over its nuclei of
# CP10k / (that nucleus's donor's whole-brain mean CP10k) = sum_d S_ud / B_d / n_u  (<file>_bnorm.bin.gz).
import sys, os, time, gzip, json, struct, numpy as np, pandas as pd
sys.path.insert(0, '.')
from cxopen import open_c_small as open_c
import tiledbsoma as soma
from scdefs import ORDER
from wb_defs import REGIONS, SCS
A = '/home/user/allen_viz/assets/'
ND_ = 4   # Siletti donors (H18.30.001 gave only a few dissections)
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
    h['S'] = np.zeros((h['U'] * ND_, G), np.float32); h['N'] = np.zeros(h['U'] * ND_)   # [unit * 3 + donor]
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
DON = {}; T = np.zeros((ND_, G)); ND = np.zeros(ND_)
done = set(open('bn_done.txt').read().split()) if os.path.exists('bn_done.txt') else set()
if os.path.exists('bn_sums.npz'):
    z = np.load('bn_sums.npz', allow_pickle=True); T = z['T']; ND = z['ND']; DON = dict(z['DON'].item())
    for name, h in F.items(): h['S'] = z[name + '_S']; h['N'] = z[name + '_N']
def save():
    np.savez('bn_tmp.npz', T=T, ND=ND, DON=np.array(DON, dtype=object), **{f'{n}_{k}': h[k] for n, h in F.items() for k in ('S', 'N')})
    os.replace('bn_tmp.npz', 'bn_sums.npz'); open('bn_done.txt', 'w').write('\n'.join(done))
if '--write' in sys.argv:
    z = np.load('bn_sums.npz', allow_pickle=True); T, ND = z['T'], z['ND']; DON = z['DON'].tolist()
    for n, h in F.items(): h['S'], h['N'] = z[f'{n}_S'], z[f'{n}_N']
    ds = ds.iloc[:0]
for k, row in ds.iterrows():
    if row.dataset_id in done: continue
    t = row.dataset_title.replace('Dissection: ', ''); t0 = time.time()
    with exp.axis_query("RNA", obs_query=soma.AxisQuery(value_filter=f"dataset_id == '{row.dataset_id}'")) as q:
        obs = q.obs(column_names=["soma_joinid", "observation_joinid", "raw_sum", "donor_id"]).concat().to_pandas()
        for dname in obs.donor_id.unique():
            if dname not in DON: DON[dname] = len(DON)
        dn = obs.donor_id.map(DON).values.astype(np.int64)
        s = obs.observation_joinid.map(sc).values; cid = obs.observation_joinid.map(cl).values
        maps = {}; n = len(obs)
        if t in dis:
            d = dis.index(t)
            maps['siletti_hipamy'] = np.array([key_region.get((d, sc_region.get(x, -9)), -1) for x in s])
            if d < 11:
                cc = code.reindex(obs.soma_joinid.values).values
                maps['siletti_custom'] = np.array([key_custom.get((d, int(x)), -1) if x == x else -1 for x in cc])
            maps['siletti_clusters'] = np.array([cl_clusters.get(int(x), -1) if x == x else -1 for x in cid])
            z6 = np.load(f'agg6/{hip_k[row.dataset_id]}.npz'); best = z6['best']; rb = z6['rb']
            if len(best) == n: maps['siletti_hpa'] = np.array([hpa_unit.get(int(b), -1) if r >= 0.15 else -1 for b, r in zip(best, rb)])
        r = reg_of.get(t.split(' - ')[0])
        if r is not None: maps['siletti_brain'] = np.array([key_brain.get((r, sc_brain.get(x, -9)), -1) for x in s])
        if t.startswith('Cerebral cortex (Cx)') and t.split(' - ')[-1] in CODES:
            a = CODES.index(t.split(' - ')[-1])
            maps['siletti_ncx'] = np.array([key_ncx.get((a, sc_ncx.get(x, -9)), -1) for x in s])
            maps['siletti_ncxcl'] = np.array([cl_ncxcl.get(int(x), -1) if x == x else -1 for x in cid])
        base = int(obs.soma_joinid.min()); span = int(obs.soma_joinid.max()) - base + 1
        dA = np.zeros(span, np.int64); dA[obs.soma_joinid.values - base] = dn
        sizeA = np.ones(span); sizeA[obs.soma_joinid.values - base] = obs.raw_sum.values
        np.add.at(ND, dn, 1)
        UA = {}
        for name, u in maps.items():
            arr = np.full(span, -1, np.int64); arr[obs.soma_joinid.values - base] = np.where(u >= 0, u * ND_ + dn, -1); UA[name] = arr
            m = u >= 0; F[name]['N'] += np.bincount((u * ND_ + dn)[m], minlength=F[name]['U'] * ND_)
        for tbl in q.X("raw").tables():
            ci = tbl["soma_dim_0"].to_numpy() - base; gj = tbl["soma_dim_1"].to_numpy(); v = tbl["soma_data"].to_numpy()
            cp = v / sizeA[ci] * 1e4; dd = dA[ci]
            T += np.bincount(dd * G + gj, weights=cp, minlength=ND_ * G).reshape(ND_, G)
            for name, arr in UA.items():
                h = F[name]; g = arr[ci]; m = g >= 0
                h['S'] += np.bincount(g[m] * G + gj[m], weights=cp[m], minlength=h['U'] * ND_ * G).reshape(h['U'] * ND_, G).astype(np.float32)
    done.add(row.dataset_id)
    if len(done) % 5 == 0: save()
    print(len(done), '/', len(ds), t[:60], n, dict(obs.donor_id.value_counts()), f'{time.time()-t0:.0f}s', flush=True)
save()
B = T / np.maximum(ND, 1)[:, None]   # per donor whole-brain mean CP10k per gene
print('donors', DON, 'nuclei', ND.tolist(), flush=True)
# <file>_bnorm.bin.gz, kept per donor so the page averages inline:
#   "BND1", uint32 D, float32 nuclei[U][D], then uint8 q[D][genes][U] with fold f_ud = S_ud / B_d / n_ud coded
#   on a log2 scale: q = 0 means 0, else log2 f = LO + (q - 1) * (HI - LO) / 254, LO = -12, HI = 8
LO, HI = -12.0, 8.0
for name, h in F.items():
    gcol = jid.reindex(h['genes']).to_numpy(); ok = ~np.isnan(gcol.astype(float)); gc = gcol[ok].astype(np.int64)
    U = h['U']; S = h['S'].reshape(U, ND_, G)[:, :, gc]; Nud = h['N'].reshape(U, ND_)
    Bg = B[:, gc]; f = np.where(Bg[None] > 0, S / np.where(Bg[None] > 0, Bg[None], 1), 0) / np.maximum(Nud, 1)[:, :, None]   # U x D x genes
    Q = np.zeros((ND_, len(h['genes']), U), np.uint8)
    with np.errstate(divide='ignore'): l2 = np.log2(f)
    q = np.where(f > 0, np.clip(np.round((np.clip(l2, LO, HI) - LO) / (HI - LO) * 254) + 1, 1, 255), 0).astype(np.uint8)
    Q[:, ok, :] = q.transpose(1, 2, 0)
    open(f'{name}_bnorm.bin.gz', 'wb').write(gzip.compress(b'BND1' + struct.pack('<I', ND_) + Nud.astype('<f4').tobytes() + Q.tobytes(), 6))
    print(name, 'written', flush=True)
print('done', flush=True)
