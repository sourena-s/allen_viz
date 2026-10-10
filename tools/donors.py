# Per-donor tables for the cells window's donor bar plots (tools/donors.py): one pass over all 105 Siletti
# dissections for a group of files, per unit, donor and gene: nuclei, share expressing, shares with >= 2..128
# counts, mean raw count and mean CP10k. Run as  donors.py A  (hipamy, custom, clusters, hpa)  then  donors.py B.
# (Header of bnorm.py, from which it derives:)
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
for name in ['siletti_hipamy', 'siletti_custom', 'siletti_clusters', 'siletti_hpa', 'siletti_brain', 'siletti_ncx', 'siletti_ncxcl']:   # all headers (for the unit maps)
    h = head(name); F[name] = h
    h['U'] = len(h['clusters']) if 'clusters' in h and name in ('siletti_clusters', 'siletti_hpa') else len(h['groups'])
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
GRP = {'A': ['siletti_hipamy', 'siletti_custom', 'siletti_clusters', 'siletti_hpa'], 'B': ['siletti_brain', 'siletti_ncx', 'siletti_ncxcl']}[sys.argv[1]]
BINS = [2, 4, 8, 16, 32, 64, 128]; NS = 3 + len(BINS)   # stats: expressing, raw sum, CP10k sum, then >= each bin
for name in GRP:
    h = F[name]; gcol = jid.reindex(h['genes']).to_numpy(); ok = ~np.isnan(gcol.astype(float))
    h['gm'] = np.full(G, -1, np.int64); h['gm'][gcol[ok].astype(np.int64)] = np.flatnonzero(ok); h['Gf'] = len(h['genes'])
    h['A'] = np.zeros((NS, h['U'] * ND_ * h['Gf']), np.float32); h['N'] = np.zeros(h['U'] * ND_)
CK = f'dn_{sys.argv[1]}'; DON = {}
done = set(open(CK + '_done.txt').read().split()) if os.path.exists(CK + '_done.txt') else set()
if done:
    for name in GRP: F[name]['A'] = np.load(f'{CK}_{name}_A.npy'); F[name]['N'] = np.load(f'{CK}_{name}_N.npy')
    DON = json.load(open(CK + '_don.json'))
def save():
    for name in GRP: np.save(f'{CK}_{name}_A.npy', F[name]['A']); np.save(f'{CK}_{name}_N.npy', F[name]['N'])
    json.dump(DON, open(CK + '_don.json', 'w')); open(CK + '_done.txt', 'w').write('\n'.join(done))
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
        for name in list(maps):
            if name not in GRP: del maps[name]
        UA = {}
        for name, u in maps.items():
            arr = np.full(span, -1, np.int64); arr[obs.soma_joinid.values - base] = np.where(u >= 0, u * ND_ + dn, -1); UA[name] = arr
            m = u >= 0; F[name]['N'] += np.bincount((u * ND_ + dn)[m], minlength=F[name]['U'] * ND_)
        for tbl in q.X("raw").tables():
            ci = tbl["soma_dim_0"].to_numpy() - base; gj = tbl["soma_dim_1"].to_numpy(); v = tbl["soma_data"].to_numpy()
            cp = v / sizeA[ci] * 1e4
            for name, arr in UA.items():
                h = F[name]; g = arr[ci]; gf = h['gm'][gj]; m = (g >= 0) & (gf >= 0); L = h['A'].shape[1]
                key = g[m] * h['Gf'] + gf[m]; vm = v[m]; uk, inv = np.unique(key, return_inverse=True); K = len(uk)
                h['A'][0, uk] += np.bincount(inv, minlength=K); h['A'][1, uk] += np.bincount(inv, weights=vm, minlength=K)
                h['A'][2, uk] += np.bincount(inv, weights=cp[m], minlength=K)
                for b, lo in enumerate(BINS):
                    mb = vm >= lo
                    if mb.any(): h['A'][3 + b, uk] += np.bincount(inv[mb], minlength=K)
    done.add(row.dataset_id)
    if len(done) % 15 == 0: save()
    print(len(done), '/', len(ds), t[:60], n, f'{time.time()-t0:.0f}s', flush=True)
save()
print('donors', DON, flush=True)
# <file>_donor.bin.gz: "DON1", uint32 D, uint32 NB, float32 nuclei[U][D], float32 rawMax[D][G], float32 cpMax[D][G],
# then uint8 [D][G][U] tables (each over its maximum, or a share of the unit-donor's nuclei, /255):
# share expressing, mean raw count, mean CP10k, then NB tables of the share with >= 2, 4, ... 128 counts
for name in GRP:
    h = F[name]; U, Gf = h['U'], h['Gf']; N = h['N'].reshape(U, ND_); A = h['A'].reshape(NS, U, ND_, Gf).transpose(0, 2, 3, 1)   # stat, D, G, U
    nn = np.maximum(N.T, 1)[:, None, :]   # D, 1, U
    q = lambda x: np.round(np.clip(x, 0, 1) * 255).astype(np.uint8)
    raw = A[1] / nn; cpm = A[2] / nn; rmx = raw.max(2); rmx[rmx == 0] = 1; cmx = cpm.max(2); cmx[cmx == 0] = 1
    out = [b'DON1', struct.pack('<II', ND_, len(BINS)), N.astype('<f4').tobytes(), rmx.astype('<f4').tobytes(), cmx.astype('<f4').tobytes(),
           q(A[0] / nn).tobytes(), q(raw / rmx[:, :, None]).tobytes(), q(cpm / cmx[:, :, None]).tobytes()] + [q(A[3 + b] / nn).tobytes() for b in range(len(BINS))]
    open(f'{name}_donor.bin.gz', 'wb').write(gzip.compress(b''.join(out), 6))
    print(name, 'written', flush=True)
print('done', flush=True)
