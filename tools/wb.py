# Whole-brain Siletti et al. 2023 (all 105 dissections, 3.37M nuclei): per (brain region, supercluster)
# whole-genome sums of log1p(cp10k), cp10k, detections, raw counts and nuclei with >= 2/4/8/16/32 counts.
# Checkpoints after every dataset (wb_sums.npz + wb_done.txt); pack with packwb.py.
import sys, os, time, numpy as np, pandas as pd
sys.path.insert(0, '.')
from cxopen import open_c_small as open_c
import tiledbsoma as soma
REGIONS = [("Cerebral cortex", ["Cerebral cortex (Cx)", "Perirhinal cortex (area 35) (A35)"]),
           ("Paleocortex and claustrum", ["Paleocortex (PalCx)", "Claustrum"]),
           ("Hippocampus", ["Head of hippocampus (HiH)", "Body of hippocampus (HiB)", "Tail of Hippocampus (HiT)"]),
           ("Amygdala", ["Amygdaloid complex (AMY)", "Extended amygdala (EXA)"]),
           ("Basal nuclei and basal forebrain", ["Basal nuclei (BN)", "Basal forebrain (BF)"]),
           ("Thalamus and epithalamus", ["Thalamus (THM)", "Epithalamus"]),
           ("Hypothalamus", ["Hypothalamus (HTH)"]),
           ("Midbrain", ["Midbrain (M)", "Midbrain (RN)"]),
           ("Pons", ["Pons (Pn)"]),
           ("Medulla", ["Myelencephalon (medulla oblongata) (Mo)"]),
           ("Cerebellum", ["Cerebellum (CB)"]),
           ("Spinal cord", ["Spinal cord"])]
reg_of = {p: i for i, (_, ps) in enumerate(REGIONS) for p in ps}
sc = pd.read_parquet('sc_map.parquet').set_index('observation_joinid').sc
SCS = sorted(sc.unique()); scix = {s: i for i, s in enumerate(SCS)}; NS = len(SCS)
U = len(REGIONS) * NS; TH = [2, 4, 8, 16, 32]
var = pd.read_parquet('var.parquet'); G = int(var.soma_joinid.max()) + 1
ds = pd.read_csv('wb_datasets.csv')
done = set(open('wb_done.txt').read().split()) if os.path.exists('wb_done.txt') else set()
if os.path.exists('wb_sums.npz'):
    z = np.load('wb_sums.npz'); A = {k: z[k] for k in ['S', 'L', 'F', 'R', 'H', 'N']}
else:
    A = dict(S=np.zeros((U, G), np.float32), L=np.zeros((U, G), np.float32), F=np.zeros((U, G), np.float32),
             R=np.zeros((U, G), np.float32), H=np.zeros((len(TH), U, G), np.float32), N=np.zeros(U))
c = open_c(); exp = c["census_data"]["homo_sapiens"]
for k, row in ds.iterrows():
    if row.dataset_id in done: continue
    title = row.dataset_title.replace('Dissection: ', ''); r = reg_of.get(title.split(' - ')[0])
    if r is None: print('skip', title, flush=True); continue
    t0 = time.time()
    with exp.axis_query("RNA", obs_query=soma.AxisQuery(value_filter=f"dataset_id == '{row.dataset_id}'")) as q:
        obs = q.obs(column_names=["soma_joinid", "observation_joinid", "raw_sum"]).concat().to_pandas()
        s = obs.observation_joinid.map(sc).map(scix)
        u = np.where(s.notna(), r * NS + s.fillna(0).astype(int), -1)
        base = int(obs.soma_joinid.min()); span = int(obs.soma_joinid.max()) - base + 1
        grpA = np.full(span, -1, np.int64); grpA[obs.soma_joinid.values - base] = u; sizeA = np.ones(span); sizeA[obs.soma_joinid.values - base] = obs.raw_sum.values
        A['N'] += np.bincount(u[u >= 0], minlength=U)
        # this dataset's region only: its NS units (small temporaries)
        acc = {key: np.zeros(NS * G) for key in ['S', 'L', 'F', 'R']}; hacc = [np.zeros(NS * G) for _ in TH]
        for tbl in q.X("raw").tables():
            ci = tbl["soma_dim_0"].to_numpy(); gj = tbl["soma_dim_1"].to_numpy(); v = tbl["soma_data"].to_numpy()
            g = grpA[ci - base]; m = g >= 0
            g, gj, v, ci = g[m].astype(np.int64), gj[m], v[m], ci[m]
            idx = (g - r * NS) * G + gj; cp = v / sizeA[ci - base] * 1e4
            acc['S'] += np.bincount(idx, weights=np.log1p(cp), minlength=NS * G); acc['L'] += np.bincount(idx, weights=cp, minlength=NS * G)
            acc['F'] += np.bincount(idx, minlength=NS * G); acc['R'] += np.bincount(idx, weights=v, minlength=NS * G)
            for i, t in enumerate(TH):
                sel = v >= t; hacc[i] += np.bincount(idx[sel], minlength=NS * G)
        for key in acc: A[key][r * NS:(r + 1) * NS] += acc[key].reshape(NS, G).astype(np.float32)
        for i in range(len(TH)): A['H'][i][r * NS:(r + 1) * NS] += hacc[i].reshape(NS, G).astype(np.float32)
        del acc, hacc
    done.add(row.dataset_id)
    if len(done) % 5 == 0 or len(done) == len(ds):
        np.savez('wb_sums_tmp.npz', **A); os.replace('wb_sums_tmp.npz', 'wb_sums.npz'); open('wb_done.txt', 'w').write('\n'.join(done))
    print(len(done), '/', len(ds), REGIONS[r][0], title[:70], len(obs), f'{time.time()-t0:.0f}s', flush=True)
np.savez('wb_sums_tmp.npz', **A); os.replace('wb_sums_tmp.npz', 'wb_sums.npz'); open('wb_done.txt', 'w').write('\n'.join(done))
print('done', flush=True)
