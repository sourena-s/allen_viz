# Per (cluster, dissection) files for the detailed view: index.json + s<k>.bin.gz, 256 genes each,
# genes in the cluster table's order (cl_sel.npy from packcl.py)
import json, gzip, os, struct, numpy as np, pandas as pd
from scdefs import ORDER
z = np.load('agg3.npz'); keys, M, Fr = z['keys'], z['means'], z['fracs']
sel = np.load('cl_sel.npy')
dsets = pd.read_csv('hipamy_datasets.csv')
dis_titles = [o[0] for o in ORDER]
pos = [dis_titles.index(t.replace('Dissection: ', '')) for t in dsets.dataset_title]
ann = pd.read_excel('cluster_annotation.xlsx').dropna(subset=['Cluster ID'])
known = set(ann['Cluster ID'].astype(int))
keep = [i for i, (c, k, n) in enumerate(keys) if c in known]
groups = [[int(keys[i][0]), pos[keys[i][1]], int(keys[i][2])] for i in keep]
M, Fr = M[keep][:, sel], Fr[keep][:, sel]
os.makedirs('siletti_cd', exist_ok=True)
open('siletti_cd/index.json.gz', 'wb').write(gzip.compress(json.dumps({"groups": groups}, separators=(',', ':')).encode(), 9))
K, PER = len(groups), 256
gmax = M.max(0); gmax[gmax == 0] = 1
qm = np.round(M / gmax * 255).astype(np.uint8); qf = np.round(Fr * 255).astype(np.uint8)
total = 0
for s0 in range(0, len(sel), PER):
    js = range(s0, min(len(sel), s0 + PER))
    parts = [b"SILCD001", struct.pack('<II', K, len(js)), np.array([gmax[j] for j in js], dtype='<f4').tobytes()]
    for j in js: parts += [qm[:, j].tobytes(), qf[:, j].tobytes()]
    data = gzip.compress(b"".join(parts), 9)
    open(f'siletti_cd/s{s0 // PER}.bin.gz', 'wb').write(data); total += len(data)
print(K, 'groups', len(sel), 'genes', (len(sel) + PER - 1) // PER, 'files', round(total / 1e6, 1), 'MB')
