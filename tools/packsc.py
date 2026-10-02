import json, gzip, struct, numpy as np, pandas as pd
from scdefs import ORDER, FIRST
var, out = pd.read_pickle('agg.pkl')
# Dissection -> (section, short label), in anatomical order
dis = [o[0] for o in ORDER]
assert set(dis) == {r['dissection'] for r in out}, set(r['dissection'] for r in out) ^ set(dis)
tot = pd.Series({}, dtype=float)
for r in out: tot[r['supercluster']] = tot.get(r['supercluster'], 0) + r['n']
scs = [s for s in FIRST if tot.get(s,0) >= 200] + [s for s in tot.sort_values(ascending=False).index if s not in FIRST and tot[s] >= 200]
keep = [r for r in out if r['supercluster'] in scs and r['n'] >= 30]
keep.sort(key=lambda r: (dis.index(r['dissection']), scs.index(r['supercluster'])))
M = np.stack([r['mean'] for r in keep]); Fr = np.stack([r['frac'] for r in keep])
names = var.feature_name.values
expressed = (Fr.max(0) >= 0.005)
# one entry per symbol (first if duplicated)
sel = [j for j in np.where(expressed)[0]]
seen=set(); sel2=[]
for j in sel:
    if names[j] in seen: continue
    seen.add(names[j]); sel2.append(j)
sel = np.array(sel2)
M = M[:, sel]; Fr = Fr[:, sel]
gmax = M.max(0); gmax[gmax == 0] = 1
qm = np.round(M / gmax * 255).astype(np.uint8); qf = np.round(Fr * 255).astype(np.uint8)
body = np.concatenate([qm.T, qf.T], axis=1)   # gene-major: means then fractions
# Linear means (counts per 10k) from agg5, appended as a third block (gene-major, one byte per group)
dsets = pd.read_csv('hipamy_datasets.csv'); titles = [t.replace('Dissection: ', '') for t in dsets.dataset_title]
parts = [np.load(f'agg5/{k}.npz') for k in range(len(dsets))]
L = np.stack([parts[titles.index(r['dissection'])]['LS'][list(parts[titles.index(r['dissection'])]['groups']).index(r['supercluster'])] / r['n'] for r in keep])[:, sel]
lmax = L.max(0); lmax[lmax == 0] = 1
ql = np.round(L / lmax * 255).astype(np.uint8)
header = dict(version=1, source="Siletti et al. 2023, Science 382:eadd7046 (Human Brain Cell Atlas v1.0, CELLxGENE Census 2025-01-30), CC BY 4.0",
  value="mean log1p(counts per 10k) per gene, scaled to the gene's maximum (max given per gene); fraction of cells with any count",
  sections=[o[1] for o in ORDER], labels=[o[2] for o in ORDER], superclusters=scs,
  groups=[[dis.index(r['dissection']), scs.index(r['supercluster']), r['n']] for r in keep],
  genes=list(names[sel]), max=[round(float(x), 4) for x in gmax],
  linear="mean counts per 10k per gene (no log), scaled to the gene's maximum (lmax); block after the means and fractions", lmax=[round(float(x), 4) for x in lmax])
hj = json.dumps(header, separators=(',',':')).encode()
pad = (-(12 + len(hj))) % 8
blob = b"SILETTI1" + struct.pack('<I', len(hj)) + hj + b"\0"*pad + body.tobytes() + np.ascontiguousarray(ql.T).tobytes()
open('siletti_hipamy.bin.gz','wb').write(gzip.compress(blob, 9))
print(len(keep), 'groups', len(sel), 'genes', len(blob)/1e6, 'MB raw')
