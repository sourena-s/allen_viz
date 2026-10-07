# Pack wbcl.py sums (neocortex per Siletti cluster) as siletti_ncxcl(.bin.gz, _hist, _raw). Each
# cluster is a "supercluster" column of the brief format; scOf gives its Siletti supercluster.
import json, gzip, struct, numpy as np, pandas as pd
z = np.load('wbcl_sums.npz'); N = z['N']
m = pd.read_parquet('cl_map.parquet').merge(pd.read_parquet('sc_map.parquet'), on='observation_joinid')
sc_of = m.groupby('cluster').sc.agg(lambda s: s.value_counts().index[0])
b = gzip.open('/home/user/allen_viz/assets/siletti_clusters.bin.gz').read(); n0 = struct.unpack('<I', b[8:12])[0]; h0 = json.loads(b[12:12 + n0])
meta = {c['id']: c for c in h0['clusters']}
var = pd.read_parquet('var.parquet'); nm = var.sort_values('soma_joinid').feature_name.to_numpy()
units = [u for u in range(len(N)) if N[u] >= 30]; n = N[units][:, None]
M = z['S'][units] / n; Lm = z['L'][units] / n; Fr = z['F'][units] / n; Rm = z['R'][units] / n; H = z['H'][:, units] / n[None]
sel, seen = [], set()
for j in np.where(Fr.max(0) >= 0.005)[0]:
    if j >= len(nm) or nm[j] in seen: continue
    seen.add(nm[j]); sel.append(j)
sel = np.array(sel); M, Lm, Fr, Rm, H = M[:, sel], Lm[:, sel], Fr[:, sel], Rm[:, sel], H[:, :, sel]
q = lambda X, mx: np.round(X / mx * 255).clip(0, 255).astype(np.uint8)
gmax = M.max(0); gmax[gmax == 0] = 1; lmax = Lm.max(0); lmax[lmax == 0] = 1; rmax = Rm.max(0); rmax[rmax == 0] = 1
names = [meta[u]['name'] if u in meta else f"c{u}" for u in units]
header = dict(version=1, source="Siletti et al. 2023, Science 382:eadd7046 (Human Brain Cell Atlas v1.0), the 25 neocortical dissections pooled, per Siletti cluster (CELLxGENE Census 2025-01-30), CC BY 4.0",
  value="mean log1p(counts per 10k) per gene, scaled to the gene's maximum; fraction of nuclei with any count",
  sections=["Neocortex"], labels=["Neocortex"], superclusters=names, groups=[[0, i, int(N[u])] for i, u in enumerate(units)],
  scOf=[str(sc_of.get(u, "")) for u in units], clId=[int(u) for u in units],
  clNt=[meta.get(u, {}).get('nt', '') for u in units], clMarkers=[meta.get(u, {}).get('markers', []) for u in units],
  genes=[str(nm[j]) for j in sel], max=[round(float(x), 4) for x in gmax], linear="mean counts per 10k (lmax)", lmax=[round(float(x), 4) for x in lmax])
body = np.concatenate([q(M, gmax).T, np.round(Fr * 255).clip(0, 255).astype(np.uint8).T], axis=1)
hj = json.dumps(header, separators=(',', ':')).encode(); pad = (-(12 + len(hj))) % 8
open('siletti_ncxcl.bin.gz', 'wb').write(gzip.compress(b"SILETTI1" + struct.pack('<I', len(hj)) + hj + b"\0" * pad + body.tobytes() + np.ascontiguousarray(q(Lm, lmax).T).tobytes(), 9))
open('siletti_ncxcl_hist.bin.gz', 'wb').write(gzip.compress(np.ascontiguousarray(np.round(np.transpose(H, (0, 2, 1)) * 255).clip(0, 255).astype(np.uint8)).tobytes(), 6))
open('siletti_ncxcl_raw.bin.gz', 'wb').write(gzip.compress(rmax.astype('<f4').tobytes() + np.ascontiguousarray(q(Rm, rmax).T).tobytes(), 6))
print(len(units), 'clusters', len(sel), 'genes', sum(1 for u in units if u in meta), 'named from the cluster table')
