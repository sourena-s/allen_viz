# Pack wbcx.py sums (per neocortical area and supercluster) as siletti_ncx(.bin.gz, _hist, _raw).
import json, gzip, struct, glob, numpy as np, pandas as pd
from wb_defs import SCS
NS = len(SCS); files = sorted(glob.glob('wbcx/*.npz'))
def label(t):
    p = t.replace('Cerebral cortex (Cx) - ', '').split(' - ')
    return f"{p[-2]} ({p[-1]})" if len(p) >= 2 else p[0]
var = pd.read_parquet('var.parquet'); nm = var.sort_values('soma_joinid').feature_name.to_numpy()
recs = []; secs = []
for d, f in enumerate(files):
    z = np.load(f); secs.append(label(str(z['title'])))
    for s in range(NS):
        n = z['N'][s]
        if n >= 30: recs.append((d, s, int(n), z['S'][s] / n, z['L'][s] / n, z['F'][s] / n, z['R'][s] / n, z['H'][:, s] / n))
M = np.stack([r[3] for r in recs]); Lm = np.stack([r[4] for r in recs]); Fr = np.stack([r[5] for r in recs]); Rm = np.stack([r[6] for r in recs]); H = np.stack([r[7] for r in recs], 1)
sel, seen = [], set()
for j in np.where(Fr.max(0) >= 0.005)[0]:
    if j >= len(nm) or nm[j] in seen: continue
    seen.add(nm[j]); sel.append(j)
sel = np.array(sel); M, Lm, Fr, Rm, H = M[:, sel], Lm[:, sel], Fr[:, sel], Rm[:, sel], H[:, :, sel]
q = lambda X, mx: np.round(X / mx * 255).clip(0, 255).astype(np.uint8)
gmax = M.max(0); gmax[gmax == 0] = 1; lmax = Lm.max(0); lmax[lmax == 0] = 1; rmax = Rm.max(0); rmax[rmax == 0] = 1
body = np.concatenate([q(M, gmax).T, np.round(Fr * 255).clip(0, 255).astype(np.uint8).T], axis=1)
header = dict(version=1, source="Siletti et al. 2023, Science 382:eadd7046 (Human Brain Cell Atlas v1.0, neocortical dissections, CELLxGENE Census 2025-01-30), CC BY 4.0; nuclei pooled per area and supercluster",
  value="mean log1p(counts per 10k) per gene, scaled to the gene's maximum; fraction of nuclei with any count",
  sections=secs, labels=secs, superclusters=SCS, groups=[[r[0], r[1], r[2]] for r in recs], genes=[str(nm[j]) for j in sel],
  max=[round(float(x), 4) for x in gmax], linear="mean counts per 10k per gene (no log), scaled to the gene's maximum (lmax)", lmax=[round(float(x), 4) for x in lmax])
hj = json.dumps(header, separators=(',', ':')).encode(); pad = (-(12 + len(hj))) % 8
open('siletti_ncx.bin.gz', 'wb').write(gzip.compress(b"SILETTI1" + struct.pack('<I', len(hj)) + hj + b"\0" * pad + body.tobytes() + np.ascontiguousarray(q(Lm, lmax).T).tobytes(), 9))
open('siletti_ncx_hist.bin.gz', 'wb').write(gzip.compress(np.ascontiguousarray(np.round(np.transpose(H, (0, 2, 1)) * 255).clip(0, 255).astype(np.uint8)).tobytes(), 6))
open('siletti_ncx_raw.bin.gz', 'wb').write(gzip.compress(rmax.astype('<f4').tobytes() + np.ascontiguousarray(q(Rm, rmax).T).tobytes(), 6))
print(len(recs), 'units', len(sel), 'genes', len(secs), 'areas')
