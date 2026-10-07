# Pack wb.py sums into the brief-file format (SILETTI1 with the linear block) as siletti_brain.bin.gz,
# with its count-histogram (_hist) and raw-mean (_raw) tables in the same gene and unit order.
import json, gzip, struct, numpy as np, pandas as pd
from wb_defs import REGIONS, SCS
z = np.load('wb_sums.npz'); N = z['N']; NS = len(SCS)
var = pd.read_parquet('var.parquet'); names = var.sort_values('soma_joinid').set_index('soma_joinid').feature_name
units = [u for u in range(len(N)) if N[u] >= 30]
n = N[units][:, None]
M = z['S'][units] / n; Lm = z['L'][units] / n; Fr = z['F'][units] / n; Rm = z['R'][units] / n; H = z['H'][:, units] / n[None]
gj_all = names.index.to_numpy(); nm = names.to_numpy()
sel, seen = [], set()
for j in np.where(Fr.max(0) >= 0.005)[0]:
    if j >= len(nm) or nm[j] in seen: continue
    seen.add(nm[j]); sel.append(j)
sel = np.array(sel)
M, Lm, Fr, Rm, H = M[:, sel], Lm[:, sel], Fr[:, sel], Rm[:, sel], H[:, :, sel]
q = lambda X, mx: np.round(X / mx * 255).clip(0, 255).astype(np.uint8)
gmax = M.max(0); gmax[gmax == 0] = 1; lmax = Lm.max(0); lmax[lmax == 0] = 1; rmax = Rm.max(0); rmax[rmax == 0] = 1
body = np.concatenate([q(M, gmax).T, np.round(Fr * 255).clip(0, 255).astype(np.uint8).T], axis=1)
header = dict(version=1, source="Siletti et al. 2023, Science 382:eadd7046 (Human Brain Cell Atlas v1.0, all 105 dissections, CELLxGENE Census 2025-01-30), CC BY 4.0; nuclei pooled per brain region and supercluster",
  value="mean log1p(counts per 10k) per gene, scaled to the gene's maximum; fraction of nuclei with any count",
  sections=[r[0] for r in REGIONS], labels=[r[0] for r in REGIONS], superclusters=SCS,
  groups=[[int(u // NS), int(u % NS), int(N[u])] for u in units], genes=[str(nm[j]) for j in sel], max=[round(float(x), 4) for x in gmax],
  linear="mean counts per 10k per gene (no log), scaled to the gene's maximum (lmax); block after the means and fractions", lmax=[round(float(x), 4) for x in lmax])
hj = json.dumps(header, separators=(',', ':')).encode(); pad = (-(12 + len(hj))) % 8
blob = b"SILETTI1" + struct.pack('<I', len(hj)) + hj + b"\0" * pad + body.tobytes() + np.ascontiguousarray(q(Lm, lmax).T).tobytes()
open('siletti_brain.bin.gz', 'wb').write(gzip.compress(blob, 9))
Hq = np.round(np.transpose(H, (0, 2, 1)) * 255).clip(0, 255).astype(np.uint8)   # [bin][gene][unit]
open('siletti_brain_hist.bin.gz', 'wb').write(gzip.compress(np.ascontiguousarray(Hq).tobytes(), 6))
open('siletti_brain_raw.bin.gz', 'wb').write(gzip.compress(rmax.astype('<f4').tobytes() + np.ascontiguousarray(q(Rm, rmax).T).tobytes(), 6))
print(len(units), 'units', len(sel), 'genes')
