# Reclustered hippocampal subfields (rc2.py labels, agg7.py sums): groups = (dissection, subfield),
# in the brief-file format (SILETTI1) with the linear block
import os, json, gzip, struct, glob, numpy as np, pandas as pd
from scdefs import ORDER
LABS = ["CA1", "CA2", "CA3", "SUB", "DG", "none"]
NAMES = ["CA1 (FIBCD1+ or FNDC1+)", "CA2 (RGS14+ and COL21A1+)", "CA3 (HS3ST4+ and NPNT+)", "Subiculum (PCP4+ or FN1+, PROX1−)", "Dentate gyrus granule (PROX1+)", "Unlabelled (none of these markers)"]
dis = [o[0] for o in ORDER]
# Rules (checked in this order; a marker is positive when >= 50% of a cluster's nuclei have any count, NPNT >= 30%)
RULES = dict(CA1="FIBCD1+ or FNDC1+ (checked 1st)", CA3="HS3ST4+ and NPNT+ (NPNT in >= 30% of nuclei; checked 2nd)", CA2="RGS14+ and COL21A1+ (checked 3rd)",
  SUB="PCP4+ or FN1+, and PROX1- (checked 4th)", DG="PROX1+ (checked 5th)", none="none of the rules above")
INFO = json.load(open(os.environ.get('INFO', 'rc3_info.json')))
var = pd.read_parquet('var.parquet'); names = var.feature_name.values
recs = []
for f in sorted(glob.glob(os.environ.get('AGG', 'agg7') + '/*.npz')):
    z = np.load(f); d = dis.index(str(z['title']))
    for i in range(len(LABS)):
        if z['n'][i] >= 30: recs.append((d, i, int(z['n'][i]), z['S'][i] / z['n'][i], z['L'][i] / z['n'][i], z['F'][i] / z['n'][i]))
recs.sort(key=lambda r: (r[0], r[1]))
M = np.stack([r[3] for r in recs]); Lm = np.stack([r[4] for r in recs]); Fr = np.stack([r[5] for r in recs])
sel, seen = [], set()
for j in np.where(Fr.max(0) >= 0.005)[0]:
    if names[j] in seen: continue
    seen.add(names[j]); sel.append(j)
sel = np.array(sel); M, Lm, Fr = M[:, sel], Lm[:, sel], Fr[:, sel]
gmax = M.max(0); gmax[gmax == 0] = 1; lmax = Lm.max(0); lmax[lmax == 0] = 1
qm = np.round(M / gmax * 255).astype(np.uint8); qf = np.round(Fr * 255).astype(np.uint8); ql = np.round(Lm / lmax * 255).astype(np.uint8)
body = np.concatenate([qm.T, qf.T], axis=1)
header = dict(version=1, source="Hippocampal principal neurons of Siletti et al. 2023 (CELLxGENE Census 2025-01-30), reclustered (tools/rc2.py) and labelled by human subfield markers",
  value="mean log1p(counts per 10k) per gene, scaled to the gene's maximum; fraction of nuclei with any count",
  sections=[o[1] for o in ORDER], labels=[o[2] for o in ORDER], superclusters=NAMES, keys=LABS, rules=RULES, rcClusters=INFO,
  groups=[[r[0], r[1], r[2]] for r in recs], genes=list(names[sel]), max=[round(float(x), 4) for x in gmax],
  linear="mean counts per 10k per gene (no log), scaled to the gene's maximum (lmax); block after the means and fractions", lmax=[round(float(x), 4) for x in lmax])
hj = json.dumps(header, separators=(',', ':')).encode(); pad = (-(12 + len(hj))) % 8
blob = b"SILETTI1" + struct.pack('<I', len(hj)) + hj + b"\0" * pad + body.tobytes() + np.ascontiguousarray(ql.T).tobytes()
open('siletti_custom.bin.gz', 'wb').write(gzip.compress(blob, 9))
print(len(recs), 'groups', len(sel), 'genes', len(blob) / 1e6, 'MB raw')
for i, l in enumerate(NAMES): print(l, sum(r[2] for r in recs if r[1] == i))
