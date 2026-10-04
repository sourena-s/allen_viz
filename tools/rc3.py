# Relabel the rc2.py Leiden clusters (resolution 1.0) with the revised human subfield markers.
# A marker is positive in a cluster when >= 50% of its nuclei have any count (NPNT: >= 30%).
# Checked in order: CA1, CA3, CA2, SUB, DG; clusters matching none are "none" (Unlabelled).
import json, numpy as np, pandas as pd, scipy.sparse as sp
X = sp.load_npz('rc_X.npz').tocsc()
g = pd.read_csv('rc_genes.csv', header=None)[0].tolist()[1:]
lab = pd.read_parquet('rc_labels_1.0.parquet')
MK = ['FIBCD1', 'FNDC1', 'NPNT', 'HS3ST4', 'RGS14', 'COL21A1', 'PCP4', 'FN1', 'PROX1']
det = {m: (X[:, g.index(m)].toarray().ravel() > 0) for m in MK}
def rule(f):
    p = lambda m, t=0.5: f[m] >= t
    if p('FIBCD1') or p('FNDC1'): return 'CA1'
    if p('HS3ST4') and p('NPNT', 0.3): return 'CA3'
    if p('RGS14') and p('COL21A1'): return 'CA2'
    if (p('PCP4') or p('FN1')) and not p('PROX1'): return 'SUB'
    if p('PROX1'): return 'DG'
    return 'none'
info, new = [], np.empty(len(lab), dtype=object)
for k in sorted(lab.cl.unique(), key=int):
    m = (lab.cl.values == k); f = {x: float(det[x][m].mean()) for x in MK}
    s = rule(f); new[m] = s
    info.append(dict(cl=int(k), n=int(m.sum()), label=s, old=lab.subfield.values[m][0], pct={x: round(100 * v) for x, v in f.items()}))
    print(k, m.sum(), lab.subfield.values[m][0], '->', s, ' '.join(f'{x}:{round(100*v)}' for x, v in f.items()))
lab['subfield'] = new; lab.to_parquet('rc_labels_v2.parquet')
json.dump(info, open('rc3_info.json', 'w'))
print(lab.subfield.value_counts())
