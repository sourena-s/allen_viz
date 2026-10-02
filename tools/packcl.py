import json, gzip, struct, numpy as np, pandas as pd, csv
z = np.load('agg2.npz'); S, F, counts, n = z['S'], z['F'], z['counts'], z['n']
var = pd.read_parquet('var.parquet'); names = var.feature_name.values
ann = pd.read_excel('cluster_annotation.xlsx').dropna(subset=['Cluster ID']).set_index('Cluster ID')
pc = {r['symbol'] for r in csv.DictReader(open('pc.txt'), delimiter='\t') if r['status'] == 'Approved'}
dsets = pd.read_csv('hipamy_datasets.csv')
from scdefs import ORDER, FIRST
dis_titles = [o[0] for o in ORDER]
col_of = [dis_titles.index(t.replace('Dissection: ', '')) for t in dsets.dataset_title]   # dataset row -> ORDER position
cnt = np.zeros_like(counts); cnt[:, col_of] = counts
keep = [c for c in range(len(n)) if n[c] >= 30 and c in ann.index]
sc_order = FIRST
def sc_rank(s): return sc_order.index(s) if s in sc_order else len(sc_order)
keep.sort(key=lambda c: (sc_rank(ann.loc[c, 'Supercluster']), -n[c]))
M = S[keep] / n[keep, None]; Fr = F[keep] / n[keep, None]
expressed = Fr.max(0) >= 0.005
sel, seen = [], set()
for j in np.where(expressed)[0]:
    if names[j] in seen: continue
    seen.add(names[j]); sel.append(j)
sel = np.array(sel)
M, Fr = M[:, sel], Fr[:, sel]; gnames = names[sel]
import re
ispc = np.array([g in pc and not re.match(r'^(RPL|RPS|MRPL|MRPS|MT-)', g) for g in gnames])   # no ribosomal / mitochondrial markers
# Markers: protein-coding genes most specific to the cluster within its supercluster (fraction
# expressing minus the siblings' mean), expressed in at least 30% of its nuclei
clusters = []
scs = [ann.loc[c, 'Supercluster'] for c in keep]
for i, c in enumerate(keep):
    sib = [k for k in range(len(keep)) if scs[k] == scs[i] and k != i]
    ref = Fr[sib].mean(0) if sib else np.zeros(Fr.shape[1])
    other = Fr[[k for k in range(len(keep)) if scs[k] != scs[i]]].mean(0)
    score = Fr[i] - (ref if sib else other)
    score[~ispc | (Fr[i] < 0.3)] = -1
    top = [gnames[j] for j in np.argsort(-score)[:4] if score[j] > 0.1]
    w = cnt[c] / max(1, cnt[c].sum())
    a = ann.loc[c]
    clusters.append(dict(id=int(c), name=str(a['Cluster name']), sc=str(a['Supercluster']), nt=str(a['Neurotransmitter auto-annotation']) if pd.notna(a['Neurotransmitter auto-annotation']) else "",
                         mtg=str(a['Transferred MTG Label']) if pd.notna(a['Transferred MTG Label']) else "", n=int(n[c]), markers=top, where=[int(round(x * 255)) for x in w]))
gmax = M.max(0); gmax[gmax == 0] = 1
qm = np.round(M / gmax * 255).astype(np.uint8); qf = np.round(Fr * 255).astype(np.uint8)
body = np.concatenate([qm.T, qf.T], axis=1)
header = dict(version=1, source="Siletti et al. 2023 clusters (Human Brain Cell Atlas v1.0, CELLxGENE Census 2025-01-30), CC BY 4.0; nuclei of the hippocampal, parahippocampal and amygdala dissections only",
  sections=[o[1] for o in ORDER], labels=[o[2] for o in ORDER], clusters=clusters, genes=list(gnames), max=[round(float(x), 4) for x in gmax])
hj = json.dumps(header, separators=(',', ':')).encode()
pad = (-(12 + len(hj))) % 8
blob = b"SILCLUS1" + struct.pack('<I', len(hj)) + hj + b"\0" * pad + body.tobytes()
open('siletti_clusters.bin.gz', 'wb').write(gzip.compress(blob, 9))
print(len(keep), 'clusters', len(sel), 'genes', len(blob) / 1e6, 'MB raw')
for cdef in [c for c in clusters if c['sc'].startswith('Hippocampal CA1') or c['sc']=='Amygdala excitatory'][:14]: print(cdef['name'], cdef['n'], cdef['markers'], np.argmax(cdef['where']))
