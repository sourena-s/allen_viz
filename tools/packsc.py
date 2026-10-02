import json, gzip, struct, numpy as np, pandas as pd
var, out = pd.read_pickle('agg.pkl')
# Dissection -> (section, short label), in anatomical order
ORDER = [
 ("Head of hippocampus (HiH) - Uncal DG-CA4", "Hippocampus · head", "DG–CA4"),
 ("Head of hippocampus (HiH) - Uncal CA2-CA3", "Hippocampus · head", "CA2–CA3"),
 ("Head of hippocampus (HiH) - Uncal CA1-CA3", "Hippocampus · head", "CA1–CA3"),
 ("Head of hippocampus (HiH) - Uncal CA1", "Hippocampus · head", "CA1"),
 ("Head of hippocampus (HiH) - Tail of Hippocampus (HiT) - Subicular cortex - Sub", "Hippocampus · head", "Subiculum"),
 ("Body of hippocampus (HiB) - Rostral DG-CA4", "Hippocampus · body", "DG–CA4"),
 ("Body of hippocampus (HiB) - Rostral CA3", "Hippocampus · body", "CA3"),
 ("Body of hippocampus (HiB) - Rostral CA1-CA3", "Hippocampus · body", "CA1–CA3"),
 ("Body of hippocampus (HiB) - Rostral CA1-2", "Hippocampus · body", "CA1–CA2"),
 ("Tail of Hippocampus (HiT) - Caudal Hippocampus - CA4-DGC", "Hippocampus · tail", "CA4–DG"),
 ("Tail of Hippocampus (HiT) - Caudal Hippocampus - CA1-CA3", "Hippocampus · tail", "CA1–CA3"),
 ("Cerebral cortex (Cx) - Anterior parahippocampal gyrus (AG) - Lateral entorhinal cortex - LEC", "Parahippocampal gyrus", "Lateral entorhinal"),
 ("Cerebral cortex (Cx) - Anterior parahippocampal gyrus, posterior part (APH) - Medial entorhinal cortex - MEC", "Parahippocampal gyrus", "Medial entorhinal"),
 ("Cerebral cortex (Cx) - Posterior parahippocampal gyrus (PPH) - TH-TL", "Parahippocampal gyrus", "Posterior PHG (TH–TL)"),
 ("Amygdaloid complex (AMY) - Basolateral nuclear group (BLN) - lateral nucleus - La", "Amygdala", "Lateral (La)"),
 ("Amygdaloid complex (AMY) - basolateral nuclear group (BLN) - basolateral nucleus (basal nucleus) - BL", "Amygdala", "Basolateral (BL)"),
 ("Amygdaloid complex (AMY) - basolateral nuclear group (BLN) - basomedial nucleus (accessory basal nucleus) - BM", "Amygdala", "Basomedial (BM)"),
 ("Amygdaloid complex (AMY) - Central nuclear group - CEN", "Amygdala", "Central (CEN)"),
 ("Amygdaloid complex (AMY) - corticomedial nuclear group - CMN", "Amygdala", "Corticomedial (CMN)"),
 ("Amygdaloid complex (AMY) - Corticomedial nuclear group (CMN) - anterior cortical nucleus - CoA", "Amygdala", "Anterior cortical (CoA)"),
 ("Extended amygdala (EXA) - Bed nucleus of stria terminalis and nearby - BNST", "Amygdala", "Bed nucleus (BNST)"),
]
dis = [o[0] for o in ORDER]
assert set(dis) == {r['dissection'] for r in out}, set(r['dissection'] for r in out) ^ set(dis)
tot = pd.Series({}, dtype=float)
for r in out: tot[r['supercluster']] = tot.get(r['supercluster'], 0) + r['n']
FIRST = ["Hippocampal dentate gyrus","Hippocampal CA4","Hippocampal CA1-3","Amygdala excitatory","Upper-layer intratelencephalic","Deep-layer intratelencephalic","Deep-layer near-projecting","Deep-layer corticothalamic and 6b","MGE interneuron","CGE interneuron","LAMP5-LHX6 and Chandelier","Medium spiny neuron","Eccentric medium spiny neuron","Splatter","Miscellaneous","Astrocyte","Oligodendrocyte precursor","Committed oligodendrocyte precursor","Oligodendrocyte","Microglia","Vascular","Fibroblast","Ependymal"]
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
header = dict(version=1, source="Siletti et al. 2023, Science 382:eadd7046 (Human Brain Cell Atlas v1.0, CELLxGENE Census 2025-01-30), CC BY 4.0",
  value="mean log1p(counts per 10k) per gene, scaled to the gene's maximum (max given per gene); fraction of cells with any count",
  sections=[o[1] for o in ORDER], labels=[o[2] for o in ORDER], superclusters=scs,
  groups=[[dis.index(r['dissection']), scs.index(r['supercluster']), r['n']] for r in keep],
  genes=list(names[sel]), max=[round(float(x), 4) for x in gmax])
hj = json.dumps(header, separators=(',',':')).encode()
pad = (-(12 + len(hj))) % 8
blob = b"SILETTI1" + struct.pack('<I', len(hj)) + hj + b"\0"*pad + body.tobytes()
open('siletti_hipamy.bin.gz','wb').write(gzip.compress(blob, 9))
print(len(keep), 'groups', len(sel), 'genes', len(blob)/1e6, 'MB raw')
