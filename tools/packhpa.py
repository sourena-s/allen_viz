# Protein Atlas brain clusters (neurons c-0..c-94, non-neurons c-0..c-31) as assigned per nucleus by
# agg6.py: per cluster pooled over the 21 dissections, in the clusters-file format (SILCLUS1)
import json, gzip, struct, re, numpy as np, pandas as pd, csv
from scdefs import ORDER
meta = pd.read_csv('hpa/meta.tsv', sep='\t'); K = len(meta)
dsets = pd.read_csv('hipamy_datasets.csv')
dis = [o[0] for o in ORDER]
col_of = [dis.index(t.replace('Dissection: ', '')) for t in dsets.dataset_title]
var = pd.read_parquet('var.parquet'); names = var.feature_name.values
S = L = F = None; cnt = np.zeros((K, len(dis)), dtype=np.int64)
for k in range(len(dsets)):
    z = np.load(f'agg6/{k}.npz')
    S = z['S'][:K].astype(np.float64) if S is None else S + z['S'][:K]
    L = z['L'][:K].astype(np.float64) if L is None else L + z['L'][:K]
    F = z['F'][:K].astype(np.float64) if F is None else F + z['F'][:K]
    cnt[:, col_of[k]] += z['n'][:K]
n = cnt.sum(1); keep = [i for i in range(K) if n[i] >= 30]
M = S[keep] / n[keep, None]; Fr = F[keep] / n[keep, None]; Lm = L[keep] / n[keep, None]
sel, seen = [], set()
for j in np.where(Fr.max(0) >= 0.005)[0]:
    if names[j] in seen: continue
    seen.add(names[j]); sel.append(j)
sel = np.array(sel); M, Fr, Lm = M[:, sel], Fr[:, sel], Lm[:, sel]; gnames = names[sel]
pc = {r['symbol'] for r in csv.DictReader(open('pc.txt'), delimiter='\t') if r['status'] == 'Approved'}
ispc = np.array([g in pc and not re.match(r'^(RPL|RPS|MRPL|MRPS|MT-)', g) for g in gnames])
cap = lambda s: s[:1].upper() + s[1:]
clusters = []
for a, i in enumerate(keep):
    r = meta.iloc[i]
    sc = Fr[a] - np.delete(Fr, a, 0).mean(0); sc[~ispc | (Fr[a] < 0.3)] = -1
    top = [gnames[j] for j in np.argsort(-sc)[:4] if sc[j] > 0.1]
    w = cnt[i] / max(1, cnt[i].sum())
    tag = 'N' if r.key.startswith('N') else 'G'
    clusters.append(dict(id=int(i), name=f"{r.cl}{'' if tag == 'N' else ' (non-neuronal)'}", sc=cap(str(r.ct)), nt='', mtg=cap(str(r.detail)), n=int(n[i]),
                         markers=top, where=[int(round(x * 255)) for x in w], canon=[], rel=str(r.rel)))
gmax = M.max(0); gmax[gmax == 0] = 1; lmax = Lm.max(0); lmax[lmax == 0] = 1
qm = np.round(M / gmax * 255).astype(np.uint8); qf = np.round(Fr * 255).astype(np.uint8); ql = np.round(Lm / lmax * 255).astype(np.uint8)
body = np.concatenate([qm.T, qf.T], axis=1)
header = dict(version=1, source="Protein Atlas brain clusters (brain neurons / non-neurons, from Siletti et al. 2023) assigned per nucleus by marker-gene correlation; nuclei of the hippocampal, parahippocampal and amygdala dissections (CELLxGENE Census 2025-01-30); HPA data CC BY-SA 4.0",
  sections=[o[1] for o in ORDER], labels=[o[2] for o in ORDER], clusters=clusters, genes=list(gnames), max=[round(float(x), 4) for x in gmax],
  linear="mean counts per 10k per gene (no log), scaled to the gene's maximum (lmax); block after the means and fractions", lmax=[round(float(x), 4) for x in lmax],
  unassigned=int(sum(np.load(f'agg6/{k}.npz')['n'][K] for k in range(len(dsets)))))
hj = json.dumps(header, separators=(',', ':')).encode(); pad = (-(12 + len(hj))) % 8
blob = b"SILCLUS1" + struct.pack('<I', len(hj)) + hj + b"\0" * pad + body.tobytes() + np.ascontiguousarray(ql.T).tobytes()
open('siletti_hpa.bin.gz', 'wb').write(gzip.compress(blob, 9))
print(len(keep), 'clusters kept of', K, len(sel), 'genes', len(blob) / 1e6, 'MB raw; unassigned', header['unassigned'])
for c in sorted(clusters, key=lambda c: -c['n'])[:25]: print(c['name'], c['sc'], '|', c['mtg'], c['n'], c['markers'])
