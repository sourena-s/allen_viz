# Brain norm reference: per gene, the whole-brain mean counts per 10k (CP10k) over all 3.37M nuclei of the
# 105 dissections (units of siletti_brain weighted by their nuclei; donors pooled, not normalised apart)
import gzip, json, struct, numpy as np
b = gzip.open('/home/user/allen_viz/assets/siletti_brain.bin.gz').read(); n = struct.unpack('<I', b[8:12])[0]; h = json.loads(b[12:12 + n])
st = ((12 + n + 7) // 8) * 8; G = len(h['genes']); U = len(h['groups']); N = np.array([g[2] for g in h['groups']], np.float64)
lin = np.frombuffer(b, np.uint8, offset=st + G * U * 2, count=G * U).reshape(G, U).astype(np.float64) / 255 * np.array(h['lmax'])[:, None]
brain = (lin * N[None]).sum(1) / N.sum()
out = dict(note="Per gene: whole-brain mean counts per 10k over all nuclei (donors pooled)", nuclei=int(N.sum()), genes=h['genes'], cp10k=[float(f'{x:.4g}') for x in brain])
open('siletti_brainref.json.gz', 'wb').write(gzip.compress(json.dumps(out, separators=(',', ':')).encode(), 9))
print(G, int(N.sum()))
