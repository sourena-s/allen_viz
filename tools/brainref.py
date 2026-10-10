# Brain norm reference: per gene, the whole-brain mean raw count per nucleus (all 3.37M nuclei of the 105
# dissections, i.e. the brain-wide raw total / nuclei), from siletti_brain_raw (units weighted by nuclei)
import gzip, json, struct, numpy as np
b = gzip.open('/home/user/allen_viz/assets/siletti_brain.bin.gz').read(); n = struct.unpack('<I', b[8:12])[0]; h = json.loads(b[12:12 + n])
G = len(h['genes']); U = len(h['groups']); N = np.array([g[2] for g in h['groups']], np.float64)
r = gzip.open('/home/user/allen_viz/assets/siletti_brain_raw.bin.gz').read()
mx = np.frombuffer(r, '<f4', count=G); Q = np.frombuffer(r, np.uint8, offset=G * 4, count=G * U).reshape(G, U)
R = Q.astype(np.float64) / 255 * mx[:, None]
brain = (R * N[None]).sum(1) / N.sum()
out = dict(note="Per gene: whole-brain mean raw count per nucleus (brain-wide raw total / nuclei, all regions)", nuclei=int(N.sum()),
           genes=h['genes'], raw=[float(f'{x:.4g}') for x in brain])
open('siletti_brainref.json.gz', 'wb').write(gzip.compress(json.dumps(out, separators=(',', ':')).encode(), 9))
print(G, int(N.sum()))
