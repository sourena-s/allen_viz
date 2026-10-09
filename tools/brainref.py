# Brain norm reference: per gene, mean and SD of the mean log1p(counts per 10k) across the whole-brain
# units of siletti_brain (12 regions x superclusters, each unit one value, as Allen z uses each sample)
import gzip, json, struct, numpy as np
b = gzip.open('/home/user/allen_viz/assets/siletti_brain.bin.gz').read(); n = struct.unpack('<I', b[8:12])[0]; h = json.loads(b[12:12 + n])
st = ((12 + n + 7) // 8) * 8; G = len(h['genes']); U = len(h['groups'])
body = np.frombuffer(b, np.uint8, offset=st, count=G * U * 2).reshape(G, 2 * U)
M = body[:, :U].astype(np.float32) / 255 * np.array(h['max'], np.float32)[:, None]
m = M.mean(1); s = M.std(1)
out = dict(note="Per gene: mean and SD of mean log1p(counts per 10k) over the whole-brain units (region x supercluster) of siletti_brain", units=U,
           genes=h['genes'], m=[round(float(x), 4) for x in m], s=[round(float(x), 4) for x in s])
open('siletti_brainref.json.gz', 'wb').write(gzip.compress(json.dumps(out, separators=(',', ':')).encode(), 9))
print(G, U, 'genes, units')
