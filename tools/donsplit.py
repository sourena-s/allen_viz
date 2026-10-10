# Split <file>_donor.bin.gz (all donors) into <file>_donor<x>.bin.gz, one donor each (GitHub's 100 MB file limit):
# "DON1", uint32 D=1, uint32 NB, float32 nuclei[U], float32 rawMax[G], float32 cpMax[G], uint8 [G][U] tables as before
import gzip, struct, json, sys, numpy as np
A = '/home/user/allen_viz/assets/'
for name in ['siletti_hipamy', 'siletti_custom', 'siletti_clusters', 'siletti_hpa', 'siletti_brain', 'siletti_ncx', 'siletti_ncxcl']:
    hb = gzip.open(A + name + '.bin.gz').read(); n = struct.unpack('<I', hb[8:12])[0]; h = json.loads(hb[12:12 + n])
    b = gzip.open(name + '_donor.bin.gz').read(); D, NB = struct.unpack('<II', b[4:12]); G = len(h['genes'])
    o = 12; nN = (len(b) - 12 - 2 * D * G * 4) // 1   # solve U from the length
    U = (len(b) - 12 - 2 * D * G * 4) // (D * 4 + (3 + NB) * D * G); assert 12 + U * D * 4 + 2 * D * G * 4 + (3 + NB) * D * G * U == len(b), name
    N = np.frombuffer(b, '<f4', U * D, o).reshape(U, D); o += U * D * 4
    rmx = np.frombuffer(b, '<f4', D * G, o).reshape(D, G); o += D * G * 4
    cmx = np.frombuffer(b, '<f4', D * G, o).reshape(D, G); o += D * G * 4
    T = np.frombuffer(b, np.uint8, (3 + NB) * D * G * U, o).reshape(3 + NB, D, G * U)
    for x in range(D):
        out = b'DON1' + struct.pack('<II', 1, NB) + N[:, x].astype('<f4').tobytes() + rmx[x].tobytes() + cmx[x].tobytes() + T[:, x].tobytes()
        open(A + f'{name}_donor{x}.bin.gz', 'wb').write(gzip.compress(out, 9))
    print(name, U, D, flush=True)
