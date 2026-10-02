import s3fs, h5py, numpy as np, pandas as pd
fs = s3fs.S3FileSystem(anon=True, client_kwargs={"region_name":"us-west-2"})
sup = pd.read_csv('superclusters.csv')
parts = []
for _, r in sup.iterrows():
    p = f"cellxgene-census-public-us-west-2/cell-census/2025-01-30/h5ads/{r.dataset_id}.h5ad"
    with fs.open(p, 'rb', block_size=2**22) as f, h5py.File(f, 'r') as h:
        o = h['obs']
        j = o['observation_joinid'][:]
        g = o['cluster_id']
        cl = g['categories'][:][g['codes'][:]] if isinstance(g, h5py.Group) else g[:]
        parts.append(pd.DataFrame({'observation_joinid': [x.decode() if isinstance(x, bytes) else x for x in j], 'cluster': cl.astype(int)}))
    print(r.dataset_title, len(j), flush=True)
pd.concat(parts).to_parquet('cl_map.parquet')
print('done')
