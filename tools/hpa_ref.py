import numpy as np, pandas as pd
d = pd.read_csv('hpa/brain_cl.tsv', sep='\t', header=None, names=['ens','sym','tis','cl','ct','reads','ncpm'])
d['key'] = np.where(d.tis=='brain neurons','N','G') + d.cl
meta = pd.read_csv('hpa/brain_clusters.tsv', sep='\t', header=None, names=['tis','cl','ct','detail','cls','n','inc','rel'])
meta['key'] = np.where(meta.tis=='brain neurons','N','G') + meta.cl
# all clusters, also those the Protein Atlas leaves out of its aggregation (few cells, e.g. c-89 CALB2+ RELN+)
M = d.pivot_table(index='ens', columns='key', values='ncpm', aggfunc='first').fillna(0)
M = M[[k for k in meta.key if k in M.columns]]
L = np.log1p(M / 100)          # log1p(counts per 10k)
# markers: per cluster, top genes by log-ratio to the mean of the others, expressed >= 10 cp10k... (nCPM>=50)
mk = set()
for k in L.columns:
    other = L.drop(columns=k).mean(1); sc = L[k] - other
    sc = sc[(M[k] >= 50)]
    mk |= set(sc.sort_values(ascending=False).index[:50])
print(len(meta), 'clusters', L.shape, 'markers', len(mk))
meta.to_csv('hpa/meta.tsv', sep='\t', index=False)
L.loc[sorted(mk)].to_parquet('hpa/ref_markers.parquet')
