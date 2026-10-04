# Recluster: counts per 10k (by each nucleus's total), log1p, scale, PCA 50, kNN 15, Leiden;
# each cluster labelled by the human markers (positive = at least half its nuclei with any count)
import numpy as np, pandas as pd, scipy.sparse as sp, scanpy as sc, anndata as ad, sys
res = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
X = sp.load_npz('rc_X.npz'); obs = pd.read_parquet('rc_obs.parquet'); genes = pd.read_csv('rc_genes.csv').iloc[:, 0].tolist()
a = ad.AnnData(X=X, obs=obs.astype({'sc': str, 'dissection': str}), var=pd.DataFrame(index=genes))
det = (a.X > 0)
a.X = sp.diags(1e4 / a.obs.raw_sum.values) @ a.X; a.X = a.X.log1p() if hasattr(a.X, 'log1p') else np.log1p(a.X)
sc.pp.scale(a, max_value=10); sc.tl.pca(a, n_comps=50)
sc.pp.neighbors(a, n_neighbors=15); sc.tl.leiden(a, resolution=res, flavor='igraph', n_iterations=2, key_added='cl')
cl = a.obs.cl.astype(int).values; K = cl.max() + 1
gi = {g: i for i, g in enumerate(genes)}
frac = {g: np.array([det[cl == k][:, gi[g]].mean() for k in range(K)]) for g in ["FIBCD1","FNDC1","RGS14","PCP4","FN1","PROX1","HS3ST4"]}
def label(k):
    pos = lambda g: frac[g][k] >= 0.5
    if pos("PROX1") and not (pos("FIBCD1") or pos("FNDC1")): return "DG"
    if pos("FIBCD1") or pos("FNDC1"): return "CA1"
    if pos("RGS14"): return "CA2"
    if (pos("PCP4") or pos("FN1")) and not pos("PROX1"): return "SUB"
    if pos("HS3ST4"): return "CA3"
    return "none"
labs = np.array([label(k) for k in range(K)])
a.obs['subfield'] = labs[cl]
tab = pd.DataFrame({'n': np.bincount(cl), 'label': labs, **{g: frac[g].round(2) for g in frac},
                    'top_dissection': [a.obs.dissection[cl == k].value_counts().index[0] for k in range(K)],
                    'siletti': [a.obs.sc[cl == k].value_counts().index[0].replace('Hippocampal ', '') for k in range(K)]})
pd.set_option('display.width', 250); print(tab.sort_values(['label', 'n'], ascending=[True, False]).to_string())
print(a.obs.groupby('subfield').size())
print(pd.crosstab(a.obs.dissection, a.obs.subfield))
a.obs[['soma_joinid', 'dissection', 'sc', 'cl', 'subfield']].to_parquet(f'rc_labels_{res}.parquet')
