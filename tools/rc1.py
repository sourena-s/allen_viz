# Hippocampal principal neurons (Siletti CA1-3, CA4, dentate gyrus superclusters) of the 11
# hippocampal dissections: raw counts of ~3,000 informative genes per nucleus, for reclustering
import sys, os, time, gzip, json, struct, numpy as np, pandas as pd, scipy.sparse as sp
sys.path.insert(0,'.')
from cxopen import open_c
import tiledbsoma as soma
from scdefs import ORDER
MARK = ["FIBCD1","FNDC1","RGS14","HS3ST4","PCP4","FN1","PROX1","WFS1","AMIGO2","CALB1","SLC17A7"]
# informative genes: variance of log means across the hippocampal principal clusters
b = gzip.open('/home/user/allen_viz/assets/siletti_clusters.bin.gz').read(); n = struct.unpack('<I', b[8:12])[0]; h = json.loads(b[12:12+n])
st = ((12+n+7)//8)*8; nU = len(h['clusters']); G = len(h['genes'])
body = np.frombuffer(b, np.uint8, offset=st, count=G*nU*2).reshape(G, 2*nU)
ks = [k for k, c in enumerate(h['clusters']) if c['sc'].startswith('Hippocampal')]
M = body[:, ks].astype(np.float32) / 255 * np.array(h['max'], np.float32)[:, None]; F = body[:, [nU + k for k in ks]] / 255
ok = (F.max(1) >= 0.1) & ~np.char.startswith(np.array(h['genes']), 'ENSG')
v = np.where(ok, M.var(1) / (M.mean(1) + 0.05), 0)
genes = [h['genes'][i] for i in np.argsort(-v)[:3000]] + MARK
var = pd.read_parquet('var.parquet'); jid = var.drop_duplicates('feature_name').set_index('feature_name').soma_joinid
genes = [g for g in dict.fromkeys(genes) if g in jid.index]; gj = jid.loc[genes].values
print(len(genes), 'genes', flush=True)
sc = pd.read_parquet('sc_map.parquet').set_index('observation_joinid').sc
dsets = pd.read_csv('hipamy_datasets.csv'); hip = {o[0] for o in ORDER[:11]}
c = open_c(); exp = c["census_data"]["homo_sapiens"]
blocks, obsall = [], []
for k, row in dsets.iterrows():
    title = row.dataset_title.replace('Dissection: ', '')
    if title not in hip: continue
    t0 = time.time()
    with exp.axis_query("RNA", obs_query=soma.AxisQuery(value_filter=f"dataset_id == '{row.dataset_id}'"), var_query=soma.AxisQuery(coords=(gj,))) as q:
        obs = q.obs(column_names=["soma_joinid","observation_joinid","raw_sum"]).concat().to_pandas()
        obs['sc'] = obs.observation_joinid.map(sc)
        obs = obs[obs.sc.isin(["Hippocampal CA1-3","Hippocampal CA4","Hippocampal dentate gyrus"])].reset_index(drop=True)
        pos = pd.Series(np.arange(len(obs)), index=obs.soma_joinid.values); col = pd.Series(np.arange(len(gj)), index=gj)
        R, C, V = [], [], []
        for tbl in q.X("raw").tables():
            ci = pos.reindex(tbl["soma_dim_0"].to_numpy()).to_numpy(); ok = ~np.isnan(ci)
            R.append(ci[ok].astype(np.int32)); C.append(col.reindex(tbl["soma_dim_1"].to_numpy()[ok]).to_numpy().astype(np.int32)); V.append(tbl["soma_data"].to_numpy()[ok].astype(np.float32))
        X = sp.csr_matrix((np.concatenate(V), (np.concatenate(R), np.concatenate(C))), shape=(len(obs), len(gj)))
    obs['dissection'] = title; blocks.append(X); obsall.append(obs)
    print(title, len(obs), 'nuclei', X.nnz, 'nnz', f'{time.time()-t0:.0f}s', flush=True)
X = sp.vstack(blocks).tocsr(); obs = pd.concat(obsall, ignore_index=True)
sp.save_npz('rc_X.npz', X); obs.to_parquet('rc_obs.parquet'); pd.Series(genes).to_csv('rc_genes.csv', index=False)
print('done', X.shape, flush=True)
