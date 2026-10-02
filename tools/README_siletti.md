# Siletti single-nucleus table (assets/siletti_hipamy.bin.gz)

Built from the Human Brain Cell Atlas v1.0 (Siletti et al. 2023, Science 382:eadd7046, CC BY 4.0)
in CELLxGENE Census release 2025-01-30, read directly from its public S3 bucket.

1. List the "Dissection: ..." datasets of the hippocampus, parahippocampal/entorhinal cortex and
   amygdala (21 datasets, ~694,000 nuclei) into `hipamy_datasets.csv`, and map every nucleus to its
   supercluster (the "Supercluster: ..." datasets, matched by `observation_joinid`) into `sc_map.parquet`.
2. `agg.py`: per dissection and supercluster, the mean of log1p(counts / total counts * 1e4) and the
   fraction of nuclei with any count, for every gene, from the raw counts.
3. `packsc.py`: keeps groups with at least 30 nuclei and genes expressed in at least 0.5% of some
   group, and writes the binary table the page reads (format described in index.html).

Requires `pip install cellxgene-census`; `cxopen.py` opens Census by its S3 path.

## Cluster table (assets/siletti_clusters.bin.gz), the "Detailed" view

4. `clmap.py`: reads `cluster_id` for every nucleus from the obs of the 31 supercluster .h5ad
   files (Census `h5ads/` folder on S3; only the obs columns are read) into `cl_map.parquet`.
5. `agg2.py`: per Siletti cluster, pooled over the 21 dissections: mean log1p(counts per 10k) and
   fraction expressing for every gene, plus the number of nuclei per dissection.
6. `packcl.py`: keeps clusters with at least 30 nuclei here, adds Siletti's names and annotations
   from `tables/cluster_annotation.xlsx` (github.com/linnarsson-lab/adult-human-brain), and picks
   four protein-coding markers per cluster (largest fraction expressing over its sibling clusters
   in the same supercluster; ribosomal and mitochondrial genes excluded). Needs HGNC's
   `gene_with_protein_product.txt` as `pc.txt`.

7. `agg5.py`: the linear means, mean(counts / total counts * 1e4) with no log, per dissection and
   supercluster and per cluster, from the same raw counts, nuclei and cell totals (about 15 minutes).
   Run it before the packers: both append these as a third block (gene-major, one byte per group or
   cluster, scaled to the gene's maximum `lmax` in the header) after the log means and fractions.
   The page uses the linear means for both the colours and co-expression.

`scdefs.py` holds the dissection order and labels and the supercluster order used by both packers.
