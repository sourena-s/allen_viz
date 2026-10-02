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
