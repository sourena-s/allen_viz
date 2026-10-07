import cellxgene_census
def open_c():
    return cellxgene_census.open_soma(uri="s3://cellxgene-census-public-us-west-2/cell-census/2025-01-30/soma/", context=cellxgene_census.get_default_soma_context(tiledb_config={"vfs.s3.region":"us-west-2","vfs.s3.no_sign_request":"true"}))
def open_c_small():
    # small read batches (whole-brain passes: the big cortical dissections otherwise use > 12 GB)
    return cellxgene_census.open_soma(uri="s3://cellxgene-census-public-us-west-2/cell-census/2025-01-30/soma/", context=cellxgene_census.get_default_soma_context(tiledb_config={"vfs.s3.region":"us-west-2","vfs.s3.no_sign_request":"true","soma.init_buffer_bytes":64*1024**2}))
