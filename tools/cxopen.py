import cellxgene_census
def open_c():
    return cellxgene_census.open_soma(uri="s3://cellxgene-census-public-us-west-2/cell-census/2025-01-30/soma/", context=cellxgene_census.get_default_soma_context(tiledb_config={"vfs.s3.region":"us-west-2","vfs.s3.no_sign_request":"true"}))
