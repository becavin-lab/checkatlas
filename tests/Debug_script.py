import scanpy

atlas_path = (
    "/Users/christophebecavin/Documents/testatlas/multiome.integrated.h5ad"
)
adata = scanpy.read_h5ad(atlas_path)
print(";".join(list(adata.obs.columns)))
print(";".join(list(adata.obsm_keys())))
print(";".join(list(adata.var_keys())))
print(";".join(list(adata.uns_keys())))
print(";".join(list(adata.varm_keys())))
print(adata)
