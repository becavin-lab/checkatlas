library(Seurat)
atlas_path = "/home/becavin/checkatlas/tests/data/pbmc3k_seurat.rds"
seurat = readRDS(atlas_path)
seurat <- RunTSNE(seurat, dims = 1:30)
gene_selected = sample(rownames(seurat), 1000)
cell_selected = sample(colnames(seurat), 1000)
seurat = seurat[gene_selected, cell_selected]
seurat = UpdateSeuratObject(seurat)
saveRDS(seurat, file = "/home/becavin/checkatlas/tests/data/pbmc3k_seurat.rds")
