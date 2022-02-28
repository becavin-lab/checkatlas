args <- commandArgs(trailingOnly = TRUE)
print(args)

workdir = args[1]
atlas_name = args[2]
#workdir = "~/Documents/testatlas"
#atlas_name = "PAH_675093_decontX"

setwd(workdir)
# Read Atlas in rds format
atlas_rds = paste0(atlas_name,".rds")
atlas = readRDS(atlas_rds)

if(class(atlas)[1] == "Seurat"){
  library(SeuratDisk)
  library(Signac)
  
  print("Convert to AnnData")
  # Parse to h5Seurat
  atlas_h5 = paste0(atlas_name,".h5Seurat")
  SaveH5Seurat(atlas, filename = atlas_h5, overwrite = T)
  
  # Convert to h5 scanpy
  Convert(atlas_h5, dest = "h5ad", overwrite = T)
  
  unlink(atlas_h5, recursive = FALSE, force = T)
}else{
  print(paste(atlas_name,"is not a Seurat or Signac object"))
}