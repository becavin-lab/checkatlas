# create personal library
dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE)
# add to the path
.libPaths(Sys.getenv("R_LIBS_USER"))
install.packages("Seurat")