# prep cellXgene

library(zellkonverter)
library(Seurat)
library(SingleCellExperiment)


data_dir <- "data"
dir.create(file.path(data_dir, "cellxgene_files"),
           showWarnings = FALSE)

so_fn <- file.path(data_dir, "so.rds")
so <- readRDS(so_fn)
sce <- as.SingleCellExperiment(so)

# need to rename projections to match scanpy format (e.g "X_pca", "X_umap")
reducedDimNames(sce) <- str_c("X_", tolower(reducedDimNames(sce)))

writeH5AD(sce,
          X_name = "logcounts",
          file = file.path(data_dir, "cellxgene_files", "so.h5ad"))


#cellxgene launch data/cellxgene_files/so.h5ad --open
