# convert data into iSEE instance

library(iSEE)
library(Seurat)
library(shiny)

# minimal code necessary for running app
data_dir <- "data"
dir.create(file.path(data_dir, "iSEE_files"),
           showWarnings = FALSE)

so_fn <- file.path(data_dir, "so.rds")
so <- readRDS(so_fn)
sce <- as.SingleCellExperiment(so)
app <- iSEE(sce,
            appTitle = "scRNA-seq analysis of lung regeneration",
            runLocal = TRUE)
runApp(app)

# to reload a saved state
state <- readRDS(file.path(data_dir, "iSEE_files", "iSEE_memory.rds"))
app <- iSEE(sce,
            appTitle = "scRNA-seq analysis of lung regeneration",
            initial = state$memory,
            runLocal = TRUE)
runApp(app)
