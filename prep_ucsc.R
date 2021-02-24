# prep UCSC cellbrowser

library(scbp)
library(presto)
library(dplyr)

data_dir <- "data"
cb_data <- file.path(data_dir, "cellbrowser_files")
dir.create(cb_data, showWarnings = FALSE)
so_fn <- file.path(data_dir, "so.rds")
so <- readRDS(so_fn)

# metadata columns to retain
cols_to_keep <- colnames(so@meta.data)
names(cols_to_keep) <- cols_to_keep

# make marker table
wilcoxauc(so, "Cell_types") %>%
  filter(logFC > 0,
         padj < 0.01) %>%
  group_by(group) %>%
  arrange(padj, logFC, .by_group = TRUE) %>%
  dplyr::slice(1:100) %>%
  write_tsv(file.path(cb_data, "mkrs.tsv"))

make_cellbrowser(so,
                 column_list = cols_to_keep,
                 project = "lung",
                 outdir = file.path(cb_data, "cb"),
                 marker_file = file.path(cb_data, "mkrs.tsv"),
                 ident = "Cell_types")

datasets <- file.path(cb_data, "cb", "lung", "cellbrowser.conf")

build_cellbrowser(datasets,
                  outdir = file.path(cb_data, "lung-cellbrowser"),
                  cbBuild_path = "/miniconda3/bin/cbBuild")

# to show in local browser
/miniconda3/bin/cbBuild \
  -i data/cellbrowser_files/cb/lung/cellbrowser.conf \
  -o data/cellbrowser_files/lung-cellbrowser \
  -p 8888
