# prep for loupe
library(Seurat)
library(tidyverse)

# example command to aggregate multiple samples
# cellranger aggr \
# --id aggregated_data \
# --csv aggregation.csv \
# --jobmode lsf.template \
# --normalize none \
# --maxjobs 12

data_dir <- "data"
dir.create(file.path(data_dir, "loupe_files"),
           showWarnings = FALSE)
so_fn <- file.path(data_dir, "so.rds")
so <- readRDS(so_fn)

# need to add 10x matching cell barcodes to metadata
# e.g.
# AAACATACGGTACT-1

# use aggregation_csv.csv to see what order libraries were added (which determines suffix integer value)

aggr_csv <- read_csv("data/aggregated_data/outs/aggregation.csv")

# no matching sample ids so need to match manually
samples <- c("ATI1expt1",
             "ATI2expt1",
             "ATII1expt1",
             "ATIIinjured2expt1",
             "ATIIinjured1expt1",
             "ATIIinjured1expt2",
             "ATIIinjured2expt2",
             "ATI1expt2",
             "ATI2expt2",
             "ATII1expt2")

names(samples) <- c(
  "3161-ATII-1",
  "3161-ATII-2",
  "3161-ATII-3",
  "3162-ATII-5",
  "3162-ATTII-4",
  "3241ATII4",
  "3241ATII5",
  "3242ATI1",
  "3242ATI2",
  "3242ATII3")

aggr_csv <- mutate(aggr_csv, sample_id = samples[library_id])

so@meta.data$cb <- str_split(colnames(so), "_", simplify = TRUE) %>% .[, 2]
so@meta.data$Barcode <- str_c(so$cb,
                                  "-",
                                  match(so$orig.ident,
                                        aggr_csv$sample_id))

mdata_to_export <- cbind(so@meta.data, so@reductions$tsne@cell.embeddings)
# export projections
# Format is csv:
# Barcode,X,Y
# AAACATACGGTACT-1,0.5,0.5
# AAACATTGCTCGCT-1,0,0
# AAACATTGGCGATT-1,-0.5,-0.5

tsne_coords <- select(mdata_to_export, Barcode, X = tSNE_1, Y = tSNE_2)
write_csv(tsne_coords, file.path(data_dir, "loupe_files", "tsne_coordinates.csv"))

# export metadata
# Only can add catagorical variables and there is an upper limit on catagories
# Format is csv:
# Barcode,Graph-based,Other
# AAACATACGGTACT-1,Cluster 1,a
# AAACATTGCTCGCT-1,Cluster 14,b
# AAACATTGGCGATT-1,Cluster 12,c

mdata_out <- select(mdata_to_export,
                      Barcode,
                      Cell_types,
                      Experiment,
                      Condition,
                      Experiment_Id)

write_csv(mdata_out,
          file.path(data_dir, "loupe_files", "custom_mdata.csv"))

