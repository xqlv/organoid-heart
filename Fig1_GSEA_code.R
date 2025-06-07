library(Seurat)
library(ggsci)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)

source("~/custom_function.R")
set.seed(717)

markers <- readRDS("~/markers_car.rds")
gsea.input <- markers

#H
gsea_res <- cat_gsea(gsea.input,
                     arrange_by = "avg_log2FC",
                     gene_name = "gene",
                     category = "H",
                     species = "human")

saveRDS(gsea_res, "~/gsea_H_car.rds")
write.csv(gsea_res, "~/gsea_H_car_table.csv",row.names = F)

#C2
gsea_res <- cat_gsea(gsea.input,
                     arrange_by = "avg_log2FC",
                     gene_name = "gene",
                     category = "C2",
                     species = "human")
					 
saveRDS(gsea_res, "~/gsea_C2_car.rds")
write.csv(gsea_res, "~/gsea_C2_car_table.csv",row.names = F)

#C5
gsea_res <- cat_gsea(gsea.input,
                     arrange_by = "avg_log2FC",
                     gene_name = "gene",
                     category = "C5",
                     species = "human")

saveRDS(gsea_res, "~/gsea_C5_car.rds")
write.csv(gsea_res, "~/gsea_C5_car_table.csv",row.names = F)

source("~/custom_plot_function.R")
gsea_C5_car |>
  cat_gseaplot(
    "GOBP_CANONICAL_WNT_SIGNALING_PATHWAY",
    title = "GOBP_CANONICAL_WNT_SIGNALING_PATHWAY",
    subplots = c(1, 2),
    pvalue_table = T
  )

gsea_C5_endo |>
  cat_gseaplot(
    "GOBP_CANONICAL_WNT_SIGNALING_PATHWAY",
    title = "GOBP_CANONICAL_WNT_SIGNALING_PATHWAY",
    subplots = c(1, 2),
    pvalue_table = T
  )

gsea_C5_ect |>
  cat_gseaplot(
    "GOBP_CANONICAL_WNT_SIGNALING_PATHWAY",
    title = "GOBP_CANONICAL_WNT_SIGNALING_PATHWAY",
    subplots = c(1, 2),
    pvalue_table = T
  )






