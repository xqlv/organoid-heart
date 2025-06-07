library(Seurat)

data_hGas_24h <- readRDS("~/GSM6341940_hGas_24h.RDS")

DimPlot(data_hGas_24h,
        reduction = "umap",
        label = T)


data_hGas_48h <- readRDS("~/GSM6341941_hGas_48h.RDS")

DimPlot(data_hGas_48h,
        reduction = "umap",
        label = T)


data_hGas_72h <- readRDS("~/GSM6341942_hGas_72h.RDS")

DimPlot(data_hGas_72h,
        reduction = "umap",
        label = T)
		
library(cowplot)
library(ggpubr)
library(RColorBrewer)
library(pheatmap)

set.seed(717)

features <- c("ADGRA2", "BAMBI", "CDH3", "COL1A1", "CTNNB1", "DACT3", "EDNRB",
              "FGFR2", "GLI3", "LGR4", "LGR5", "LYPD6", "RBMS3", 
              "SMAD3", "SOX4", "TBL1X", "TMEM88", "WWTR1")

# Function to process data and plot heatmap
process_and_plot <- function(path, time_point) {
  data <- readRDS(path)
  data@meta.data[["celltype"]] <- data@active.ident
  seurat.obj <- data
  DefaultAssay(seurat.obj)
  seurat.obj <- NormalizeData(seurat.obj)
  data <- AverageExpression(seurat.obj, assays = 'RNA', slot = "data", group.by = "celltype") %>% as.data.frame()
  colnames(data) <- gsub('RNA.', '', colnames(data))
  data2 <- data[features, ]
  
  pheatmap(data2, 
           scale = 'row',
           cluster_rows = T, cluster_cols = F,
           cellwidth = 30, cellheight = 10,
           treeheight_col = 2, treeheight_row = 6,
           clustering_method = "ward.D2",
           color = colorRampPalette(rev(c(brewer.pal(11, "RdBu"))))(100),
           angle_col = c("90"),
           border_color = "white",
           fontsize = 10,
           main = paste("intersection_genes_in", time_point, "hGas"))
}

# Process and plot for each time point
process_and_plot("~/GSM6341940_hGas_24h.RDS", "24h")
process_and_plot("~/GSM6341941_hGas_48h.RDS", "48h")
process_and_plot("~/GSM6341942_hGas_72h.RDS", "72h")
		