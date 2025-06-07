library(cowplot)
library(clusterProfiler)
library(ggpubr)
library(Seurat)

markers <- readRDS("~/markers_car.rds")
car_MESO.obj <- readRDS("~/car_MESO.rds")
endo_ENDO.obj <- readRDS("~/endo_ENDO.rds")
ect_ECT.obj <- readRDS("~/ect_ECT.rds")

merged_obj <- merge(car_MESO.obj, y = c(endo_ENDO.obj, ect_ECT.obj))
merged_obj$cell_type3 <- factor(merged_obj$cell_type3, levels = c('car', 'endo', 'ect'))
seurat.obj <- merged_obj

C5_gene_sets <- msigdbr::msigdbr(species = "human",
                                 category = "C5") %>%
  dplyr::select(gs_name, gene_symbol)

selected_gene_sets <- C5_gene_sets %>%
  filter(gs_name %in% c('GOBP_CANONICAL_WNT_SIGNALING_PATHWAY')) %>% pull(gene_symbol)

features <- intersect(x=selected_gene_sets, y = markers_car[markers_car$avg_log2FC > 0 & markers_car$p_val < 0.05,]$gene)

DefaultAssay(seurat.obj)
seurat.obj <- NormalizeData(seurat.obj)
data <- AverageExpression(seurat.obj,
                          assays = 'RNA',slot = "data",group.by = "cell_type3") %>% as.data.frame()
colnames(data) <- gsub('RNA.', '', colnames(data))

data2 <- data[features,]
pheatmap::pheatmap(data2, 
                   scale = 'row',
                   cluster_rows=T, cluster_cols=F,
                   cellwidth = 30, cellheight = 10,
                   treeheight_col = 2, treeheight_row = 6,
                   clustering_method="ward.D2",
                   color = colorRampPalette(rev(c(brewer.pal(11, "RdBu"))))(100),
                   angle_col = c("0"),
                   border_color = "white",
                   fontsize=10)

features_24h <- c("ADGRA2","AMER2","AMER3","ANKRD6","APC","AXIN2","BAMBI","BCL9","BCL9L","BMP2","CAPRIN2","CDH2","CDH3","CHD8","COL1A1","CSNK1A1","CSNK1E","CTNNB1",
                  "CTNNBIP1","CTNND1","DAAM2","DACT1","DACT3","DKK1","DKK3","DVL2","DVL3","EDNRB","EGR1","EXT1","FGFR2","FZD4","FZD7","FZD8","GLI3","ILK",
                  "ISL1","JUP","KREMEN1","LATS2","LGR4","LGR5","LMBR1L","LRP6","LYPD6","MESP1","NKD1","NOTUM","NRARP","PLPP3","RBMS3","ROR2","RSPO1","SCYL2",
                  "SEMA5A","SFRP1","SFRP5","SMAD3","SNAI2","SOX4","SRC","TBL1X","TBL1XR1","TMEM88","TNKS","TPBG","UBR5","USP34","USP47","USP8","VCP","WNT11",
                  "WNT3","WNT4","WNT6","WNT9B","WWTR1","XIAP","YAP1","ZNF703")
features_48h <- c("ADGRA2","AMER3","ANKRD6","APOE","AXIN2","BAMBI","BCL9","BCL9L","BMP2","CCAR2","CDH3","COL1A1","CSNK1E","CSNK1G3","CTHRC1","CTNNB1","CYLD","DAB2",
                  "DAB2IP","DACT3","EDNRA","EDNRB","FERMT1","FGF10","FGFR2","FOXO1","FOXO3","FZD4","FZD6","GATA3","GLI3","ILK","ISL1","JADE1","JUP","KANK1",
                  "LGR4","LGR5","LMBR1L","LYPD6","LZTS2","MAPK14","MCC","MDK","MED12","NKX2-5","PFDN5","RBMS3","SCYL2","SFRP5","SLC9A3R1","SMAD3","SOX13","SOX4",
                  "STK3","TBL1X","TBL1XR1","TCF7L2","TMEM88","VPS35","WLS","WNT11","WNT9B","WWTR1","XIAP","YAP1","ZBED3", "ZEB2")
features_72h <- c("ADGRA2","AMER3","APOE","BAMBI","BCL9L","BMP2","CDH3","COL1A1","CTHRC1","CTNNB1","DAB2","DACT3","DDIT3","DVL1","DVL3","EDNRA","EDNRB","FERMT1",
                  "FGF10","FGFR2","FOXO1","FOXO3","FZD4","FZD6","GATA3","GLI3","GPC5","ISL1","JADE1","JUP","KANK1","LATS2","LGR4","LGR5","LMBR1L","LYPD6",
                  "LZTS2","MAPK14","MCC","MDK","MED12","NKX2-5","PFDN5","PIN1","PTEN","PTK7","RBMS3","RYK","SCYL2","SDHAF2","SEMA5A","SFRP1","SFRP5","SLC9A3R1",
                  "SMAD3","SOX4","STK4","TBL1X","TBL1XR1","TCF7L1","TCF7L2","TGFB1","TLE6","TMEM198","TMEM88","TNKS","UBAC2","UBR5","USP34","USP47","VCP","WNT11",
                  "WWTR1","YAP1","ZNF703")


features_intersect <- Reduce(intersect, list(features, features_24h, features_48h, features_72h))
data3 <- data[features_intersect,]
pheatmap::pheatmap(data3, 
                   scale = 'row',
                   cluster_rows=T, cluster_cols=F,
                   cellwidth = 30, cellheight = 10,
                   treeheight_col = 2, treeheight_row = 6,
                   clustering_method="ward.D2",
                   color = colorRampPalette(rev(c(brewer.pal(11, "RdBu"))))(100),
                   angle_col = c("0"),
                   border_color = "white",
                   fontsize=10)


library(VennDiagram)
venn.plot <- venn.diagram(
  x = list(features = features, features_24h = features_24h, features_48h = features_48h, features_72h = features_72h),
  category.names = c("features", "features_24h",'features_48h','features_72h'),
  filename = NULL,
  output = TRUE
)
grid.draw(venn.plot)


