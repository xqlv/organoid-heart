#### Load packages ----
library(Seurat)
library(DropletUtils)
library(tidyverse)
library(ggsci)

##### Load data CCO_day8 ----
CCO_day8_counts <- Read10X(data.dir = "~/CCO_day8/1.Basic_analysis/1.2.filtered_feature_bc_matrix/")
dim(CCO_day8_counts)
CCO_day8 <-
  CreateSeuratObject(
    CCO_day8_counts,
    project = "CCO_day8",
    min.cells = 3,
    min.features = 200
  )
CCO_day8[["group"]] <- "CCO_day8"

#### QC ---
is.mt <- grep("^MT-",rownames(CCO_day8))
rownames(CCO_day8)[is.mt]
CCO_day8 <- PercentageFeatureSet(CCO_day8,
                              pattern = "^MT-",
                              col.name = "percent.mt")

is.rb <- grep("^RP[SL]",rownames(CCO_day8))
rownames(CCO_day8)[is.rb]
CCO_day8 <- PercentageFeatureSet(CCO_day8,
                              pattern = "^RP[SL]",
                              col.name = "percent.rb")

CCO_day8 <- CellCycleScoring(object = CCO_day8,
                          g2m.features = cc.genes$g2m.genes,
                          s.features = cc.genes$s.genes)


VlnPlot(
  CCO_day8,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt",
               "G2M.Score","S.Score", "percent.rb"),
  pt.size = 0.01
)

CCO_day8 <- subset(CCO_day8, subset =
                  nFeature_RNA >= 200 &
                  nFeature_RNA <= 5000 &
                  nCount_RNA <= 10000 &
                  percent.mt <= 5)

saveRDS(CCO_day8, '~/CCO_day8_QC.rds')