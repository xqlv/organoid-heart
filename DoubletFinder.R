library(DoubletFinder)
library(tidyverse)
library(Seurat)
library(ggsci)

#### CCO_day8 ----
CCO_day8 <- readRDS('~/CCO_day8_QC.rds')
CCO_day8 <- NormalizeData(CCO_day8, normalization.method = "LogNormalize", scale.factor = 10000)
CCO_day8 <- FindVariableFeatures(CCO_day8, selection.method = "vst", nfeatures = 2000)
CCO_day8 <- ScaleData(CCO_day8, 
                   # features = all.genes, 
                   vars.to.regress = c("percent.mt","nFeature_RNA", "nCount_RNA",
                                       "percent.rb","G2M.Score","S.Score"))
CCO_day8 <- RunPCA(CCO_day8, features = VariableFeatures(object = CCO_day8))
ElbowPlot(CCO_day8, ndims = 50)
CCO_day8 <- RunUMAP(CCO_day8, dims = 1:30)
CCO_day8 <- RunTSNE(CCO_day8, dims = 1:30)
DimPlot(CCO_day8,reduction = "umap")+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+
  xlab('UMAP1')+ylab('UMAP2')+scale_color_npg()

CCO_day8 <- FindNeighbors(CCO_day8,
                       dims = 1:30)
CCO_day8 <- FindClusters(CCO_day8,
                      resolution = 1) 
sweep.res.list <- CCO_day8 %>%
  paramSweep_v3(PCs = 1:30,
                sct = FALSE,
                num.cores = 10)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <-
  bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()

doublet_rate <- 0.076 
homotypic.prop <-
  modelHomotypic(CCO_day8$seurat_clusters)
nExp_poi <- round(doublet_rate * ncol(CCO_day8))
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

CCO_day8 <-
  doubletFinder_v3(
    CCO_day8,
    PCs = 1:30,
    pN = 0.25,
    pK = pK_bcmvn,
    nExp = nExp_poi.adj,
    reuse.pANN = F,
    sct = F
  )

table(CCO_day8$DF.classifications_0.25_0.26_1144)
DimPlot(CCO_day8, reduction = "umap", 
        group.by = "DF.classifications_0.25_0.26_1144")

sub_CCO_day8 <- CCO_day8 %>% 
  subset(DF.classifications_0.25_0.26_1144 == "Singlet")
DimPlot(sub_CCO_day8, reduction = "umap", 
        group.by = "DF.classifications_0.25_0.26_1144")
saveRDS(sub_CCO_day8, "~/CCO_day8_Singlet.rds")