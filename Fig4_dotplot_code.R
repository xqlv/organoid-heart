library(clusterProfiler)
library(enrichplot)
library(ggplot2)

deg_CV_C <- readRDS('~/D15CV_vs_D15C-deg_Upregulated_enrichgo.rds')
library(RColorBrewer)
cols <- c(brewer.pal(11, "RdBu"))
dotplot(deg_CV_C,
        x="GeneRatio",color="p.adjust",size="Count",
        showCategory=c('positive regulation of cytosolic calcium ion concentration',
                       'skeletal system morphogenesis',
                       'embryonic organ development',
                       'wound healing',
                       'toll-like receptor 4 signaling pathway',
                       'anterior/posterior pattern specification',
                       'embryonic skeletal system development',
                       'regulation of membrane invagination',
                       'tissue remodeling',
                       'calcium-mediated signaling',
                       'potassium ion transport',
                       'embryonic skeletal system morphogenesis',
                       'pattern specification process',
                       'cell-substrate adhesion',
                       'regulation of actin cytoskeleton organization',
                       'response to metal ion',
                       'regulation of actin filament-based process',
                       'positive regulation of angiogenesis',
                       'positive regulation of vasculature development',
                       'calcium ion transport'
        ),
        split=NULL,font.size=12,title="D15CV_vs_D15C_up")+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(), axis.line = element_line())+
  scale_fill_gradientn(colours = cols)+ scale_size (range=c (3, 7))

