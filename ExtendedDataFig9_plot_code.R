library(clusterProfiler)
library(enrichplot)
library(ggplot2)

deg_V_C_df <- read.csv('~/D15V_vs_D15C-deg_Upregulated_enrichgo.csv')
deg_V_C <- readRDS('~/D15V_vs_D15C-deg_Upregulated_enrichgo.rds')
library(RColorBrewer)
cols <- c(brewer.pal(11, "RdBu"))
dotplot(deg_V_C,
        x="GeneRatio",color="p.adjust",size="Count",
        showCategory=c('immune response-regulating signaling pathway',
                       'cell activation involved in immune response',
                       'positive regulation of cell adhesion',
                       'regulation of angiogenesis',
                      ' wound healing',
                       'regulation of vasculature development',
                       'cell-substrate adhesion',
                       'endothelial cell migration',
                       'positive regulation of angiogenesis',
                       'positive regulation of vasculature development',
                       'extracellular matrix organization',
                       'extracellular structure organization',
                       'skeletal system morphogenesis',
                       'cytokine production involved in immune response',
                       'cell-matrix adhesion',
                       'vascular process in circulatory system',
                       'cell adhesion mediated by integrin',
                       'endothelium development',
                       'regulation of wound healing',
                       'blood vessel endothelial cell migration'
        ),
        split=NULL,font.size=12,title="D15V_vs_D15C_up")+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(), axis.line = element_line())+
  scale_fill_gradientn(colours = cols)+ scale_size (range=c (3, 7))
  
deg_V_C_df <- read.csv('~/D15V_vs_D15C-deg_Downregulated_enrichgo.csv')
deg_V_C <- readRDS('~/D15V_vs_D15C-deg_Downregulated_enrichgo.rds')
library(RColorBrewer)
cols <- c(brewer.pal(11, "RdBu"))
dotplot(deg_V_C,
        x="GeneRatio",color="p.adjust",size="Count",
        showCategory=c('muscle tissue development',
                       'striated muscle tissue development',
                       'cardiac muscle tissue development',
                       'muscle cell development',
                       'myofibril assembly',
                       'striated muscle cell development',
                       'striated muscle cell differentiation',
                       'muscle cell differentiation',
                       'heart contraction',
                       'heart process',
                       'muscle contraction',
                       'regulation of heart contraction',
                       'cardiac cell development',
                       'muscle system process',
                       'regulation of blood circulation',
                       'cardiac muscle cell development',
                       'sarcomere organization',
                       'striated muscle contraction',
                       'cellular component assembly involved in morphogenesis',
                       'muscle organ development'
        ),
        split=NULL,font.size=12,title="D15V_vs_D15C_down")+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(), axis.line = element_line())+
  scale_fill_gradientn(colours = cols)+ scale_size (range=c (3, 7))  
  
  


