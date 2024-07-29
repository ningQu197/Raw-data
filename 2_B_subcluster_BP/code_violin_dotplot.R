library(ggplot2)
library(Seurat)
library(SeuratObject)


mydata = readRDS("../1_preprocessing/mydata_cluster.rds")
mydata = subset(mydata, cell_type %in% c("B cells 1", "B cells 2", "B cells 3", "B cells 4", "B cells 5"))
colors = c("#4E79A7", "#A0CBE8", "#B6992D", "#F1CE63", "#DD585B")

### 
### 
gene_set1 = c("IGHM", "JCHAIN", "IGKC", "IGHD")
p = VlnPlot(mydata, features=gene_set1, pt.size=0, cols=colors, ncol=2)+NoLegend()+theme(axis.title.x=element_blank())
ggsave("violin_Plasma_B.pdf", p, width=6, height=6)

gene_set2 = c("BIRC5", "KIF22", "CDT1", "H2AFX", "MAD2L1", "ZWINT", "PTTG1", "KPNA2", "UBE2S")
p = DotPlot(mydata, features=gene_set2)+coord_flip()+scale_color_distiller(palette="RdYlBu")+
theme(text=element_text(family="Times"), axis.text.x=element_text(angle=30, hjust=1, size=12, face="bold"), axis.text.y=element_text(face="bold", size=12), axis.title.x=element_blank(), axis.title.y=element_blank())
ggsave("dotplot_B_cells_1.pdf", p, width=6, height=6)

### 
### 












