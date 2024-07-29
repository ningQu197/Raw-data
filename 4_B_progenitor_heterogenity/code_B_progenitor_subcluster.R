library(Seurat)
library(harmony)
library(tidyverse)
library(sctransform)
library(glmGamPoi)

mydata = readRDS("../1_preprocessing/mydata_cluster.rds")
mydata = SCTransform(mydata, method="glmGamPoi", vars.to.regress=c("nCount_RNA", "percent.mito", "S.Score", "G2M.Score"), verbose=FALSE)
mydata = RunPCA(mydata, verbose=FALSE, assay="SCT")
mydata = RunHarmony(mydata, group.by.vars="Sample", max.iter.harmony=50, lambda=0.5, assay.use="SCT")
ElbowPlot(mydata)
mydata <- FindNeighbors(mydata, dims=1:20, reduction="harmony")
mydata <- RunUMAP(mydata, dims=1:20, reduction="harmony")
colors = c4a("tableau.classic20", 12)
colors = c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C", "#98DF8A", "#D62728")
######################################  
DimPlot(mydata, reduction="umap", group.by="Sample", pt.size=1)+
theme(legend.position="right", plot.title=element_blank())+
scale_color_manual(values=colors)
#################################################################################
#############################  
mydata <- FindClusters(mydata, resolution=0.4)
UMAPPlot(mydata, pt.size=2, label=T, cols=colors, label.size=6)+NoLegend()
#################################################################################
# 
cell_label = c("B progenitor cells 1", "B progenitor cells 2", "B progenitor cells 3", "B progenitor cells 4",
"B progenitor cells 5", "B progenitor cells 6", "B progenitor cells 7")
#################################################################################
## 
names(cell_label) <- levels(mydata)
mydata <- RenameIdents(mydata, cell_label)
mydata[["cell_type"]] = Idents(mydata)
p = UMAPPlot(mydata, pt.size=2, label=T, label.size=5)+NoLegend()+scale_color_manual(values=colors)+theme(text=element_text(family="Times"))
ggsave("UMAP_B_progenitor_subcluster.pdf", p, width=6, height=6)
saveRDS(mydata, "./B_progenitor_subcluster.rds")
markers <- FindAllMarkers(mydata, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
write.table(markers, "B_progenitor_subcluster_DEG.txt", col.names=T, row.names=T, quote=F, sep="\t")
###################################################################################
########################### 
library(AUCell)
set_df = read.table("./chr6p_amp_genes.txt", header=T, sep="\t")
geneset = split(set_df$Gene, set_df$Term)
CellRank <- AUCell_buildRankings(as.matrix(mydata@assays$RNA@data))
cells_AUC <- AUCell_calcAUC(geneset, CellRank, nCores=5, aucMaxRank=nrow(CellRank)*0.05)
select_set <- "chr6p_amp"
aucs <- getAUC(cells_AUC)[select_set, ]
########################### 
celltype_df = subset(mydata@meta.data, select=c("cell_type"))
df = merge(x=celltype_df, y=aucs, by=intersect(names(celltype_df), names(aucs)))
write.table(df, "./B_progenitor_sub_chr6p_amp_score.txt", col.names=T, row.names=T, sep="\t", quote=F)
p = ggplot(df, aes(x=cell_type, y=chr6p_amp, color=cell_type, fill=cell_type))+
scale_color_manual(values=colors)+
scale_fill_manual(values=colors)+
geom_violin(alpha=0.7)+
geom_boxplot(fill="white")+
geom_jitter(width=0.2, shape=21, size=0.5, aes(fill=cell_type))+
geom_signif(comparisons=list(cell_label), textsize=5, test=kruskal.test, step_increase=0.2, map_signif_level=F)+
ylab("Score of AUCell")+
theme_bw()+
theme(text=element_text(family="Times"), axis.text=element_text(size=15, face="bold"), axis.title.y=element_text(size=15, face="bold"), axis.title.x=element_blank(), legend.position="none")
ggsave("B_progenitor_sub_chr6p_amp_score.pdf", p, width=6, height=6)
###################################################################################
# 画B progenitor cells 5
genes = read.table("./B_progenitor_cells_5_proliferative_genes.txt", header=T, sep="\t")
genes = genes$Gene
p = DotPlot(mydata, features=genes)+coord_flip()+scale_color_distiller(palette="RdYlBu")+
theme(text=element_text(family="Times"), axis.text.x=element_text(angle=30, hjust=1, size=12, face="bold"), axis.text.y=element_text(face="bold", size=12), axis.title.x=element_blank(), axis.title.y=element_blank())
ggsave("dotplot_proliferation_genes.pdf", p, width=6, height=6)










