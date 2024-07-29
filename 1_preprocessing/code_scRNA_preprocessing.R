library(Seurat)
library(tidyverse)
library(sctransform)
library(harmony)
###############################################################################
assays <- dir("./data/")
dir <- paste0("./data/", assays)
# 
samples_name = assays

# 
scRNAlist <- list()
for(i in 1:length(dir)){
  counts <- Read10X(data.dir=dir[i])
  scRNAlist[[i]] <- CreateSeuratObject(counts, project=samples_name[i], min.cells=3, min.features=200)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]], add.cell.id=samples_name[i])   
  if(T){    
    scRNAlist[[i]][["percent.mito"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern="^MT-") 
  }
}
### 
names(scRNAlist) <- samples_name
# 
scRNA <- merge(scRNAlist[[1]], scRNAlist[2:length(scRNAlist)])
scRNA = subset(scRNA, subset=nFeature_RNA>200&nFeature_RNA<5000&percent.mito<10)  #
p = VlnPlot(scRNA, features=c("nFeature_RNA"), pt.size=0, cols=colors)
ggsave("nFeature_RNA_violin.pdf", p, width=6, height=6)
colors = c4a("tableau.20", 12)
#########################################################################################################
s.genes=Seurat::cc.genes.updated.2019$s.genes
g2m.genes=Seurat::cc.genes.updated.2019$g2m.genes
scRNA <- CellCycleScoring(scRNA, s.features=s.genes, g2m.features=g2m.genes, set.ident=TRUE)
scRNA = SCTransform(scRNA, method="glmGamPoi", vars.to.regress=c("nCount_RNA", "percent.mito", "S.Score", "G2M.Score"), verbose=FALSE)
scRNA = RunPCA(scRNA, verbose=FALSE)
scRNA = RunHarmony(scRNA, group.by.vars="Sample", max.iter.harmony=50, lambda=0.5, assay.use="SCT")
ElbowPlot(scRNA)
scRNA <- FindNeighbors(scRNA, dims=1:20, reduction="harmony")
scRNA <- RunUMAP(scRNA, dims=1:20, reduction="harmony")
p = DimPlot(scRNA, reduction="umap", group.by="Sample", pt.size=1)+
theme(legend.position="right", plot.title=element_blank())+scale_color_manual(values=colors)
ggsave("UMAP_Sample.pdf", p, height=6, width=6)
######################################################### 
mydata <- FindClusters(scRNA, resolution=0.2)
UMAPPlot(mydata, pt.size=1, label=T, cols=colors, label.size=5)+NoLegend()
markers <- FindAllMarkers(mydata, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
write.table(markers, "./All_celltype_DEG.txt", col.names=T, row.names=T, quote=F, sep="\t")
############## 
VlnPlot(mydata, features=c("S100A9"), pt.size=0, cols=colors)+NoLegend()+theme(axis.title.x=element_blank())
# 0,1,6,7,10: B cells:              CD19, BANK1, CD24
# 2:          Naive T cells:        CCR7, CD3D
# 3, 4, 8:    Erythrocytes:         HBB, HBA1, HBA2
# 5:          Monocytes:            CST3, CD14
# 9:          Cytotoxic NK/T cells: NKG7, CCL5, GZMA

# 
cell_label = c(
"B cells 1", "B cells 2", "Naive T cells", "Erythrocytes 1", "Erythrocytes 2", "Monocytes", "B cells 3",
"B cells 4", "Erythrocytes 3", "Cytotoxic NK/T cells", "B cells 5"
)
names(cell_label) <- levels(mydata)
mydata <- RenameIdents(mydata, cell_label)
mydata[["cell_type"]] = Idents(mydata)
p = UMAPPlot(mydata, pt.size=1, label=T, label.size=5)+NoLegend()+scale_color_manual(values=colors)
ggsave("UMAP_subcluster.pdf", p, width=6, height=6)

genes = c("CD19", "BANK1", "CD24", "CD3D", "CCR7", "HBB", "HBA1", "CST3", "CD14", "NKG7", "CCL5")
p = DotPlot(mydata, features=genes, cols=c("snow", "chartreuse4"))+coord_flip()+
theme(text=element_text(family="Times"), axis.text.x=element_text(angle=30, hjust=1, size=12, face="bold"), axis.text.y=element_text(face="bold", size=12), axis.title.x=element_blank(), axis.title.y=element_blank())

# 
set = c("CD24", "CD19", "CCR7", "CD3D", "HBB", "CST3", "CCL5", "NKG7")
FeaturePlot(mydata, features=set, cols=c("snow", "red"), ncol=4)

genes = c("AKR1C3", "WNT7A", "FAM72B", "RERG", "IDO1", "HEY1")
VlnPlot(mydata, features=genes, pt.size=0, cols=colors, ncol=2)+NoLegend()+theme(axis.title.x=element_blank())

#####################################################################  
bar = mydata@meta.data %>% group_by(Sample, cell_type) %>% count()
ggplot(data=bar, aes(x=Sample, y=n, fill=cell_type))+ 
geom_bar(stat="identity", position=position_fill())+
scale_fill_manual(values=colors)+theme_classic()+
theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(face="bold", size=12, angle=15, hjust=1), axis.text.y=element_text(face="bold", size=12), legend.text=element_text(face="bold", size=12), legend.title=element_blank(), legend.position="right")

#####################################################################  
Type_label = c("Control", "24h Treated", "72h Treated")
bar$Type = factor(bar$Type, levels=Type_label)
bar = bar %>% group_by(Type) %>% mutate(percent=100*n/sum(n))

ggplot(data=bar, aes(x=cell_type, y=percent, fill=Type))+
geom_bar(stat="identity", position=position_dodge())+
scale_fill_manual(values=colors)+theme_classic()+
ggtitle("Percent(%)")+
theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(angle=30, hjust=1, size=12, face="bold"), axis.text.y=element_text(face="bold", size=12), legend.text=element_text(face="bold", size=12), legend.title=element_blank())
ggsave("barplot_pair_number.pdf", p, width=6, height=6)


########################## 
library(ggplot2)
library(ggpubr)

df = read.table("./celltype_number_percent.txt", header=T, sep="\t")
ggplot(df, aes(x=reorder(cell_type, -percent, sum), y=percent, fill=Type))+
scale_fill_manual(values=c("#9EB9F3", "#FC8D62"))+
geom_boxplot(outlier.size=0.1, width=0.3)+
theme_bw()+
stat_compare_means(aes(group=Type), label="p.signif", method="t.test")+
theme(axis.text.x=element_text(angle=15, hjust=1, face="bold", size=10), axis.text.y=element_text(face="bold", size=10), axis.title.x=element_blank())





