# https://zhuanlan.zhihu.com/p/601431156
# http://www.360doc.com/content/22/1206/08/65403234_1059126221.shtml
# https://blog.csdn.net/qq_42090739/article/details/126264415
# 
# 
############################################ 
$ conda activate cellphoneDB

# 
$ pip install CellphoneDB==3.0.0     # 
$ conda install -c conda-forge rpy2
# 
$ cellphonedb database download
# 
# 
# 
cellphonedb method statistical_analysis meta.txt counts.txt --counts-data=gene_name
# 

# 
# deconvoluted.txt:      
# mean.txt:              
# pvalues.txt            
# signigicant_means.txt: 

# 
cellphonedb plot dot_plot
cellphonedb plot heatmap_plot meta.txt
# 
# count_network.txt       
# means.txt               
# pvalues.txt             

# 
pip install markupsafe==2.0.1
# 
conda install -c conda-forge r-base=3.6.0
conda install -c conda-forge  r-ggplot2
conda install -c conda-forge r-pheatmap
conda install -c conda-forge r-Seurat

####################################################################################################
# 
library(CellChat)
library(tidyr)
count_inter = read.delim("count_network.txt", check.names=FALSE)
count_inter = spread(count_inter, TARGET, count)
rownames(count_inter) = count_inter$SOURCE
count_inter = count_inter[, -1]
count_inter = as.matrix(count_inter)
netVisual_circle(count_inter, weight.scale=T, color.use=c("red", "green", "blue"), label.edge=T, arrow.size=0.8, vertex.weight=1)

# 
par(mfrow=c(2, 2), xpd=TRUE)
colors = c("orange", "deepskyblue", "orangered")
for (i in 1:nrow(count_inter)) {
	mat2 = matrix(0, nrow=nrow(count_inter), ncol=ncol(count_inter), dimnames=dimnames(count_inter))
	mat2[i, ] = count_inter[i, ]
	netVisual_circle(mat2,
	weight.scale=T,
	label.edge=T,
	color.use=colors,
	edge.weight.max=max(count_inter),
	title.name=rownames(count_inter)[i],
	arrow.size=0.6)
}
# 
####################################################################################################
# 
# 
devtools::install_github("zktuong/ktplots", dependencies=TRUE)
library(ktplots)
pvals = read.delim("pvalues.txt", check.names=FALSE)
means = read.delim("means.txt", check.names=FALSE)
# 
# 
# 
# keep_significant_only=T
scRNA = readRDS("./scRNA.rds")
p = plot_cpdb(cell_type1="B progenitor cells 5", cell_type2="", scdata=scRNA, idents="cell_type",
means=means, pvals=pvals,, highlight_size=0, genes=c("KLRC", "CD47", "LILRB2", "FAM3C"),
keep_significant_only=T,
col_option=colorRampPalette(c("snow", "dodgerblue1", "peachpuff1","hotpink"))(50)
)+
theme(axis.text=element_text(size=12, color="black", face="bold"), axis.text.x=element_text(angle=15, hjust=0))
ggsave("cellphoneDB_dotplot.pdf", p, width=6, height=6)
# c("snow", "springgreen", "peachpuff1", "tan1")
# c("snow", "chartreuse", "gold1", "brown1")



