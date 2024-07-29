# https://cloud.tencent.com/developer/article/1692240
# https://zhuanlan.zhihu.com/p/358986392
library(Seurat)
library(doParallel)
library(SCopeLoomR)

mydata = readRDS("./B_progenitor_cells_5.rds")
exprMat <- as.matrix(mydata@assays$RNA@counts)
cellInfo = mydata@meta.data[, c(12, 3, 2)]
colnames(cellInfo)=c('CellType', 'nGene' ,'nUMI')

### Initialize settings
library(SCENIC)
# 
mydbDIR <- "./database"
mydbs <- c("hg19-500bp-upstream-7species.mc9nr.feather", "hg19-tss-centered-10kb-7species.mc9nr.feather")
names(mydbs) <- c("500bp", "10kb")
scenicOptions <- initializeScenic(org="hgnc", dbDir=mydbDIR, dbs=mydbs, nCores=4, datasetTitle="Epithelial cells") 
saveRDS(scenicOptions, "int/scenicOptions.rds")

### Co-expression network
genesKept <- geneFiltering(exprMat, scenicOptions, minCountsPerGene=3*0.01*ncol(exprMat), minSamples=ncol(exprMat)*0.01)
exprMat_filtered <- exprMat[genesKept, ]
## 
runCorrelation(exprMat_filtered, scenicOptions)
## 
exprMat_filtered_log <- log2(exprMat_filtered+1)
runGenie3(exprMat_filtered_log, scenicOptions, nParts=20)
## 

### Build and score the GRN
exprMat_log <- log2(exprMat+1)
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) # Toy run settings
runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType="AUC")        # choose settings

####################################################################################################################
# 
# 
##
AUCmatrix <- readRDS("int/3.4_regulonAUC.Rds")
AUCmatrix <- AUCmatrix@assays@data@listData$AUC
AUCmatrix <- data.frame(t(AUCmatrix), check.names=F)
RegulonName_AUC <- colnames(AUCmatrix)
RegulonName_AUC <- gsub(' \\(','_',RegulonName_AUC)
RegulonName_AUC <- gsub('\\)','',RegulonName_AUC)
colnames(AUCmatrix) <- RegulonName_AUC
mydata <- AddMetaData(mydata, AUCmatrix) #
##################################################################
## 

## 
##
##

##

df = read.table("./TF_AUC_heatmap.txt", header=T, row.names=1, sep="\t")
pheatmap(df, scale="row", angle_col="45", color=colorRampPalette(c("#0da9ce", "white", "#e74a32"))(500), clustering_method="median", border_color="white")




