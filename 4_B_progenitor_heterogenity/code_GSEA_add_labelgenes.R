# 
# 
devtools::install_github("junjunlab/GseaVis")
library(stringr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(GseaVis)

mydata = readRDS("./B_progenitor_subcluster.rds")
DEG = FindMarkers(mydata, group.by="cell_type", ident.1="B progenitor cells 5", ident.2=".", logfc.threshold=0, min.pct=0)
DEG = arrange(DEG, -log2FoldChange)

symbol <- rownames(DEG)
entrez <- bitr(symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
genelist <- DEG$avg_log2FC
names(genelist) <- rownames(DEG)
#
genelist <- genelist[names(genelist) %in% entrez[,1]]
names(genelist) <- entrez[match(names(genelist), entrez[,1]), 2]
genelist = sort(genelist, decreasing=T)

# 
R.utils::setOption("clusterProfiler.download.method", "auto") 

KEGG_ges <- gseKEGG(
geneList=genelist,
organism="hsa",
minGSSize=10,
maxGSSize=500,
pvalueCutoff=0.05,
pAdjustMethod="BH",
verbose=FALSE,
eps=0
)

#################################################################################
#
KEGG_ges2 <- setReadable(KEGG_ges, OrgDb=org.Hs.eg.db, keyType="SYMBOL")
core <- result$core_enrichment[52]
core_genes <- str_split(core ,'/')[[1]]
p = gseaNb(
    object = KEGG_ges2,
    geneSetID = KEGG_ges2@result$ID[52],
    addPval = T,
    pvalX = 0.95,
    pvalY = 0.8,
    newGsea = T,
    addPoint = F,
    newCurveCol = c("red","lightgray", "blue"),
    newHtCol = c("blue", "lightgray", "red"),
    addGene = core_genes, #
    kegg = F,
    geneCol = 'black',
    max.overlaps = 30
)
ggsave("GSEA_KEGG_cellcycle.pdf", p, width=6, height=6)


