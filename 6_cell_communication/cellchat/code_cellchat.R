library(CellChat)
#################### 
cellchat = createCellChat(object=scRNA, group.by="cell_type", meta=scRNA@meta.data)
groupSize = as.numeric(table(cellchat@idents))
CellChatDB = CellChatDB.human
#####  
[1] "Secreted Signaling"	"ECM-Receptor"	"Cell-Cell Contact"  #
CellChatDB.use <- subsetDB(CellChatDB, search="Secreted Signaling")
cellchat@DB <- CellChatDB.use
#####  
cellchat = subsetData(cellchat)
future::plan("multisession", workers=4)
## 
cellchat = identifyOverExpressedInteractions(cellchat)  #
## 
cellchat = projectData(cellchat, PPI.human)
## 
#####################################################################################
### 
# 
cellchat = computeCommunProb(cellchat, raw.use=FALSE, population.size=TRUE)
# 
cellchat = filterCommunication(cellchat, min.cells=10)
#####  
cellchat = computeCommunProbPathway(cellchat)
cellchat = aggregateNet(cellchat)
#####################################################################################
#######  

netVisual_bubble(cellchat, sources.use=c(1), targets.use=c(2,3), signaling=c("TGFb"), color.heatmap="Spectral")+
theme(text=element_text(family="Times"), axis.text.x=element_text(angle=15, hjust=0.5, face="bold", size=10), axis.text.y=element_text(face="bold", size=10))

# 
