setwd("D:/R/001/workspace")
library(Seurat)
str<-readRDS("str.Rds")

str$new.cluster.ids<-str$seurat_clusters
library(tidyverse)
Idents(str)<-'seurat_clusters'
str$new.cluster.ids<-recode(str$new.cluster.ids,'0'="Mu-MCAM-NDUFA4L2",'1'="EC-ACKR1",'2'="Fb-FAP-MMP11",'3'="EC-KDR-ESM1",'4'="Fb-MMP11"
                            ,'5'="Fb-KLF4",'6'="Fb-FAP",'7'="EC-KDR-IGFBP3",'8'="Fb-CXCL12",'9'="EC-STMN1",
                            '10'="SMC",'11'="immune like",'12'="undefined-2",'13'="undefined-1",'14'="EC - Mu",'15'="EC-TFF3"
                            ,'16'="Mu-NDUFA4L2",'17'="Fb-MPZ",'18'="epithelial like")

Idents(str)<-'new.cluster.ids'


str1=subset(x = str ,idents=c("EC-ACKR1","EC-KDR-IGFBP3"))
str2=subset(x = str ,idents=c("EC-KDR-ESM1","EC-STMN1","EC-TFF3"))


pbmc.combined <- merge(str1, y = str2,merge.data = TRUE)

Idents(pbmc.combined)<-'new.cluster.ids'

table(pbmc.combined$new.cluster.ids)

rm("str1","str2","str")



#############cellchat安装
#library(devtools)
#BiocManager::install("ComplexHeatmap")
#devtools::install_github( "sqjin/CellChat" )

library(CellChat)

##提取表达矩阵和细胞类别创建cellchat对象
data.input <- GetAssayData(pbmc.combined, assay = "RNA", slot = "data")
identity <- subset(pbmc.combined@meta.data, select = "new.cluster.ids")
cellchat <- createCellChat(object = data.input, meta = identity,  group.by = "new.cluster.ids")

#导入配受体数据库（人或者鼠）
CellChatDB <- CellChatDB.human

#查看可以选择的侧面（选择特定的信息描述细胞间的作用）
unique(CellChatDB$interaction$annotation)
CellChatDB.use <- subsetDB(CellChatDB, search = "Cell-Cell Contact")
cellchat@DB <- CellChatDB.use
#对数据进行子集化，节省计算成本
cellchat <- subsetData(cellchat)
# 识别过表达基因
cellchat <- identifyOverExpressedGenes(cellchat)
# 识别配体-受体对
cellchat <- identifyOverExpressedInteractions(cellchat)
# 将配体、受体投射到PPI网络
cellchat <- projectData(cellchat, PPI.human)

#计算通信概率推断细胞互作的通信网络

cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
#过滤掉低质量的细胞间的通信
cellchat <- filterCommunication(cellchat, min.cells = 3)

#提取所有推断的配体/受体级别的细胞-细胞通信
df.net <- subsetCommunication(cellchat)

#在信号通路水平上推断细胞间的通讯
cellchat <- computeCommunProbPathway(cellchat)
##汇总细胞间的通讯
cellchat <- aggregateNet(cellchat)

#计算聚合细胞互作通信网络
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
#互作的数量
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")


###
##每个细胞如何跟别的细胞互作（number+of+interaction图）
mat <- cellchat@net$count
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}



##气泡图EC-KDR-IGFBP3
levels(cellchat@idents)
netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1:5), remove.isolate = FALSE)



setwd("D:/R/001/picture")
pdf(file="9A.pdf")
netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1:5), remove.isolate = FALSE)
dev.off()



##气泡图ACKR1
levels(cellchat@idents)
netVisual_bubble(cellchat, sources.use = 1, targets.use = c(1:5), remove.isolate = FALSE)


pdf(file="9B.pdf")
netVisual_bubble(cellchat, sources.use = 1, targets.use = c(1:5), remove.isolate = FALSE)
dev.off()
#查看信号通路
cellchat@netP$pathways




#计算和可视化网络中心性得分
pathways.show <- c("APP","CD99","JAM","NOTCH","CDH5","ESAM",  
                   "CD34",   "CD46" ,  "EPHB" ,  "PECAM1", "PTPRM",  "SEMA6", 
                   "SELE",   "MPZ"    ,"EPHA" ,  "VCAM" ,  "CEACAM" ,"SEMA4") 
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("APP","CD99","JAM","NOTCH","CDH5","ESAM",  
                                                                 "CD34",   "CD46" ,  "EPHB" ,  "PECAM1", "PTPRM",  "SEMA6", 
                                                                 "SELE",   "MPZ"    ,"EPHA" ,  "VCAM" ,  "CEACAM" ,"SEMA4"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg2



pdf(file="9C.pdf",width=9,height=7)
gg2
dev.off()