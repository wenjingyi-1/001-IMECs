library(dplyr)
library(Seurat)
library(patchwork) #加载数据集 


###########################################
setwd("D:/R/001/workspace")
str<-readRDS("str.Rds")
str$new.cluster.ids<-str$seurat_clusters
library(tidyverse)
Idents(str)<-'seurat_clusters'
str$new.cluster.ids<-recode(str$new.cluster.ids,'0'="Mu",'1'="ECACKR1",'2'="Fb",'3'="ECKDRESM1",'4'="Fb"
                            ,'5'="Fb",'6'="Fb",'7'="ECKDRIGFBP3",'8'="Fb",'9'="ECSTMN1",
                            '10'="SMC",'11'="strundefined",'12'="strundefined",'13'="strundefined",'14'="ECMu",'15'="ECTFF3",'16'="Mu"
                            ,'17'="Fb",'18'="strundefined")
Idents(str)<-'new.cluster.ids'
str=subset(x = str ,idents=c("ECACKR1","ECKDRESM1","ECKDRIGFBP3","ECSTMN1","Mu","Fb","SMC","ECMu","ECTFF3"))

setwd("D:/R/001/txt")
all1.markers<-read.table("firstcluster.markers.txt",header = T)

pbmc<-str
rm("str")

########################################

df <- as.data.frame(GetAssayData(pbmc, assay = "RNA", slot = "counts"))#获取标准化后表达矩阵数据 
pbmc_idents <- as.data.frame(t(as.data.frame(Idents(object = pbmc))))#获取细胞ID和亚群的对应关系 
identical(colnames(df),colnames(pbmc_idents))#判断顺序 
colnames(df) <- pbmc_idents[1,]#把细胞ID替换成亚群名称
x <- data.frame(GeneSymbol = rownames(df)) 
rownames(df) <- NULL 
final <- cbind(x,df) #得到没有观测名，第一列是Gene symbol后面是细胞亚群对应的表达矩阵 


#这个就是进行CIBERSORTX的单细胞参考表达数据


index2<-which(all1.markers$cluster=='12')
index3<-which(all1.markers$cluster=='13')
index4<-which(all1.markers$cluster=='15')
all1.markers=all1.markers[c(index2,index3,index4),]



all1.markers<-all1.markers[,7]
all1.markers<-as.data.frame(all1.markers)
library(VennDiagram)

venn_list <- list(group1 = all1.markers$all1.markers, group2 = final$GeneSymbol)
inter <- get.venn.partitions(venn_list)
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')

interset_mRNA<-as.data.frame(inter$..values..[1])

interset_mRNA<-interset_mRNA[,1]

genshin=final
rownames(genshin)<-genshin$GeneSymbol
genshin=genshin[interset_mRNA,]



write.table( genshin,file = "cibersort.markers.txt",row.names = F,sep = "\t")








###补充，Signature matrix file的文件格式转换
setwd("D:/R/001/txt")
LMEC<-read.table("LMEC.txt",header = T)
setwd("D:/R/001/csv")
write.csv(LMEC,'LMEC.csv',row.names = F)

