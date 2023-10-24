library(Seurat)
library(ggplot2)
library(clustree)
library(dplyr)
library(cowplot)
options(scipen = 200)

{setwd("D:/R/001/data/132465")
  library(data.table)
  CRC61<-fread(file="GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt",header=T,sep="\t")
  
  
  b61<-read.table("GSE132465_GEO_processed_CRC_10X_cell_annotation.txt",header=T,sep="\t")
  
  CRC62<-CRC61[,-1]
  
  row.names(CRC62)<-CRC61$Index
  
  identical(colnames(CRC62),b61$Index)
  #将两个表的探针建立对应关系
  
  index2<-which(b61$Class=='Normal')
  tum62.group=b61[index2,]
  head(tum62.group)
  #b61表中分出癌旁组织组
  
  CRC62.tum=CRC62[,index2,with=F]
  row.names(CRC62.tum)<-CRC61$Index
  #CRC62表中也分出癌旁组织组
  
  sce_meta=data.frame(PatientTypeID=tum62.group$Sample,row.names = tum62.group$Index)
  PBMC62<-CreateSeuratObject(counts = CRC62.tum,meta.data = sce_meta,min.cells = 3,min.features = 50)
  
  
  PBMC62[["percent.mt"]]<-PercentageFeatureSet(PBMC62,pattern = "^MT-")
  
  
  
  table(PBMC62@meta.data$percent.mt)
  VlnPlot(object = PBMC62,features = c("nFeature_RNA","nCount_RNA","percent.mt"),group.by = 'PatientTypeID')
  VlnPlot(object = PBMC62,features = c("nFeature_RNA","nCount_RNA","percent.mt"),group.by = 'orig.ident')
  
  PBMC621<-subset(PBMC62,subset=nFeature_RNA>500&nFeature_RNA<6000&percent.mt<20)
  
  plot1<-FeatureScatter(PBMC621,feature1 = "nCount_RNA",feature2 = "percent.mt",group.by = 'PatientTypeID',pt.size = 1.5)
  plot1
  plot2<-FeatureScatter(PBMC621,feature1 = "nCount_RNA",feature2 = "nFeature_RNA",group.by = 'PatientTypeID',pt.size = 1.5)
  plot2
  
  rm(b61,CRC61,CRC62,CRC62.tum,plot1,plot2,sce_meta,tum62.group,index2,PBMC62)
}#GSE132465

{setwd("D:/R/001/data/144735")
  CRC63<-read.table("GSE144735_processed_KUL3_CRC_10X_raw_UMI_count_matrix.txt",header=T,sep="\t")
  
  b63<-read.table("GSE144735_processed_KUL3_CRC_10X_annotation.txt",header=T,sep="\t")
  
  CRC64<-CRC63[,-1]
  
  row.names(CRC64)<-CRC63$Index
  
  colnames(CRC64)=b63$Index
  
  identical(colnames(CRC64),b63$Index)
  #将两个表的探针建立对应关系
  
  index3<-which(b63$Class=='Normal')
  tum64.group=b63[index3,]
  head(tum64.group)
  #b63表中分出肿瘤组织组
  
  CRC64.tum=CRC64[,index3]
  row.names(CRC64.tum)<-CRC63$Index
  #CRC64表中也分出肿瘤组织组
  
  sce_meta=data.frame(PatientTypeID=tum64.group$Sample,row.names = tum64.group$Index)
  PBMC64<-CreateSeuratObject(counts = CRC64.tum,meta.data = sce_meta,min.cells = 3,min.features = 50)
  
  PBMC64[["percent.mt"]]<-PercentageFeatureSet(PBMC64,pattern = "^MT-")
  
  table(PBMC64@meta.data$percent.mt)
  VlnPlot(object = PBMC64,features = c("nFeature_RNA","nCount_RNA","percent.mt"),group.by = 'PatientTypeID')
  VlnPlot(object = PBMC64,features = c("nFeature_RNA","nCount_RNA","percent.mt"),group.by = 'orig.ident')
  PBMC64<-subset(PBMC64,subset=nFeature_RNA>500&nFeature_RNA<6000&percent.mt<20)
  
  plot1<-FeatureScatter(PBMC64,feature1 = "nCount_RNA",feature2 = "percent.mt",group.by = 'PatientTypeID',pt.size = 1.5)
  plot1
  plot2<-FeatureScatter(PBMC64,feature1 = "nCount_RNA",feature2 = "nFeature_RNA",group.by = 'PatientTypeID',pt.size = 1.5)
  plot2
  
  rm(b63,CRC63,CRC64,CRC64.tum,plot1,plot2,sce_meta,tum64.group,index3)
}#GSE144735

{setwd("D:/R/001/data/178341")
  library(hdf5r)
  CRC20<-Seurat::Read10X_h5(filename = "GSE178341_crc10x_full_c295v4_submit.h5",use.names = T)
  
  c1<-read.csv("GSE178341_crc10x_full_c295v4_submit_metatables.csv",header=T,sep=",")
  
  CRC20<-CreateSeuratObject(counts = CRC20,min.cells = 3,min.features = 200)
  
  row.names(c1)<-c1$cellID
  CRC20<-AddMetaData(CRC20,metadata =c1 )
  
  
  
  head(CRC20@meta.data, 8)
  tail(CRC20@meta.data, 8)
  table(CRC20@meta.data$orig.ident) #查看细胞数？
  
  CRC20=subset(x = CRC20, subset = SPECIMEN_TYPE == "N")
  
  
  
  CRC20[["percent.mt"]]<-PercentageFeatureSet(CRC20,pattern = "^MT")
  table(CRC20@meta.data$percent.mt)
  
  
  VlnPlot(object = CRC20,features = c("nFeature_RNA","nCount_RNA","percent.mt"))
  
  
  CRC20<-subset(CRC20,subset=nFeature_RNA>500&nFeature_RNA<6000&percent.mt<20)
  
  plot1<-FeatureScatter(CRC20,feature1 = "nCount_RNA",feature2 = "percent.mt",group.by = 'PatientTypeID',pt.size = 1.5)
  plot1
  plot2<-FeatureScatter(CRC20,feature1 = "nCount_RNA",feature2 = "nFeature_RNA",group.by = 'PatientTypeID',pt.size = 1.5)
  plot2
  
  rm(c1,plot1,plot2)
}#GSE178341


#合并PBMC64,PBMC621
ALL1<-merge(PBMC64,PBMC621,add.cell.ids = c("PBMC64","PBMC621"),project = "ALL")
rm(PBMC64,PBMC621)
table(ALL1$PatientTypeID)

#合并FF1
ALL2<-merge(ALL1,CRC20,add.cell.ids = c("ALL1","CRC20"),project = "ALL")
rm(CRC20,ALL1)
table(ALL2$PatientTypeID)




ALL1 <- ALL2
ALL1 <-  ALL1%>%
  Seurat::NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2500) %>% 
  ScaleData()
ALL1 <- RunPCA(ALL1, npcs = 50)
print(ALL1[["pca"]],dims=1:5,nfeatures=5)
#PCA


if(!require(harmony))devtools::install_github("immunogenomics/harmony")

ALL1=ALL1%>%RunHarmony("PatientTypeID",plot_convergence=TRUE)

ALL1<-ALL1%>%
  RunUMAP(reduction="harmony",dims=1:30)%>%
  FindNeighbors(reduction="harmony",dims=1:30)%>%
  FindClusters(resolution=0.5)%>%
  identity()
#去批次


##
DimPlot(ALL1,reduction = "umap",label = TRUE,pt.size = 0.01)

setwd("D:/R/001/picture")
pdf(file="S3A.pdf")
DimPlot(ALL1,reduction = "umap",label = TRUE,pt.size = 0.01)
dev.off()


p3<-DimPlot(ALL1,reduction = "umap",group.by = "PatientTypeID",pt.size = 0.01)+theme(
  axis.line = element_blank(),
  axis.ticks = element_blank(),axis.text = element_blank()
)

p4<-DimPlot(ALL1,reduction = "umap",group.by = "ident",pt.size = 0.01,label = TRUE,repel = TRUE)+theme(
  axis.line = element_blank(),
  axis.ticks = element_blank(),axis.text = element_blank()
)
p3|p4




DotPlot(object=ALL1,features=c("PECAM1","COL1A1","VWF"),group.by="seurat_clusters")

pdf(file="S3B.pdf")
DotPlot(object=ALL1,features=c("PECAM1","COL1A1","VWF"),group.by="seurat_clusters")
dev.off()

#提取分群

Idents(ALL1)<-ALL1$seurat_clusters

asstromalcell=subset(x = ALL1, idents=c("22","20","18","17","14","12","10","6"))#基质

setwd("D:/R/001/workspace")
saveRDS(asstromalcell,file = "normalasstromalcell.Rds")
asstromalcell<-readRDS("normalasstromalcell.Rds")


#二次聚类
{str <- asstromalcell
  str <-  str%>%
    Seurat::NormalizeData() %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    ScaleData()
  str <- RunPCA(str, npcs = 50)
  print(str[["pca"]],dims=1:5,nfeatures=5)
  #PCA
  
  
  if(!require(harmony))devtools::install_github("immunogenomics/harmony")
  
  str=str%>%RunHarmony("PatientTypeID",plot_convergence=TRUE)
  
  str<-str%>%
    RunUMAP(reduction="harmony",dims=1:20)%>%
    FindNeighbors(reduction="harmony",dims=1:20)%>%
    FindClusters(resolution=0.4)%>%
    identity()
  #去批次
  
  DimPlot(str,reduction = "umap",label = TRUE,pt.size = 0.01)
  
  
  setwd("D:/R/001/picture")
  pdf(file="S3C.pdf")
  DimPlot(str,reduction = "umap",label = TRUE,pt.size = 0.01)
  dev.off()
  
  
  p3<-DimPlot(str,reduction = "umap",group.by = "PatientTypeID",pt.size = 0.01)+theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),axis.text = element_blank()
  )
}#基质


#手动注释
{#Marker2手动注释
  str.markers<-FindAllMarkers(str,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
  str.markers%>%group_by(cluster)%>%top_n(n=10,wt=avg_log2FC)
  
  top10<-str.markers%>%group_by(cluster)%>%top_n(n=10,wt=avg_log2FC)
  
  
  
}#基质





{#(10,6,8,2/9)
  pdf(file="S3D.pdf",width = 8,height = 10)
  DotPlot(object=str,
          features=c("SELE","ACKR1","VWF","CD74","PECAM1","ESM1"
                     ,"HES1","CD34","KDR","INSR","FLT1","PODXL"
                     ,"IGFBP3","CLDN5","PTPRB","TM4SF1","CXCL12","GJA5","ID1","PLCG2",
                     "MCAM","PCLAF","UBE2C","TK1","STMN1","CENPK"
                     ,"FABP4","PROX1","MMRN1","TFF3","CCL21","LYVE1","COLEC12"),group.by="seurat_clusters")+
    coord_flip()
  
  dev.off()
}#基质

setwd("D:/R/001/workspace")
saveRDS(str,file = "normalstr.Rds")
str<-readRDS("normalstr.Rds")


str$new.cluster.ids<-str$seurat_clusters
library(tidyverse)
Idents(str)<-'seurat_clusters'
str$new.cluster.ids<-recode(str$new.cluster.ids,'2'="EC-KDR-1",
                            '6'="EC-ACKR1",'8'="EC-IGFBP3",'9'="EC-KDR-2",
                            '10'="EC-TFF3")
Idents(str)<-'new.cluster.ids'


ecc=subset(x = str, idents=c("EC-ACKR1","EC-KDR-1","EC-KDR-2","EC-IGFBP3","EC-TFF3"))
DimPlot(ecc,reduction = "umap",label = TRUE,pt.size = 0.01)
###气泡热图绘制


#图1E

col <- c("#a27bb1","#f6b878","#96cd86","#e4a0ca","#6b7682")
#96cd86  green

ecc=subset(x = str, idents=c("EC-ACKR1","EC-KDR-1","EC-IGFBP3","EC-KDR-2","EC-TFF3"))
DimPlot(ecc,reduction = "umap",label = TRUE,pt.size = 0.5,cols=col)

setwd("D:/R/001/picture")
pdf(file="1E.pdf",width = 8,height = 8)
DimPlot(ecc,reduction = "umap",label = TRUE,pt.size = 0.5,cols=col)
dev.off()

#图1F
#pdf(file="1F.pdf",height = 10,width = 9)
#DotPlot(object=ecc,features=c("SELE","ACKR1","VWF","CD74","PECAM1","ESM1"
#                              ,"HES1","CD34","KDR","INSR","FLT1","PODXL"
#                              ,"IGFBP3","CLDN5","PTPRB","TM4SF1","CXCL12","GJA5","ID1","PLCG2",
#                              "MCAM","PCLAF","UBE2C","TK1","STMN1","CENPK"
#                              ,"FABP4","PROX1","MMRN1","TFF3","CCL21","LYVE1","COLEC12"),group.by="new.cluster.ids")+
#  scale_color_gradientn(colours = c('navy',"white",'firebrick3'))+
#  coord_flip()
#dev.off()



pdf(file="1F.pdf",height = 4,width = 27)
DotPlot(object=ecc,features=c("SELE","ACKR1","VWF","CD74","PECAM1","ESM1"
                              ,"HES1","CD34","KDR","INSR","FLT1","PODXL"
                              ,"IGFBP3","CLDN5","PTPRB","TM4SF1","CXCL12","GJA5","ID1","PLCG2",
                              "MCAM","PCLAF","UBE2C","TK1","STMN1","CENPK"
                              ,"FABP4","PROX1","MMRN1","TFF3","CCL21","LYVE1","COLEC12"),group.by="new.cluster.ids")+
  scale_color_gradientn(colours = c('navy',"white",'firebrick3'))
dev.off()