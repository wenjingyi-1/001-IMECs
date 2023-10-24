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
  
  index2<-which(b61$Class=='Tumor')
  tum62.group=b61[index2,]
  head(tum62.group)
  #b61表中分出肿瘤组织组
  
  CRC62.tum=CRC62[,index2,with=F]
  row.names(CRC62.tum)<-CRC61$Index
  #CRC62表中也分出肿瘤组织组
  
  sce_meta=data.frame(PatientTypeID=tum62.group$Sample,row.names = tum62.group$Index)
  P62<-CreateSeuratObject(counts = CRC62.tum,meta.data = sce_meta,min.cells = 3,min.features = 50)
  
  
  P62[["percent.mt"]]<-PercentageFeatureSet(P62,pattern = "^MT-")
  
  
  
  table(P62@meta.data$percent.mt)
  VlnPlot(object = P62,features = c("nFeature_RNA","nCount_RNA","percent.mt"),group.by = 'PatientTypeID')
  VlnPlot(object = P62,features = c("nFeature_RNA","nCount_RNA","percent.mt"),group.by = 'orig.ident')
  
  P621<-subset(P62,subset=nFeature_RNA>500&nFeature_RNA<6000&percent.mt<20)
  
  plot1<-FeatureScatter(P621,feature1 = "nCount_RNA",feature2 = "percent.mt",group.by = 'PatientTypeID',pt.size = 1.5)
  plot1
  plot2<-FeatureScatter(P621,feature1 = "nCount_RNA",feature2 = "nFeature_RNA",group.by = 'PatientTypeID',pt.size = 1.5)
  plot2
  
  rm(b61,CRC61,CRC62,CRC62.tum,plot1,plot2,sce_meta,tum62.group,index2,P62)
}#GSE132465

{setwd("D:/R/001/data/144735")
  
  CRC63<-read.table("GSE144735_processed_KUL3_CRC_10X_raw_UMI_count_matrix.txt",header=T,sep="\t")
  
  b63<-read.table("GSE144735_processed_KUL3_CRC_10X_annotation.txt",header=T,sep="\t")
  
  CRC64<-CRC63[,-1]
  
  row.names(CRC64)<-CRC63$Index
  
  colnames(CRC64)=b63$Index
  
  identical(colnames(CRC64),b63$Index)
  #将两个表的探针建立对应关系
  
  index3<-which(b63$Class=='Tumor')
  tum64.group=b63[index3,]
  head(tum64.group)
  #b63表中分出肿瘤组织组
  
  CRC64.tum=CRC64[,index3]
  row.names(CRC64.tum)<-CRC63$Index
  #CRC64表中也分出肿瘤组织组
  
  sce_meta=data.frame(PatientTypeID=tum64.group$Sample,row.names = tum64.group$Index)
  P64<-CreateSeuratObject(counts = CRC64.tum,meta.data = sce_meta,min.cells = 3,min.features = 50)
  
  P64[["percent.mt"]]<-PercentageFeatureSet(P64,pattern = "^MT-")
  
  table(P64@meta.data$percent.mt)
  VlnPlot(object = P64,features = c("nFeature_RNA","nCount_RNA","percent.mt"),group.by = 'PatientTypeID')
  VlnPlot(object = P64,features = c("nFeature_RNA","nCount_RNA","percent.mt"),group.by = 'orig.ident')
  P64<-subset(P64,subset=nFeature_RNA>500&nFeature_RNA<6000&percent.mt<20)
  
  plot1<-FeatureScatter(P64,feature1 = "nCount_RNA",feature2 = "percent.mt",group.by = 'PatientTypeID',pt.size = 1.5)
  plot1
  plot2<-FeatureScatter(P64,feature1 = "nCount_RNA",feature2 = "nFeature_RNA",group.by = 'PatientTypeID',pt.size = 1.5)
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
  
  CRC20=subset(x = CRC20, subset = SPECIMEN_TYPE == "T")
  
  
  
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

{setwd("D:/R/001/data/188711")
  #06
  CRC06<-Read10X(data.dir = '5688706')
  P06<-CreateSeuratObject(counts = CRC06,project = "P06",min.cells = 3,min.features = 200)
  sce_meta06=data.frame(PatientTypeID=P06$orig.ident)
  P06<-AddMetaData(P06,metadata =sce_meta06 )
  rm(CRC06,sce_meta06)
  
  
  #07
  CRC07<-Read10X(data.dir = '5688707')
  P07<-CreateSeuratObject(counts = CRC07,project = "P07",min.cells = 3,min.features = 200)
  sce_meta07=data.frame(PatientTypeID=P07$orig.ident)
  P07<-AddMetaData(P07,metadata =sce_meta07 )
  rm(CRC07,sce_meta07)
  
  #08
  CRC08<-Read10X(data.dir = '5688708')
  P08<-CreateSeuratObject(counts = CRC08,project = "P08",min.cells = 3,min.features = 200)
  sce_meta08=data.frame(PatientTypeID=P08$orig.ident)
  P08<-AddMetaData(P08,metadata =sce_meta08 )
  rm(CRC08,sce_meta08)
  
  #09
  CRC09<-Read10X(data.dir = '5688709')
  P09<-CreateSeuratObject(counts = CRC09,project = "P09",min.cells = 3,min.features = 200)
  sce_meta09=data.frame(PatientTypeID=P09$orig.ident)
  P09<-AddMetaData(P09,metadata =sce_meta09 )
  rm(CRC09,sce_meta09)
  
  #10
  CRC10<-Read10X(data.dir = '5688710')
  P10<-CreateSeuratObject(counts = CRC10,project = "P10",min.cells = 3,min.features = 200)
  sce_meta10=data.frame(PatientTypeID=P10$orig.ident)
  P10<-AddMetaData(P10,metadata =sce_meta10 )
  rm(CRC10,sce_meta10)
  
  #11
  CRC11<-Read10X(data.dir = '5688711')
  P11<-CreateSeuratObject(counts = CRC11,project = "P11",min.cells = 3,min.features = 200)
  sce_meta11=data.frame(PatientTypeID=P11$orig.ident)
  P11<-AddMetaData(P11,metadata =sce_meta11 )
  rm(CRC11,sce_meta11)
  
  #合并06,07
  FF1<-merge(P06,P07,add.cell.ids = c("P06","P07"),project = "ALL")
  rm(P06,P07)
  table(FF1$PatientTypeID)
  
  #合并08
  FF2<-merge(FF1,P08,add.cell.ids = c("FF1","P08"),project = "ALL")
  rm(P08)
  table(FF2$PatientTypeID)
  
  #合并09
  FF1<-merge(FF2,P09,add.cell.ids = c("FF2","P09"),project = "ALL")
  rm(P09)
  table(FF1$PatientTypeID)
  
  #合并10
  FF2<-merge(FF1,P10,add.cell.ids = c("FF1","P10"),project = "ALL")
  rm(P10)
  table(FF2$PatientTypeID)
  
  #合并11
  FF1<-merge(FF2,P11,add.cell.ids = c("FF2","P11"),project = "ALL")
  rm(P11)
  table(FF1$PatientTypeID)
  
  rm(FF2)
  #########
  FF1[["percent.mt"]]<-PercentageFeatureSet(FF1,pattern = "^MT")
  table(FF1@meta.data$percent.mt)
  
  
  VlnPlot(object = FF1,features = c("nFeature_RNA","nCount_RNA","percent.mt"))
  
  
  FF1<-subset(FF1,subset=nFeature_RNA>500&nFeature_RNA<6000&percent.mt<20)
  
  plot1<-FeatureScatter(FF1,feature1 = "nCount_RNA",feature2 = "percent.mt",group.by = 'PatientTypeID',pt.size = 1.5)
  plot1
  plot2<-FeatureScatter(FF1,feature1 = "nCount_RNA",feature2 = "nFeature_RNA",group.by = 'PatientTypeID',pt.size = 1.5)
  plot2
  
  rm(plot1,plot2)
}#GSE188711


{#合并PBMC64,PBMC621
  ALL1<-merge(P64,P621,add.cell.ids = c("P64","P621"),project = "ALL")
  rm(P64,P621)
  table(ALL1$PatientTypeID)
  
  #合并FF1
  ALL2<-merge(ALL1,FF1,add.cell.ids = c("ALL1","FF1"),project = "ALL")
  rm(FF1,ALL1)
  table(ALL2$PatientTypeID)
  
  #合并CRC20
  ALL1<-merge(ALL2,CRC20,add.cell.ids = c("ALL2","CRC20"),project = "ALL")
  rm(CRC20,ALL2)
  table(ALL1$PatientTypeID)
}#合并

ALL1 <- ALL1
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
pdf(file="S1A.pdf")
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

DotPlot(object=ALL1,features=c("CD3D","CD8A","CD4","FOXP3","TRDC","NKG7","CD79A","MS4A1","IGHG4"),group.by="seurat_clusters")
DotPlot(object=ALL1,features=c("CD14","FCGR3A","CD68","CD163","CD1C","LAMP3","TPSAB1","CSF3R","S100A8"),group.by="seurat_clusters")
DotPlot(object=ALL1,features=c("PECAM1","COL1A1","VWF"),group.by="seurat_clusters")
DotPlot(object=ALL1,features=c("EPCAM"),group.by="seurat_clusters")


pdf(file="S1B.pdf",width = 10,height = 5)
DotPlot(object=ALL1,features=c("CD3D","CD8A","CD4","FOXP3","TRDC","NKG7","CD79A","MS4A1","IGHG4"),group.by="seurat_clusters")
dev.off()


pdf(file="S1C.pdf",width = 10,height = 5)
DotPlot(object=ALL1,features=c("CD14","FCGR3A","CD68","CD163","CD1C","LAMP3","TPSAB1","CSF3R","S100A8"),group.by="seurat_clusters")
dev.off()

pdf(file="S1D.pdf")
DotPlot(object=ALL1,features=c("PECAM1","COL1A1","VWF"),group.by="seurat_clusters")
dev.off()

pdf(file="S1E.pdf")
DotPlot(object=ALL1,features=c("EPCAM"),group.by="seurat_clusters")
dev.off()

#Marker2手动注释
all1.markers<-FindAllMarkers(ALL1,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
all1.markers%>%group_by(cluster)%>%top_n(n=10,wt=avg_log2FC)

top10<-all1.markers%>%group_by(cluster)%>%top_n(n=10,wt=avg_log2FC)

setwd("D:/R/001/txt")
write.table( all1.markers,file = "firstcluster.markers.txt",row.names = F,quote = F)
all1.markers<-read.table("firstcluster.markers.txt",header = T)

setwd("D:/R/001/csv")
write.csv(all1.markers,'firstcluster.markers.csv',row.names = F)
all1.markers<-read.csv("firstcluster.markers.csv",header = T)


#提取分群

Idents(ALL1)<-ALL1$seurat_clusters

asMyeloidcell=subset(x = ALL1, idents=c("3","5","14","16","18","19"))#髓系

aslymphocytecell=subset(x = ALL1, idents=c("1","2","4","6","7","10","20","21","23","22"))#淋巴

asepithelialcell=subset(x = ALL1, idents=c("0","8","9","11","17"))#上皮

asstromalcell=subset(x = ALL1, idents=c("12","13","15"))#基质

setwd("D:/R/001/workspace")
saveRDS(aslymphocytecell,file = "aslymphocytecell.Rds")
aslymphocytecell<-readRDS("aslymphocytecell.Rds")


saveRDS(asstromalcell,file = "asstromalcell.Rds")
asstromalcell<-readRDS("asstromalcell.Rds")



ALL1$new.cluster.ids<-ALL1$seurat_clusters
library(tidyverse)
Idents(ALL1)<-'seurat_clusters'
ALL1$new.cluster.ids<-recode(ALL1$new.cluster.ids,'0'="epithelialcell",'1'="lymphocytecell",'2'="lymphocytecell",'3'="Myeloidcell",'4'="lymphocytecell"
                             ,'5'="Myeloidcell",'6'="lymphocytecell",'7'="lymphocytecell",'8'="epithelialcell",'9'="epithelialcell",
                             '10'="lymphocytecell",'11'="epithelialcell",'12'="stromalcell",'13'="stromalcell",'14'="Myeloidcell",'15'="stromalcell"
                             ,'16'="Myeloidcell",'17'="epithelialcell",'18'="Myeloidcell",'19'="Myeloidcell",'20'="lymphocytecell",'21'="lymphocytecell",'22'="lymphocytecell",'23'="lymphocytecell")
Idents(ALL1)<-'new.cluster.ids'

DimPlot(ALL1,reduction = "umap",label = TRUE,pt.size = 0.01)


########图1A
col <- c("#eb7b1d","#32a9df","#3fa64c","#dca818")
DimPlot(ALL1,reduction = "umap",label = TRUE,pt.size = 0.01,cols=col)

setwd("D:/R/001/picture")
pdf(file="1A.pdf")
col <- c("#eb7b1d","#32a9df","#3fa64c","#dca818")
DimPlot(ALL1,reduction = "umap",label = TRUE,pt.size = 0.01,cols=col)
dev.off()


#二次聚类
setwd("D:/R/001/workspace")
aslymphocytecell<-readRDS("aslymphocytecell.Rds")

{lym <- aslymphocytecell
  lym <-  lym%>%
    Seurat::NormalizeData() %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    ScaleData()
  lym <- RunPCA(lym, npcs = 50)
  print(lym[["pca"]],dims=1:5,nfeatures=5)
  #PCA
  
  
  if(!require(harmony))devtools::install_github("immunogenomics/harmony")
  
  lym=lym%>%RunHarmony("PatientTypeID",plot_convergence=TRUE)
  
  lym<-lym%>%
    RunUMAP(reduction="harmony",dims=1:20)%>%
    FindNeighbors(reduction="harmony",dims=1:20)%>%
    FindClusters(resolution=1.5)%>%
    identity()
  #去批次
  

  DimPlot(lym,reduction = "umap",label = TRUE,pt.size = 0.01)
  
  setwd("D:/R/001/picture")
  pdf(file="S5A.pdf",width = 10,height = 10)
  DimPlot(lym,reduction = "umap",label = TRUE,pt.size = 0.01)
  dev.off()
  
  
  p3<-DimPlot(lym,reduction = "umap",group.by = "PatientTypeID",pt.size = 0.01)+theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),axis.text = element_blank()
  )
}#lym


#手动注释
{#Marker2手动注释
  lym.markers<-FindAllMarkers(lym,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
  lym.markers%>%group_by(cluster)%>%top_n(n=10,wt=avg_log2FC)
  
  top10<-lym.markers%>%group_by(cluster)%>%top_n(n=10,wt=avg_log2FC)
  
}#lym


###气泡热图绘制
{   #19(CD4T-CXCL13-TCF7)  #7(CD4T-CXCL13)
  pdf(file="S5B.pdf",height = 8,width = 8)
  DotPlot(object=lym,
          features=c("TCF7","MAL","CCR7","WDR74","CXCL13","PDCD1","CD4","CD8A"
          ),group.by="seurat_clusters")
  dev.off()
  

  
  
  #3(CD4T-FOXP3),13(CD4T-FOXP3-CTLA4),29(CD4T-FOXP3-STMN1)
  pdf(file="S5C.pdf",height = 8,width = 9)
  DotPlot(object=lym,
          features=c("ENTPD1","RTKN2","FOXP3","TNFRSF18" ,"CTLA4","IL2RA","CD4","CD8A" ),group.by="seurat_clusters")
  dev.off()
  
  #2(CD4T-CCR7),1(CD4T-CD69)
  pdf(file="S5D.pdf",height = 8,width = 9)
  DotPlot(object=lym,
          features=c("TCF7","MAL","CCR7","SELL","WDR74" ,"GPR183","CCR6","CD69","CD4","CD8A" ),group.by="seurat_clusters")
  dev.off()
  
  
}#淋巴

lym$new.cluster.ids<-lym$seurat_clusters
library(tidyverse)
Idents(lym)<-'seurat_clusters'
lym$new.cluster.ids<-recode(lym$new.cluster.ids,'0'="memoryB",'1'="CD4T-CD69",'2'="CD4T-CCR7",'3'="CD4T-FOXP3",'4'="CD8T"
                            ,'5'="CD8T-PDCD1",'6'="CD8T",'7'="CD4T-CXCL13",'8'="IgA-PC",'9'="IgA-PC",
                            '10'="NK",'11'="naiveB",'12'="T-STMN1",'13'="CD4T-FOXP3-CTLA4",'14'="CD8T",'15'="IgG-PC",'16'="T-STMN1"
                            ,'17'="NK",'18'="IgG-PC",'19'="CD4T-CXCL13-TCF7",'20'="memoryB",'21'="unuse",
                            '22'="unuse",'23'="B-STMN1",'24'="IgM-PC",'25'="MAIT",'26'="memoryB",'27'="B-STMN1",
                            '28'="memoryB",'29'="CD4T-FOXP3-STMN1",'30'="IgA-PC",'31'="IgA-PC")



Idents(lym)<-'new.cluster.ids'

lym=subset(x = lym ,idents=c("CD4T-CXCL13-TCF7","CD4T-CXCL13","CD4T-FOXP3","CD4T-FOXP3-CTLA4","CD4T-FOXP3-STMN1",
                             "CD4T-CCR7","CD4T-CD69"))

DimPlot(lym,reduction = "umap",label = TRUE,pt.size = 0.01)

setwd("D:/R/001/workspace")
saveRDS(lym,file = "lym.Rds")
lym<-readRDS("lym.Rds")


#二次聚类
asstromalcell<-readRDS("asstromalcell.Rds")

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
  pdf(file="S2A.pdf")
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
  
  setwd("D:/R/001/txt")
  write.table( str.markers,file = "str.markers.txt",row.names = F,quote = F)
  setwd("D:/R/001/csv")
  write.csv(str.markers,'str.markers.csv',row.names = F)
  
  
}#基质

setwd("D:/R/001/workspace")
saveRDS(str,file = "str.Rds")
str<-readRDS("str.Rds")



###气泡热图绘制

setwd("D:/R/001/picture")
{#1，3，7，9，15内皮细胞
  pdf(file="S2B.pdf",width = 8,height = 10)
 DotPlot(object=str,
                  features=c("SELE","ACKR1","VWF","CD74","PECAM1","ESM1"
                             ,"HES1","CD34","KDR","INSR","FLT1","PODXL"
                             ,"IGFBP3","CLDN5","PTPRB","TM4SF1","GJA5","CXCL12","ID1","PLCG2",
                             "MCAM","PCLAF","UBE2C","TK1","STMN1","CENPK"
                             ,"FABP4","PROX1","MMRN1","TFF3","CCL21","LYVE1","COLEC12"),group.by="seurat_clusters")+
   coord_flip()
          
  
  dev.off()
  
  #0,14,10,16周细胞、SMC、内皮-周皮标记基因共表达细胞群
  pdf(file="S2C.pdf",height = 7,width = 8)
  DotPlot(object=str,
          features=c("MYL9","NDUFA4L2","MCAM","NOTCH3","SOD3","MYH11"
                     ,"ABCC9","KCNJ8","CNN1","IGFBP5","STMN1","ACTA2"
                     ,"S1PR3","PDGFRB","RGS5"),group.by="seurat_clusters")+
    coord_flip()
  dev.off()
  
  
  #11,18，12,13舍弃的细胞
  pdf(file="S2D.pdf",height = 4,width = 8)
  DotPlot(object=str,
          features=c("EPCAM","IGHA1","JCHAIN","LYZ","TYROBP","HLA-DRB1","ACTA2" ),group.by="seurat_clusters")+
  coord_flip()
  dev.off()
  
  #2,4,6成纤维
  pdf(file="S2E.pdf",height = 7,width = 8)
  DotPlot(object=str,
          features=c("C1R","MMP11","FAP","CTHRC1","COL1A1","MMP2"
                     ,"VCAN","COL3A1","PDGFRA","DCN","LUM","COL11A1"
                     ,"TGFB1","POSTN","VEGFA","VEGFB","CTGF","PDPN"),group.by="seurat_clusters")+
    coord_flip()
  dev.off()
  #5,8,17成纤维
  pdf(file="S2F.pdf",height = 7,width = 8)
  DotPlot(object=str,
          features=c("SERPING1","IGFBP6","C3","C1R","C7","CCL2"
                     ,"DCN","LUM","PDGFRA","MMP2","MMP11"
                     ,"CXCL12","CFD","LGI4","MPZ","S100B","CD74"
                     ,"PLP1","SOX2","SOX10"),group.by="seurat_clusters")+
  coord_flip()
  dev.off()
}#基质

#图1B
str$new.cluster.ids<-str$seurat_clusters
library(tidyverse)
Idents(str)<-'seurat_clusters'
str$new.cluster.ids<-recode(str$new.cluster.ids,'0'="Mu-MCAM-NDUFA4L2",'1'="EC-ACKR1",'2'="Fb-FAP-MMP11",'3'="EC-KDR-ESM1",'4'="Fb-MMP11"
                            ,'5'="Fb-KLF4",'6'="Fb-FAP",'7'="EC-KDR-IGFBP3",'8'="Fb-CXCL12",'9'="EC-STMN1",
                            '10'="SMC",'11'="immune like",'12'="undefined-2",'13'="undefined-1",'14'="EC - Mu",'15'="EC-TFF3"
                            ,'16'="Mu-NDUFA4L2",'17'="Fb-MPZ",'18'="epithelial like")
Idents(str)<-'new.cluster.ids'
DimPlot(str,reduction = "umap",label = TRUE,pt.size = 0.5)

pdf(file="1B.pdf",width = 9,height = 9)
DimPlot(str,reduction = "umap",label = TRUE,pt.size = 0.5)
dev.off()


#图1C

col <- c("#f6b878","#d27bb1","#96cd86","#16b4c2","#6b7682")
ecc=subset(x = str, idents=c("EC-ACKR1","EC-KDR-ESM1","EC-KDR-IGFBP3","EC-STMN1","EC-TFF3"))
DimPlot(ecc,reduction = "umap",label = TRUE,pt.size = 0.5,cols=col)

pdf(file="1C.pdf",width = 8,height = 8)
DimPlot(ecc,reduction = "umap",label = TRUE,pt.size = 0.5,cols=col)
dev.off()



###图1D
#pdf(file="1D.pdf",height = 10,width = 11)
#DotPlot(object=ecc,features=c("SELE","ACKR1","VWF","CD74","PECAM1","ESM1"
#                              ,"HES1","CD34","KDR","INSR","FLT1","PODXL"
#                              ,"IGFBP3","CLDN5","PTPRB","TM4SF1","CXCL12","GJA5","ID1","PLCG2",
#                              "MCAM","PCLAF","UBE2C","TK1","STMN1","CENPK"
#                              ,"FABP4","PROX1","MMRN1","TFF3","CCL21","LYVE1","COLEC12"),group.by="new.cluster.ids")+
#  scale_color_gradientn(colours = c('navy',"white",'firebrick3'))+
#coord_flip()
#dev.off()



pdf(file="1D.pdf",height = 4,width = 27)
DotPlot(object=ecc,features=c("SELE","ACKR1","VWF","CD74","PECAM1","ESM1"
                              ,"HES1","CD34","KDR","INSR","FLT1","PODXL"
                              ,"IGFBP3","CLDN5","PTPRB","TM4SF1","CXCL12","GJA5","ID1","PLCG2",
                              "MCAM","PCLAF","UBE2C","TK1","STMN1","CENPK"
                              ,"FABP4","PROX1","MMRN1","TFF3","CCL21","LYVE1","COLEC12"),group.by="new.cluster.ids")+
  scale_color_gradientn(colours = c('navy',"white",'firebrick3'))
dev.off()