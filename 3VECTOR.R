install.packages('circlize')   # 0.4.11
install.packages('gatepoints') # 0.1.3
install.packages('stringr')    # 1.4.0
install.packages('igraph')     # 1.2.6
install.packages('gmodels')    # 2.18.1


library(stringr)
library(Seurat)
library(ggplot2)
library(clustree)
library(dplyr)
library(cowplot)
options(scipen = 200)
setwd("D:/R/001/workspace")
asstromalcell<-readRDS("asstromalcell.Rds")

#二次聚类
{str <- asstromalcell
  str <-  str%>%
    Seurat::NormalizeData() %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    ScaleData()
  str <- RunPCA(str, npcs = 150)
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
  pdf(file="S3E.pdf")
  DimPlot(str,reduction = "umap",label = TRUE,pt.size = 0.01)
  dev.off()
  
  
  p3<-DimPlot(str,reduction = "umap",group.by = "PatientTypeID",pt.size = 0.01)+theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),axis.text = element_blank()
  )
}#基质


pdf(file="S3F.pdf",height = 10,width = 9)
DotPlot(object=str,
        features=c("SELE","ACKR1","VWF","CD74","PECAM1","ESM1"
                   ,"HES1","CD34","KDR","INSR","FLT1","PODXL"
                   ,"IGFBP3","CLDN5","PTPRB","TM4SF1","CXCL12","GJA5","ID1","PLCG2",
                   "MCAM","PCLAF","UBE2C","TK1","STMN1","CENPK"
                   ,"FABP4","PROX1","MMRN1","TFF3","CCL21","LYVE1","COLEC12"),group.by="seurat_clusters")+
coord_flip()
dev.off()




str$new.cluster.ids<-str$seurat_clusters
library(tidyverse)
Idents(str)<-'seurat_clusters'
#这一部分只对前文注释的内皮细胞聚类进行了严谨的判断，其余聚类名与亚群序号之间并无严格的对应关系。
str$new.cluster.ids<-recode(str$new.cluster.ids,'0'="Mu-MCAM-NDUFA4L2",'1'="EC-ACKR1",'2'="EC-KDR-ESM1",'3'="later",'4'="Fb-MMP11"
                            ,'5'="Fb-KLF4",'6'="Fb-FAP",'7'="EC-KDR-IGFBP3",'8'="Fb-CXCL12",'9'="EC-STMN1",
                            '10'="SMC",'11'="undefined-1",'12'="undefined-2",'13'="undefined-3",'14'="EC - Mu",'15'="later"
                            ,'16'="EC-TFF3",'17'="Fb-MPZ",'18'="epi")
Idents(str)<-'new.cluster.ids'

str=subset(x = str, idents=c("EC-ACKR1","EC-KDR-ESM1","EC-KDR-IGFBP3","EC-STMN1","EC-TFF3"))

DimPlot(str,reduction = "umap",label = TRUE,pt.size = 0.1)



#图2
col <- c("#64a7eb","#fee052","#c9c995","#ee9ca7","#6b7682")
ecc=subset(x = str, idents=c("EC-ACKR1","EC-KDR-ESM1","EC-KDR-IGFBP3","EC-STMN1","EC-TFF3"))
DimPlot(ecc,reduction = "umap",label = TRUE,pt.size = 0.5,cols=col)

pdf(file="2A.pdf",width = 8,height = 8)
DimPlot(ecc,reduction = "umap",label = TRUE,pt.size = 0.5,cols=col)
dev.off()

VEC = str@reductions$umap@cell.embeddings
rownames(VEC) = colnames(str)
PCA = str@reductions$pca@cell.embeddings

unloadNamespace("Seurat")
unloadNamespace("sctransform")
unloadNamespace("devtools")
unloadNamespace("profvis")
unloadNamespace("reshape2")
unloadNamespace("tidyverse")

setwd("D:/R/001/work")
source('vector read.R')

# Remove quantile-based colinearity among PCs (new feature in VECTOR 0.0.3):   
PCA=vector.rankPCA(PCA)


source('vector read.R')

# Define pixel
OUT=vector.buildGrid(VEC, N=30,SHOW=TRUE)

# Build network
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)

# Calculate Quantile Polarization (QP) score
OUT=vector.getValue(OUT, PCA, SHOW=TRUE)

# Get pixel's QP score
OUT=vector.gridValue(OUT,SHOW=TRUE)

# Find starting point
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)

# Infer vector
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL, SHOW.SUMMIT=TRUE)

# OUT$P.PS : Peseudotime Score (PS) of each cell



setwd("D:/R/001/picture")
pdf(file="2B.pdf",width = 8,height = 8)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL, SHOW.SUMMIT=TRUE)
dev.off()