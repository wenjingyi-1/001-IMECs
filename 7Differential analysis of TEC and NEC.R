setwd("D:/R/001/workspace")
library(Seurat)
#差异基因分析
######EC-ACKR1-normal
str<-readRDS("normalstr.Rds")

str1<-str

Idents(str1)<-'seurat_clusters'

normal=subset(x = str1, idents=c("6","8"))

{normal$new.cluster.ids<-normal$seurat_clusters
  library(tidyverse)
  Idents(normal)<-'seurat_clusters'
  normal$new.cluster.ids<-recode(normal$new.cluster.ids,'6'="EC-ACKR1-normal",'8'="EC-IGFBP3-normal")
  
  Idents(normal)<-'new.cluster.ids'
  
  DimPlot(normal,reduction = "umap",label = TRUE,pt.size = 0.01)
}#细胞分组




#打开肿瘤组的RDS
######EC-ACKR1-tumor
str<-readRDS("str.Rds")

Idents(str)<-'seurat_clusters'

tumor=subset(x = str, idents=c("1","7"))

{tumor$new.cluster.ids<-tumor$seurat_clusters
  library(tidyverse)
  Idents(tumor)<-'seurat_clusters'
  tumor$new.cluster.ids<-recode(tumor$new.cluster.ids,'1'="EC-ACKR1-tumor",'7'="EC-IGFBP3-tumor")
  
  Idents(tumor)<-'new.cluster.ids'
  
  DimPlot(tumor,reduction = "umap",label = TRUE,pt.size = 0.01)
}#细胞分组


##合并
pbmc.combined <- merge(tumor, y = normal, add.cell.ids = c("tumor","normal" ), project = "PBMC12K",merge.data = TRUE)






##
VlnPlot(pbmc.combined, features = "HLA-DRA",pt.size=0)+ 
  #width控制箱体宽度，col控制边框颜色，fill控制填充颜色  
  geom_boxplot(width=.2,col="black",fill="white")+  
  NoLegend()


##
VlnPlot(pbmc.combined, features = "HLA-DRB1",pt.size=0)+ 
  #width控制箱体宽度，col控制边框颜色，fill控制填充颜色  
  geom_boxplot(width=.2,col="black",fill="white")+  
  NoLegend()

########
setwd("D:/R/001/picture")
pdf(file="6C.pdf")
VlnPlot(pbmc.combined, features = "HLA-DPA1",pt.size=0)+ 
  #width控制箱体宽度，col控制边框颜色，fill控制填充颜色  
  geom_boxplot(width=.2,col="black",fill="white")+  
  NoLegend()
dev.off()
##
VlnPlot(pbmc.combined, features = "HLA-DPB1",pt.size=0)+ 
  #width控制箱体宽度，col控制边框颜色，fill控制填充颜色  
  geom_boxplot(width=.2,col="black",fill="white")+  
  NoLegend()






