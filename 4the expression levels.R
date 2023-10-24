library(Seurat)
setwd("D:/R/001/workspace")
str<-readRDS("str.Rds")

str$new.cluster.ids<-str$seurat_clusters
library(tidyverse)
Idents(str)<-'seurat_clusters'
str$new.cluster.ids<-recode(str$new.cluster.ids,'0'="Mu-MCAM-NDUFA4L2",'1'="EC-ACKR1",'2'="Fb-FAP-MMP11",'3'="EC-KDR-ESM1",'4'="Fb-MMP11"
                            ,'5'="Fb-KLF4",'6'="Fb-FAP",'7'="EC-KDR-IGFBP3",'8'="Fb-CXCL12",'9'="EC-STMN1",
                            '10'="SMC",'11'="immune like",'12'="undefined-2",'13'="undefined-1",'14'="EC - Mu",'15'="EC-TFF3"
                            ,'16'="Mu-NDUFA4L2",'17'="Fb-MPZ",'18'="epithelial like")

Idents(str)<-'new.cluster.ids'
str=subset(x = str, idents=c("EC-ACKR1","EC-KDR-ESM1","EC-KDR-IGFBP3","EC-STMN1","EC-TFF3"))

DimPlot(str,reduction = "umap",label = TRUE,pt.size = 0.1)


#基因
setwd("D:/R/001/picture")
color <- c('#000000','#400f70','#b23676','#b23676','#f88f6c','#feedb1')#设置颜色  

pdf(file="3A.pdf")
FeaturePlot(str, features = 'HLA-DRA',cols = color, pt.size = 1.2)
dev.off()


pdf(file="3B.pdf")
 FeaturePlot(str, features = 'HLA-DRB1',cols = color, pt.size = 1.2)
 dev.off()

 pdf(file="3C.pdf")
 FeaturePlot(str, features = 'HLA-DPA1',cols = color, pt.size = 1.2)
 dev.off()

 pdf(file="3D.pdf")
 FeaturePlot(str, features = 'HLA-DPB1',cols = color, pt.size = 1.2)
 dev.off()


 pdf(file="3E.pdf")
 FeaturePlot(str, features = 'HLA-DQB1',cols = color, pt.size = 1.2)
 dev.off()

 pdf(file="3F.pdf")
 FeaturePlot(str, features = 'HLA-DMA',cols = color, pt.size = 1.2)
 dev.off()



#######共刺激
 pdf(file="3G.pdf")
DotPlot(object=str,features=c("CD80","CD86","ICAM1","ICAM2"),group.by="new.cluster.ids")+
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))
dev.off()

#######################VEGF受体靶点
#药物能够覆盖主要的血管生成亚群，但部分增殖亚群、EC-ACKR1亚群无法充分覆盖

TIKcol<-c('#000000','#400f70','#b23676','#b23676','#f88f6c')
TIKS<-c("FGFR1","FGFR2","FLT1","KDR","FLT4","PDGFRA","PDGFRB","KIT","MET")
TIKS<-as.data.frame(TIKS)
VlnPlot(str, features = TIKS$TIKS,
        cols = TIKcol,  split.by =  "new.cluster.ids",
        stack = T,sort = F,flip = F,pt.size=0)+
  NoLegend()

setwd("D:/R/001/picture")
pdf(file="7A.pdf",width = 7,height = 5)
TIKcol<-c('#000000','#400f70','#b23676','#b23676','#f88f6c')
TIKS<-c("FGFR1","FGFR2","FLT1","KDR","FLT4","PDGFRA","PDGFRB","KIT","MET")
TIKS<-as.data.frame(TIKS)
VlnPlot(str, features = TIKS$TIKS,
        cols = TIKcol,  split.by =  "new.cluster.ids",
        stack = T,sort = F,flip = F,pt.size=0)+
  NoLegend()
dev.off()



pdf(file="7B.pdf",width = 7,height = 5)
TIKcol<-c('#000000','#400f70','#af2934','#b23676','#f88f6c')
TIKS<-c("MAPK1","KRAS","SHC1")
TIKS<-as.data.frame(TIKS)
VlnPlot(str, features = TIKS$TIKS,
        cols = TIKcol,  split.by =  "new.cluster.ids",
        stack = T,sort = F,flip = T,pt.size=0)+ 
  NoLegend()
dev.off()