library(Seurat)
setwd("D:/R/001/workspace")
#差异基因分析
######EC-ACKR1-normal
str<-readRDS("normalstr.Rds")

str1<-str

Idents(str1)<-'seurat_clusters'

normal=subset(x = str1, idents=c("6"))

{normal$new.cluster.ids<-normal$seurat_clusters
  library(tidyverse)
  Idents(normal)<-'seurat_clusters'
  normal$new.cluster.ids<-recode(normal$new.cluster.ids,'6'="EC-ACKR1-normal")
  
  Idents(normal)<-'new.cluster.ids'
  
  DimPlot(normal,reduction = "umap",label = TRUE,pt.size = 0.01)
}#细胞分组




#打开肿瘤组的RDS
######EC-ACKR1-tumor
str<-readRDS("str.Rds")

Idents(str)<-'seurat_clusters'

tumor=subset(x = str, idents=c("1"))

{tumor$new.cluster.ids<-tumor$seurat_clusters
  library(tidyverse)
  Idents(tumor)<-'seurat_clusters'
  tumor$new.cluster.ids<-recode(tumor$new.cluster.ids,'1'="EC-ACKR1-tumor")
  
  Idents(tumor)<-'new.cluster.ids'
  
  DimPlot(tumor,reduction = "umap",label = TRUE,pt.size = 0.01)
}#细胞分组


##合并
pbmc.combined <- merge(tumor, y = normal, add.cell.ids = c("tumor","normal" ), project = "PBMC12K",merge.data = TRUE)


diff_hot1 <- FindMarkers(pbmc.combined, min.pct = 0.25, 
                         logfc.threshold = 0.25,
                         group.by = "new.cluster.ids",
                         ident.1 ="EC-ACKR1-normal",
                         ident.2="EC-ACKR1-tumor")  

VlnPlot(pbmc.combined,features = c("CLU"))

#BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)
p1<-EnhancedVolcano(diff_hot1,
                    lab = rownames(diff_hot1),
                    x = 'avg_log2FC',
                    y = 'p_val_adj',
                    pCutoff = 0.05,
                    FCcutoff = 0.25,
                    pointSize = 3.0,
                    labSize = 6.0,
                    title = 'EC-ACKR1',
                    selectLab = c("HLA-DRA","HLA-DRB1","HLA-DPA1","HLA-DPB1","HLA-DQB1","HLA-DMA" ),
                    drawConnectors = TRUE,
                    widthConnectors = 0.75)

p1

setwd("D:/R/001/picture")
pdf(file="6A.pdf",width = 15,height = 15)
p1
dev.off()
########


if(!require(rvcheck))devtools::install_version("rvcheck",version="0.1.8",repos="http://cran.us.r-project.org")
if(!require(clusterProfiler))BiocManager::install("clusterProfiler")
if(!require(org.Mm.eg.db))BiocManager::install("org.Mm.eg.db")
if(!require(org.Hs.eg.db))BiocManager::install("org.Hs.eg.db")
library(patchwork)

#######基因注释


##

######################
setwd("D:/R/001/csv/to DAVID")
typhon<-diff_hot1
write.csv(typhon,'all-EC-ACKR1.csv',row.names = TRUE)
###拆分出上调和下调
###DAVID线上分析

setwd("D:/R/001/txt")
typhon<-read.table("EC-ACKR1-tumorup.txt",header = T,sep = "\t")


setwd("D:/R/001/csv")
write.csv(typhon,'EC-ACKR1-tumorup.csv',row.names = TRUE)
#基因富集分析EC-ACKR1

######################
{#1
  ggoMF<-enrichGO(typhon$To,org.Hs.eg.db,keyType = "ENTREZID",ont = "MF",
                  pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.2,minGSSize = 10,
                  maxGSSize = 500,readable = FALSE,pool = FALSE)
  
  ggoMF <- simplify(ggoMF, cutoff=0.7, by="p.adjust", select_fun=min)
  
  
  ggoCC<-enrichGO(typhon$To,org.Hs.eg.db,keyType = "ENTREZID",ont = "CC",
                  pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.2,minGSSize = 10,
                  maxGSSize = 500,readable = FALSE,pool = FALSE)
  
  ggoCC <- simplify(ggoCC, cutoff=0.7, by="p.adjust", select_fun=min)
  
  
  ggoBP<-enrichGO(typhon$To,org.Hs.eg.db,keyType = "ENTREZID",ont = "BP",
                  pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.2,minGSSize = 10,
                  maxGSSize = 500,readable = FALSE,pool = FALSE)
  
  ggoBP <- simplify(ggoBP, cutoff=0.7, by="p.adjust", select_fun=min)
  
  
  erich2plot<-function(data4plot){
    library(ggplot2)
    data4plot<-data4plot[order(data4plot$qvalue,decreasing = F)[1:20],]
    data4plot$BgRatio<-
      apply(data4plot,1,function(x){
        as.numeric(strsplit(x[3],'/')[[1]][1])
      })/apply(data4plot,1,function(x){
        as.numeric(strsplit(x[4],'/')[[1]][1])
      })
    
    p<-ggplot(data4plot,aes(BgRatio,Description))
    p<-p+geom_point()
    
    pbubble<-p+geom_point(aes(size=Count,color=-1*log10(qvalue)))
    
    pr<-pbubble+scale_colour_gradient(low="#90EE90",high = "red")+
      labs(color=expression(-log[10](qvalue)),size="observed.gene.count",
           x="Richfactor",y="term.description")
    pr<-pr + theme_bw()
    pr
  }
  
  
  
  p1<-erich2plot(ggoMF@result)
  p2<-erich2plot(ggoCC@result)
  p3<-erich2plot(ggoBP@result)
  
  
  p1
  p2
  p3
  
   setwd("D:/R/001/picture")
  
  pdf(file = "S6A.pdf", height = 7, width = 25)
  p1+p2+p3
  dev.off()
}#基因富集分析EC-ACKR1上调











setwd("D:/R/001/txt")
typhon<-read.table("EC-ACKR1-tumordown.txt",header = T,sep = "\t")


setwd("D:/R/001/csv")
write.csv(typhon,'EC-ACKR1-tumordown.csv',row.names = TRUE)

#基因富集分析EC-ACKR1

######################
{#1
  ggoMF<-enrichGO(typhon$To,org.Hs.eg.db,keyType = "ENTREZID",ont = "MF",
                  pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.2,minGSSize = 10,
                  maxGSSize = 500,readable = FALSE,pool = FALSE)
  
  ggoMF <- simplify(ggoMF, cutoff=0.7, by="p.adjust", select_fun=min)
  
  
  ggoCC<-enrichGO(typhon$To,org.Hs.eg.db,keyType = "ENTREZID",ont = "CC",
                  pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.2,minGSSize = 10,
                  maxGSSize = 500,readable = FALSE,pool = FALSE)
  
  ggoCC <- simplify(ggoCC, cutoff=0.7, by="p.adjust", select_fun=min)
  
  
  ggoBP<-enrichGO(typhon$To,org.Hs.eg.db,keyType = "ENTREZID",ont = "BP",
                  pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.2,minGSSize = 10,
                  maxGSSize = 500,readable = FALSE,pool = FALSE)
  
  ggoBP <- simplify(ggoBP, cutoff=0.7, by="p.adjust", select_fun=min)
  
  
  erich2plot<-function(data4plot){
    library(ggplot2)
    data4plot<-data4plot[order(data4plot$qvalue,decreasing = F)[1:20],]
    data4plot$BgRatio<-
      apply(data4plot,1,function(x){
        as.numeric(strsplit(x[3],'/')[[1]][1])
      })/apply(data4plot,1,function(x){
        as.numeric(strsplit(x[4],'/')[[1]][1])
      })
    
    p<-ggplot(data4plot,aes(BgRatio,Description))
    p<-p+geom_point()
    
    pbubble<-p+geom_point(aes(size=Count,color=-1*log10(qvalue)))
    
    pr<-pbubble+scale_colour_gradient(low="#90EE90",high = "red")+
      labs(color=expression(-log[10](qvalue)),size="observed.gene.count",
           x="Richfactor",y="term.description")
    pr<-pr + theme_bw()
    pr
  }
  
  erichACKR1p1plot<-function(data4plot){
    library(ggplot2)
    data4plot<-data4plot[order(data4plot$qvalue,decreasing = F)[1:15],]
    data4plot$BgRatio<-
      apply(data4plot,1,function(x){
        as.numeric(strsplit(x[3],'/')[[1]][1])
      })/apply(data4plot,1,function(x){
        as.numeric(strsplit(x[4],'/')[[1]][1])
      })
    
    p<-ggplot(data4plot,aes(BgRatio,Description))
    p<-p+geom_point()
    
    pbubble<-p+geom_point(aes(size=Count,color=-1*log10(qvalue)))
    
    pr<-pbubble+scale_colour_gradient(low="#90EE90",high = "red")+
      labs(color=expression(-log[10](qvalue)),size="observed.gene.count",
           x="Richfactor",y="term.description")
    pr<-pr + theme_bw()
    pr
  }
  
  erichACKR1p2plot<-function(data4plot){
    library(ggplot2)
    data4plot<-data4plot[order(data4plot$qvalue,decreasing = F)[1:19],]
    data4plot$BgRatio<-
      apply(data4plot,1,function(x){
        as.numeric(strsplit(x[3],'/')[[1]][1])
      })/apply(data4plot,1,function(x){
        as.numeric(strsplit(x[4],'/')[[1]][1])
      })
    
    p<-ggplot(data4plot,aes(BgRatio,Description))
    p<-p+geom_point()
    
    pbubble<-p+geom_point(aes(size=Count,color=-1*log10(qvalue)))
    
    pr<-pbubble+scale_colour_gradient(low="#90EE90",high = "red")+
      labs(color=expression(-log[10](qvalue)),size="observed.gene.count",
           x="Richfactor",y="term.description")
    pr<-pr + theme_bw()
    pr
  }
  
  p1<-erichACKR1p1plot(ggoMF@result)
  p2<-erichACKR1p2plot(ggoCC@result)
  p3<-erich2plot(ggoBP@result)
  
  p1
  p2
  p3
  
  
  setwd("D:/R/001/picture")
  
  pdf(file = "S6B.pdf", height = 7, width = 25)
  p1+p2+p3
  dev.off()
}#基因富集分析EC-ACKR1下调




setwd("D:/R/001/workspace")
#差异基因分析
######EC-IGFBP3-normal
str<-readRDS("normalstr.Rds")

Idents(str1)<-'seurat_clusters'

normal=subset(x = str1, idents=c("8"))

{normal$new.cluster.ids<-normal$seurat_clusters
  library(tidyverse)
  Idents(normal)<-'seurat_clusters'
  normal$new.cluster.ids<-recode(normal$new.cluster.ids,'8'="EC-IGFBP3-normal")
  
  Idents(normal)<-'new.cluster.ids'
  
  DimPlot(normal,reduction = "umap",label = TRUE,pt.size = 0.01)
}#细胞分组




#打开肿瘤组
######EC-IGFBP3-tumor
str<-readRDS("str.Rds")

Idents(str)<-'seurat_clusters'

tumor=subset(x = str, idents=c("7"))

{tumor$new.cluster.ids<-tumor$seurat_clusters
  library(tidyverse)
  Idents(tumor)<-'seurat_clusters'
  tumor$new.cluster.ids<-recode(tumor$new.cluster.ids,'7'="EC-KDR-IGFBP3-tumor")
  
  Idents(tumor)<-'new.cluster.ids'
  
  DimPlot(tumor,reduction = "umap",label = TRUE,pt.size = 0.01)
}#细胞分组


##合并
pbmc.combined <- merge(tumor, y = normal, add.cell.ids = c("tumor","normal" ), project = "PBMC12K",merge.data = TRUE)


diff_hot1 <- FindMarkers(pbmc.combined, min.pct = 0.25, 
                         logfc.threshold = 0.25,
                         group.by = "new.cluster.ids",
                         ident.1 ="EC-IGFBP3-normal",
                         ident.2="EC-KDR-IGFBP3-tumor")  

VlnPlot(pbmc.combined,features = c("CLU"))

#BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)
p1<-EnhancedVolcano(diff_hot1,
                    lab = rownames(diff_hot1),
                    x = 'avg_log2FC',
                    y = 'p_val_adj',
                    pCutoff = 0.05,
                    FCcutoff = 0.25,
                    pointSize = 3.0,
                    labSize = 6.0,
                    title = 'EC-KDR-IGFBP3',
                    selectLab = c("HLA-DRA","HLA-DRB1","HLA-DPA1","HLA-DPB1","HLA-DQB1","HLA-DMA" ),
                    drawConnectors = TRUE,
                    widthConnectors = 0.75)

p1

setwd("D:/R/001/picture")
pdf(file="6B.pdf",width = 15,height = 15)
p1
dev.off()


########


if(!require(rvcheck))devtools::install_version("rvcheck",version="0.1.8",repos="http://cran.us.r-project.org")
if(!require(clusterProfiler))BiocManager::install("clusterProfiler")
if(!require(org.Mm.eg.db))BiocManager::install("org.Mm.eg.db")
if(!require(org.Hs.eg.db))BiocManager::install("org.Hs.eg.db")
library(patchwork)

#######基因注释

setwd("D:/R/001/csv/to DAVID")
typhon<-diff_hot1
write.csv(typhon,'all-EC-IGFBP3.csv',row.names = TRUE)
###拆分出上调和下调
###DAVID线上分析

setwd("D:/R/001/txt")
typhon<-read.table("EC-IGFBP3-tumordown.txt",header = T,sep = "\t")


setwd("D:/R/001/csv")
write.csv(typhon,'EC-IGFBP3-tumordown.csv',row.names = TRUE)

##




######################
{#1
  ggoMF<-enrichGO(typhon$To,org.Hs.eg.db,keyType = "ENTREZID",ont = "MF",
                  pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.2,minGSSize = 10,
                  maxGSSize = 500,readable = FALSE,pool = FALSE)
  
  ggoMF <- simplify(ggoMF, cutoff=0.7, by="p.adjust", select_fun=min)
  
  
  ggoCC<-enrichGO(typhon$To,org.Hs.eg.db,keyType = "ENTREZID",ont = "CC",
                  pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.2,minGSSize = 10,
                  maxGSSize = 500,readable = FALSE,pool = FALSE)
  
  ggoCC <- simplify(ggoCC, cutoff=0.7, by="p.adjust", select_fun=min)
  
  
  ggoBP<-enrichGO(typhon$To,org.Hs.eg.db,keyType = "ENTREZID",ont = "BP",
                  pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.2,minGSSize = 10,
                  maxGSSize = 500,readable = FALSE,pool = FALSE)
  
  ggoBP <- simplify(ggoBP, cutoff=0.7, by="p.adjust", select_fun=min)
  
  
  erich2plot<-function(data4plot){
    library(ggplot2)
    data4plot<-data4plot[order(data4plot$qvalue,decreasing = F)[1:20],]
    data4plot$BgRatio<-
      apply(data4plot,1,function(x){
        as.numeric(strsplit(x[3],'/')[[1]][1])
      })/apply(data4plot,1,function(x){
        as.numeric(strsplit(x[4],'/')[[1]][1])
      })
    
    p<-ggplot(data4plot,aes(BgRatio,Description))
    p<-p+geom_point()
    
    pbubble<-p+geom_point(aes(size=Count,color=-1*log10(qvalue)))
    
    pr<-pbubble+scale_colour_gradient(low="#90EE90",high = "red")+
      labs(color=expression(-log[10](qvalue)),size="observed.gene.count",
           x="Richfactor",y="term.description")
    pr<-pr + theme_bw()
    pr
  }
  
  erichIGFBP3p1plot<-function(data4plot){
    library(ggplot2)
    data4plot<-data4plot[order(data4plot$qvalue,decreasing = F)[1:11],]
    data4plot$BgRatio<-
      apply(data4plot,1,function(x){
        as.numeric(strsplit(x[3],'/')[[1]][1])
      })/apply(data4plot,1,function(x){
        as.numeric(strsplit(x[4],'/')[[1]][1])
      })
    
    p<-ggplot(data4plot,aes(BgRatio,Description))
    p<-p+geom_point()
    
    pbubble<-p+geom_point(aes(size=Count,color=-1*log10(qvalue)))
    
    pr<-pbubble+scale_colour_gradient(low="#90EE90",high = "red")+
      labs(color=expression(-log[10](qvalue)),size="observed.gene.count",
           x="Richfactor",y="term.description")
    pr<-pr + theme_bw()
    pr
  }
  
  erichIGFBP3p2plot<-function(data4plot){
    library(ggplot2)
    data4plot<-data4plot[order(data4plot$qvalue,decreasing = F)[1:15],]
    data4plot$BgRatio<-
      apply(data4plot,1,function(x){
        as.numeric(strsplit(x[3],'/')[[1]][1])
      })/apply(data4plot,1,function(x){
        as.numeric(strsplit(x[4],'/')[[1]][1])
      })
    
    p<-ggplot(data4plot,aes(BgRatio,Description))
    p<-p+geom_point()
    
    pbubble<-p+geom_point(aes(size=Count,color=-1*log10(qvalue)))
    
    pr<-pbubble+scale_colour_gradient(low="#90EE90",high = "red")+
      labs(color=expression(-log[10](qvalue)),size="observed.gene.count",
           x="Richfactor",y="term.description")
    pr<-pr + theme_bw()
    pr
  }
  
  p1<-erichIGFBP3p1plot(ggoMF@result)
  p2<-erichIGFBP3p2plot(ggoCC@result)
  p3<-erich2plot(ggoBP@result)

  
  p1
  p2
  p3
  
  setwd("D:/R/001/picture")
  
  pdf(file = "S7B.pdf", height = 7, width = 25)
  p1+p2+p3
  dev.off()
  
}#基因富集分析EC-IGFBP3下调









setwd("D:/R/001/txt")
typhon<-read.table("EC-IGFBP3-tumorup.txt",header = T,sep = "\t")


setwd("D:/R/001/csv")
write.csv(typhon,'EC-IGFBP3-tumorup.csv',row.names = TRUE)


######################
{#1
  ggoMF<-enrichGO(typhon$To,org.Hs.eg.db,keyType = "ENTREZID",ont = "MF",
                  pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.2,minGSSize = 10,
                  maxGSSize = 500,readable = FALSE,pool = FALSE)
  
  ggoMF <- simplify(ggoMF, cutoff=0.7, by="p.adjust", select_fun=min)
  
  
  ggoCC<-enrichGO(typhon$To,org.Hs.eg.db,keyType = "ENTREZID",ont = "CC",
                  pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.2,minGSSize = 10,
                  maxGSSize = 500,readable = FALSE,pool = FALSE)
  
  ggoCC <- simplify(ggoCC, cutoff=0.7, by="p.adjust", select_fun=min)
  
  
  ggoBP<-enrichGO(typhon$To,org.Hs.eg.db,keyType = "ENTREZID",ont = "BP",
                  pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.2,minGSSize = 10,
                  maxGSSize = 500,readable = FALSE,pool = FALSE)
  
  ggoBP <- simplify(ggoBP, cutoff=0.7, by="p.adjust", select_fun=min)
  
  
  erich2plot<-function(data4plot){
    library(ggplot2)
    data4plot<-data4plot[order(data4plot$qvalue,decreasing = F)[1:20],]
    data4plot$BgRatio<-
      apply(data4plot,1,function(x){
        as.numeric(strsplit(x[3],'/')[[1]][1])
      })/apply(data4plot,1,function(x){
        as.numeric(strsplit(x[4],'/')[[1]][1])
      })
    
    p<-ggplot(data4plot,aes(BgRatio,Description))
    p<-p+geom_point()
    
    pbubble<-p+geom_point(aes(size=Count,color=-1*log10(qvalue)))
    
    pr<-pbubble+scale_colour_gradient(low="#90EE90",high = "red")+
      labs(color=expression(-log[10](qvalue)),size="observed.gene.count",
           x="Richfactor",y="term.description")
    pr<-pr + theme_bw()
    pr
  }
  
  
  
  p1<-erich2plot(ggoMF@result)
  p2<-erich2plot(ggoCC@result)
  p3<-erich2plot(ggoBP@result)
  
  p1
  p2
  p3
 

  setwd("D:/R/001/picture")
  
  pdf(file = "S7A.pdf", height = 7, width = 25)
  p1+p2+p3
  dev.off()
}#基因富集分析EC-IGFBP3上调