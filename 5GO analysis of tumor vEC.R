library(Seurat)
setwd("D:/R/001/workspace")
str<-readRDS("str.Rds")
#共有部分
str$new.cluster.ids<-str$seurat_clusters
library(tidyverse)
Idents(str)<-'seurat_clusters'
str$new.cluster.ids<-recode(str$new.cluster.ids,'0'="Mu-MCAM-NDUFA4L2",'1'="EC-ACKR1",'2'="Fb-FAP-MMP11",'3'="EC-KDR-ESM1",'4'="Fb-MMP11"
                            ,'5'="Fb-KLF4",'6'="Fb-FAP",'7'="EC-KDR-IGFBP3",'8'="Fb-CXCL12",'9'="EC-STMN1",
                            '10'="SMC",'11'="immune like",'12'="undefined-2",'13'="undefined-1",'14'="EC - Mu",'15'="EC-TFF3"
                            ,'16'="Mu-NDUFA4L2",'17'="Fb-MPZ",'18'="epithelial like")

Idents(str)<-'new.cluster.ids'

DimPlot(str,reduction = "umap",label = TRUE,pt.size = 0.01)

setwd("D:/R/001/txt")
str.markers<-read.table("str.markers.txt",header = T)


{sig.markers<-str.markers%>%select(gene,everything())%>%
    subset(p_val_adj<0.05&abs(str.markers$avg_log2FC)>1)
  dim(sig.markers)
  
  setwd("D:/R/001/csv/to DAVID")
  
  marker1<-sig.markers[c(sig.markers$cluster=='1'),]
  write.csv(marker1,'marker1.csv',row.names = TRUE)
  
  
  marker3<-sig.markers[c(sig.markers$cluster=='3'),]
  write.csv(marker3,'marker3.csv',row.names = TRUE)
  
  
  marker7<-sig.markers[c(sig.markers$cluster=='7'),]
  write.csv(marker7,'marker7.csv',row.names = TRUE)
  
  
  marker9<-sig.markers[c(sig.markers$cluster=='9'),]
  write.csv(marker9,'marker9.csv',row.names = TRUE)
}#markers转换的处理

#DAVID线上分析（csv拿去分析，分析出一个txt）
setwd("D:/R/001/txt")
marker1<-read.table("marker1.txt",header = T,sep = "\t")

marker3<-read.table("marker3.txt",header = T,sep = "\t")

marker7<-read.table("marker7.txt",header = T,sep = "\t")

marker9<-read.table("marker9.txt",header = T,sep = "\t")

########
setwd("D:/R/001/csv")
write.csv(marker1,'marker1.csv',row.names = F)
write.csv(marker3,'marker3.csv',row.names = F)
write.csv(marker7,'marker7.csv',row.names = F)
write.csv(marker9,'marker9.csv',row.names = F)

if(!require(rvcheck))devtools::install_version("rvcheck",version="0.1.8",repos="http://cran.us.r-project.org")
if(!require(clusterProfiler))BiocManager::install("clusterProfiler")
if(!require(org.Mm.eg.db))BiocManager::install("org.Mm.eg.db")
if(!require(org.Hs.eg.db))BiocManager::install("org.Hs.eg.db")
library(patchwork)

#######基因注释


##

######################
{#1
  ggoMF<-enrichGO(marker1$To,org.Hs.eg.db,keyType = "ENTREZID",ont = "MF",
                  pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.2,minGSSize = 10,
                  maxGSSize = 500,readable = FALSE,pool = FALSE)
  
  ggoMF <- simplify(ggoMF, cutoff=0.7, by="p.adjust", select_fun=min)
  
  
  ggoCC<-enrichGO(marker1$To,org.Hs.eg.db,keyType = "ENTREZID",ont = "CC",
                  pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.2,minGSSize = 10,
                  maxGSSize = 500,readable = FALSE,pool = FALSE)
  
  ggoCC <- simplify(ggoCC, cutoff=0.7, by="p.adjust", select_fun=min)
  
  
  ggoBP<-enrichGO(marker1$To,org.Hs.eg.db,keyType = "ENTREZID",ont = "BP",
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
           x="Richfactor",y="term.description",title = "EC-ACKR1 Enrichment Process")
    pr<-pr + theme_bw()
    pr
  }
  
  erichACKR1p1plot<-function(data4plot){
    library(ggplot2)
    data4plot<-data4plot[order(data4plot$qvalue,decreasing = F)[1:17],]
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
           x="Richfactor",y="term.description",title = "EC-ACKR1 Enrichment Process")
    pr<-pr + theme_bw()
    pr
  }
  
  p1<-erichACKR1p1plot(ggoMF@result)
  p2<-erich2plot(ggoCC@result)
  p3<-erich2plot(ggoBP@result)
  
  p1
  p2
  p3
  
  
  setwd("D:/R/001/picture")
  pdf(file = "4A.pdf", height = 7, width = 25)
  p1+p2+p3
  dev.off()
}#基因富集分析EC-ACKR1













######################
{#3
  
  ggoMF<-enrichGO(marker3$To,org.Hs.eg.db,keyType = "ENTREZID",ont = "MF",
                  pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.2,minGSSize = 10,
                  maxGSSize = 500,readable = FALSE,pool = FALSE)
  
  ggoMF <- simplify(ggoMF, cutoff=0.7, by="p.adjust", select_fun=min)
  
  ggoCC<-enrichGO(marker3$To,org.Hs.eg.db,keyType = "ENTREZID",ont = "CC",
                  pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.2,minGSSize = 10,
                  maxGSSize = 500,readable = FALSE,pool = FALSE)
  
  ggoCC <- simplify(ggoCC, cutoff=0.7, by="p.adjust", select_fun=min)
  
  ggoBP<-enrichGO(marker3$To,org.Hs.eg.db,keyType = "ENTREZID",ont = "BP",
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
           x="Richfactor",y="term.description",title = "EC-KDR-ESM1 Enrichment Process")
    pr<-pr + theme_bw()
    pr
  }
  
  erich2longplot<-function(data4plot){
    library(ggplot2)
    data4plot<-data4plot[order(data4plot$qvalue,decreasing = F)[1:16],]
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
           x="Richfactor",y="term.description",title = "EC-KDR-ESM1 Enrichment Process")
    pr<-pr + theme_bw()
    pr
  }
  
  p1<-erich2longplot(ggoMF@result)
  p2<-erich2longplot(ggoCC@result)
  p3<-erich2plot(ggoBP@result)
  
  
  p1
  p2
  p3
  
  

  pdf(file = "S4A.pdf", height = 7, width = 25)
  p1+p2+p3
  dev.off()
}#基因富集分析EC-KDR-ESM1





######################
{#7
  #
  ggoMF<-enrichGO(marker7$To,org.Hs.eg.db,keyType = "ENTREZID",ont = "MF",
                  pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.2,minGSSize = 10,
                  maxGSSize = 500,readable = FALSE,pool = FALSE)
  
  ggoMF <- simplify(ggoMF, cutoff=0.7, by="p.adjust", select_fun=min)
  
  
  ggoCC<-enrichGO(marker7$To,org.Hs.eg.db,keyType = "ENTREZID",ont = "CC",
                  pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.2,minGSSize = 10,
                  maxGSSize = 500,readable = FALSE,pool = FALSE)
  
  ggoCC <- simplify(ggoCC, cutoff=0.7, by="p.adjust", select_fun=min)
  
  
  ggoBP<-enrichGO(marker7$To,org.Hs.eg.db,keyType = "ENTREZID",ont = "BP",
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
           x="Richfactor",y="term.description",title = "EC-KDR-IGFBP3 Enrichment Process")
    pr<-pr + theme_bw()
    pr
  }
  
  erichIGFBP3p1plot<-function(data4plot){
    library(ggplot2)
    data4plot<-data4plot[order(data4plot$qvalue,decreasing = F)[1:13],]
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
           x="Richfactor",y="term.description",title = "EC-KDR-IGFBP3 Enrichment Process")
    pr<-pr + theme_bw()
    pr
  }
  
  p1<-erichIGFBP3p1plot(ggoMF@result)
  p2<-erich2plot(ggoCC@result)
  p3<-erich2plot(ggoBP@result)
  
  p1
  p2
  p3
  
  
  pdf(file = "4B.pdf", height = 7, width = 25)
  p1+p2+p3
  dev.off()
  
}#基因富集分析EC-KDR-IGFBP3



######################
{#9
  
  ggoMF<-enrichGO(marker9$To,org.Hs.eg.db,keyType = "ENTREZID",ont = "MF",
                  pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.2,minGSSize = 10,
                  maxGSSize = 500,readable = FALSE,pool = FALSE)
  
  ggoMF <- simplify(ggoMF, cutoff=0.7, by="p.adjust", select_fun=min)
  
  
  ggoCC<-enrichGO(marker9$To,org.Hs.eg.db,keyType = "ENTREZID",ont = "CC",
                  pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.2,minGSSize = 10,
                  maxGSSize = 500,readable = FALSE,pool = FALSE)
  
  ggoCC <- simplify(ggoCC, cutoff=0.7, by="p.adjust", select_fun=min)
  
  
  ggoBP<-enrichGO(marker9$To,org.Hs.eg.db,keyType = "ENTREZID",ont = "BP",
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
           x="Richfactor",y="term.description",title = "EC-STMN1 Enrichment Process")
    pr<-pr + theme_bw()
    pr
  }
  erichSTMN1p2plot<-function(data4plot){
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
           x="Richfactor",y="term.description",title = "EC-STMN1 Enrichment Process")
    pr<-pr + theme_bw()
    pr
  }
  
  p1<-erich2plot(ggoMF@result)
  p2<-erichSTMN1p2plot(ggoCC@result)
  p3<-erich2plot(ggoBP@result)
  # enrichplot::cnetplot(ggoMF,circular=TRUE,colorEdge=TRUE)
  
  p1
  p2
  p3
  

  pdf(file = "S4B.pdf", height = 7, width = 25)
  p1+p2+p3
  dev.off()
}#基因富集分析EC-STMN1