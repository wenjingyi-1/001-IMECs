setwd("D:/R/001/csv")
gene_clinical<-read.csv("GSE19862andGSE19860.csv",header = T)

library(ggpubr)
# 采用ggboxplot绘图
p <- ggboxplot(gene_clinical,x = "respond", y= "EC.KDR.IGFBP3", color="respond", 
               #palette=c("#00AFBB","#E7B800"),
               add="jitter", shape="respond")+
  stat_compare_means(method = "wilcox.test")#t.test

p


p1 <- ggboxplot(gene_clinical,x = "respond", y= "EC.KDR.ESM1", color="respond", 
                #palette=c("#00AFBB","#E7B800"),
                add="jitter", shape="respond")+
  stat_compare_means(method = "wilcox.test")#wilcox.test

p1


p2 <- ggboxplot(gene_clinical,x = "respond", y= "EC.ACKR1", color="respond", 
                #palette=c("#00AFBB","#E7B800"),
                add="jitter", shape="respond")+
  stat_compare_means(method = "wilcox.test")#wilcox.test

p2



setwd("D:/R/001/picture")

pdf(file="10A.pdf")
p1
dev.off()


pdf(file="10B.pdf")
p
dev.off()


pdf(file="10C.pdf")
p2
dev.off()