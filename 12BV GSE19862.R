#GSE19862
setwd("D:/R/001/data/19862")
a<-read.table(file = 'GSE19862_series_matrix.txt.gz',sep = '\t',
              comment.char = "!",header = T)

#这里的改指删去表格!Sample_title以上的注释文件（只有一列的内容）
b<-read.table(file = 'GSE19862_series_matrix(改).txt',sep = '\t',
              fill = T)
###基因symbol对应
{#BiocManager::install("hgu133plus2.db")
  library(hgu133plus2.db)
  ids=toTable(hgu133plus2SYMBOL)
  head(ids)
  
  a2<-a
  a2$probe_id<-a$ID_REF
  a2<-merge(ids,a2,by='probe_id')
  
  A<-a2[,-c(1,3)]
  A1<-A
  
  matrix1 <- aggregate(.~symbol,data=A1,mean) 
  A1<-matrix1
  
  rownames(A1)<-A1$symbol
  A1<-A1[,-c(1)]
  
  
}
library(tidyverse)
A1<-rownames_to_column(A1)

setwd("D:/R/001/txt")
write.table( A1,file = "GSE19862.txt",row.names = F,sep = "\t")
#拿这个表去cibersortx推断


A1<-read.table("GSE19862CIBERSORTx_Job40_Results.txt",header = T,sep = "\t")


###临床信息对应
c<-b[c(1:34),]
f<-b[34,]
colnames(c)<-f

c<-t(c)
c<-as.data.frame(c)
c1<-c[-1,]
c1$Mixture<-rownames(c1)



###基因对临床
identical(A1$Mixture,c1$Mixture)

#去冗余
#rm("a","b","f","a1","GPL570","A","a2","c","matrix1")



library(glmnet)
library(readxl)
library(plyr)
library(corrplot)
library(ggplot2)
library(Hmisc)
library(openxlsx)
library(survival)
library(survminer)

library(tidyverse)


data8<-merge(A1,c1,by="Mixture")




#mark
CTD<-c1
CTD$respond<-ifelse(CTD$`11`=="treatment response: BV_Non-responder",'BV_Non-responder','BV_Responder')

gene_clinical<-CTD


gene_clinical<-merge(gene_clinical,A1,by="Mixture")




setwd("D:/R/001/csv")
write.csv(gene_clinical,'GSE19862Clinical.csv',row.names = TRUE)
