
#load the required packages

library(tibble)
library(dplyr)
library(DESeq2)
library(EnhancedVolcano)
library(ggplot2)

#set your working directory

setwd('/Users/davidmunglah/Desktop/')

# load the count and metadata 

count <- t(read.csv("/Users/davidmunglah/Desktop/EC/",row.names = 1,check.names = F))

metadata<- read.csv("/Users/davidmunglah/Desktop/EC/")

#take the control as reference

metadata$Treatments <- relevel(factor(metadata$Treatment), ref = "Citristim")

#filter the samples that are only present within the metdata into the count data 

count<-count[,colnames(count)%in%metadata$Sample.ID]

#filter the samples which are present in the count data into the metadata

metadata<-metadata[metadata$Sample.ID%in%colnames(count),]

#check if the order of the samples match in the count and metadata

all(metadata$Sample.ID  == colnames(count))

#order the dataframe if the previous line is false

reorder_idx <- match(metadata$Sample.ID,colnames(count))

count <- count[,reorder_idx]

EC_count<-t(count)

#Remove N/A

EC_count[is.na(EC_count)]<- 0


#remove the zero column 

columns_sum<-colSums(EC_count[,1:ncol(EC_count)])==0
EC_count <- EC_count[,!columns_sum]
count<-t(EC_count)



#Perform DESeq2

count_data <- count+1

dds_count <- DESeqDataSetFromMatrix(countData=count, 
                                          colData=metadata, 
                                          design=~Treatments, tidy = FALSE)


dds_count <- DESeq(dds_count)
res_count <- results(dds_count)
FC_count <- (results(dds_count , tidy=TRUE))


#plot 


pdf('Neg Control_Citristim_21_EC.pdf',width=20,height=10) #pdf

EnhancedVolcano(res_count,lab = rownames(res_count),x = 'log2FoldChange',y = 'pvalue',
                title = 'Neg Control_Citristim_21_EC',pCutoff = 0.05,FCcutoff = 2,
                pointSize = 3.0,labSize = 5.0) 


dev.off()

#significance table output 

write.csv(FC_count,"Neg Control_Citristim_21_EC.csv", row.names = FALSE)

