library(Seurat)
library(future)
plan("multiprocess", workers = 8)
options(future.globals.maxSize = 10*1024^3)

dir <- "./"

## take Neu_CCL3
Neutro <- readRDS("Neutrophil.rds")
G1 <- subset(Neutro,idents=1)
G1_count <- as.matrix(G1@assays$RNA@data)
G1_meta <- data.frame(Cell=rownames(G1@meta.data), celltype = "Neu_CCL3")

## take Macro_SPP1
MF <- readRDS("MF.rds")
C1 <- subset(MF,idents=1)
C1_count <- as.matrix(C1@assays$RNA@data)
C1_meta <- data.frame(Cell=rownames(C1@meta.data), celltype = "Macro_SPP1")

count <- merge(G1_count,C1_count,by=0)
rownames(count) <- count[,1]
count <- count[,-1]
meta_data <- rbind(G1_meta,C1_meta)

write.table(count, 'cellphonedb/G1_C1/G1_C1_count.txt', sep='\t', quote=F, row.names=T)
write.table(meta_data, 'cellphonedb/G1_C1/G1_C1_meta.txt', sep='\t', quote=F, row.names=F)

################# cmd code  #####################
cellphonedb method statistical_analysis --counts-data=gene_name --threads=20 --output-path=cellphonedb/G1_C1/ --pvalue 0.01 cellphonedb/G1_C1/G1_C1_meta.txt  cellphonedb/G1_C1/G1_C1_count.txt
cd cellphonedb/G1_C1/ 
cellphonedb plot dot_plot --means-path ./means.txt --pvalues-path ./pvalues.txt --output-path ./
cellphonedb plot heatmap_plot G1_C1_meta.txt --pvalues-path ./pvalues.txt --output-path ./

## pick some L-R pair to plot
library(dplyr)
library(ggplot2)
plot_dir <- paste0(dir,"cellphonedb/G1_C1/")
mypvals <- read.delim(paste0(plot_dir,"pvalues.txt"), check.names = FALSE)
mymeans <- read.delim(paste0(plot_dir,"significant_means.txt"), check.names = FALSE)
# select cell pair
keep <- c(which(mypvals$`Neu_CCL3|Macro_SPP1` < 0.01),which(mypvals$`Macro_SPP1|Neu_CCL3` < 0.01))
keep <- unique(keep)
mypvals_keep <- mypvals[keep,]
mypair <- intersect(mypvals_keep$interacting_pair,mymeans$interacting_pair)
# pick out chemokines and so on
chemokines <- grep("^CXC|CCL|CCR|CX3|XCL|XCR|^IL[0-9]", mypair,value = T)
write.table(as.data.frame(chemokines),paste0(plot_dir,"G1_C1_chemokines.txt"),col.names=F,quote=F,row.names=F)
costimulatory <- grep("CD86|CD80|CD48|LILRB2|LILRB4|TNF|CD2|ICAM|SLAM|LT[AB]|NECTIN2|CD40|CD70|CD27|CD28|CD58|TSLP|PVR|CD44|CD55|CD[1-9]", 
              mypair,value = T)
write.table(as.data.frame(costimulatory),paste0(plot_dir,"G1_C1_costimulatory.txt"),col.names=F,quote=F,row.names=F)
coinhibitory <- grep("SIRP|CD47|ICOS|TIGIT|CTLA4|PC1D1|CD274|LAG3|HAVCR|VSIR|ENTPD1", 
                     mypair,value = T)
write.table(as.data.frame(coinhibitory),paste0(plot_dir,"G1_C1_coinhibitory.txt"),col.names=F,quote=F,row.names=F)
SPP1 <- grep("SPP1", mypair,value = T)
write.table(as.data.frame(SPP1),paste0(plot_dir,"G1_C1_SPP1.txt"),col.names=F,quote=F,row.names=F)
All_LR <- c(chemokines,costimulatory,coinhibitory,SPP1)
write.table(as.data.frame(All_LR),paste0(plot_dir,"G1_C1_All_LR.txt"),col.names=F,quote=F,row.names=F)

#### plot the selected L-R
cd cellphonedb/G1_C1/
cellphonedb plot dot_plot --means-path ./means.txt --pvalues-path ./pvalues.txt --output-path ./ --rows G1_C1_chemokines.txt --columns pair_plot.txt --output-name G1_C1_chemokines.pdf
cellphonedb plot dot_plot --means-path ./means.txt --pvalues-path ./pvalues.txt --output-path ./ --rows G1_C1_costimulatory.txt --columns pair_plot.txt --output-name G1_C1_costimulatory.pdf
cellphonedb plot dot_plot --means-path ./means.txt --pvalues-path ./pvalues.txt --output-path ./ --rows G1_C1_coinhibitory.txt --columns pair_plot.txt --output-name G1_C1_coinhibitory.pdf
cellphonedb plot dot_plot --means-path ./means.txt --pvalues-path ./pvalues.txt --output-path ./ --rows G1_C1_SPP1.txt --columns pair_plot.txt --output-name G1_C1_SPP1.pdf
cellphonedb plot dot_plot --means-path ./means.txt --pvalues-path ./pvalues.txt --output-path ./ --rows G1_C1_All_LR.txt --columns pair_plot.txt --output-name G1_C1_All_LR.pdf
cellphonedb plot dot_plot --means-path ./means.txt --pvalues-path ./pvalues.txt --output-path ./ --rows G1_C1_plot.txt --columns pair_plot.txt --output-name G1_C1_plot.pdf
