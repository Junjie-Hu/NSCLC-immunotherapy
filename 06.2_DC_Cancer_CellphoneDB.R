
library(Seurat)

dir <- "./"

## take cDC2 cell
MF <- readRDS("MF.rds")
DC <- subset(MF,idents=c(5))
DC_count <- as.matrix(DC@assays$RNA@data)
DC_meta <- data.frame(Cell=rownames(DC@meta.data), celltype = "cDC2")

## take Cancerthelium cell
RNA_mag <- readRDS("Epi_mag.rds")
Cancer <- subset(RNA_mag,Group=="MPR")
Cancer_count <- as.matrix(Cancer@assays$RNA@data)
Cancer_meta <- data.frame(Cell=rownames(Cancer@meta.data), celltype = "Tumor")

count <- merge(DC_count,Cancer_count,by=0)
meta_data <- rbind(DC_meta,Cancer_meta)

write.table(count, 'cellphonedb/DC_Cancer/DC_Cancer_count.txt', sep='\t', quote=F, row.names=F)
write.table(meta_data, 'cellphonedb/DC_Cancer/DC_Cancer_meta.txt', sep='\t', quote=F, row.names=F)

################# cmd code  #####################
cellphonedb method statistical_analysis --counts-data=gene_name --threads=20 --output-path=cellphonedb/DC_Cancer/ --pvalue 0.01 cellphonedb/DC_Cancer/DC_Cancer_meta.txt  cellphonedb/DC_Cancer/DC_Cancer_count.txt
cd cellphonedb/DC_Cancer/ 
cellphonedb plot dot_plot --means-path ./means.txt --pvalues-path ./pvalues.txt --output-path ./
cellphonedb plot heatmap_plot DC_Cancer_meta.txt --pvalues-path ./pvalues.txt --output-path ./

## pick some L-R pair to plot
library(dplyr)
library(ggplot2)
plot_dir <- paste0(dir,"cellphonedb/DC_Cancer/")
mypvals <- read.delim(paste0(plot_dir,"pvalues.txt"), check.names = FALSE)
mymeans <- read.delim(paste0(plot_dir,"significant_means.txt"), check.names = FALSE)
# select cell pair
keep <- c(which(mypvals$`cDC2|Tumor` < 0.01),which(mypvals$`Tumor|cDC2` < 0.01))
mypvals_keep <- mypvals[keep,]
mypair <- intersect(mypvals_keep$interacting_pair,mymeans$interacting_pair)
# pick out chemokines and so on
chemokines <- grep("^CXC|CCL|CCR|CX3|XCL|XCR", mypair,value = T)
write.table(as.data.frame(chemokines),paste0(plot_dir,"DC_Cancer_chemokines.txt"),col.names=F,quote=F,row.names=F)
costimulatory <- grep("CD86|CD80|CD48|LILRB2|LILRB4|TNF|CD2|ICAM|SLAM|LT[AB]|NECTIN2|CD40|CD70|CD27|CD28|CD58|TSLP|PVR|CD44|CD55|CD[1-9]", 
              mypair,value = T)
write.table(as.data.frame(costimulatory),paste0(plot_dir,"DC_Cancer_costimulatory.txt"),col.names=F,quote=F,row.names=F)
coinhibitory <- grep("SIRP|CD47|ICOS|TIGIT|CTLA4|PDCD1|CD274|LAG3|HAVCR|VSIR|ENTPD1", 
                     mypair,value = T)
write.table(as.data.frame(coinhibitory),paste0(plot_dir,"DC_Cancer_coinhibitory.txt"),col.names=F,quote=F,row.names=F)
All_LR <- c(chemokines,costimulatory,coinhibitory)
write.table(as.data.frame(All_LR),paste0(plot_dir,"DC_Cancer_All_LR.txt"),col.names=F,quote=F,row.names=F)

#### plot the selected L-R
cd cellphonedb/DC_Cancer/
cellphonedb plot dot_plot --means-path ./means.txt --pvalues-path ./pvalues.txt --output-path ./ --rows DC_Cancer_chemokines.txt --columns pair_plot.txt --output-name DC_Cancer_chemokines.pdf
cellphonedb plot dot_plot --means-path ./means.txt --pvalues-path ./pvalues.txt --output-path ./ --rows DC_Cancer_costimulatory.txt --columns pair_plot.txt --output-name DC_Cancer_costimulatory.pdf
cellphonedb plot dot_plot --means-path ./means.txt --pvalues-path ./pvalues.txt --output-path ./ --rows DC_Cancer_coinhibitory.txt --columns pair_plot.txt --output-name DC_Cancer_coinhibitory.pdf
cellphonedb plot dot_plot --means-path ./means.txt --pvalues-path ./pvalues.txt --output-path ./ --rows DC_Cancer_All_LR.txt --columns pair_plot.txt --output-name DC_Cancer_All_LR.pdf
cellphonedb plot dot_plot --means-path ./means.txt --pvalues-path ./pvalues.txt --output-path ./ --rows DC_Cancer_plot.txt --columns pair_plot.txt --output-name DC_Cancer_plot.pdf

