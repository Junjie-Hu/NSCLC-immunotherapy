library(Seurat)

dir <- "./"

## take B4_FCRL4
Bcell <- readRDS("Bcell.rds")
B4 <- subset(Bcell,idents=4)
B4_count <- as.matrix(B4@assays$RNA@data)
B4_meta <- data.frame(Cell=rownames(B4@meta.data), celltype = "B4_FCRL4")

## Tfh
Tcell <- readRDS("Tcell.rds")
T4 <- subset(Tcell,T_cluster_annotation %in% c("CD4_MAF","CD4_CXCL13"))
T4_count <- as.matrix(T4@assays$RNA@data)
T4_meta <- data.frame(Cell=rownames(T4@meta.data), celltype = "Tfh")

## take CD8 cell
T8 <- subset(Tcell,T_cluster_annotation %in% c("CD8_IL7R","CD8_GZMK","CD8_GZMB","CD8_HAVCR2","CD8_STMN1"))
T8_count <- as.matrix(T8@assays$RNA@data)
T8_meta <- data.frame(Cell=rownames(T8@meta.data), celltype = "CD8 T")

count <- cbind(B4_count,T4_count,T8_count)
meta_data <- rbind(B4_meta,T4_meta,T8_meta)

write.table(count, 'cellphonedb/B4_Tfh_T8/B4_Tfh_T8_count.txt', sep='\t', quote=F, row.names=T)
write.table(meta_data, 'cellphonedb/B4_Tfh_T8/B4_Tfh_T8_meta.txt', sep='\t', quote=F, row.names=F)

################# cmd code  #####################
cellphonedb method statistical_analysis --counts-data=gene_name --threads=20 --output-path=cellphonedb/B4_Tfh_T8/ --pvalue 0.01 cellphonedb/B4_Tfh_T8/B4_Tfh_T8_meta.txt  cellphonedb/B4_Tfh_T8/B4_Tfh_T8_count.txt
cd cellphonedb/B4_Tfh_T8/ 
cellphonedb plot dot_plot --means-path ./means.txt --pvalues-path ./pvalues.txt --output-path ./
cellphonedb plot heatmap_plot B4_Tfh_T8_meta.txt --pvalues-path ./pvalues.txt --output-path ./

## pick some L-R pair to plot
library(dplyr)
library(ggplot2)
plot_dir <- paste0(dir,"cellphonedb/B4_Tfh_T8/")
mypvals <- read.delim(paste0(plot_dir,"pvalues.txt"), check.names = FALSE)
mymeans <- read.delim(paste0(plot_dir,"significant_means.txt"), check.names = FALSE)
# select cell pair
keep <- c(which(mypvals$`B4_FCRL4|Tfh` < 0.01),which(mypvals$`Tfh|B4_FCRL4` < 0.01),
	which(mypvals$`B4_FCRL4|CD8 T` < 0.01),which(mypvals$`CD8 T|B4_FCRL4` < 0.01),
	which(mypvals$`CD8 T|Tfh` < 0.01),which(mypvals$`Tfh|CD8 T` < 0.01))
keep <- unique(keep)
mypvals_keep <- mypvals[keep,]
mypair <- intersect(mypvals_keep$interacting_pair,mymeans$interacting_pair)
# pick out chemokines and so on
chemokines <- grep("^CXC|CCL|CCR|CX3|XCL|XCR|^IL[1-9]", mypair,value = T)
write.table(as.data.frame(chemokines),paste0(plot_dir,"B4_Tfh_T8_chemokines.txt"),col.names=F,quote=F,row.names=F)
costimulatory <- grep("CD86|CD80|CD48|LILRB2|LILRB4|TNF|CD2|ICAM|SLAM|LT[AB]|NECTIN2|CD40|CD70|CD27|CD28|CD58|TSLP|PVR|CD44|CD55|CD[1-9]", 
              mypair,value = T)
write.table(as.data.frame(costimulatory),paste0(plot_dir,"B4_Tfh_T8_costimulatory.txt"),col.names=F,quote=F,row.names=F)
coinhibitory <- grep("SIRP|CD47|ICOS|TIGIT|CTLA4|PDCD1|CD274|LAG3|HAVCR|VSIR|ENTPD1", 
                     mypair,value = T)
write.table(as.data.frame(coinhibitory),paste0(plot_dir,"B4_Tfh_T8_coinhibitory.txt"),col.names=F,quote=F,row.names=F)
All_LR <- c(chemokines,costimulatory,coinhibitory)
write.table(as.data.frame(All_LR),paste0(plot_dir,"B4_Tfh_T8_All_LR.txt"),col.names=F,quote=F,row.names=F)

#### plot the selected L-R
cd cellphonedb/B4_Tfh_T8/
cellphonedb plot dot_plot --means-path ./means.txt --pvalues-path ./pvalues.txt --output-path ./ --rows B4_Tfh_T8_chemokines.txt --columns pair_plot.txt --output-name B4_Tfh_T8_chemokines.pdf
cellphonedb plot dot_plot --means-path ./means.txt --pvalues-path ./pvalues.txt --output-path ./ --rows B4_Tfh_T8_costimulatory.txt --columns pair_plot.txt --output-name B4_Tfh_T8_costimulatory.pdf
cellphonedb plot dot_plot --means-path ./means.txt --pvalues-path ./pvalues.txt --output-path ./ --rows B4_Tfh_T8_coinhibitory.txt --columns pair_plot.txt --output-name B4_Tfh_T8_coinhibitory.pdf
cellphonedb plot dot_plot --means-path ./means.txt --pvalues-path ./pvalues.txt --output-path ./ --rows B4_Tfh_T8_All_LR.txt --columns pair_plot.txt --output-name B4_Tfh_T8_All_LR.pdf
cellphonedb plot dot_plot --means-path ./means.txt --pvalues-path ./pvalues.txt --output-path ./ --rows B4_Tfh_T8_plot.txt --columns pair_plot.txt --output-name B4_Tfh_T8_plot.pdf
