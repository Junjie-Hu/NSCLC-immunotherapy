library(Seurat)

dir <- "./"

## take Cancer cell
RNA_mag <- readRDS("Epi_mag.rds")
Cancer <- subset(RNA_mag,Group=="MPR")
Cancer_count <- as.matrix(Cancer@assays$RNA@data)
Cancer_meta <- data.frame(Cell=rownames(Cancer@meta.data), celltype = "Cancer")

## take NK_FCGR3A cell
RNA <- readRDS("Tcell.rds")
NK <- subset(RNA,idents=6)
NK_MPR <- subset(NK,Group=="MPR")
NK_count <- as.matrix(NK_MPR@assays$RNA@data)
NK_meta <- data.frame(Cell=rownames(NK_MPR@meta.data), celltype = "NK_FCGR3A")

count <- cbind(NK_count,Cancer_count)
meta_data <- rbind(Cancer_meta,NK_meta)

write.table(count, 'cellphonedb/Cancer_NK/Cancer_NK_count.txt', sep='\t', quote=F, row.names=T)
write.table(meta_data, 'cellphonedb/Cancer_NK/Cancer_NK_meta.txt', sep='\t', quote=F, row.names=F)

################# cmd code  #####################
cellphonedb method statistical_analysis --counts-data=gene_name --threads=20 --output-path=cellphonedb/Cancer_NK/ --pvalue 0.01 cellphonedb/Cancer_NK/Cancer_NK_meta.txt  cellphonedb/Cancer_NK/Cancer_NK_count.txt
cd cellphonedb/Cancer_NK/ 
cellphonedb plot dot_plot --means-path ./means.txt --pvalues-path ./pvalues.txt --output-path ./
cellphonedb plot heatmap_plot Cancer_NK_meta.txt --pvalues-path ./pvalues.txt --output-path ./

## pick some L-R pair to plot
library(dplyr)
library(ggplot2)
plot_dir <- paste0(dir,"cellphonedb/Cancer_NK/")
mypvals <- read.delim(paste0(plot_dir,"pvalues.txt"), check.names = FALSE)
mymeans <- read.delim(paste0(plot_dir,"significant_means.txt"), check.names = FALSE)
# select cell pair
keep <- c(which(mypvals$`NK_FCGR3A|Cancer` < 0.01),which(mypvals$`Cancer|NK_FCGR3A` < 0.01))
keep <- unique(keep)
mypvals_keep <- mypvals[keep,]
mypair <- intersect(mypvals_keep$interacting_pair,mymeans$interacting_pair)
# pick out chemokines and so on
chemokines <- grep("^CXC|CCL|CCR|CX3|XCL|XCR|^IL[0-9]", mypair,value = T)
write.table(as.data.frame(chemokines),paste0(plot_dir,"Cancer_NK_chemokines.txt"),col.names=F,quote=F,row.names=F)
costimulatory <- grep("CD86|CD80|CD48|LILRB2|LILRB4|TNF|CD2|ICAM|SLAM|LT[AB]|NECTIN2|CD40|CD70|CD27|CD28|CD58|TSLP|PVR|CD44|CD55|CD[1-9]", 
              mypair,value = T)
write.table(as.data.frame(costimulatory),paste0(plot_dir,"Cancer_NK_costimulatory.txt"),col.names=F,quote=F,row.names=F)
coinhibitory <- grep("SIRP|CD47|ICOS|TIGIT|CTLA4|PDCD1|CD274|LAG3|HAVCR|VSIR|ENTPD1", 
                     mypair,value = T)
write.table(as.data.frame(coinhibitory),paste0(plot_dir,"Cancer_NK_coinhibitory.txt"),col.names=F,quote=F,row.names=F)
All_LR <- c(chemokines,costimulatory,coinhibitory)
write.table(as.data.frame(All_LR),paste0(plot_dir,"Cancer_NK_All_LR.txt"),col.names=F,quote=F,row.names=F)

#### plot the selected L-R
cd cellphonedb/Cancer_NK/
cellphonedb plot dot_plot --means-path ./means.txt --pvalues-path ./pvalues.txt --output-path ./ --rows Cancer_NK_chemokines.txt --columns pair_plot.txt --output-name Cancer_NK_chemokines.pdf
cellphonedb plot dot_plot --means-path ./means.txt --pvalues-path ./pvalues.txt --output-path ./ --rows Cancer_NK_costimulatory.txt --columns pair_plot.txt --output-name Cancer_NK_costimulatory.pdf
cellphonedb plot dot_plot --means-path ./means.txt --pvalues-path ./pvalues.txt --output-path ./ --rows Cancer_NK_coinhibitory.txt --columns pair_plot.txt --output-name Cancer_NK_coinhibitory.pdf
cellphonedb plot dot_plot --means-path ./means.txt --pvalues-path ./pvalues.txt --output-path ./ --rows Cancer_NK_All_LR.txt --columns pair_plot.txt --output-name Cancer_NK_All_LR.pdf
cellphonedb plot dot_plot --means-path ./means.txt --pvalues-path ./pvalues.txt --output-path ./ --rows Cancer_NK_plot.txt --columns pair_plot.txt --output-name Cancer_NK_plot.pdf
