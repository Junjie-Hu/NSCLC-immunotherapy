
library(Seurat)
library(future)
plan("multiprocess", workers = 8)
options(future.globals.maxSize = 10*1024^3)

dir <- "./"

## take PMo
RNA <- readRDS("MF.rds")
PMo <- subset(RNA,idents=7)
PMo_MPR <- subset(RNA,Group=="MPR")
PMo_count <- as.matrix(PMo_MPR@assays$RNA@data)
PMo_meta <- data.frame(Cell=rownames(PMo_MPR@meta.data), celltype = "PMo")

## take Cancerthelium cell
RNA <- readRDS("epithelium.rds")
RNA_mag <- subset(RNA,TC.cluster %in% c(0,1,3,4,7))
Cancer <- subset(RNA_mag,Group=="MPR")
Cancer_count <- as.matrix(Cancer@assays$RNA@data)
Cancer_meta <- data.frame(Cell=rownames(Cancer@meta.data), celltype = "Cancer")

count <- merge(PMo_count,Cancer_count,by=0)
meta_data <- rbind(PMo_meta,Cancer_meta)

write.table(count, 'cellphonedb/PMo_Cancer/PMo_Cancer_count.txt', sep='\t', quote=F, row.names=F)
write.table(meta_data, 'cellphonedb/PMo_Cancer/PMo_Cancer_meta.txt', sep='\t', quote=F, row.names=F)

################# cmd code  #####################
cellphonedb method statistical_analysis --counts-data=gene_name --threads=20 --output-path=cellphonedb/PMo_Cancer/ --pvalue 0.01 cellphonedb/PMo_Cancer/PMo_Cancer_meta.txt  cellphonedb/PMo_Cancer/PMo_Cancer_count.txt
cd cellphonedb/PMo_Cancer/ 
cellphonedb plot dot_plot --means-path ./means.txt --pvalues-path ./pvalues.txt --output-path ./
cellphonedb plot heatmap_plot PMo_Cancer_meta.txt --pvalues-path ./pvalues.txt --output-path ./

## pick some L-R pair to plot
library(dplyr)
library(ggplot2)
plot_dir <- paste0(dir,"cellphonedb/PMo_Cancer/")
mypvals <- read.delim(paste0(plot_dir,"pvalues.txt"), check.names = FALSE)
mymeans <- read.delim(paste0(plot_dir,"significant_means.txt"), check.names = FALSE)
# select cell pair
keep <- c(which(mypvals$`PMo|Cancer` < 0.01),which(mypvals$`Cancer|PMo` < 0.01))
mypvals_keep <- mypvals[keep,]
mypair <- intersect(mypvals_keep$interacting_pair,mymeans$interacting_pair)
# pick out chemokines and so on
chemokines <- grep("^CXC|CCL|CCR|CX3|XCL|XCR", mypair,value = T)
write.table(as.data.frame(chemokines),paste0(plot_dir,"PMo_Cancer_chemokines.txt"),col.names=F,quote=F,row.names=F)
costimulatory <- grep("CD86|CD80|CD48|LILRB2|LILRB4|TNF|CD2|ICAM|SLAM|LT[AB]|NECTIN2|CD40|CD70|CD27|CD28|CD58|TSLP|PVR|CD44|CD55|CD[1-9]", 
              mypair,value = T)
write.table(as.data.frame(costimulatory),paste0(plot_dir,"PMo_Cancer_costimulatory.txt"),col.names=F,quote=F,row.names=F)
coinhibitory <- grep("SIRP|CD47|ICOS|TIGIT|CTLA4|PDCD1|CD274|LAG3|HAVCR|VSIR|ENTPD1", 
                     mypair,value = T)
write.table(as.data.frame(coinhibitory),paste0(plot_dir,"PMo_Cancer_coinhibitory.txt"),col.names=F,quote=F,row.names=F)
All_LR <- c(chemokines,costimulatory,coinhibitory)
write.table(as.data.frame(All_LR),paste0(plot_dir,"PMo_Cancer_All_LR.txt"),col.names=F,quote=F,row.names=F)

#### plot the selected L-R
cd cellphonedb/PMo_Cancer/
cellphonedb plot dot_plot --means-path ./means.txt --pvalues-path ./pvalues.txt --output-path ./ --rows PMo_Cancer_chemokines.txt --columns pair_plot.txt --output-name PMo_Cancer_chemokines.pdf
cellphonedb plot dot_plot --means-path ./means.txt --pvalues-path ./pvalues.txt --output-path ./ --rows PMo_Cancer_costimulatory.txt --columns pair_plot.txt --output-name PMo_Cancer_costimulatory.pdf
cellphonedb plot dot_plot --means-path ./means.txt --pvalues-path ./pvalues.txt --output-path ./ --rows PMo_Cancer_coinhibitory.txt --columns pair_plot.txt --output-name PMo_Cancer_coinhibitory.pdf
cellphonedb plot dot_plot --means-path ./means.txt --pvalues-path ./pvalues.txt --output-path ./ --rows PMo_Cancer_All_LR.txt --columns pair_plot.txt --output-name PMo_Cancer_All_LR.pdf


# plot 
pldf %>% filter(means >1) %>% 
  ggplot(aes(CC.x,interacting_pair.x) )+ 
  geom_point(aes(color=means,size=-log10(pvals)) ) +
  scale_size_continuous(range = c(1,3))+
  scale_color_gradient2(high="red",mid = "yellow",low ="darkblue",midpoint = 2)+ theme_bw()+ 
  theme(axis.text.x = element_text(angle = -45,hjust = -0.1,vjust = 0.8)) 
