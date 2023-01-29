# Load libraries
library(data.table)
library(Seurat)
library(future)
plan("multiprocess", workers = 8)
options(future.globals.maxSize = 10*1024^3)

# Set the full path to the github directory 
rm(list=ls())
dir <- "./"

# read in each matrix
txt_list <- list.files(path=paste0(dir,"count_matrix_trans/"),pattern=".txt")
txt_num <- length(txt_list)
txt_id <- gsub(".txt","",txt_list)
for(i in 1:txt_num){
	data <- read.table(paste0(dir,"count_matrix_trans/",txt_list[i]),sep="\t",header=T,stringsAsFactors=F,row.names=1)
	assign(txt_id[i],data)
}

###### QC for each sample one by one
sample_name <- BD_immune01
## create Seurat object
RNA <- CreateSeuratObject(sample_name)
## remove low quality cells
# Find  mitochondrial genes, compute the mitochondrial rate
mito.genes <- grep("^MT-", rownames(sample_name), value = TRUE)
percent.mito <- Matrix::colSums(sample_name[mito.genes, ])/Matrix::colSums(sample_name)
RNA <- AddMetaData(object = RNA, metadata = percent.mito, col.name = "percent.mito")
# Find ribosomal genes, compute the mitochondrial rate
ribo.genes <- grep(pattern = "^RP[SL][[:digit:]]", x = rownames(x = sample_name), value = TRUE)
percent.ribo <- Matrix::colSums(sample_name[ribo.genes, ])/Matrix::colSums(sample_name)
RNA <- AddMetaData(object = RNA, metadata = percent.ribo, col.name = "percent.ribo")
# Calculate housekeeping gene score(UMI sum)
housekeeking_marker <- c("ACTB","GAPDH","MALAT1")
hw_data <- sample_name[housekeeking_marker,]
housekeeping_UMI <- colSums(hw_data)
RNA <- AddMetaData(object = RNA, metadata = housekeeping_UMI, col.name = "housekeeping_score")
# summary
dim(RNA);median(RNA$nFeature_RNA);median(RNA$nCount_RNA);median(RNA@meta.data$percent.mito)
table(RNA$nFeature_RNA < 500);table(RNA@meta.data$percent.mito > 0.2);table(RNA@meta.data$percent.ribo > 0.5);table(RNA@meta.data$housekeeping_score < 1)
# filter
RNA <- subset(RNA, nFeature_RNA >= 500) 
RNA <- subset(RNA, percent.mito <= 0.2)
RNA <- subset(RNA, percent.ribo <= 0.5) 
RNA <- subset(RNA, housekeeping_score >= 1)
dim(RNA) 

## save filtered reslts 
write.table(RNA@assays$RNA@counts,paste0(dir,"count_matrix_filter/","BD_immune01","_filter.txt"),quote=F,sep="\t")
