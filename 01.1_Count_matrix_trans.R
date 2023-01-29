# Load libraries
library(data.table)
library(future)
plan("multiprocess", workers = 8)
options(future.globals.maxSize = 10*1024^3)

# Set the full path to the github directory 
rm(list=ls())
dir <- "./"

# Load raw count matrix from BD Rhapsody, transpose the matrix and save it
counts_list <- list.files(path=paste0(dir,"count_matrix_raw/"),pattern="_RSEC_MolsPerCell.csv")
sample_num <- length(counts_list)
sample_id <- gsub("_RSEC_MolsPerCell.csv","",counts_list)
for(i in 1:sample_num){
   data <- fread(paste0(dir,"count_matrix_raw/",counts_list[i]),sep=",",header= F,stringsAsFactors=F)
   dim(data)
   head(data[,1:5])
   data <- t(data)
   dim(data)
   data[1,2:ncol(data)] <- paste(sample_id[i],data[1,2:ncol(data)],sep = "_")
   head(data[,1:5])
   fwrite(data,file=paste0(dir,"count_matrix_trans/",sample_id[i],".txt"),row.names=F,col.names=F,quote=F,sep="\t")
}
