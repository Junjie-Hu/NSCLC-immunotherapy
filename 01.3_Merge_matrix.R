# Load libraries
library(data.table)
library(future)
plan("multiprocess", workers = 8)
options(future.globals.maxSize = 10*1024^3)

# Set the full path to the github directory 
rm(list=ls())
dir <- "./"

# Merge all count matrix and save it
txt_list <- list.files(path=paste0(dir,"count_matrix_scrublet/"),pattern=".csv")
txt_num <- length(txt_list)
txt_id <- gsub(".csv","",txt_list)
for(i in 1:txt_num){
	data <- read.csv(paste0(dir,"count_matrix_scrublet/",txt_list[i]),header=T,stringsAsFactors=F,row.names=1)
	assign(txt_id[i],data)
}
ls()

count_all <- merge(get(txt_id[1]),get(txt_id[2]),by=0)
rownames(count_all) <- count_all[,1]
count_all <- count_all[,-1]
for(i in 3:txt_num){
	count_all <- merge(count_all,get(txt_id[i]),by=0)
        rownames(count_all) <- count_all[,1]
        count_all <- count_all[,-1]
}
dim(count_all)
write.table(count_all,file=paste0(dir,"count_all_scrublet.txt"),quote=F,sep="\t")

# prepare metadata and save it
clin_info <- read.csv("clin_info.csv",header=T,stringsAsFactors=F)
metadata <- data.frame(Cell = colnames(count_all),
	                   Sample_ID = gsub("_[0-9]{1,}","",colnames(count_all)))
metadata <- merge(metadata,clin_info,by="Sample_ID")
metadata <- metadata[,c(2,1,3:ncol(metadata))]
write.csv(metadata,paste0(dir,"metadata_scrublet.csv"),row.names=F,quote=F)


