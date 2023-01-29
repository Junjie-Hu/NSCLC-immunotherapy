# load libraries
library(Seurat)
library(ggplot2)
library(future)
library(tidyverse)
library(gridExtra)
library(ggridges)
library(ggplot2)
library(ggExtra)
plan("multiprocess", workers = 8)
options(future.globals.maxSize = 10*1024^3)

# set work directory
rm(list=ls())
dir <- "./"

### read raw count matrix
raw.data <- read.table(paste0(dir,"count_all_scrublet.txt"),header=T,row.names=1,stringsAsFactors=F)

### read metadata
metadata <- read.csv("metadata_scrublet.csv", row.names=1, header=T,stringsAsFactors=F)

# Create the Seurat object with all the data (filtered)
RNA <- CreateSeuratObject(counts = raw.data)
# add metadata to Seurat object 
RNA <- AddMetaData(object = RNA, metadata = metadata)
# Head to check
head(RNA@meta.data)

# Save Seurat object 
saveRDS(RNA, file=paste0(dir,"all_cell.rds"))

## check gene and count number distribution
RNA <- readRDS(file=paste0(dir,"all_cell.rds"))
library(ggplot2)
library(cowplot)
p <- ggplot(data=RNA@meta.data,aes(x=nFeature_RNA))+geom_histogram(bins=50,fill="#219ebc",color="black")+labs(x="Genes",y="Cell counts",title="")+theme_cowplot()+
scale_x_continuous(limits=c(400,8500),breaks=seq(500,8500,1000))+
theme(axis.text.x=element_text(size=15,angle = 45, hjust=1, vjust=1),axis.text.y=element_text(size = 15),axis.title.y=element_text(size = 18),axis.title.x=element_text(size = 18),axis.ticks.length = unit(0.3,"cm"))
ggsave(paste0(dir,"/plot_out/all_cell_gene_his.pdf"),p,width=8,heigh=6)

p <- ggplot(data=RNA@meta.data,aes(x=nCount_RNA))+geom_histogram(bins=50,fill="#219ebc",color="black")+labs(x="UMIs",y="Cell counts",title="")+theme_cowplot()+
scale_x_continuous(limits=c(500,12500),breaks=seq(500,12500,1500))+
theme(axis.text.x=element_text(size=15,angle = 45, hjust=1, vjust=1),axis.text.y=element_text(size = 15),axis.title.y=element_text(size = 18),axis.title.x=element_text(size = 18),axis.ticks.length = unit(0.3,"cm"))
ggsave(paste0(dir,"/plot_out/all_cell_count_his.pdf"),p,width=8,heigh=6)

### first cluster
RNA <- NormalizeData(object = RNA)
RNA <- FindVariableFeatures(object = RNA)
RNA <- ScaleData(object = RNA)
RNA <- RunPCA(object = RNA, do.print = FALSE)
ElbowPlot(RNA)
RNA <- FindNeighbors(object = RNA, verbose = T, dims = 1:20)
RNA <- FindClusters(object = RNA, verbose = T, resolution = 0.6)
RNA <- RunUMAP(RNA, dims = 1:20)
# Visualize UMAP colored by cluster
p=DimPlot(RNA, reduction = "umap", label = T)
ggsave(paste0(dir,"plot_out/all_cell_UMAP.pdf"),p,width=6.5, height=6)
# Check batch 
pdf(paste(dir,"plot_out/all_cell_meta.pdf", sep=""),7,6)
DimPlot(RNA, reduction = "umap", label = F, group.by = "Patient_ID")
DimPlot(RNA, reduction = "umap", label = F, group.by = "Group")
DimPlot(RNA, reduction = "umap", label = F, group.by = "PD1_name")
dev.off()
## there were probabaly batch effects on Patients

### Remove batch effect by CCA
data.list <- SplitObject(RNA, split.by = "Sample_ID")
for (i in 1:length(data.list)) {
    data.list[[i]] <- NormalizeData(data.list[[i]], verbose = FALSE)
    data.list[[i]] <- FindVariableFeatures(data.list[[i]],
        selection.method = "vst", nfeatures = 3000,
        verbose = FALSE)
}
anchors <- FindIntegrationAnchors(object.list = data.list,
    dims = 1:20, anchor.features = 3000)
RNA.integrated <- IntegrateData(anchorset = anchors, dims = 1:20)
RNA.integrated@project.name <- paste0(RNA@project.name, "_CCA")

### cluster
RNA.integrated <- ScaleData(object = RNA.integrated)
RNA.integrated <- RunPCA(object = RNA.integrated, do.print = FALSE)
p <- ElbowPlot(RNA.integrated)
ggsave(paste0(dir,"/plot_out/all_cell_elbow_CCA.png"),p,width=5,heigh=5)

RNA.integrated <- FindNeighbors(object = RNA.integrated, verbose = T, dims = 1:20)
RNA.integrated <- FindClusters(object = RNA.integrated, verbose = T, resolution = 0.6)
RNA.integrated <- RunUMAP(RNA.integrated, dims = 1:20)
# plot on UMAP
p=DimPlot(RNA.integrated, reduction = "umap", label = T)
ggsave(paste0(dir,"plot_out/all_cell_UMAP_CCA.pdf"),p,width=7, height=6)
# recheck batch
pdf(paste(dir,"plot_out/all_cell_meta_umap_CCA.pdf", sep=""),7,6)
DimPlot(RNA.integrated, reduction = "umap", label = F, group.by = "Patient_ID")
DimPlot(RNA.integrated, reduction = "umap", label = F, group.by = "Group")
DimPlot(RNA.integrated, reduction = "umap", label = F, group.by = "PD1_name")
dev.off()

#### integrate CCA cluster res into pre-CCA seurat object
RNA[["umap"]] <- RNA.integrated[["umap"]]
RNA[["tsne"]] <- RNA.integrated[["tsne"]]
Idents(RNA) <- Idents(RNA.integrated)
RNA@meta.data$integrated_snn_res.0.6 <- RNA.integrated@meta.data$integrated_snn_res.0.6
# plot on UMAP
p=DimPlot(RNA, reduction = "umap", label = T)
ggsave(paste0(dir,"plot_out/all_cell_UMAP.pdf"),p,width=7, height=6)
## find marker for each cluster
markers <- FindAllMarkers(object = RNA, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.csv(markers, file=paste0(dir,"data_out/all_cell_clusters.csv"),quote=F)

# Visualize nFeature and nCount on UMAP
data <- Embeddings(object = RNA[["umap"]])[, c(1,2)]
data <- as.data.frame(data)
data$nFeature <- RNA$nFeature_RNA
data$nCount <- RNA$nCount_RNA
library(cowplot)
# plot nFeature
p <- ggplot(data = data,aes(x=data[,1],y=data[,2]))+
     geom_point(aes_string(color=data$nFeature))+
     theme(axis.title.x = element_blank(),axis.title.y = element_blank())+theme_cowplot()+
     scale_color_gradient(low = "yellow", high = "#DE2D26")+
     guides(color = guide_colorbar(title = "Genes"))+
     labs(x="UMAP_1",y="UMAP_2")
ggsave(paste0(dir,"plot_out/all_cell_nFeature_CCA.pdf"),p,width =7,height =6)
# plot nCount
p <- ggplot(data = data,aes(x=data[,1],y=data[,2]))+
     geom_point(aes_string(color=data$nCount))+
     theme(axis.title.x = element_blank(),axis.title.y = element_blank())+theme_cowplot()+
     scale_color_gradient(low = "lightblue", high = "blue")+
     guides(color = guide_colorbar(title = "UMIs"))+
     labs(x="UMAP_1",y="UMAP_2")
ggsave(paste0(dir,"plot_out/all_cell_nCount_CCA.pdf"),p,width =7,height =6)
# plot PD1_name
data$PD1_name <- factor(RNA$PD1_name,levels=c("Toripalimab","Sintilimab","Camrelizumab"))
p <- ggplot(data = data,aes(x=data[,1],y=data[,2]))+
     geom_point(aes_string(color=data$PD1_name),size=0.3)+
     theme(axis.title.x = element_blank(),axis.title.y = element_blank())+theme_cowplot()+
     scale_color_manual(name = "",values = c("Toripalimab"="#619cff","Sintilimab"="#00ba38","Camrelizumab"="#f8766d"))+
     labs(x="UMAP_1",y="UMAP_2")
ggsave(paste0(dir,"plot_out/all_cell_PD1.pdf"),p,width =7,height =6)

# Barplot of clusters per patient
library(cowplot)
tab1 <- cbind(as.data.frame(RNA@meta.data$Patient_ID),as.data.frame(RNA@meta.data$integrated_snn_res.0.6))
colnames(tab1) <- c("Patient", "Clusters")
pdf(paste0(dir,"plot_out/all_cell_proportion_per_patient.pdf"),width=12,height=8)
ggplot(data=tab1,aes(x = Patient, fill = factor(Clusters))) + theme_bw()+theme_cowplot()+
  geom_bar(position = "fill",width=0.6) + theme(axis.text.x = element_text(size = 15,hjust = 1, vjust = 1,angle = 45),axis.text.y = element_text(size = 15),axis.ticks.length = unit(0.3,"cm"))+
  labs(x="",y="Cellular fraction",fill="Cluster")+theme(axis.title.y=element_text(size = 16),legend.title = element_text(size =15))
dev.off()
# Barplot of patients per cluster
pdf(paste(dir,"plot_out/all_cell_proportion_per_cluster.pdf", sep=""),width=10,height=15)
ggplot(tab1) +
  aes(x = factor(Clusters,levels=c(0:25)), fill = Patient) + theme_bw()+theme_cowplot()+
  geom_bar(position = "fill",width=0.6) + labs(x="Clusters",y="")+
  theme(axis.text.y = element_text(size=15),axis.text.x = element_text(size=15),legend.title=element_blank(),axis.ticks.length = unit(0.3,"cm"))
dev.off()

# Manually ananotate immmune markers to identify major cell clusters
cell_marker=c("PTPRC","CD3E","LYZ","CD79A","MS4A1","IGHG1","CSF3R","FGFBP2","KIT","LILRA4","VWF","COL1A1","EPCAM","MKI67")
library(ggExtra)
marker_plot=function(SeuratObj,marker){
    pn=length(marker)
    pp=list()
    for(i in 1:pn){
    	pg=FeaturePlot(object = SeuratObj, features = marker[i], cols = c("grey", "red"), reduction = "umap")+NoLegend()+labs(x="",y="")+theme(plot.title = element_text(hjust = 0.5))
    	pp[[marker[i]]]=pg
    }
    return(pp)
}
p=marker_plot(RNA,cell_marker)
for(i in 1:14){
ggsave(paste0(dir,"plot_out/all_cell_",cell_marker[i],".png"),p[[cell_marker[i]]],width =6,height =6)
}

###### Manually pick and ananotate each cluster to immmune, stromal and epithelial cells
immune_cluster <- c(0,1,2,3,5,6,7,8,9,10,11,12,14,15,18,19,21,23,24,25)
stromal_cluster <- c(20)
epi_cluster <- c(4,13,16,17,22)
RNA@meta.data$major_cluster_annotation <- ifelse(RNA@meta.data$integrated_snn_res.0.6 %in% immune_cluster, "Immune cell",ifelse(RNA@meta.data$integrated_snn_res.0.6 %in% stromal_cluster,"Stromal cell","Epithelium"))
# plot tsne
pdf(paste(dir,"plot_out/major_cluster_annotation.pdf", sep=""),width = 7,height = 6)
DimPlot(RNA, reduction = "tsne", label = T, group.by = 'major_cluster_annotation')
dev.off()
# plot umap
pdf(paste(dir,"plot_out/major_cluster_annotation_umap.pdf", sep=""),width = 7,height = 6)
DimPlot(RNA, reduction = "umap", label = T, group.by = 'major_cluster_annotation')
dev.off()

# Manually ananotate cell markers
# stash current cluster IDs
RNA[["all.cluster"]] <- Idents(object = RNA)
# enumerate current cluster IDs and the labels for them
cluster.ids <- 0:(length(unique(RNA@meta.data$all.cluster))-1)
# Annotate each of the clusters 
free_annotation <- c("T cell","B cell","T cell","Neutrophil","Epithelium","T cell","T cell","Myeloid cell","Myeloid cell","Myeloid cell",
  "Neutrophil","Plasma cell","T cell","Epithelium","T cell","Cycling immune cell","Epithelium","Epithelium","T cell","NK cell","Fibroblast/Endothelium",
  "Mast cell","Epithelium","pDC","Plasma cell","Myeloid cell")         
# Map free annotation to cluster numbers and store as all_subtype_annotation
RNA@meta.data[,'all_cluster_annotation'] <- plyr::mapvalues(x = RNA@meta.data$all.cluster, from = cluster.ids, to = free_annotation)
# plot umap
pdf(paste(dir,"plot_out/all_cluster_annotation_umap.pdf", sep=""),width = 8,height = 6)
DimPlot(RNA, reduction = "umap", label = F, group.by = 'all_cluster_annotation')
dev.off()

### save annotation
saveRDS(RNA,paste0(dir,"all_cell.rds"))

RNA <- readRDS(paste0(dir,"all_cell.rds"))
### UMAP plot in TN, MPR and NMPR, respectively
for(i in c("TN","MPR","NMPR")){
    RNA_sub <- subset(RNA,Group==i)
    pdf(paste0(dir,"plot_out/all_cluster_",i,"_umap.pdf"),width = 8,height = 6)
    DimPlot(RNA_sub, reduction = "umap", label = F, group.by = 'all_cluster_annotation')
    dev.off()
}

###Calculate cell fractions in different response groups
source(paste0(dir,"code/Cal_fraction.R"))
tab.1 <- cal_fraction(RNA,"all_cluster_annotation")
# plot
p1<- ggplot(tab.1, aes(x=Group, y=Estimate)) +
    geom_bar(stat = "identity", aes(fill=cell)) + facet_grid(cols =  vars(cell)) + 
    theme_bw() +  
    theme(axis.text.x = element_text(size=13,angle = 45, hjust=1, vjust=1), axis.text.y = element_text(size=13),axis.ticks.length = unit(0.3,"cm"),legend.position= "none") + 
    xlab("") + ylab("Cellular fraction (%)") +
    theme(axis.title.y=element_text(size=16),strip.text = element_text(size = 13))+
    geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,position=position_dodge(0.05))
# Save plot 
ggsave(p1, filename = paste(dir,"plot_out/all_cell_across_treatment.pdf", sep=""),width = 15, height = 5)

##take out the immune cell
RNA_immune <- subset(RNA, integrated_snn_res.0.6 %in% immune_cluster)
saveRDS(RNA_immune, paste0(dir,"immune.rds"))
##take out the non-immune cell
RNA_non_immune <- subset(RNA, integrated_snn_res.0.6 %in% c(stromal_cluster,epi_cluster))
saveRDS(RNA_non_immune, paste0(dir,"non_immune.rds"))
##take out the stromal cell
RNA_stromal <- subset(RNA, integrated_snn_res.0.6 %in% stromal_cluster)
saveRDS(RNA_stromal, paste0(dir, "stromal.rds"))
##take out the epithelium
RNA_epi <- subset(RNA, integrated_snn_res.0.6 %in% epi_cluster)
saveRDS(RNA_epi, paste0(dir, "epithelium.rds"))


###################
###################  new plot
###Calculate cell fractions in different response groups

source(paste0(dir,"code/Cal_fraction_patients.R"))
tab.2 <- cal_fraction_patients(RNA,"all_cluster_annotation")
colnames(tab.2)[4] <- "Patient_ID"

clin_info <- read.csv("clin_info.csv",header=T,stringsAsFactors=F)
sample_info <- clin_info[,c("Patient_ID","Group")]
table.plot <- merge(tab.2,sample_info,by="Patient_ID")
table.plot$Group <- factor(table.plot$Group,levels=c("TN","MPR","NMPR"))

Group_color_panel <- c("TN"="#228B22","MPR"="#FF4500","NMPR"="#757bc8")

p1<- ggplot(tab.1, aes(x=Group, y=Estimate)) +
    geom_bar(stat = "identity", aes(fill=Group)) +
    scale_fill_manual(name = "Group", values = Group_color_panel) + facet_grid(~ cell) + 
    theme_bw() + 
    theme(axis.text.x = element_text(size=13,angle = 45, hjust=1, vjust=1), axis.text.y = element_text(size=13),axis.ticks.length = unit(0.2,"cm")) + 
    xlab("") + ylab("Cellular fraction (%) of all cells") +
    theme(axis.title.y=element_text(size=13),strip.text = element_text(size = 13))+
    geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,position=position_dodge(0.05))

pdf(paste(dir,"plot_out/all_cell_across_treatment3.pdf", sep=""),width = 16,height = 5)
    p1+geom_point(data=table.plot,aes(x=Group, y=Estimate,color=Patient_ID),size = 1.5)
dev.off()

ggplot(table.plot, aes(x=cell, y=Estimate,fill = Group)) +
stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, 
    geom = 'bar', width = 0.3, size = 0.3,fill = 'transparents') + 
stat_summary(fun.data = function(x) median_hilow(x, 0.5), 
    geom = 'errorbar', width = 0.25, size = 0.2) +  
geom_jitter(size = 1.5, width = 0.2,stat = "identity") +  
scale_fill_manual(name = "Group", values = Group_color_panel) + facet_grid(~ cell) +
theme_bw() + 
theme(axis.text.x = element_text(size=13,angle = 45, hjust=1, vjust=1), axis.text.y = element_text(size=13),axis.ticks.length = unit(0.2,"cm")) + 
xlab("") + ylab("Cellular fraction (%) of all cells") +
theme(axis.title.y=element_text(size=13),strip.text = element_text(size = 13))

# plot
p2<- ggplot(table.plot) +
    geom_boxplot(aes(x=cell, y=Estimate,color=Group,outlier.colour = NA),width = 0.6,position=position_dodge(0.75))+
    geom_jitter(aes(x=cell, y=Estimate,color=Group),size = 1.5, position = position_jitterdodge())+
    scale_color_manual(name="Group",values = Group_color_panel)+
    theme_bw() +  
    theme(axis.text.x = element_text(size=13), axis.text.y = element_text(size=13),axis.ticks.length = unit(0.3,"cm")) + 
    xlab("") + ylab("Cellular fraction (%) of all cells") +
    theme(axis.title.y=element_text(size=15),strip.text = element_text(size = 13))
# Save plot 
ggsave(p2, filename = paste(dir,"plot_out/all_cell_across_treatment2.pdf", sep=""),width = 20, height = 6)


library(ggsignif)
compaired <- list(c("TN","MPR"),c("TN","NMPR"),c("MPR","NMPR"))
ps <- list()
for(i in unique(table.plot$cell)){
ps[[i]]<- ggplot(subset(table.plot,cell==i),aes(x=Group, y=Estimate,color=Group)) +
    geom_boxplot(width = 0.5,position=position_dodge(0.75),outlier.colour = NA)+
    geom_jitter(size = 1.5, position = position_jitterdodge())+
    scale_color_manual(name="Group",values = Group_color_panel)+
    theme_bw() +  
    theme(axis.text.x = element_text(size=13,angle = 45, hjust=1, vjust=1), axis.text.y = element_text(size=13),axis.ticks.length = unit(0.3,"cm")) + 
    labs(x="",y="Cellular fraction (%) of all cells",title=i) +
    theme(axis.title.y=element_text(size=15),strip.text = element_text(size = 13))+
    geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test)
}

pdf(paste0(dir,"plot_out/all_cell_fraction_each_cell.pdf"),height=5,width=4)
for(i in unique(table.plot$cell)){print(ps[[i]])}
dev.off()


# plot
p3<- ggplot(table.plot) +
    geom_boxplot(aes(x=cell, y=Estimate,color=Group,outlier.colour = NA),width = 0.6,position=position_dodge(0.75))+
    geom_jitter(aes(x=cell, y=Estimate,color=Group),size = 1.5, position = position_jitterdodge())+
    scale_color_manual(name="Group",values = Group_color_panel)+
    theme_bw() +  
    theme(axis.text.x = element_text(size=13), axis.text.y = element_text(size=13),axis.ticks.length = unit(0.3,"cm")) + 
    xlab("") + ylab("Cellular fraction (%) of all cells") +
    theme(axis.title.y=element_text(size=15),strip.text = element_text(size = 13))
# Save plot 
ggsave(p3, filename = paste(dir,"plot_out/all_cell_across_treatment2.pdf", sep=""),width = 14, height = 6)

# split plot
RNA@meta.data$Group <- factor(RNA@meta.data$Group,levels=c("TN","MPR","NMPR"))
p <- DimPlot(RNA, reduction = "umap", label = F, group.by = "all_cluster_annotation",split.by="Group") + NoLegend()
ggsave("plot_out/all_cell_group_split.png",p,width=9,height=4)


