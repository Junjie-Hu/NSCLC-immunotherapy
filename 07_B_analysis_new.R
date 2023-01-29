library(Seurat)
library(ggplot2)
library(tidyverse)

dir <- "./"

RNA <- readRDS("Bcell.rds")


# Manually ananotate immmune markers to identify major cell clusters
B_marker=c("MS4A1","CD27","GPR183","IGHD","RGS13","IGHM")
marker_plot=function(SeuratObj,marker){
    pn=length(marker)
    pp=list()
    for(i in 1:pn){
      pg=FeaturePlot(object = SeuratObj, features = marker[i], cols = c("grey", "red"), reduction = "umap")+NoLegend()+labs(x="",y="")+theme(plot.title = element_text(hjust = 0.5))
      pp[[marker[i]]]=pg
    }
    return(pp)
}
p=marker_plot(RNA,B_marker)
for(i in 1:6){
ggsave(paste0(dir,"plot_out/B_",B_marker[i],".png"),p[[B_marker[i]]],width =6,height =6)
}

# Dotplot of top N DE expressed genes 
markers<- read.csv(paste0(dir,"/data_out/B_clusters.csv"),row.names=1,stringsAsFactors=F)
markers.small  <- markers %>% group_by(cluster) %>% top_n(5, avg_logFC)
genes_to_check <- unique(markers.small$gene)
# Create Dotplot 
pdf(paste(dir,"plot_out/B_dotplot.pdf", sep=""),6,8)
DotPlot(RNA, features = genes_to_check) +coord_flip()+ theme(axis.title.y = element_text(size=15),axis.title.x = element_text(size=15))
dev.off()

# Manually ananotate B markers
# stash current cluster IDs
RNA[["B.cluster"]] <- Idents(object = RNA)
# enumerate current cluster IDs and the labels for them
cluster.ids <- 0:(length(unique(RNA@meta.data$B.cluster))-1)
# Annotate each of the clusters 
free_annotation <- c("Memory B", "Memory B", "Memory B","Naive B","Memory B","Memory B","GC B")              
# Map free annotation to cluster numbers and store as all_subtype_annotation
RNA@meta.data[,'B_cluster_annotation'] <- plyr::mapvalues(x = RNA@meta.data$B.cluster, from = cluster.ids, to = free_annotation)

# Annotate each of the clusters 
free_annotation2 <- c("B0-MS4A1", "B1-IGHM", "B2-HSPA1A","B3-IGHD","B4-FCRL4","B5-CD83","B6-RGS13")              
# Map free annotation to cluster numbers and store as all_subtype_annotation
RNA@meta.data[,'B_cluster_annotation2'] <- plyr::mapvalues(x = RNA@meta.data$B.cluster, from = cluster.ids, to = free_annotation2)
# plot umap
pdf(paste(dir,"plot_out/B_cluster_annotation_umap.pdf", sep=""),width = 7,height = 6)
DimPlot(RNA, reduction = "umap", label = F, group.by = 'B_cluster_annotation')
dev.off()
pdf(paste(dir,"plot_out/B_cluster_annotation_umap2.pdf", sep=""),width = 7,height = 6)
DimPlot(RNA, reduction = "umap", label = F, group.by = 'B_cluster_annotation2')
dev.off()

saveRDS(RNA,"Bcell.rds")


###Calculate cell fractions in different response groups
source(paste0(dir,"code/Cal_fraction.R"))
tab.1 <- cal_fraction(RNA,"B_cluster_annotation2")
tab.1$cell <- factor(tab.1$cell,levels =c("B3 IGHD","B0 MS4A1", "B1 IGHM", "B2 HSPA1A","B4 FCRL4","B5 CD83","B6 RGS13"))
Group_color_panel <- c("TN"="#228B22","MPR"="#FF4500","NMPR"="#757bc8")
# plot
p1<- ggplot(tab.1, aes(x=Group, y=Estimate)) +
    geom_bar(stat = "identity", aes(fill=Group)) +
    scale_fill_manual(name = "Group", values = Group_color_panel) + facet_grid(~ cell) +
    theme_bw() +  
    theme(axis.text.x = element_text(size=15,angle = 45, hjust=1, vjust=1), axis.text.y = element_text(size=15),axis.ticks.length = unit(0.2,"cm"),legend.position= "bottom") + 
    xlab("") + ylab("Cellular fraction (%) of B cell") +
    theme(axis.title.y=element_text(size=16),strip.text = element_text(size = 13))+
    geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,position=position_dodge(0.05))
ggsave(p1, filename = paste0(dir,"plot_out/B_across_treatment.pdf"),width = 9, height = 5.5)


# boxplot
p2<- ggplot(table.plot) +
    geom_boxplot(aes(x=cell, y=Estimate,color=Group,outlier.colour = NA),width = 0.6,position=position_dodge(0.75))+
    geom_jitter(aes(x=cell, y=Estimate,color=Group),size = 1.5, position = position_jitterdodge())+
    scale_color_manual(name="Group",values = Group_color_panel)+
    theme_bw() +  
    theme(axis.text.x = element_text(size=13,angle = 45, hjust=1, vjust=1), axis.text.y = element_text(size=13),axis.ticks.length = unit(0.2,"cm"),legend.position= "bottom") + 
    xlab("") + ylab("Cellular fraction (%) of B cells") +
    theme(axis.title.y=element_text(size=15),strip.text = element_text(size = 13))
# Save plot 
ggsave(p2, filename = paste0(dir,"plot_out/B_across_treatment2.pdf"),width = 9, height = 5.5)

library(ggsignif)
compaired <- list(c("TN","MPR"),c("TN","NMPR"),c("MPR","NMPR"))
ps <- list()
for(i in unique(table.plot$cell)){
ps[[i]]<- ggplot(subset(table.plot,cell==i),aes(x=Group, y=Estimate,color=Group)) +
    geom_boxplot(width = 0.5,position=position_dodge(0.75),outlier.colour = NA)+
    geom_jitter(size = 1.5, position = position_jitterdodge())+
    scale_color_manual(name="Group",values = Group_color_panel)+
    theme_bw() +  
    theme(axis.text.x = element_text(size=13,angle = 45, hjust=1, vjust=1), axis.text.y = element_text(size=13),axis.ticks.length = unit(0.2,"cm")) + 
    labs(x="",y="Cellular fraction (%) of B cells",title=i) +
    theme(axis.title.y=element_text(size=15),strip.text = element_text(size = 13))+
    geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test)
}

pdf(paste0(dir,"plot_out/B_fraction_each_cell.pdf"),height=5,width=4)
for(i in unique(table.plot$cell)){print(ps[[i]])}
dev.off()

###Calculate cell fractions in major clusters
source(paste0(dir,"code/Cal_fraction.R"))
tab.1 <- cal_fraction(RNA,"B_cluster_annotation")
tab.1$cell <- factor(tab.1$cell,levels =c("Memory B","Naive B","GC B"))
Group_color_panel <- c("TN"="#228B22","MPR"="#FF4500","NMPR"="#757bc8")
# plot
p1<- ggplot(tab.1, aes(x=Group, y=Estimate)) +
    geom_bar(stat = "identity", aes(fill=Group)) +
    scale_fill_manual(name = "Group", values = Group_color_panel) + facet_grid(~ cell) +
    theme_bw() +  
    theme(axis.text.x = element_text(size=15,angle = 45, hjust=1, vjust=1), axis.text.y = element_text(size=15),axis.ticks.length = unit(0.2,"cm"),legend.position= "bottom") + 
    xlab("") + ylab("Cellular fraction (%) of B cell") +
    theme(axis.title.y=element_text(size=16),strip.text = element_text(size = 13))+
    geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,position=position_dodge(0.05))

source(paste0(dir,"code_clean/Cal_fraction_patients.R"))
tab.2 <- cal_fraction_patients(RNA_sub,"B_cluster_annotation")
tab.2$cell <- gsub("\\.[0-9]{1,}","",rownames(tab.2))
tab.2$cell <- gsub("\\."," ",tab.2$cell)
tail(tab.2)
tab.2$cell <- factor(tab.2$cell, levels =c("Memory B","Naive B","GC B"))
#colnames(tab.2)[4] <- "Patient_ID"
clin_info <- read.csv("clin_info.csv",header=T,stringsAsFactors=F)
sample_info <- clin_info[,c("Patient_ID","Group")]
table.plot <- merge(tab.2,sample_info,by="Patient_ID")
table.plot$Group <- factor(table.plot$Group,levels=c("TN","MPR","NMPR"))

pdf(paste0(dir,"plot_out/B_major_across_treatment.pdf"),width = 4,height = 5)
    p1+geom_point(data=table.plot,aes(x=Group, y=Estimate,color=Patient_ID),size = 1.5)
dev.off()

# boxplot
p2<- ggplot(table.plot) +
    geom_boxplot(aes(x=cell, y=Estimate,color=Group,outlier.colour = NA),width = 0.6,position=position_dodge(0.75))+
    geom_jitter(aes(x=cell, y=Estimate,color=Group),size = 1.5, position = position_jitterdodge())+
    scale_color_manual(name="Group",values = Group_color_panel)+
    theme_bw() +  
    theme(axis.text.x = element_text(size=13,angle = 45, hjust=1, vjust=1), axis.text.y = element_text(size=13),axis.ticks.length = unit(0.2,"cm"),legend.position= "bottom") + 
    xlab("") + ylab("Cellular fraction (%) of B cells") +
    theme(axis.title.y=element_text(size=15),strip.text = element_text(size = 13))
# Save plot 
ggsave(p2, filename = paste(dir,"plot_out/B_across_treatment2.pdf", sep=""),width = 4, height = 5)

library(ggsignif)
compaired <- list(c("TN","MPR"),c("TN","NMPR"),c("MPR","NMPR"))
ps <- list()
for(i in unique(table.plot$cell)){
ps[[i]]<- ggplot(subset(table.plot,cell==i),aes(x=Group, y=Estimate,color=Group)) +
    geom_boxplot(width = 0.5,position=position_dodge(0.75),outlier.colour = NA)+
    geom_jitter(size = 1.5, position = position_jitterdodge())+
    scale_color_manual(name="Group",values = Group_color_panel)+
    theme_bw() +  
    theme(axis.text.x = element_text(size=13,angle = 45, hjust=1, vjust=1), axis.text.y = element_text(size=13),axis.ticks.length = unit(0.2,"cm")) + 
    labs(x="",y="Cellular fraction (%) of B cells",title=i) +
    theme(axis.title.y=element_text(size=15),strip.text = element_text(size = 13))+
    geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test)
}

pdf(paste0(dir,"plot_out/B_fraction_each_major.pdf"),height=5,width=4)
for(i in unique(table.plot$cell)){print(ps[[i]])}
dev.off()

###### B Trajectory
library(monocle)
# select B
B_tj<-subset(RNA, idents = c(0:5))
# Extract data, phenotype data, and feature data from the SeuratObject
data <- as(as.matrix(B_tj@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = B_tj@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
# Construct monocle cds
monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())

monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)
##
B_DEG_genes <- differentialGeneTest(monocle_cds, fullModelFormulaStr = '~B_cluster_annotation', cores = 30)
ordering_genes <- row.names (subset(B_DEG_genes, qval < 0.01))
gene_check <- row.names(B_DEG_genes)[order(B_DEG_genes$qval)][1:400]
#gene_id <- c()
monocle_cds <- setOrderingFilter(monocle_cds,gene_check)
monocle_cds <- reduceDimension( monocle_cds,
  max_components = 2,
  method = 'DDRTree')
monocle_cds <- orderCells(monocle_cds)
ClusterName_color_panel <- c("B0" = "#f28482","B1-IGHM" = "#118ab2","B2-HSPA1A" = "#84a59d","B3-IGHD"="#06d6a0","B4-FCRL4"="#e07a5f","B5-CD83"="#cdb4db")
#"B6-RGS13"="#006d77"
library(cowplot)
pdf(paste0(dir,"plot_out/B_Pseudotime.pdf"),width = 8,height = 6)
plot_cell_trajectory(monocle_cds, color_by = "Pseudotime",show_branch_points=F)+theme_cowplot()
plot_cell_trajectory(monocle_cds, color_by = "B_cluster_annotation2",cell_size = 1,show_branch_points=F)+theme_cowplot()+
   scale_color_manual(name = "",values = ClusterName_color_panel)+guides(color = guide_legend(override.aes = list(size=5)))+
   theme(axis.text.y = element_text(size=15),axis.text.x = element_text(size=15),axis.title.y = element_text(size=18),axis.title.x = element_text(size=18),axis.ticks.length = unit(0.3,"cm"))
dev.off()

saveRDS(monocle_cds,"./data_out/B_monocle.rds")



