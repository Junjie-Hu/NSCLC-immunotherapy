library(Seurat)
library(ggplot2)
library(future)
library(tidyverse)
plan("multiprocess", workers = 8)
options(future.globals.maxSize = 10*1024^3)

dir <- "./"

RNA <- readRDS("Tcell.rds")

# Manually ananotate immmune markers to identify major cell clusters
T_marker=c("CD8A","CD4","FCGR3A","CTLA4","PDCD1","TIGIT","HAVCR2","GZMA","GZMB","GZMK","NKG7","SELL","TCF7","FOXP3","ZNF683","ITGAE")
marker_plot=function(SeuratObj,marker){
    pn=length(marker)
    pp=list()
    for(i in 1:pn){
      pg=FeaturePlot(object = SeuratObj, features = marker[i], cols = c("grey", "red"), reduction = "umap")+NoLegend()+labs(x="",y="")+theme(plot.title = element_text(hjust = 0.5))
      pp[[marker[i]]]=pg
    }
    return(pp)
}
p=marker_plot(RNA,T_marker)
for(i in 1:16){
ggsave(paste0(dir,"plot_out/T_",T_marker[i],".png"),p[[T_marker[i]]],width =6,height =6)
}

# Manually ananotate T markers
# stash current cluster IDs
RNA[["T.cluster"]] <- Idents(object = RNA)
# enumerate current cluster IDs and the labels for them
cluster.ids <- 0:(length(unique(RNA@meta.data$T.cluster))-1)
# Annotate each of the clusters 
free_annotation <- c("CD8_GZMK","CD8_GZMB", "T_IL7R", "Treg_SELL","CD4_CCR7","CD4_MAF" ,"NK_FCGR3A", "CD8_HAVCR2", "CD4_CXCL13", "NK_KLRD1","T_MKI67","CD8_STMN1","Treg_CTLA4")
# Map free annotation to cluster numbers and store as all_subtype_annotation
RNA@meta.data[,'T_cluster_annotation'] <- plyr::mapvalues(x = RNA@meta.data$T.cluster, from = cluster.ids, to = free_annotation)
# plot umap
pdf(paste(dir,"plot_out/T_cluster_annotation_umap.pdf", sep=""),width = 7.4,height = 6)
DimPlot(RNA, reduction = "umap", label = F, group.by = 'T_cluster_annotation')
dev.off()

### read Tmem
Tmem <- readRDS("Tmem.rds")
meta <- Tmem@meta.data
CD8Tmem <- subset(meta,Tmem_cluster == "CD8_IL7R")
CD4Tmem <- subset(meta,Tmem_cluster == "CD4_IL7R")
CD8_ord <- which(rownames(RNA@meta.data) %in% rownames(CD8Tmem))
CD4_ord <- which(rownames(RNA@meta.data) %in% rownames(CD4Tmem))
head(CD8_ord);head(CD4_ord)
RNA$T_cluster_annotation <- as.character(RNA$T_cluster_annotation)
RNA$T_cluster_annotation[CD8_ord] <- "CD8_IL7R"
RNA$T_cluster_annotation[CD4_ord] <- "CD4_IL7R"
# plot umap
RNA$T_cluster_annotation <- factor(RNA$T_cluster_annotation,
    levels=c("CD8_GZMK","CD8_GZMB", "CD4_IL7R", "Treg_SELL","CD4_CCR7","CD4_MAF","CD8_IL7R","NK_FCGR3A","CD8_HAVCR2", "CD4_CXCL13", "NK_KLRD1","T_MKI67","CD8_STMN1","Treg_CTLA4"))
pdf(paste(dir,"plot_out/T_cluster_annotation_umap2.pdf", sep=""),width = 7.4,height = 6)
DimPlot(RNA, reduction = "umap", label = F, group.by = 'T_cluster_annotation')
dev.off()

saveRDS(RNA,"Tcell.rds")

###### Heatmap plot
features=c("CD8A","CD8B","CD4","FCGR3A","KLRD1","FOXP3","IL2RA","IKZF2","TCF7","SELL","LEF1","CCR7","LAG3","TIGIT","PDCD1","HAVCR2","CTLA4","LAYN","ENTPD1","GZMA","GZMB","GZMK","GNLY","IFNG","PRF1","NKG7",
  "CD28","ICOS","CD40LG","TNFRSF4","TNFRSF9","TNFRSF18","ZNF683","ITGAE","CCL3","CCL5","CXCL13","IL21","EOMES","MAF","TOX2","ID2","TBX21","HOPX","MKI67","STMN1")
row_split = c(rep("CD8/CD4",3),rep("NK",2),rep("Treg",3),rep("Naive",4),rep("Exhausted",7),rep("Cytotoxic",7),rep("Co-stimulatory",6),rep("TRM",2),rep("Chemokines",4),rep("TFs",6),rep("Proliferating",2))
row_split = factor(row_split,levels = c("CD8/CD4","NK","Treg","Naive","Exhausted","Cytotoxic","Co-stimulatory","TRM","Chemokines","TFs","Proliferating"))
### plot heatmap
source(paste0(dir,"code_clean/Heat_Dot_data.R"))
### set colnames order
plot_ord <- c("NK_FCGR3A","NK_KLRD1","CD8_IL7R","CD8_GZMK","CD8_GZMB","CD8_HAVCR2","CD8_STMN1","T_MKI67","CD4_IL7R","CD4_CCR7","CD4_MAF", "CD4_CXCL13","Treg_SELL","Treg_CTLA4")

data.plot <- Heat_Dot_data(object=RNA,features=features,group.by="T_cluster_annotation")
exp.mat <- data.plot %>% select(features.plot,id,avg.exp.scaled) %>% spread(id,avg.exp.scaled)
rownames(exp.mat) <- exp.mat$features.plot
exp.mat$features.plot <- NULL
exp.mat <- exp.mat[,plot_ord]
per.mat <- data.plot %>% select(features.plot,id,pct.exp) %>% spread(id,pct.exp)
rownames(per.mat) <- per.mat$features.plot
per.mat$features.plot <- NULL
per.mat <- per.mat[,plot_ord]/100

### plot heatmap
library(ComplexHeatmap)
library(circlize) ## color
# set color gradient
col_fun <- colorRamp2(c(-1.5, 0, 2.5), c("#118ab2", "#fdffb6", "#e63946"))
# split heatmap
col_split = c(rep("NK",2),rep("CD8 T",5),"Prol",rep("CD4 Tconv",4),rep("CD4 Treg",2))
col_split =factor(col_split,levels = c("NK","CD8 T","Prol","CD4 Tconv","CD4 Treg"))
# left annotation
annot = c("CD8/CD4","NK","Treg","Naive","Exhausted","Cytotoxic","Co-stimulatory","TRM","Chemokines","TFs","Proliferating")
row_color=c("#e76f51","#ffafcc","#0077b6","#ddbea9","#00b4d8","#dc2f02","#2a9d8f","#57cc99","#b5838d","#8a5a44","#023047")

ha = HeatmapAnnotation(df = data.frame(Marker=row_split),which = "row",
                       col = list(Marker = c("CD8/CD4" = "#e76f51", "NK" = "#ffafcc","Treg"="#0077b6","Naive"="#ddbea9",
                                           "Exhausted"="#00b4d8","Cytotoxic"="#dc2f02","Co-stimulatory"="#fca311","TRM"="#57cc99","Chemokines"="#b5838d","TFs"="#8a5a44","Proliferating"="#023047")))
pdf(paste0(dir,"plot_out/T_maker_heat_new.pdf"),width = 9,height = 13)
Heatmap(exp.mat, col = col_fun,cluster_columns = F,cluster_rows = F,
        show_column_names = T,show_row_names = T,rect_gp=gpar(col="grey"),
        column_names_side = "top",row_names_side = "right",
        row_split = row_split,column_split = col_split,
        row_gap = unit(3, "mm"),column_gap =unit(3, "mm"), 
        left_annotation = ha,
        heatmap_legend_param=list(title = "Expression",legend_height=unit(3, "cm")))
dev.off()

pdf(paste0(dir,"plot_out/T_maker_heat_new2.pdf"),width = 7,height = 15)
Heatmap(exp.mat, col = col_fun,cluster_columns = F,cluster_rows = F,
        show_column_names = T,show_row_names = T,rect_gp=gpar(type = "none"),
        cell_fun = function(j, i, x, y, width, height, fill){
          grid.rect(x = x, y = y, width = width, height = height,gp = gpar(col = "grey", fill = NA))
          grid.circle(x = x, y = y,r=per.mat[i,j]/2 * max(unit.c(width, height)),gp = gpar(fill = col_fun(exp.mat[i, j]), col = NA))},
        column_names_side = "top",row_names_side = "right",
        row_split = row_split,column_split = col_split,
        row_gap = unit(3, "mm"),column_gap =unit(3, "mm"),
        left_annotation = ha,
        heatmap_legend_param=list(title = "Expression",legend_height=unit(3, "cm")))
dev.off()

###### heatmap cluster each group
source(paste0(dir,"code_clean/Heat_Dot_split_group.R"))
plot_ord2 <- c()
for(i in plot_ord){
  for(j in c("TN","MPR","NMPR")){
    plot_ord2 <- c(plot_ord2,paste(i,j,sep="_"))
  }
}
data.plot2 <- Heat_Dot_split_group(object=RNA,features=features,group.by="T_cluster_annotation",split.by="Group")
exp.mat2 <- data.plot2 %>% select(features.plot,id,avg.exp.scaled) %>% spread(id,avg.exp.scaled)
rownames(exp.mat2) <- exp.mat2$features.plot
exp.mat2$features.plot <- NULL
exp.mat2 <- exp.mat2[,plot_ord2]
per.mat2 <- data.plot2 %>% select(features.plot,id,pct.exp) %>% spread(id,pct.exp)
rownames(per.mat2) <- per.mat2$features.plot
per.mat2$features.plot <- NULL
per.mat2 <- per.mat2[,plot_ord2]/100

library(ComplexHeatmap)
library(circlize) ## color 
# set color gradient
col_fun <- colorRamp2(c(-1.5, 0, 2.5), c("#118ab2", "#fdffb6", "#e63946"))
# left annotation
annot = c("CD8/CD4","NK","Treg","Naive","Exhausted","Cytotoxic","Co-stimulatory","TRM","Chemokines","TFs","Proliferating")
row_color=c("#e76f51","#ffafcc","#0077b6","#ddbea9","#00b4d8","#dc2f02","#2a9d8f","#57cc99","#b5838d","#8a5a44","#023047")

ha = HeatmapAnnotation(df = data.frame(Marker=row_split),which = "row",
                       col = list(Marker = c("CD8/CD4" = "#e76f51", "NK" = "#ffafcc","Treg"="#0077b6","Naive"="#ddbea9",
                                           "Exhausted"="#00b4d8","Cytotoxic"="#dc2f02","Co-stimulatory"="#fca311","TRM"="#57cc99","Chemokines"="#b5838d","TFs"="#8a5a44","Proliferating"="#023047")))
# top annotation
ha_top = HeatmapAnnotation(df = data.frame(Group = rep(c("TN","MPR","NMPR"),length(unique(RNA$T_cluster_annotation)))),
                       which = "column",
                       col = list(Group=c("TN"="#228B22","MPR"="#FF4500","NMPR"="#757bc8")))
pdf(paste0(dir,"plot_out/T_maker_heat_group.pdf"),width = 16,height = 15)
Heatmap(exp.mat2, col = col_fun,cluster_columns = F,cluster_rows = F,
        show_column_names = T,show_row_names = T,rect_gp=gpar(col="grey"),
        column_names_side = "top",row_names_side = "right",
        row_split = row_split,
        row_gap = unit(3, "mm"),
        left_annotation = ha,top_annotation =ha_top,
        heatmap_legend_param=list(title = "Expression",legend_height=unit(3, "cm")))
dev.off()

pdf(paste0(dir,"plot_out/T_maker_heat2_group.pdf"),width = 16,height = 15)
Heatmap(exp.mat2, col = col_fun,cluster_columns = F,cluster_rows = F,
        show_column_names = T,show_row_names = T,rect_gp=gpar(type = "none"),
        cell_fun = function(j, i, x, y, width, height, fil){
          grid.rect(x = x, y = y, width = width, height = height,gp = gpar(col = "grey", fill = NA))
          grid.circle(x = x, y = y,r=per.mat2[i,j]/2 * max(unit.c(width, height)),gp = gpar(fill = col_fun(exp.mat2[i, j]), col = NA))},
        column_names_side = "top",row_names_side = "right",
        row_split = row_split,
        row_gap = unit(3, "mm"),
        left_annotation = ha,top_annotation =ha_top,
        heatmap_legend_param=list(title = "Expression",legend_height=unit(3, "cm")))
dev.off()

###Calculate cell fractions in different response groups
source(paste0(dir,"code_clean/Cal_fraction.R"))
tab.1 <- cal_fraction(RNA,"T_cluster_annotation")
tab.1$cell <- factor(tab.1$cell, levels =c("NK_FCGR3A","NK_KLRD1","CD8_IL7R","CD8_GZMK","CD8_GZMB","CD8_HAVCR2","CD8_STMN1","T_MKI67","CD4_IL7R","CD4_CCR7","CD4_MAF", "CD4_CXCL13","Treg_SELL","Treg_CTLA4"))
Group_color_panel <- c("TN"="#228B22","MPR"="#FF4500","NMPR"="#757bc8")
# plot
p1<- ggplot(tab.1, aes(x=Group, y=Estimate)) +
    geom_bar(stat = "identity", aes(fill=Group)) +
    scale_fill_manual(name = "Group", values = Group_color_panel) + facet_grid(~ cell) +
    theme_bw() +  
    theme(axis.text.x = element_text(size=15,angle = 45, hjust=1, vjust=1), axis.text.y = element_text(size=15),axis.ticks.length = unit(0.2,"cm"),legend.position= "bottom") + 
    xlab("") + ylab("Cellular fraction (%) of all T/NK cell") +
    theme(axis.title.y=element_text(size=16),strip.text = element_text(size = 13))+
    geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,position=position_dodge(0.05))
# Save plot 
ggsave(p1, filename = paste0(dir,"plot_out/T_across_treatment.pdf"),width = 18, height = 6)

source(paste0(dir,"code_clean/Cal_fraction_patients.R"))
tab.2 <- cal_fraction_patients(RNA,"T_cluster_annotation")
tab.2$cell <- gsub("\\.[0-9]{1,}","",rownames(tab.2))
tab.2$cell <- factor(tab.2$cell, levels =c("NK_FCGR3A","NK_KLRD1","CD8_IL7R","CD8_GZMK","CD8_GZMB","CD8_HAVCR2","CD8_STMN1","T_MKI67","CD4_IL7R","CD4_CCR7","CD4_MAF","CD4_CXCL13","Treg_SELL","Treg_CTLA4"))
#colnames(tab.2)[4] <- "Patient_ID"
clin_info <- read.csv("clin_info.csv",header=T,stringsAsFactors=F)
sample_info <- clin_info[,c("Patient_ID","Group")]
table.plot <- merge(tab.2,sample_info,by="Patient_ID")
table.plot$Group <- factor(table.plot$Group,levels=c("TN","MPR","NMPR"))

pdf(paste0(dir,"plot_out/T_across_treatment3.pdf"),width = 18,height = 7)
    p1+geom_point(data=table.plot,aes(x=Group, y=Estimate,color=Patient_ID),size = 1.5)
dev.off()

# boxplot
p2<- ggplot(table.plot) +
    geom_boxplot(aes(x=cell, y=Estimate,color=Group,outlier.colour = NA),width = 0.6,position=position_dodge(0.75))+
    geom_jitter(aes(x=cell, y=Estimate,color=Group),size = 1.5, position = position_jitterdodge())+
    scale_color_manual(name="Group",values = Group_color_panel)+
    theme_bw() +  
    theme(axis.text.x = element_text(size=13,angle = 45, hjust=1, vjust=1), axis.text.y = element_text(size=13),axis.ticks.length = unit(0.2,"cm"),legend.position= "bottom") + 
    xlab("") + ylab("Cellular fraction (%) of T/NK cells") +
    theme(axis.title.y=element_text(size=15),strip.text = element_text(size = 13))
# Save plot 
ggsave(p2, filename = paste(dir,"plot_out/T_across_treatment2.pdf", sep=""),width = 18, height = 6)

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
    labs(x="",y="Cellular fraction (%) of all cells",title=i) +
    theme(axis.title.y=element_text(size=15),strip.text = element_text(size = 13))+
    geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test)
}

pdf(paste0(dir,"plot_out/T_fraction_each_cell.pdf"),height=5,width=4)
for(i in unique(table.plot$cell)){print(ps[[i]])}
dev.off()


saveRDS(RNA,"Tcell.rds")

####### CD8 T trajactory analysis ######
library(monocle)
# select CD8 T
CD8_tj<-subset(RNA, T_cluster_annotation %in% c("CD8_IL7R","CD8_GZMK","CD8_GZMB","CD8_HAVCR2"))
# Extract data, phenotype data, and feature data from the SeuratObject
data <- as(as.matrix(CD8_tj@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = CD8_tj@meta.data)
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
clustering_DEG_genes <- differentialGeneTest(monocle_cds, fullModelFormulaStr = '~T_cluster_annotation', cores = 20)
write.csv(clustering_DEG_genes,"data_out/T8_monocle_DEG.csv",quote=F)
ordering_genes <- row.names (subset(clustering_DEG_genes, qval < 0.01))
gene_id <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]
##
monocle_cds <- setOrderingFilter(monocle_cds, gene_id)
monocle_cds <- reduceDimension(
  monocle_cds,
  max_components = 2,
  method = 'DDRTree')
monocle_cds <- orderCells(monocle_cds)
ClusterName_color_panel <- c("CD8_GZMK" = "#f28482", "CD8_GZMB" = "#06d6a0", 
  "CD8_IL7R" = "#84a59d","CD8_HAVCR2"="#118ab2")
library(cowplot)
pdf(paste0(dir,"plot_out/T8_Pseudotime2.pdf"),width = 8,height = 5)
plot_cell_trajectory(monocle_cds, color_by = "Pseudotime",show_branch_points=F)+theme_cowplot()
plot_cell_trajectory(monocle_cds, color_by = "T_cluster_annotation",cell_size = 1,show_branch_points=F)+theme_cowplot()+
   scale_color_manual(name = "", values = ClusterName_color_panel)+guides(color = guide_legend(override.aes = list(size=5)))+
   theme(axis.text.y = element_text(size=15),axis.text.x = element_text(size=15),axis.title.y = element_text(size=18),axis.title.x = element_text(size=18),axis.ticks.length = unit(0.2,"cm"))
dev.off()

saveRDS(monocle_cds,"./data_out/T8_monocle2.rds")

######## CD8 Activating_Exhausted plot
# select CD8 T
CD8_T <-subset(RNA, T_cluster_annotation %in% c("CD8_IL7R","CD8_GZMK","CD8_GZMB","CD8_HAVCR2","CD8_STMN1"))
# define feature
act_features <- list(c("GZMA","GZMB","GZMK","GNLY","IFNG","PRF1","NKG7"))
exh_features <- list(c("LAG3","TIGIT","PDCD1","HAVCR2","CTLA4","LAYN","ENTPD1"))
# calculate score
CD8_T <- AddModuleScore(object=CD8_T,features=act_features,name="Cytotoxic_score")
CD8_T <- AddModuleScore(object=CD8_T,features=exh_features,name="Exhausted_score")
# plot
data.plot <- CD8_T@meta.data[,c("Group","T_cluster_annotation","Cytotoxic_score1","Exhausted_score1")]
library(ggsignif)
library(cowplot)
Group_color_panel <- c("TN"="#228B22","MPR"="#FF4500","NMPR"="#757bc8")
compaired <- list(c("TN","MPR"),c("TN","NMPR"),c("MPR","NMPR"))
p <- ggplot(data = data.plot,aes(x=factor(Group,levels=c("TN","MPR","NMPR")),y=Cytotoxic_score1,fill=Group))+
     geom_violin(color="white")+
     geom_boxplot(width=0.1,position=position_dodge(0.9),fill="white")+theme_cowplot()+
     scale_fill_manual(name="Group",values = Group_color_panel)+
     theme(axis.text.x = element_text(size=10,angle = 45, hjust=1, vjust=1), axis.text.y = element_text(size=10),axis.ticks.length = unit(0.2,"cm"))+
     NoLegend()+labs(x="",y="Cytotoxic score")+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test)
ggsave(paste0(dir,"plot_out/T8_cyto_across_treatment.pdf"),p,width =3,height =5)
p <- ggplot(data = data.plot,aes(x=factor(Group,levels=c("TN","MPR","NMPR")),y=Exhausted_score1,fill=Group))+
     geom_violin(color="white")+
     geom_boxplot(width=0.1,position=position_dodge(0.9),fill="white")+theme_cowplot()+
     scale_fill_manual(name="Group",values = Group_color_panel)+
     theme(axis.text.x = element_text(size=10,angle = 45, hjust=1, vjust=1), axis.text.y = element_text(size=10),axis.ticks.length = unit(0.2,"cm"))+
     NoLegend()+labs(x="",y="Exhausted score")+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test)
ggsave(paste0(dir,"plot_out/T8_exhau_across_treatment.pdf"),p,width =3,height =5)
#### each cluster
pc <- ggplot(data = data.plot,aes(x=factor(Group,levels=c("TN","MPR","NMPR")),y=Cytotoxic_score1,fill=Group))+
     geom_violin(color="white")+
     geom_boxplot(width=0.1,position=position_dodge(0.9),fill="white")+facet_grid(~T_cluster_annotation) +theme_bw()+
     scale_fill_manual(name="Group",values = Group_color_panel)+
     theme(axis.text.x = element_text(size=10,angle = 45, hjust=1, vjust=1), axis.text.y = element_text(size=10),axis.ticks.length = unit(0.2,"cm"),legend.position="bottom")+
     labs(x="",y="Cytotoxic score")+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test)
ggsave(paste0(dir,"plot_out/T8_cyto_cluster_across_treatment.pdf"),pc,width =10,height =5)
pe <- ggplot(data = data.plot,aes(x=factor(Group,levels=c("TN","MPR","NMPR")),y=Exhausted_score1,fill=Group))+
     geom_violin(color="white")+
     geom_boxplot(width=0.1,position=position_dodge(0.9),fill="white")+facet_grid(~T_cluster_annotation) +theme_bw()+
     scale_fill_manual(name="Group",values = Group_color_panel)+
     theme(axis.text.x = element_text(size=10,angle = 45, hjust=1, vjust=1), axis.text.y = element_text(size=10),axis.ticks.length = unit(0.2,"cm"),legend.position="bottom")+
     labs(x="",y="Exhausted score")+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test)
ggsave(paste0(dir,"plot_out/T8_exhau_cluster_across_treatment.pdf"),pe,width =10,height =5)

######## NK Activating_Exhausted plot
# select NK T
NK_T <-subset(RNA, idents = c(6,9))
# calculate score
NK_T <- AddModuleScore(object=NK_T,features=act_features,name="Cytotoxic_score")
NK_T <- AddModuleScore(object=NK_T,features=exh_features,name="Exhausted_score")
# plot
data.plot <- NK_T@meta.data[,c("Group","T_cluster_annotation","Cytotoxic_score1","Exhausted_score1")]
library(ggsignif)
compaired <- list(c("TN","MPR"),c("TN","NMPR"),c("MPR","NMPR"))
p <- ggplot(data = data.plot,aes(x=factor(Group,levels=c("TN","MPR","NMPR")),y=Cytotoxic_score1,fill=Group))+
     geom_violin(color="white")+
     geom_boxplot(width=0.1,position=position_dodge(0.9),fill="white")+theme_cowplot()+
     scale_fill_manual(name="Group",values = Group_color_panel)+
     theme(axis.text.x = element_text(size=10,angle = 45, hjust=1, vjust=1), axis.text.y = element_text(size=10),axis.ticks.length = unit(0.2,"cm"))+
     NoLegend()+labs(x="",y="Cytotoxic score")+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test)
ggsave(paste0(dir,"plot_out/NK_cyto_across_treatment.pdf"),p,width =3,height =5)
p <- ggplot(data = data.plot,aes(x=factor(Group,levels=c("TN","MPR","NMPR")),y=Exhausted_score1,fill=Group))+
     geom_violin(color="white")+
     geom_boxplot(width=0.1,position=position_dodge(0.9),fill="white")+theme_cowplot()+
     scale_fill_manual(name="Group",values = Group_color_panel)+
     theme(axis.text.x = element_text(size=10,angle = 45, hjust=1, vjust=1), axis.text.y = element_text(size=10),axis.ticks.length = unit(0.2,"cm"))+
     NoLegend()+labs(x="",y="Exhausted score")+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test)
ggsave(paste0(dir,"plot_out/NK_exhau_across_treatment.pdf"),p,width =3,height =5)
#### each cluster
pc <- ggplot(data = data.plot,aes(x=factor(Group,levels=c("TN","MPR","NMPR")),y=Cytotoxic_score1,fill=Group))+
     geom_violin(color="white")+
     geom_boxplot(width=0.1,position=position_dodge(0.9),fill="white")+facet_grid(~T_cluster_annotation) +theme_bw()+
     scale_fill_manual(name="Group",values = Group_color_panel)+
     theme(axis.text.x = element_text(size=10,angle = 45, hjust=1, vjust=1), axis.text.y = element_text(size=10),axis.ticks.length = unit(0.2,"cm"),legend.position="bottom")+
     labs(x="",y="Cytotoxic score")+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test)
ggsave(paste0(dir,"plot_out/NK_cyto_cluster_across_treatment.pdf"),pc,width =4.34,height =5)
pe <- ggplot(data = data.plot,aes(x=factor(Group,levels=c("TN","MPR","NMPR")),y=Exhausted_score1,fill=Group))+
     geom_violin(color="white")+
     geom_boxplot(width=0.1,position=position_dodge(0.9),fill="white")+facet_grid(~T_cluster_annotation) +theme_bw()+
     scale_fill_manual(name="Group",values = Group_color_panel)+
     theme(axis.text.x = element_text(size=10,angle = 45, hjust=1, vjust=1), axis.text.y = element_text(size=10),axis.ticks.length = unit(0.2,"cm"),legend.position="bottom")+
     labs(x="",y="Exhausted score")+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test)
ggsave(paste0(dir,"plot_out/NK_exhau_cluster_across_treatment.pdf"),pe,width =4.34,height =5)

#### memory T analysis
Tm <- subset(RNA,T_cluster_annotation=="CD8_IL7R")
## DEG after therapy
Idents(Tm) <- Tm@meta.data$Group
DEG<- FindAllMarkers(object = Tm, only.pos = F, logfc.threshold = 0,min.pct = 0.25)
DEG <- subset(DEG,cluster=="TN")
write.csv(DEG, file=paste0(dir,"data_out/T8_mem_DEG.csv"),quote=F)
#vocano plot
library(cowplot)
library(ggrepel)
label_gene <- c("NR4A1","NR4A2","NR4A3","TENT5C","CD8A","CCL5","GZMK","GZMA","CD74","HLA-DRA")
DEG$avg_logFC <- 0-DEG$avg_logFC
DEG$group <- ifelse(abs(DEG$avg_logFC) >= 0.3 & DEG$p_val_adj < 0.05,ifelse(DEG$avg_logFC >= 0.3,"Up","Down"),"Not")
DEG$label=NA
ord=match(label_gene,DEG$gene)
for(i in 1:length(ord)){DEG$label[ord[i]]=label_gene[i]}
gp <- ggplot(DEG) + 
  geom_point(aes(x = avg_logFC, y = -log10(p_val_adj),color = factor(group)), size = 1, alpha = 0.8, na.rm = TRUE) + # add gene points
  xlab(expression(log[2]("FC"))) + # x-axis label
  ylab(expression(-log[10]("FDR"))) + # y-axis label
  scale_color_manual(name = "", values = c("red", "blue", "grey"), limits = c("Up", "Down", "Not"))+
  geom_text_repel(aes(x = avg_logFC, y = -log10(p_val_adj),label=label),
                   arrow = arrow(length=unit(0.01, "npc")),
                   box.padding=unit(0.5, "lines"), point.padding=unit(0.5, "lines"),
                   segment.color = "black",show.legend = FALSE)+
  theme_cowplot() + 
  NoLegend()
ggsave(paste0(dir,"plot_out/T8_mem_DEG_vocano.pdf"),gp,width =5,height =6)

### Exhausted score in Treg
# select T3, T12
RNA_Treg <-subset(RNA, idents = c(3,12))
# calculate score
RNA_Treg <- AddModuleScore(object=RNA_Treg,features=act_features,name="Cytotoxic_score")
RNA_Treg <- AddModuleScore(object=RNA_Treg,features=exh_features,name="Exhausted_score")
# plot
data.plot <- RNA_Treg@meta.data[,c("Group","T_cluster_annotation","Cytotoxic_score1","Exhausted_score1")]
library(ggsignif)
compaired <- list(c("TN","MPR"),c("TN","NMPR"),c("MPR","NMPR"))
p <- ggplot(data = data.plot,aes(x=factor(Group,levels=c("TN","MPR","NMPR")),y=Exhausted_score1,fill=Group))+
     geom_violin(color="white")+
     geom_boxplot(width=0.1,position=position_dodge(0.9),fill="white")+theme_cowplot()+
     scale_fill_manual(name="Group",values = Group_color_panel)+
     theme(axis.text.x = element_text(size=10,angle = 45, hjust=1, vjust=1), axis.text.y = element_text(size=10),axis.ticks.length = unit(0.2,"cm"))+
     NoLegend()+labs(x="",y="Exhausted score")+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test)
ggsave(paste0(dir,"plot_out/Treg_exhau_across_treatment.pdf"),p,width =3,height =5)
pe <- ggplot(data = data.plot,aes(x=factor(Group,levels=c("TN","MPR","NMPR")),y=Exhausted_score1,fill=Group))+
     geom_violin(color="white")+
     geom_boxplot(width=0.1,position=position_dodge(0.9),fill="white")+facet_grid(~T_cluster_annotation) +theme_bw()+
     scale_fill_manual(name="Group",values = Group_color_panel)+
     theme(axis.text.x = element_text(size=10,angle = 45, hjust=1, vjust=1), axis.text.y = element_text(size=10),axis.ticks.length = unit(0.2,"cm"),legend.position="bottom")+
     labs(x="",y="Exhausted score")+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test)
ggsave(paste0(dir,"plot_out/Treg_exhau_cluster_across_treatment.pdf"),pe,width =4.34,height =5)

#### CD8_HAVCR2 analysis
T7 <- subset(RNA,idents=7)
p <- VlnPlot(T7,features=c("HAVCR2","TIGIT","GZMA","GZMB","NKG7","GNLY"), pt.size = 0,group.by="Group")
ggsave(paste0(dir,"plot_out/Tm_across_treatment.pdf"),p,width =5,height =3)
## DEG after therapy
Idents(T7) <- T7@meta.data$Group
DEG<- FindAllMarkers(object = T7, only.pos = F, logfc.threshold = 0,min.pct = 0.25)
DEG <- subset(DEG,cluster=="TN")
write.csv(DEG, file=paste0(dir,"data_out/T7_DEG.csv"),quote=F)
#vocano plot
library(cowplot)
library(ggrepel)
label_gene <- c("TNF","NR4A2","NR4A3","IL7R","LAG3","TIGIT","GZMA","GZMB","NKG7","GNLY","CXCL13","GZMH","PRF1")
DEG$avg_logFC <- 0-DEG$avg_logFC
DEG$group <- ifelse(abs(DEG$avg_logFC) >= 0.3 & DEG$p_val_adj < 0.05,ifelse(DEG$avg_logFC >= 0.3,"Up","Down"),"Not")
DEG$label=NA
ord=match(label_gene,DEG$gene)
for(i in 1:length(ord)){DEG$label[ord[i]]=label_gene[i]}
gp <- ggplot(DEG) + 
  geom_point(aes(x = avg_logFC, y = -log10(p_val_adj),color = factor(group)), size = 1, alpha = 0.8, na.rm = TRUE) + # add gene points
  xlab(expression(log[2]("FC"))) + # x-axis label
  ylab(expression(-log[10]("FDR"))) + # y-axis label
  scale_color_manual(name = "", values = c("red", "blue", "grey"), limits = c("Up", "Down", "Not"))+
  geom_text_repel(aes(x = avg_logFC, y = -log10(p_val_adj),label=label),
                   arrow = arrow(length=unit(0.01, "npc")),
                   box.padding=unit(0.5, "lines"), point.padding=unit(0.5, "lines"),
                   segment.color = "black",show.legend = FALSE)+
  theme_cowplot() + 
  NoLegend()
ggsave(paste0(dir,"plot_out/T7_DEG_vocano.pdf"),gp,width =5,height =6)

### CD8_IL7R density
T8_cds <- readRDS("./data_out/T8_monocle2.rds")
T8mem <- subset(RNA,T_cluster_annotation=="CD8_IL7R")
dense_cds <- T8_cds[,colnames(T8mem)]
plotdf=pData(dense_cds)
plotdf$Pseudotime <- 15 - plotdf$Pseudotime
plotdf$Pseudotime[plotdf$Pseudotime > 10] <- 8
plotdf$Group <- factor(plotdf$Group,levels=c("NMPR","MPR","TN"))
library(ggridges)
ggplot(plotdf, aes(x=Pseudotime,y=Group,fill=Group))+
  geom_density_ridges(scale=1) +
  geom_vline(xintercept = c(3.2,6),linetype=2)+
  scale_y_discrete("")+
  theme_minimal()+
  theme(
    panel.grid = element_blank()
  )
ggsave("plot_out/T8mem_pseudo_dense.pdf",width = 8,height = 3)
## fraction
plotdf$phase <- ifelse(plotdf$Pseudotime > 3.2,
    ifelse(plotdf$Pseudotime > 6,"phase1","phase2"),"phase0")
table(plotdf$phase,plotdf$Group)

### heatmap
plotdf$DEG <- ifelse(plotdf$Pseudotime > 3.2,"trans","naive")
T8mem@meta.data$DEG <- plotdf$DEG
Idents(T8mem) <- T8mem@meta.data$DEG
markers <- FindAllMarkers(object = T8mem, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.csv(markers, file=paste0(dir,"data_out/T8mem_pseudo_DEG.csv"),quote=F)

plot_gene <- c("NR4A1","NR4A2","NR4A3","FOSB","JUND","TENT5C","CCL5","GZMA","NKG7","GNLY","CD74","HLA-DRA","HLA-DPB1")
heat_cds <- dense_cds[plot_gene,]
plot_pseudotime_heatmap(heat_cds,
                        num_clusters = 2,
                        cores = 1,
                        show_rownames = T)


