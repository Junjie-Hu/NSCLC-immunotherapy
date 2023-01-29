library(Seurat)
library(ggplot2)
library(tidyverse)

dir <- "./"
RNA <- readRDS("Neutrophil2_new.rds")

###### Heatmap plot
features=c("FCGR3B","SELL","CXCR4","CXCR2","CXCR1","CXCL8","IL1B","CCL3","CCL4","CCL4L2","S100A12","S100A9","S100A8","CYBB","ELANE","MMP9","PADI4","HMGB1","TNF","CXCL9","CXCL10","OSM","PTGS2",
  "ARG1","TGFB1","VEGFA","PROK2","IFIT1","IFIT2","IFIT3","RSAD2","CD274","IDO1")
row_split = c(rep("Receptors",5),rep("Chemokines",5),rep("Granules",6),rep("NETs",2),rep("Pro-inflammatory",5),rep("Anti-inflammatory",4),rep("ISGs",4),rep("Checkpoints",2))
row_split = factor(row_split,levels = c("Receptors","Chemokines","Granules","NETs","Pro-inflammatory","Anti-inflammatory","ISGs","Checkpoints"))
### plot heatmap
source(paste0(dir,"code_clean/Heat_Dot_data.R"))
### set colnames order
plot_ord <- c("Neu_S100A12","Neu_OSM","Neu_IFIT3","Neu_CCL3")

data.plot <- Heat_Dot_data(object=RNA,features=features,group.by="N_cluster_annotation")
exp.mat <- data.plot %>% select(features.plot,id,avg.exp.scaled) %>% spread(id,avg.exp.scaled)
rownames(exp.mat) <- exp.mat$features.plot
exp.mat$features.plot <- NULL
exp.mat <- exp.mat[,plot_ord]
per.mat <- data.plot %>% select(features.plot,id,pct.exp) %>% spread(id,pct.exp)
rownames(per.mat) <- per.mat$features.plot
per.mat$features.plot <- NULL
per.mat <- per.mat[,plot_ord]/100
min(exp.mat);max(exp.mat)

### plot heatmap
library(ComplexHeatmap)
library(circlize) ## color 
# set color gradient
col_fun <- colorRamp2(c(-1.5, 0, 1.5), c("#118ab2", "#fdffb6", "#e63946"))
# left annotation
annot = c("Receptors","Chemokines","Granules","NETs","Pro-inflammatory","Anti-inflammatory","ISGs","Checkpoints")

ha = HeatmapAnnotation(df = data.frame(Marker=row_split),which = "row",
                       col = list(Marker = c("Receptors" = "#e76f51", "Chemokines"="#0077b6","Granules"="#ddbea9",
                                           "NETs"="#00b4d8","Pro-inflammatory"="#dc2f02","Anti-inflammatory"="#f4a261",
                                           "ISGs"="#57cc99","Checkpoints"="#b5838d")))
pdf(paste0(dir,"plot_out/N_maker_heat_new_new.pdf"),width = 5,height = 11)
Heatmap(exp.mat, col = col_fun,cluster_columns = F,cluster_rows = F,
        show_column_names = T,show_row_names = T,rect_gp=gpar(col="grey"),
        column_names_side = "top",row_names_side = "right",
        row_split = row_split,
        row_gap = unit(3, "mm"),
        left_annotation = ha,
        heatmap_legend_param=list(title = "Expression",legend_height=unit(3, "cm")))
dev.off()

pdf(paste0(dir,"plot_out/N_maker_heat_new2_new.pdf"),width = 4.3,height = 14)
Heatmap(exp.mat, col = col_fun,cluster_columns = F,cluster_rows = F,
        show_column_names = T,show_row_names = T,rect_gp=gpar(type = "none"),
        cell_fun = function(j, i, x, y, width, height, fill){
          grid.rect(x = x, y = y, width = width, height = height,gp = gpar(col = "grey", fill = NA))
          grid.circle(x = x, y = y,r=per.mat[i,j]/2 * max(unit.c(width, height)),gp = gpar(fill = col_fun(exp.mat[i, j]), col = NA))},
        column_names_side = "top",row_names_side = "right",
        row_split = row_split,
        row_gap = unit(3, "mm"),
        left_annotation = ha,
        heatmap_legend_param=list(title = "Expression",legend_height=unit(3, "cm")))
dev.off()

###Calculate cell fractions in different response groups
source(paste0(dir,"code_clean/Cal_fraction.R"))
tab.1 <- cal_fraction(RNA,"N_cluster_annotation")
tab.1$cell <- factor(tab.1$cell, levels=c("Neu_OSM","Neu_S100A12","Neu_CCL3","Neu_IFIT3"))
Group_color_panel <- c("TN"="#228B22","MPR"="#FF4500","NMPR"="#757bc8")
# plot
p1<- ggplot(tab.1, aes(x=Group, y=Estimate)) +
    geom_bar(stat = "identity", aes(fill=Group)) +
    scale_fill_manual(name = "Group", values = Group_color_panel) + facet_grid(~ cell) +
    theme_bw() +  
    theme(axis.text.x = element_text(size=15,angle = 45, hjust=1, vjust=1), axis.text.y = element_text(size=15),axis.ticks.length = unit(0.2,"cm"),legend.position= "bottom") + 
    xlab("") + ylab("Cellular fraction (%) of neutrophil") +
    theme(axis.title.y=element_text(size=16),strip.text = element_text(size = 13))+
    geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,position=position_dodge(0.05))
# Save plot 
ggsave(p1, filename = paste0(dir,"plot_out/N_across_treatment_new.pdf"),width = 6, height = 5)

source(paste0(dir,"code_clean/Cal_fraction_patients.R"))
tab.2 <- cal_fraction_patients(RNA,"N_cluster_annotation")
tab.2$cell <- gsub("\\.[0-9]{1,}","",rownames(tab.2))
tab.2$cell <- gsub("\\."," ",tab.2$cell)
tab.2$cell <- factor(tab.2$cell, levels =c("Neu_OSM","Neu_S100A12","Neu_CCL3","Neu_IFIT3"))
#colnames(tab.2)[4] <- "Patient_ID"
clin_info <- read.csv("clin_info.csv",header=T,stringsAsFactors=F)
sample_info <- clin_info[,c("Patient_ID","Group")]
table.plot <- merge(tab.2,sample_info,by="Patient_ID")
table.plot$Group <- factor(table.plot$Group,levels=c("TN","MPR","NMPR"))

pdf(paste0(dir,"plot_out/N_across_treatment3_new.pdf"),width =6,height = 5)
    p1+geom_point(data=table.plot,aes(x=Group, y=Estimate,color=Patient_ID),size = 1.5)
dev.off()

# boxplot
p2<- ggplot(table.plot) +
    geom_boxplot(aes(x=cell, y=Estimate,color=Group,outlier.colour = NA),width = 0.6,position=position_dodge(0.75))+
    geom_jitter(aes(x=cell, y=Estimate,color=Group),size = 1.5, position = position_jitterdodge())+
    scale_color_manual(name="Group",values = Group_color_panel)+
    theme_bw() +  
    theme(axis.text.x = element_text(size=13,angle = 45, hjust=1, vjust=1), axis.text.y = element_text(size=13),axis.ticks.length = unit(0.2,"cm"),legend.position= "bottom") + 
    xlab("") + ylab("Cellular fraction (%) of neutrophil") +
    theme(axis.title.y=element_text(size=15),strip.text = element_text(size = 13))
# Save plot 
ggsave(p2, filename = paste0(dir,"plot_out/N_across_treatment2_new.pdf"),width = 6, height = 5)

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

pdf(paste0(dir,"plot_out/N_fraction_each_cell_new.pdf"),height=5,width=4)
for(i in unique(table.plot$cell)){print(ps[[i]])}
dev.off()

####### trajactory analysis ######
library(monocle)
data <- as(as.matrix(RNA@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = RNA@meta.data)
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
monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(monocle_cds),num_cells_expressed >= 10))
Neu_DEG_genes <- differentialGeneTest(monocle_cds[expressed_genes,], fullModelFormulaStr = '~N_cluster_annotation', cores = 20)
write.csv(Neu_DEG_genes,"data_out/Neu_monocle_DEG2.csv",quote=F)
ordering_genes <- row.names (subset(Neu_DEG_genes, qval < 0.001))
length(ordering_genes)
gene_id <- row.names(Neu_DEG_genes)[order(Neu_DEG_genes$qval)][1:1000]
##
monocle_cds <- setOrderingFilter(monocle_cds, gene_id)
monocle_cds <- reduceDimension(
  monocle_cds,
  max_components = 2,
  method = 'DDRTree')
monocle_cds <- orderCells(monocle_cds)
ClusterName_color_panel <- c("Neu_OSM" = "#f8766d", "Neu_S100A12" = "#7cae00", 
  "Neu_CCL3" = "#00bfc4","Neu_IFIT3"="#c77cff")
library(cowplot)
pdf(paste0(dir,"plot_out/Neu_Pseudotime_new.pdf"),width = 8,height = 5)
plot_cell_trajectory(monocle_cds, color_by = "Pseudotime",show_branch_points=F)+theme_cowplot()
plot_cell_trajectory(monocle_cds, color_by = "N_cluster_annotation",cell_size = 1,show_branch_points=F)+theme_cowplot()+
   scale_color_manual(name = "", values = ClusterName_color_panel)+guides(color = guide_legend(override.aes = list(size=5)))+
   theme(axis.text.y = element_text(size=15),axis.text.x = element_text(size=15),axis.title.y = element_text(size=18),axis.title.x = element_text(size=18),axis.ticks.length = unit(0.2,"cm"))
dev.off()

saveRDS(monocle_cds,"./data_out/Neu_monocle_new.rds")

plot_genes <- c("CXCR4","CXCR2")
cds_subset <- monocle_cds[plot_genes,]
plot_genes_in_pseudotime(cds_subset, color_by = "Group",ncol=2)
source("code_clean/Monocle_plot_gene.R")
monocle_theme_opts <- function () {
  theme(strip.background = element_rect(colour = "white", 
    fill = "white")) + theme(panel.border = element_blank()) + 
    theme(axis.line.x = element_line(size = 0.25, color = "black")) + 
    theme(axis.line.y = element_line(size = 0.25, color = "black")) + 
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) + 
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + 
    theme(panel.background = element_rect(fill = "white")) + 
    theme(legend.key = element_blank())
}
pdf(paste0(dir,"plot_out/N_pseudo_CXCR4_new.pdf"),width = 4,height = 8)
My_plot_gene_pseudotime(cds_subset, color_by = "Group",ncol=1)
dev.off()

### density
plotdf=pData(monocle_cds)
library(ggridges)
ggplot(plotdf, aes(x=Pseudotime,y=N_cluster_annotation,fill=N_cluster_annotation))+
  geom_density_ridges(scale=1) +
  geom_vline(xintercept = c(5,10),linetype=2)+
  scale_y_discrete("")+
  theme_minimal()+
  theme(
    panel.grid = element_blank()
  )
ggsave("plot_out/N_monocle_cluster_density.pdf",width = 6,height = 4)

ggplot(plotdf, aes(x=Pseudotime,y=Group,fill=Group))+
  geom_density_ridges(scale=1) +
  geom_vline(xintercept = c(5,10),linetype=2)+
  scale_y_discrete("")+
  theme_minimal()+
  theme(
    panel.grid = element_blank()
  )
ggsave("plot_out/N_monocle_Group_density.pdf",width =5,height = 4)

# CCL3 and CCL4 
N1 <- subset(RNA,N_cluster_annotation=="Neu_CCL3")
N1.gene <- FetchData(object = N1, vars = c("CCL3","CCL4","CXCL8","VEGFA"))
N1.gene$Group <- N1$Group
library(ggsignif)
Group_color_panel <- c("TN"="#228B22","MPR"="#FF4500","NMPR"="#757bc8")
compaired <- list(c("TN","MPR"),c("TN","NMPR"),c("MPR","NMPR"))
p <- ggplot(data=N1.gene,aes(x=factor(Group,levels=c("TN","MPR","NMPR")),y=CXCL8,fill=Group)) + 
     geom_violin(color="white")+
     geom_boxplot(width=0.1,position=position_dodge(0.9),fill="white") + theme_cowplot()+
     scale_fill_manual(name = "Group", values = Group_color_panel)+
     theme(axis.text.x = element_text(size=10,angle = 45, hjust=1, vjust=1), axis.text.y = element_text(size=10),axis.ticks.length = unit(0.2,"cm"))+
     NoLegend()+labs(x="")+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test)
ggsave(paste0(dir,"plot_out/N1_CXCL8_treatment_new.pdf"),p,width =3,height =4)

#### branch heatmap
BEAM_res <- BEAM(monocle_cds, branch_point = 1, cores = 1)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
tmp1=plot_genes_branched_heatmap(monocle_cds[row.names(subset(BEAM_res,
                                                  qval < 1e-6)),],
                            branch_point = 1,
                            num_clusters = 3,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = F,
                           return_heatmap = T)
pdf(paste0(dir,"plot_out/N_monocle_brach_heatmap.pdf"),width=5,height=7)
tmp1$ph_res
dev.off()
## GO-BP analysis 
gene_group <- tmp1$annotation_row
gene_group$gene <- rownames(gene_group)
library(clusterProfiler)
library(org.Hs.eg.db)
allcluster_go <- data.frame()
for(i in unique(gene_group$Cluster)){
    samll_gene_group= filter(gene_group,gene_group$Cluster==i)
    df_name = bitr(samll_gene_group$gene,fromType="SYMBOL",toType=c("ENTREZID"),OrgDb="org.Hs.eg.db")
    go <- enrichGO(gene=unique(df_name$ENTREZID),
        OrgDb=org.Hs.eg.db,
        keyType="ENTREZID",
        ont="BP",
        pAdjustMethod="BH",
        pvalueCutoff=0.05,
        qvalueCutoff=0.05,
        readable=T)
    go_res=go@result
    if(dim(go_res)[1] !=0){
        go_res$cluster=i
        allcluster_go=rbind(allcluster_go,go_res)
    }
}
write.csv(allcluster_go,"data_out/N_monocle_brach_BP.csv")

###### pseudo-time heatmap
heat_cds <- monocle_cds[gene_id ,]
tmp2=plot_pseudotime_heatmap(heat_cds,
                        num_clusters = 4,
                        cores = 1,
                        show_rownames = F,
                       return_heatmap = T)
pdf(paste0(dir,"plot_out/N_monocle_pseudotime_heatmap.pdf"),width=5,height=7)
print(tmp2)
dev.off()
## GO-BP analysis 
gene_group <- data.frame(Cluster = factor(cutree(tmp2$tree_row,4)))
gene_group$gene <- rownames(gene_group)
head(gene_group)
library(clusterProfiler)
library(org.Hs.eg.db)
allcluster_go <- data.frame()
for(i in unique(gene_group$Cluster)){
    samll_gene_group= filter(gene_group,gene_group$Cluster==i)
    df_name = bitr(samll_gene_group$gene,fromType="SYMBOL",toType=c("ENTREZID"),OrgDb="org.Hs.eg.db")
    go <- enrichGO(gene=unique(df_name$ENTREZID),
        OrgDb=org.Hs.eg.db,
        keyType="ENTREZID",
        ont="BP",
        pAdjustMethod="BH",
        pvalueCutoff=0.05,
        qvalueCutoff=0.05,
        readable=T)
    go_res=go@result
    if(dim(go_res)[1] !=0){
        go_res$cluster=i
        allcluster_go=rbind(allcluster_go,go_res)
    }
}
write.csv(allcluster_go,"data_out/N_monocle_pseudotime_BP.csv")

###Calculate fraction of major population in different response groups
source(paste0(dir,"code_clean/Cal_fraction.R"))
tab.1 <- cal_fraction(RNA,"N_major_cluster")
tab.1$cell <- factor(tab.1$cell, levels=c("mature","aged"))
Group_color_panel <- c("TN"="#228B22","MPR"="#FF4500","NMPR"="#757bc8")
# plot
p1<- ggplot(tab.1, aes(x=Group, y=Estimate)) +
    geom_bar(stat = "identity", aes(fill=Group)) +
    scale_fill_manual(name = "Group", values = Group_color_panel) + facet_grid(~ cell) +
    theme_bw() +  
    theme(axis.text.x = element_text(size=15,angle = 45, hjust=1, vjust=1), axis.text.y = element_text(size=15),axis.ticks.length = unit(0.2,"cm"),legend.position= "bottom") + 
    xlab("") + ylab("Cellular fraction (%) of neutrophil") +
    theme(axis.title.y=element_text(size=16),strip.text = element_text(size = 13))+
    geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,position=position_dodge(0.05))
# Save plot 
ggsave(p1, filename = paste0(dir,"plot_out/N_major_across_treatment.pdf"),width = 3.2, height = 5)

source(paste0(dir,"code_clean/Cal_fraction_patients.R"))
tab.2 <- cal_fraction_patients(RNA,"N_major_cluster")
tab.2$cell <- gsub("\\.[0-9]{1,}","",rownames(tab.2))
tab.2$cell <- gsub("\\."," ",tab.2$cell)
tab.2$cell <- factor(tab.2$cell, levels =c("mature","aged"))
#colnames(tab.2)[4] <- "Patient_ID"
clin_info <- read.csv("clin_info.csv",header=T,stringsAsFactors=F)
sample_info <- clin_info[,c("Patient_ID","Group")]
table.plot <- merge(tab.2,sample_info,by="Patient_ID")
table.plot$Group <- factor(table.plot$Group,levels=c("TN","MPR","NMPR"))

pdf(paste0(dir,"plot_out/N_major_across_treatment3.pdf"),width = 3.2,height = 5)
    p1+geom_point(data=table.plot,aes(x=Group, y=Estimate,color=Patient_ID),size = 1.5)
dev.off()

# boxplot
p2<- ggplot(table.plot) +
    geom_boxplot(aes(x=cell, y=Estimate,color=Group,outlier.colour = NA),width = 0.6,position=position_dodge(0.75))+
    geom_jitter(aes(x=cell, y=Estimate,color=Group),size = 1.5, position = position_jitterdodge())+
    scale_color_manual(name="Group",values = Group_color_panel)+
    theme_bw() +  
    theme(axis.text.x = element_text(size=13,angle = 45, hjust=1, vjust=1), axis.text.y = element_text(size=13),axis.ticks.length = unit(0.2,"cm"),legend.position= "bottom") + 
    xlab("") + ylab("Cellular fraction (%) of neutrophil") +
    theme(axis.title.y=element_text(size=15),strip.text = element_text(size = 13))
# Save plot 
ggsave(p2, filename = paste0(dir,"plot_out/N_major_across_treatment2.pdf"),width = 3.2, height = 5)

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

pdf(paste0(dir,"plot_out/N_major_fraction_each_cell.pdf"),height=5,width=4)
for(i in unique(table.plot$cell)){print(ps[[i]])}
dev.off()



