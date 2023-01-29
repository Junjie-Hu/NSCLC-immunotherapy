library(Seurat)
library(ggplot2)
library(tidyverse)

dir <- "./"

RNA <- readRDS("MF.rds")

###### Heatmap plot
## select MF
MF <-subset(RNA, idents = c(0:4,6:8))
features=c("VCAN","FCN1","S100A8","S100A9","FABP4","MCEMP1","MARCO","C1QA","C1QB","GPNMB","APOE","SPP1","SELENOP","MRC1","TGFB1","CD163","CCL18","MSR1","VEGFA","TNF","CXCL9","CXCL10","CXCL11","HLA-DRA","HLA-DQA1","HLA-DPA1","CD74","MKI67","TOP2A")
row_split = c(rep("Monocyte",4),rep("Alveolar macrophage",3),rep("Immunosuppressive",4),rep("M2 signature",8),rep("M1 signature",4),rep("MHC II",4),rep("Proliferating",2))
row_split = factor(row_split,levels = c("Monocyte","Alveolar macrophage","Immunosuppressive","M2 signature","M1 signature","MHC II","Proliferating"))
### plot heatmap
source(paste0(dir,"code_clean/Heat_Dot_data.R"))
### set colnames order
plot_ord <- c("Mono_CX3CR1","Mono_VEGFA","Macro_FABP4","Macro_SPP1", "Macro_CXCL9","Macro_SELENOP","Macro_MKI67","Macro_C1QA")

data.plot <- Heat_Dot_data(object=MF,features=features,group.by="MF_cluster_annotation")
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
col_split = c(rep("Monocyte",2),rep("Macrophage",6))
col_split =factor(col_split,levels = c("Monocyte","Macrophage"))
# left annotation
annot = c("Monocyte","Alveolar macrophage","Immunosuppressive","M2 signature","M1 signature","MHC II","Proliferating")

ha = HeatmapAnnotation(df = data.frame(Marker=row_split),which = "row",
                       col = list(Marker = c("Monocyte" = "#e76f51", "Alveolar macrophage"="#57cc99","Immunosuppressive"="#0077b6","M2 signature"="#ddbea9",
                                           "M1 signature"="#00b4d8","MHC II"="#dc2f02","Proliferating"="#8a5a44")))
pdf(paste0(dir,"plot_out/MF_maker_heat_new.pdf"),width = 6,height = 8)
Heatmap(exp.mat, col = col_fun,cluster_columns = F,cluster_rows = F,
        show_column_names = T,show_row_names = T,rect_gp=gpar(col="grey"),
        column_names_side = "top",row_names_side = "right",
        row_split = row_split,column_split = col_split,
        row_gap = unit(3, "mm"),column_gap =unit(3, "mm"), 
        left_annotation = ha,
        heatmap_legend_param=list(title = "Expression",legend_height=unit(3, "cm")))
dev.off()

pdf(paste0(dir,"plot_out/MF_maker_heat_new2.pdf"),width = 6,height = 13)
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

###### Heatmap plot of Mono
## select Mono
Mono <-subset(RNA, idents = c(0,7))
Mono_features=c("CD14","FCGR3A","CX3CR1","NR4A1","CD81","VEGFA","OLR1","FPR3","GPNMB","MRC1","CD163","MSR1","FN1","FCGR2B",
    "SELL","VNN2","S100A8","S100A9","S100A12","RIPOR2","POU2F2","CFP","MEGF9","MNDA")
# data.features
data.features <- FetchData(object = Mono, vars = Mono_features)
data.features$id <- Mono$MF_cluster_annotation
if (!is.factor(x = data.features$id)) {
  data.features$id <- factor(x = data.features$id)
}
id.levels <- levels(x = data.features$id)
data.features$id <- as.vector(x = data.features$id)
## prepare data.plot
data.exp <- NULL
data.plot <- sapply(X = unique(x = data.features$id), FUN = function(ident) {
  data.use <- data.features[data.features$id == ident, 
                            1:(ncol(x = data.features) - 1), drop = FALSE]
  avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
    return(mean(x = expm1(x = x)))
  })
  if(!is.null(data.exp)){
    data.exp <- cbind(data.exp,avg.exp)
  }else{data.exp <- avg.exp}
  return(data.exp)
})
min(data.plot);max(data.plot)
## log
data.plot <- log2(data.plot+1)
min(data.plot);max(data.plot)

## per.exp
source(paste0(dir,"code_clean/Heat_Dot_data.R"))
exp.plot <- Heat_Dot_data(object=Mono,features=Mono_features,group.by="MF_cluster_annotation")
exp.mat <- exp.plot %>% select(features.plot,id,avg.exp.scaled) %>% spread(id,avg.exp.scaled)
rownames(exp.mat) <- exp.mat$features.plot
exp.mat$features.plot <- NULL
per.mat <- exp.plot %>% select(features.plot,id,pct.exp) %>% spread(id,pct.exp)
rownames(per.mat) <- per.mat$features.plot
per.mat$features.plot <- NULL
per.mat <- per.mat/100

exp.mat <- data.plot

### plot heatmap
library(ComplexHeatmap)
library(circlize) ## color 
# set color gradient
col_fun2 <- colorRamp2(c(0, 2, 4), c("#118ab2", "#fdffb6", "#e63946"))
pdf(paste0(dir,"plot_out/Mono_maker_heat.pdf"),width = 2.5,height = 5)
Heatmap(data.plot, col = col_fun2,cluster_columns = F,cluster_rows = F,
        show_column_names = T,show_row_names = T,rect_gp=gpar(col="grey"),
        column_names_side = "top",row_names_side = "right",
        heatmap_legend_param=list(title = "Expression",legend_height=unit(3, "cm")))
dev.off()

pdf(paste0(dir,"plot_out/Mono_maker_heat2.pdf"),width = 2.5,height = 9.5)
Heatmap(exp.mat, col = col_fun2,cluster_columns = F,cluster_rows = F,
        show_column_names = T,show_row_names = T,rect_gp=gpar(type = "none"),
        cell_fun = function(j, i, x, y, width, height, fill){
          grid.rect(x = x, y = y, width = width, height = height,gp = gpar(col = "grey", fill = NA))
          grid.circle(x = x, y = y,r=per.mat[i,j]/2 * max(unit.c(width, height)),gp = gpar(fill = col_fun(exp.mat[i, j]), col = NA))},
        column_names_side = "top",row_names_side = "right",
        heatmap_legend_param=list(title = "Expression",legend_height=unit(3, "cm")))
dev.off()


###### MF Trajectory
library(monocle)
# select MF
MF_tj <- subset(RNA, idents = c(0:3,7:8))
# Extract data, phenotype data, and feature data from the SeuratObject
data <- as(as.matrix(MF_tj@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = MF_tj@meta.data[,c("MF.cluster","MF_cluster_annotation")])
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
MF_DEG_genes <- differentialGeneTest(monocle_cds, fullModelFormulaStr = '~MF_cluster_annotation', cores = 20)
ordering_genes <- row.names (subset(MF_DEG_genes, qval < 0.001))
gene_check <- row.names(MF_DEG_genes)[order(MF_DEG_genes$qval)][1:1000]
write.csv(MF_DEG_genes,"data_out/MF_monocle_DEG.csv",quote=F)
## 
monocle_cds <- setOrderingFilter(monocle_cds,gene_check)
monocle_cds <- reduceDimension(
  monocle_cds,
  max_components = 2,
  method = 'DDRTree')
monocle_cds <- orderCells(monocle_cds)
ClusterName_color_panel <- c("Mono_VEGFA" = "#F8766D", "Macro_SPP1" = "#DB8E00","Macro_CXCL9" = "#AEA200","Macro_SELENOP"="#64B200","Mono_CX3CR1"="#00A6FF","Macro_C1QA"="#B385FF")
library(cowplot)
pdf(paste0(dir,"plot_out/MF_Pseudotime_new.pdf"),width = 8,height = 5)
plot_cell_trajectory(monocle_cds, color_by = "Pseudotime",show_branch_points=F)+theme_cowplot()
plot_cell_trajectory(monocle_cds, color_by = "MF_cluster_annotation",cell_size = 1,show_branch_points=F)+theme_cowplot()+
   scale_color_manual(name = "", values = ClusterName_color_panel)+guides(color = guide_legend(override.aes = list(size=5)))+
   theme(axis.text.y = element_text(size=13),axis.text.x = element_text(size=13),axis.title.y = element_text(size=16),axis.title.x = element_text(size=15),axis.ticks.length = unit(0.2,"cm"))
dev.off()

saveRDS(monocle_cds,"data_out/MF_monocle_new.rds")

## CFP in PMo across treatment
PMo <- subset(RNA, idents = 7)
data.pmo <- FetchData(object = PMo, vars = c("CFP","CX3CR1"))
data.pmo$Group <- PMo$Group
library(ggsignif)
compaired <- list(c("TN","MPR"),c("TN","NMPR"),c("MPR","NMPR"))
Group_color_panel <- c("TN"="#228B22","MPR"="#FF4500","NMPR"="#757bc8")
# CFP
p <- ggplot(data=data.pmo,aes(x=factor(Group,levels=c("TN","MPR","NMPR")),y=CFP,fill=Group)) + 
     geom_violin(color="white")+
     geom_boxplot(width=0.1,position=position_dodge(0.9),fill="white") + theme_cowplot()+
     scale_fill_manual(name = "Group", values = Group_color_panel)+
     theme(axis.text.x = element_text(size=10,angle = 45, hjust=1, vjust=1), axis.text.y = element_text(size=10),axis.ticks.length = unit(0.2,"cm"))+
     NoLegend()+labs(x="")+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test)
ggsave(paste0(dir,"plot_out/PMo_CFP_treatment.pdf"),p,width =3,height =4)

## VEGFA in Mono_VEGFA across treatment
Mono <- subset(RNA, idents = 0)
data.mono <- FetchData(object = Mono, vars = c("VEGFA"))
data.mono$Group <- Mono$Group
library(ggsignif)
compaired <- list(c("TN","MPR"),c("TN","NMPR"),c("MPR","NMPR"))
Group_color_panel <- c("TN"="#228B22","MPR"="#FF4500","NMPR"="#757bc8")
# CFP
p <- ggplot(data=data.mono,aes(x=factor(Group,levels=c("TN","MPR","NMPR")),y=VEGFA,fill=Group)) + 
     geom_violin(color="white")+
     geom_boxplot(width=0.1,position=position_dodge(0.9),fill="white") + theme_cowplot()+
     scale_fill_manual(name = "Group", values = Group_color_panel)+
     theme(axis.text.x = element_text(size=10,angle = 45, hjust=1, vjust=1), axis.text.y = element_text(size=10),axis.ticks.length = unit(0.2,"cm"))+
     NoLegend()+labs(x="")+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test)
ggsave(paste0(dir,"plot_out/Mono_VEGFA_treatment.pdf"),p,width =3,height =4)

###Calculate cell fractions in different response groups
source(paste0(dir,"code/Cal_fraction.R"))
tab.1 <- cal_fraction(RNA,"MF_cluster_annotation")
tab.1$cell <- factor(tab.1$cell, levels=c("Mono_CX3CR1","Mono_VEGFA","Macro_CXCL9","Macro_SPP1","Macro_SELENOP","Macro_C1QA","Macro_FABP4","Macro_MKI67","DC_CD1C","DC_LAMP3","DC_XCR1"))
Group_color_panel <- c("TN"="#228B22","MPR"="#FF4500","NMPR"="#757bc8")
# plot
p1<- ggplot(tab.1, aes(x=Group, y=Estimate)) +
    geom_bar(stat = "identity", aes(fill=Group)) +
    scale_fill_manual(name = "Group", values = Group_color_panel) + facet_grid(~ cell) +
    theme_bw() +  
    theme(axis.text.x = element_text(size=15,angle = 45, hjust=1, vjust=1), axis.text.y = element_text(size=15),axis.ticks.length = unit(0.2,"cm"),legend.position= "right") + 
    xlab("") + ylab("Cellular fraction (%) of myeloid cell") +
    theme(axis.title.y=element_text(size=16),strip.text = element_text(size = 13))+
    geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,position=position_dodge(0.05))
## theme(strip.background = element_rect(fill="white"))
# Save plot 
ggsave(p1, filename = paste0(dir,"plot_out/MF_across_treatment.pdf"),width = 16, height = 5)

source(paste0(dir,"code_clean/Cal_fraction_patients.R"))
tab.2 <- cal_fraction_patients(RNA,"MF_cluster_annotation")
tab.2$cell <- gsub("\\.[0-9]{1,}","",rownames(tab.2))
tab.2$cell <- gsub("\\."," ",tab.2$cell)
tab.2$cell <- factor(tab.2$cell, levels=c("Mono_CX3CR1","Mono_VEGFA","Macro_CXCL9","Macro_SPP1","Macro_SELENOP","Macro_C1QA","Macro_FABP4","Macro_MKI67","DC_CD1C","DC_LAMP3","DC_XCR1"))
#colnames(tab.2)[4] <- "Patient_ID"
clin_info <- read.csv("clin_info.csv",header=T,stringsAsFactors=F)
sample_info <- clin_info[,c("Patient_ID","Group")]
table.plot <- merge(tab.2,sample_info,by="Patient_ID")
table.plot$Group <- factor(table.plot$Group,levels=c("TN","MPR","NMPR"))

pdf(paste0(dir,"plot_out/MF_across_treatment3.pdf"),width = 16,height = 5)
    p1+geom_point(data=table.plot,aes(x=Group, y=Estimate,color=Patient_ID),size = 1.5)
dev.off()

# boxplot
p2<- ggplot(table.plot) +
    geom_boxplot(aes(x=cell, y=Estimate,color=Group,outlier.colour = NA),width = 0.6,position=position_dodge(0.75))+
    geom_jitter(aes(x=cell, y=Estimate,color=Group),size = 1.5, position = position_jitterdodge())+
    scale_color_manual(name="Group",values = Group_color_panel)+
    theme_bw() +  
    theme(axis.text.x = element_text(size=13,angle = 45, hjust=1, vjust=1), axis.text.y = element_text(size=13),axis.ticks.length = unit(0.2,"cm"),legend.position= "right") + 
    xlab("") + ylab("Cellular fraction (%) of myeloid cell") +
    theme(axis.title.y=element_text(size=15),strip.text = element_text(size = 13))
# Save plot 
ggsave(p2, filename = paste0(dir,"plot_out/MF_across_treatment2.pdf"),width = 16, height = 5)

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

pdf(paste0(dir,"plot_out/MF_fraction_each_cell.pdf"),height=5,width=4)
for(i in unique(table.plot$cell)){print(ps[[i]])}
dev.off()

saveRDS(RNA,"MF.rds")

######## MF M1_M2 plot
# select MF
MF <-subset(RNA, idents = c(1,2,3,4,6,8))
# feature from CIBERSORT LM22
M0_features <- list(c("ACP5","BHLHE41","C5AR1","CCDC102B","CCL22","CCL7","COL8A2","CSF1","CXCL3","CXCL5","CYP27A1","DCSTAMP","GPC4","HK3","IGSF6","MARCO","MMP9","NCF2","PLA2G7","PPBP","QPCT","SLAMF8","SLC12A8","TNFSF14","VNN1"))
M1_features <- list(c("ACHE","ADAMDEC1","APOBEC3A","APOL3","APOL6","AQP9","ARRB1","CCL19","CCL5","CCL8","CCR7","CD38","CD40","CHI3L1","CLIC2","CXCL10","CXCL11","CXCL13","CXCL9","CYP27B1","DHX58","EBI3","GGT5","HESX1","IDO1","IFI44L","IL2RA","KIAA0754","KYNU","LAG3","LAMP3","LILRA3","LILRB2","NOD2","PLA1A","PTGIR","RASSF4","RSAD2","SIGLEC1","SLAMF1","SLC15A3","SLC2A6","SOCS1","TLR7","TLR8","TNFAIP6","TNIP3","TRPM4"))
M2_features <- list(c("ADAMDEC1","AIF1","ALOX15","CCL13","CCL14","CCL18","CCL23","CCL8","CD209","CD4","CD68","CFP","CHI3L1","CLEC10A","CLEC4A","CLIC2","CRYBB1","EBI3","FAM198B","FES","FRMD4A","FZD2","GGT5","GSTT1","HRH1","HTR2B","MS4A6A","NME8","NPL","P2RY13","PDCD1LG2","RENBP","SIGLEC1","SLC15A3","TLR8","TREM2","WNT5B"))
# calculate score
MF <- AddModuleScore(object=MF,features=M1_features,name="M1_score")
MF <- AddModuleScore(object=MF,features=M2_features,name="M2_score")
MF <- AddModuleScore(object=MF,features=M0_features,name="M0_score")
# plot
Group_color_panel <- c("TN"="#228B22","MPR"="#FF4500","NMPR"="#757bc8")
# violin plot
data.plot <- MF@meta.data[,c("Group","MF_cluster_annotation","M0_score1","M1_score1","M2_score1")]
library(ggsignif)
library(cowplot)
compaired <- list(c("TN","MPR"),c("TN","NMPR"),c("MPR","NMPR"))
p3 <- ggplot(data = data.plot,aes(x=factor(Group,levels=c("TN","MPR","NMPR")),y=M1_score1,fill=Group))+
     geom_violin(color="white")+
     geom_boxplot(width=0.1,position=position_dodge(0.9),fill="white")+theme_cowplot()+
     scale_fill_manual(name="Group",values = Group_color_panel)+
     theme(axis.text.x = element_text(size=10,angle = 45, hjust=1, vjust=1), axis.text.y = element_text(size=10),axis.ticks.length = unit(0.2,"cm"))+
     NoLegend()+labs(x="",y="M1 signature")+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test)
ggsave(paste0(dir,"plot_out/M1_across_treatment.pdf"),p3,width =3,height =5)
p4 <- ggplot(data = data.plot,aes(x=factor(Group,levels=c("TN","MPR","NMPR")),y=M2_score1,fill=Group))+
     geom_violin(color="white")+
     geom_boxplot(width=0.1,position=position_dodge(0.9),fill="white")+theme_cowplot()+
     scale_fill_manual(name="Group",values = Group_color_panel)+
     theme(axis.text.x = element_text(size=10,angle = 45, hjust=1, vjust=1), axis.text.y = element_text(size=10),axis.ticks.length = unit(0.2,"cm"))+
     NoLegend()+labs(x="",y="M2 signature")+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test)
ggsave(paste0(dir,"plot_out/M2_across_treatment.pdf"),p4,width =3,height =5)
p5 <- ggplot(data = data.plot,aes(x=factor(Group,levels=c("TN","MPR","NMPR")),y=M0_score1,fill=Group))+
     geom_violin(color="white")+
     geom_boxplot(width=0.1,position=position_dodge(0.9),fill="white")+theme_cowplot()+
     scale_fill_manual(name="Group",values = Group_color_panel)+
     theme(axis.text.x = element_text(size=10,angle = 45, hjust=1, vjust=1), axis.text.y = element_text(size=10),axis.ticks.length = unit(0.2,"cm"))+
     NoLegend()+labs(x="",y="M0 signature")+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test)
ggsave(paste0(dir,"plot_out/M0_across_treatment.pdf"),p5,width =3,height =5)
#### each cluster
pm0 <- ggplot(data = data.plot,aes(x=factor(Group,levels=c("TN","MPR","NMPR")),y=M0_score1,fill=Group))+
     geom_violin(color="white")+
     geom_boxplot(width=0.1,position=position_dodge(0.9),fill="white")+facet_grid(~MF_cluster_annotation) +theme_bw()+
     scale_fill_manual(name="Group",values = Group_color_panel)+
     theme(axis.text.x = element_text(size=10,angle = 45, hjust=1, vjust=1), axis.text.y = element_text(size=10),axis.ticks.length = unit(0.2,"cm"),legend.position="bottom")+
     labs(x="",y="M0 signature")+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test)
ggsave(paste0(dir,"plot_out/M0_cluster_across_treatment.pdf"),pm0,width =12,height =5)
pm1 <- ggplot(data = data.plot,aes(x=factor(Group,levels=c("TN","MPR","NMPR")),y=M1_score1,fill=Group))+
     geom_violin(color="white")+
     geom_boxplot(width=0.1,position=position_dodge(0.9),fill="white")+facet_grid(~MF_cluster_annotation) +theme_bw()+
     scale_fill_manual(name="Group",values = Group_color_panel)+
     theme(axis.text.x = element_text(size=10,angle = 45, hjust=1, vjust=1), axis.text.y = element_text(size=10),axis.ticks.length = unit(0.2,"cm"),legend.position="bottom")+
     labs(x="",y="M1 signature")+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test)
ggsave(paste0(dir,"plot_out/M1_cluster_across_treatment.pdf"),pm1,width =12,height =5)
pm2 <- ggplot(data = data.plot,aes(x=factor(Group,levels=c("TN","MPR","NMPR")),y=M2_score1,fill=Group))+
     geom_violin(color="white")+
     geom_boxplot(width=0.1,position=position_dodge(0.9),fill="white")+facet_grid(~MF_cluster_annotation) +theme_bw()+
     scale_fill_manual(name="Group",values = Group_color_panel)+
     theme(axis.text.x = element_text(size=10,angle = 45, hjust=1, vjust=1), axis.text.y = element_text(size=10),axis.ticks.length = unit(0.2,"cm"),legend.position="bottom")+
     labs(x="",y="M2 signature")+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test)
ggsave(paste0(dir,"plot_out/M2_cluster_across_treatment.pdf"),pm2,width =12,height =5)

##### DC 
DC <- subset(RNA,idents=c(5,9,10))
## heapmap plot
features=c("XCR1","CLEC9A","CADM1","CLNK","CD1C","CD1E","FCER1A","CLEC10A","LAMP3","FSCN1","CCR7","LAD1","CCL17","CCL19","CCL22","CX3CR1","CD86","CD83","CD80","CD40","ICOSLG",
"IDO1","CD274","PDCD1LG2","CD200","TLR1","TLR2","TLR4","TLR7","TLR10")
### plot heatmap
source(paste0(dir,"code_clean/Heat_Dot_data.R"))
### set colnames order
plot_ord <- c("DC_XCR1","DC_CD1C","DC_LAMP3")

data.plot <- Heat_Dot_data(object=DC,features=features,group.by="MF_cluster_annotation")
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
col_fun <- colorRamp2(c(-1.5, 0, 1.5), c("#118ab2", "#fdffb6", "#e63946"))
# split heatmap
row_split = c(rep("cDC1",4),rep("cDC2",4),rep("Activation and migration",4),rep("Chemotaxis",4),rep("Co-stimulatory",5),rep("Inhibitory",4),rep("TLRs",5))
row_split = factor(row_split,levels = c("cDC1","cDC2","Activation and migration","Chemotaxis","Co-stimulatory","Inhibitory","TLRs"))
# left annotation
ha = HeatmapAnnotation(df = data.frame(Marker=row_split),which = "row",
                       col = list(Marker = c("cDC1" = "#e76f51", "cDC2"="#0077b6","Activation and migration"="#ddbea9",
                                  "Chemotaxis"="#00b4d8","Co-stimulatory"="#dc2f02","Inhibitory"="#fca311","TLRs"="#8a5a44")))
pdf(paste0(dir,"plot_out/DC_maker_heat_new.pdf"),width = 5,height = 9)
Heatmap(exp.mat, col = col_fun,cluster_columns = F,cluster_rows = F,
        show_column_names = T,show_row_names = T,rect_gp=gpar(col="grey"),
        column_names_side = "top",row_names_side = "right",
        row_split = row_split,
        row_gap = unit(3, "mm"),
        left_annotation = ha,
        heatmap_legend_param=list(title = "Expression",legend_height=unit(3, "cm")))
dev.off()

pdf(paste0(dir,"plot_out/DC_maker_heat_new2.pdf"),width = 4.5,height = 12)
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

# DC activating signature
DC_score_gene <- read.csv(paste0(dir,"data_out/DC_score_gene.csv"),header=T,stringsAsFactors=F)
DC_ap <- list(DC_score_gene$AP)
DC_is <- list(DC_score_gene$IS)
# calculate score
DC <- AddModuleScore(object=DC,features=DC_ap,name="DC_ap")
DC <- AddModuleScore(object=DC,features=DC_is,name="DC_is")
library(ggsignif)
data.plot <- DC@meta.data[,c("Group","MF_cluster_annotation","DC_ap1","DC_is1")]
compaired <- list(c("TN","MPR"),c("TN","NMPR"),c("MPR","NMPR"))
p <- ggplot(data = data.plot,aes(x=factor(Group,levels=c("TN","MPR","NMPR")),y=DC_ap1,fill=Group))+
     geom_violin(color="white")+
     geom_boxplot(width=0.1,position=position_dodge(0.9),fill="white")+theme_cowplot()+
     scale_fill_manual(name="Group",values = Group_color_panel)+
     theme(axis.text.x = element_text(size=10,angle = 45, hjust=1, vjust=1), axis.text.y = element_text(size=10),axis.ticks.length = unit(0.2,"cm"))+
     NoLegend()+labs(x="",y="Antigen-presenting signature")+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test)
ggsave(paste0(dir,"plot_out/DC_ap_across_treatment.pdf"),p,width =3,height =5)
p <- ggplot(data = data.plot,aes(x=factor(Group,levels=c("TN","MPR","NMPR")),y=DC_is1,fill=Group))+
     geom_violin(color="white")+
     geom_boxplot(width=0.1,position=position_dodge(0.9),fill="white")+theme_cowplot()+
     scale_fill_manual(name="Group",values = Group_color_panel)+
     theme(axis.text.x = element_text(size=10,angle = 45, hjust=1, vjust=1), axis.text.y = element_text(size=10),axis.ticks.length = unit(0.2,"cm"))+
     NoLegend()+labs(x="",y="Immunosuppressive signature")+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test)
ggsave(paste0(dir,"plot_out/DC_is_across_treatment.pdf"),p,width =3,height =5)
## cluster
pa <- ggplot(data = data.plot,aes(x=factor(Group,levels=c("TN","MPR","NMPR")),y=DC_ap1,fill=Group))+
     geom_violin(color="white")+
     geom_boxplot(width=0.1,position=position_dodge(0.9),fill="white")+facet_grid(~MF_cluster_annotation) +theme_bw()+
     scale_fill_manual(name="Group",values = Group_color_panel)+
     theme(axis.text.x = element_text(size=10,angle = 45, hjust=1, vjust=1), axis.text.y = element_text(size=10),axis.ticks.length = unit(0.2,"cm"),legend.position="bottom")+
     labs(x="",y="Antigen-presenting signature")+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test)
ggsave(paste0(dir,"plot_out/DC_cluster_ap_across_treatment.pdf"),pa,width =6.5,height =5)
pi <- ggplot(data = data.plot,aes(x=factor(Group,levels=c("TN","MPR","NMPR")),y=DC_is1,fill=Group))+
     geom_violin(color="white")+
     geom_boxplot(width=0.1,position=position_dodge(0.9),fill="white")+facet_grid(~MF_cluster_annotation) +theme_bw()+
     scale_fill_manual(name="Group",values = Group_color_panel)+
     theme(axis.text.x = element_text(size=10,angle = 45, hjust=1, vjust=1), axis.text.y = element_text(size=10),axis.ticks.length = unit(0.2,"cm"),legend.position="bottom")+
     labs(x="",y="Immunosuppressive signature")+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test)
ggsave(paste0(dir,"plot_out/DC_cluster_is_across_treatment.pdf"),pi,width =6.5,height =5)




