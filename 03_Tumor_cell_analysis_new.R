library(Seurat)
library(ggplot2)
library(tidyverse)
dir <- "./"

RNA <- readRDS("epithelium.rds")

### CNV from copyKAT
copykat.test <- readRDS(paste0(dir,"data_out/copykat_res.rds"))
pred.test <- data.frame(copykat.test$prediction)

# some cells were removed when running copykat
RNA@meta.data$cell.names <- rownames(RNA@meta.data)
RNA_sub <- subset(RNA,cell.names %in% pred.test$cell.names)
keep <- which(pred.test$cell.names %in% RNA_sub$cell.names)
pred.test <- pred.test[keep,]
RNA_sub@meta.data$copykat <- pred.test$copykat.pred
# plot
p=DimPlot(RNA_sub, reduction = "umap", label = T, group.by = "copykat")
ggsave(paste0(dir,"plot_out/TC_copykat_umap.pdf"),p,width=6.5, height=6)

## pick out E1-KRT17 to tell normal and tumor cell apart
RNA_E1 <- subset(RNA_sub,idents=1)
# plot
p=DimPlot(RNA_E1, reduction = "umap", label = F, group.by = "copykat")
ggsave(paste0(dir,"plot_out/TC_E1_umap.pdf"),p,width=4, height=4)
# annotate E1-KRT17 normal and malignant cells
cell_mag <- subset(RNA_E1@meta.data,copykat=="aneuploid")
cell_norm <- subset(RNA_E1@meta.data,copykat=="diploid")
keep1 <- RNA_sub@meta.data$cell.names %in% cell_mag$cell.names
E1_mag <- which(keep1==T)
keep2 <- RNA_sub@meta.data$cell.names %in% cell_norm$cell.names
E1_norm <- which(keep2==T)
RNA_sub$TC_cluster_annotation <- as.vector(RNA_sub$TC_cluster_annotation)
RNA_sub$TC_cluster_annotation[E1_mag] <- "E1-KRT17-malig"
RNA_sub$TC_cluster_annotation[E1_norm] <- "E1-KRT17-norm"
# plot umap
pdf(paste(dir,"plot_out/TC_cluster_annotation_E1sub.pdf", sep=""),width = 7,height = 6)
DimPlot(RNA_sub, reduction = "umap", label = F, group.by = 'TC_cluster_annotation')
dev.off()

saveRDS(RNA_sub,"Epithelium_sub.rds")

RNA_sub <- readRDS("Epithelium_sub.rds")
###Calculate cell fractions in different response groups
source(paste0(dir,"code/Cal_fraction.R"))
RNA_sub$TC_cluster_annotation <- factor(RNA_sub$TC_cluster_annotation)
tab.1 <- cal_fraction(RNA_sub,"TC_cluster_annotation")
tab.1$cell <- gsub("\\.[0-9]{1,}$","",rownames(tab.1))
tab.1$cell <- gsub("\\.","_",tab.1$cell)
tab.1$cell <- factor(tab.1$cell,levels=c("E1_KRT17_norm","E2_S100P","E5_SFTPA2","E6_SPARCL1","E8_SCGB3A1","E9_TPPP3",
        "E0_DST","E1_KRT17_malig","E3_PCNA","E4_TOP2A","E7_SERPINB9"))
Group_color_panel <- c("TN"="#228B22","MPR"="#FF4500","NMPR"="#757bc8")
# plot
p1<- ggplot(tab.1, aes(x=Group, y=Estimate)) +
    geom_bar(stat = "identity", aes(fill=Group)) +
    scale_fill_manual(name = "Group", values = Group_color_panel) + facet_grid(~ cell) +
    theme_bw() +  
    theme(axis.text.x = element_text(size=13,angle = 45, hjust=1, vjust=1), axis.text.y = element_text(size=13),axis.ticks.length = unit(0.2,"cm"),legend.position= "bottom") + 
    xlab("") + ylab("Cellular fraction (%) of epithelial cells") +theme(axis.title.y=element_text(size=15),strip.text = element_text(size = 11))+
    geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,position=position_dodge(0.05))
# Save plot 
ggsave(p1, filename = paste0(dir,"plot_out/TC_across_treatment.pdf"),width = 14.8, height = 5.6)

source(paste0(dir,"code_clean/Cal_fraction_patients.R"))
tab.2 <- cal_fraction_patients(RNA_sub,"TC_cluster_annotation")
tab.2$cell <- gsub("\\.[0-9]{1,}$","",rownames(tab.2))
tab.2$cell <- gsub("\\.","_",tab.2$cell)
tab.2$cell <- factor(tab.2$cell, levels=c("E1_KRT17_norm","E2_S100P","E5_SFTPA2","E6_SPARCL1","E8_SCGB3A1","E9_TPPP3",
        "E0_DST","E1_KRT17_malig","E3_PCNA","E4_TOP2A","E7_SERPINB9"))
#colnames(tab.2)[4] <- "Patient_ID"
clin_info <- read.csv("clin_info.csv",header=T,stringsAsFactors=F)
sample_info <- clin_info[,c("Patient_ID","Group")]
table.plot <- merge(tab.2,sample_info,by="Patient_ID")
table.plot$Group <- factor(table.plot$Group,levels=c("TN","MPR","NMPR"))

pdf(paste0(dir,"plot_out/TC_across_treatment3.pdf"),width = 14.8,height = 6.3)
    p1+geom_point(data=table.plot,aes(x=Group, y=Estimate,color=Patient_ID),size = 1.5)
dev.off()

# boxplot
p2<- ggplot(table.plot) +
    geom_boxplot(aes(x=cell, y=Estimate,color=Group,outlier.colour = NA),width = 0.6,position=position_dodge(0.75))+
    geom_jitter(aes(x=cell, y=Estimate,color=Group),size = 1.5, position = position_jitterdodge())+
    scale_color_manual(name="Group",values = Group_color_panel)+
    theme_bw() +  
    theme(axis.text.x = element_text(size=13,angle = 45, hjust=1, vjust=1), axis.text.y = element_text(size=13),axis.ticks.length = unit(0.2,"cm"),legend.position= "bottom") + 
    xlab("") + ylab("Cellular fraction (%) of epithelial cells") +
    theme(axis.title.y=element_text(size=15),strip.text = element_text(size = 13))
# Save plot 
ggsave(p2, filename = paste(dir,"plot_out/TC_across_treatment2.pdf", sep=""),width = 14.8, height = 6.3)

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
    labs(x="",y="Cellular fraction (%) of epithelial cells",title=i) +
    theme(axis.title.y=element_text(size=15),strip.text = element_text(size = 13))+
    geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test)
}

pdf(paste0(dir,"plot_out/TC_fraction_each_cell.pdf"),height=5,width=4)
for(i in unique(table.plot$cell)){print(ps[[i]])}
dev.off()


RNA_mag <- readRDS("Epi_mag.rds")

#### MPR genes check
RNA_mag$Group <- factor(RNA_mag$Group,levels=c("TN","MPR","NMPR"))
p <- VlnPlot(RNA_mag,features=c("AKR1B1","AKR1B10","AKR1C1","AKR1C2","AKR1C3"),cols=c("#228B22","#FF4500","#757bc8"),pt.size = 0,group.by="Group",ncol=3)
ggsave(paste0(dir,"plot_out/TC_NMPR_marker.pdf"),p,width =7,height =6)

# NMPR genes check
p <- VlnPlot(RNA_mag,features=c("SPARCL1","CX3CL1","CD74","HLA-DRA","HLA-DPA1","HLA-DQA1","SERPINB9","HLA-DRB1","HLA-DPB1"),cols=c("#228B22","#FF4500","#757bc8"),pt.size = 0,group.by="Group",ncol=3)
ggsave(paste0(dir,"plot_out/TC_MPR_marker.pdf"),p,width =7,height =9)

# gene expr distribution
gene.plot <- c("AKR1B1","AKR1B10","AKR1C1","AKR1C2","AKR1C3",
    "SPARCL1","CX3CL1","CD74","HLA-DRA","HLA-DPA1","HLA-DQA1","SERPINB9","HLA-DRB1","HLA-DPB1")
data.exp <- FetchData(object = RNA_mag, vars = gene.plot)
library(cowplot)
library(ggsignif)
library(gridExtra)
data.exp$Group <- RNA_mag$Group
compaired <- list(c("TN","MPR"),c("TN","NMPR"),c("MPR","NMPR"))
# plot
pdf(paste0(dir,"plot_out/TC_marker_compare.pdf"),width =3,height = 4)
for(i in 1:length(gene.plot)){
pg <- ggplot(data=data.exp,aes(x=factor(Group,levels=c("TN","MPR","NMPR")),y=data.exp[,i],fill=Group)) + 
     geom_violin(color="white")+ theme_cowplot()+
     scale_fill_manual(name="Group",values = Group_color_panel)+
     NoLegend()+labs(x="",y="Expression level",title=gene.plot[i])+theme(plot.title = element_text(hjust = 0.5))+
     geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test)
print(pg)
}
dev.off()

#### DEG
## TN vs MPR
RNA1 <- subset(RNA_mag,Group != "NMPR")
Idents(RNA1) <- RNA1@meta.data$Group
DEG<- FindAllMarkers(object = RNA1, only.pos = T, logfc.threshold = 0.25,min.pct = 0.25)
write.csv(DEG, file=paste0(dir,"data_out/TC_TN-MPR.csv"),quote=F)

## TN vs NMPR
RNA2 <- subset(RNA_mag,Group != "MPR")
Idents(RNA2) <- RNA2@meta.data$Group
DEG<- FindAllMarkers(object = RNA2, only.pos = T, logfc.threshold = 0.25,min.pct = 0.25)
write.csv(DEG, file=paste0(dir,"data_out/TC_TN-NMPR.csv"),quote=F)

## MPR vs NMPR
RNA3 <- subset(RNA_mag,Group != "TN")
Idents(RNA3) <- RNA3@meta.data$Group
DEG<- FindAllMarkers(object = RNA3, only.pos = T, logfc.threshold = 0.25,min.pct = 0.25)
write.csv(DEG, file=paste0(dir,"data_out/TC_MPR-NMPR.csv"),quote=F)


## GSVA TN vs MPR
source(paste0(dir,"code/GSVA.R"))
RNA1 <- subset(RNA_mag,Group != "NMPR")
hallmark <- "/media/inspur/AS2150G2/HJJ/scrna/h.all.v7.2.symbols.gmt"
kegg <- "/media/inspur/AS2150G2/HJJ/scrna/c2.cp.kegg.v7.2.symbols.gmt"
reactome <- "/media/inspur/AS2150G2/HJJ/scrna/c2.cp.reactome.v7.2.symbols.gmt"
res <- GSVA_run(SeuratObj=RNA1,cluster_name="Group",cluster1="TN",cluster2="MPR",
                gmtFile=hallmark,kcdf="Gaussian")
res$pathway <- gsub("HALLMARK_","",rownames(res))
write.csv(res,paste0(dir,"data_out/TC_TN-MPR_gsva.csv"),row.names=F,quote=F)
res_plot <- res %>% dplyr::arrange(t) %>% group_by(cluster) %>% top_n(10,abs(t))
res_plot$pathway <- factor(res_plot$pathway,levels = res_plot$pathway)
# plot
p <- ggplot(res_plot,aes(x=pathway,y=t,fill=cluster))+geom_bar(stat ='identity')+
  geom_hline(yintercept = 0)+theme_bw() + theme(panel.grid =element_blank())+
  scale_fill_manual(name = "", values = c("C1"="#0077b6", "C2"="#f08080"))+
  labs(x="",y="t value")+
  theme(axis.text.x = element_text(size = 14),axis.title.x = element_text(size = 15),axis.text.y = element_text(size = 14),axis.ticks.length.y = unit(0,"mm"),panel.border = element_blank())+
  coord_flip()
ggsave(paste0(dir,"plot_out/TC_TN-MPR_gsva.pdf"),p,width =8,height =6)

## GSVA TN vs NMPR
source(paste0(dir,"code/GSVA.R"))
RNA2 <- subset(RNA_mag,Group != "MPR")
hallmark <- "/media/inspur/AS2150G2/HJJ/scrna/h.all.v7.2.symbols.gmt"
kegg <- "/media/inspur/AS2150G2/HJJ/scrna/c2.cp.kegg.v7.2.symbols.gmt"
reactome <- "/media/inspur/AS2150G2/HJJ/scrna/c2.cp.reactome.v7.2.symbols.gmt"
res <- GSVA_run(SeuratObj=RNA2,cluster_name="Group",cluster1="TN",cluster2="NMPR",
                gmtFile=hallmark,kcdf="Gaussian")
res$pathway <- gsub("HALLMARK_","",rownames(res))
write.csv(res,paste0(dir,"data_out/TC_TN-NMPR_gsva.csv"),row.names=F,quote=F)
res_plot <- res %>% dplyr::arrange(t) %>% group_by(cluster) %>% top_n(10,abs(t))
res_plot$pathway <- factor(res_plot$pathway,levels = res_plot$pathway)
# plot
p <- ggplot(res_plot,aes(x=pathway,y=t,fill=cluster))+geom_bar(stat ='identity')+
  geom_hline(yintercept = 0)+theme_bw() + theme(panel.grid =element_blank())+
  scale_fill_manual(name = "", values = c("C1"="#0077b6", "C2"="#7b2cbf"))+
  labs(x="",y="t value")+
  theme(axis.text.x = element_text(size = 14),axis.title.x = element_text(size = 15),axis.text.y = element_text(size = 14),axis.ticks.length.y = unit(0,"mm"),panel.border = element_blank())+
  coord_flip()
ggsave(paste0(dir,"plot_out/TC_TN-NMPR_gsva.pdf"),p,width =8,height =6)

## GSVA MPR vs NMPR
source(paste0(dir,"code/GSVA.R"))
RNA3 <- subset(RNA_mag,Group != "TN")
hallmark <- "/media/inspur/AS2150G2/HJJ/scrna/h.all.v7.2.symbols.gmt"
kegg <- "/media/inspur/AS2150G2/HJJ/scrna/c2.cp.kegg.v7.2.symbols.gmt"
reactome <- "/media/inspur/AS2150G2/HJJ/scrna/c2.cp.reactome.v7.2.symbols.gmt"
res <- GSVA_run(SeuratObj=RNA3,cluster_name="Group",cluster1="MPR",cluster2="NMPR",
                gmtFile=hallmark,kcdf="Gaussian")
res$pathway <- gsub("HALLMARK_","",rownames(res))
write.csv(res,paste0(dir,"data_out/TC_MPR-NMPR_gsva.csv"),row.names=F,quote=F)
res_plot <- res %>% dplyr::arrange(t) %>% group_by(cluster) %>% top_n(10,abs(t))
res_plot$pathway <- factor(res_plot$pathway,levels = res_plot$pathway)
# plot
p <- ggplot(res_plot,aes(x=pathway,y=t,fill=cluster))+geom_bar(stat ='identity')+
  geom_hline(yintercept = 0)+theme_bw() + theme(panel.grid =element_blank())+
  scale_fill_manual(name = "", values = c("C1"="#f08080", "C2"="#7b2cbf"))+
  labs(x="",y="t value")+
  theme(axis.text.x = element_text(size = 14),axis.title.x = element_text(size = 15),axis.text.y = element_text(size = 14),axis.ticks.length.y = unit(0,"mm"),panel.border = element_blank())+
  coord_flip()
ggsave(paste0(dir,"plot_out/TC_MPR-NMPR_gsva.pdf"),p,width =8,height =6)
