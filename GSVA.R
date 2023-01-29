##GSVA
require(GSVA)
require(GSEABase)
require(GSVAdata)
require(clusterProfiler)
data(c2BroadSets)
library(limma)
## kcdf: kcdf="Gaussian" for continuous and 'Poisson for integer counts'

GSVA_run <- function(SeuratObj,cluster_name,cluster1,cluster2,gmtFile,kcdf){
	hallgmt <- read.gmt(gmtFile)
	hall_list = split(hallgmt$gene, hallgmt$ont)
	expr=as.matrix(SeuratObj@assays$RNA@data)
	hall <- gsva(expr, hall_list, parallel.sz=10,kcdf=kcdf) 
	## get DE pathways
	group <- factor(SeuratObj@meta.data[,cluster_name],levels = c(cluster1,cluster2),ordered = F)
	design <- model.matrix(~0+group)
	colnames(design) <- c("C1","C2")
	rownames(design) <- colnames(hall)

	fit <- lmFit(hall,design)
	cont.matrix=makeContrasts('C1-C2',levels = design)
	fit2=contrasts.fit(fit,cont.matrix)
	fit2=eBayes(fit2)

	gs <- topTable(fit2,adjust='BH', number=Inf, p.value=0.05)
	gs$cluster <- ifelse(gs$t > 0 , "C1", "C2")
	return(gs)
}

GSVA_run2 <- function(SeuratObj,cluster_name,cluster1,cluster2,kcdf){
	expr=as.matrix(SeuratObj@assays$RNA@data)
	## change gene symbol to geneid
	gene_entrezid <- bitr(rownames(expr), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
	expression_filt <- expr[gene_entrezid$SYMBOL,]
	rownames(expression_filt) <- gene_entrezid$ENTREZID
	expression_filt <- as.matrix(expression_filt)

	hall <- gsva(expression_filt, c2BroadSets, parallel.sz=10,kcdf=kcdf) 
	## get DE pathways
	group <- factor(SeuratObj@meta.data[,cluster_name],levels = c(cluster1,cluster2),ordered = F)
	design <- model.matrix(~0+group)
	colnames(design) <- c("C1","C2")
	rownames(design) <- colnames(hall)

	fit <- lmFit(hall,design)
	cont.matrix=makeContrasts('C1-C2',levels = design)
	fit2=contrasts.fit(fit,cont.matrix)
	fit2=eBayes(fit2)

	gs <- topTable(fit2,adjust='BH', number=Inf, p.value=0.05)
	gs$cluster <- ifelse(gs$t > 0 , "C1", "C2")
	return(gs)
}
