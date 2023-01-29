library(Seurat)

dir <- "./"

RNA <- readRDS("Neutrophil2_new.rds")

### input data prepare
cellInfo <- RNA@meta.data[,c("N_cluster_annotation","nFeature_RNA","nCount_RNA")]
colnames(cellInfo) <- c("CellType","nGene","nUMI")
cellInfo$CellType <- as.character(cellInfo$CellType)
dir.create("SCENIC")
saveRDS(cellInfo, file=paste0(dir,"SCENIC/N_cellInfo.Rds"))
colVars <- list(CellType=c("Neu_OSM" = "#f8766d","Neu_S100A12" = "#7cae00", 
  "Neu_CCL3" = "#00bfc4","Neu_IFIT3"="#c77cff"))
colVars$CellType <- colVars$CellType[intersect(names(colVars$CellType), cellInfo$CellType)]
saveRDS(colVars, file=paste0(dir,"SCENIC/N_colVars.Rds"))

exprMat <- as.matrix(RNA@assays$RNA@counts)

### Initialize SCENIC settings
library(SCENIC)
org <- "hgnc" # or hgnc, or dmel
dbDir <- "/media/inspur/AS2150G2/HJJ/scrna/cisTarget_databases" # RcisTarget databases location
myDatasetTitle <- "Fibro" # choose a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=20)

scenicOptions@inputDatasetInfo$cellInfo <- "SCENIC/N_cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "SCENIC/N_colVars.Rds"

saveRDS(scenicOptions, file=paste0(dir,"SCENIC/N_scenicOptions.Rds")) 

### Co-expression network
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)

##  check whether any known relevant genes are filtered-out
interestingGenes <- c("NR4A1", "ID2", "CEBPB")
# any missing?
interestingGenes[which(!interestingGenes %in% genesKept)]

exprMat_filtered <- exprMat[genesKept, ]
dim(exprMat_filtered)

## Correlation
runCorrelation(exprMat_filtered, scenicOptions)

## GENIE3
# Optional: add log (if it is not logged/normalized already)
exprMat_filtered_log <- log2(exprMat_filtered+1)

# Run GENIE3 for 3K-5K cells, long time
runGenie3(exprMat_filtered_log, scenicOptions)
# for large cells, run GRNBoost (in Python) instead of GENIE3
# exportsForArboreto(exprMat,scenicOption,dir = "./")

### Build and score the GRN (runSCENIC_…)
exprMat_log <- log2(exprMat+1)

scenicOptions <- readRDS("SCENIC/N_scenicOptions.Rds")
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 1
scenicOptions@settings$seed <- 123

# For a very quick run: 
# coexMethod=c("top5perTarget")
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # For toy run
# save...

scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) #** Only for toy run!!
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)

## Regulators for known cell types or clusters
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$CellType),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
pdf("SCENIC/N_SCENIC.pdf",width =5,height =6)
ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity",
	 column_names_side = "top",row_names_side = "right")
dev.off()

topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
viewTable(topRegulators)

saveRDS(scenicOptions, file=paste0(dir,"SCENIC/N_scenicOptions.Rds"))
saveRDS(regulonActivity_byCellType_Scaled, file=paste0(dir,"SCENIC/N_regulonActivity_byCellType_Scaled.Rds"))

## Projection the AUC and TF expression onto t-SNEs
aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat)
savedSelections <- shiny::runApp(aucellApp)
print(tsneFileName(scenicOptions))
tSNE_scenic <- readRDS(tsneFileName(scenicOptions))
aucell_regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
# Show TF expression:
TF_check <- c("SPI1", "Sox10", "Sox9","Irf1", "Stat6")
par(mfrow=c(2,3))
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, exprMat, aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[TF_check],], plots="Expression")
# Save AUC as PDF:
Cairo::CairoPDF("output/Step4_BinaryRegulonActivity_tSNE_colByAUC.pdf", width=20, height=15)
par(mfrow=c(4,6))
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, cellsAUC=aucell_regulonAUC, plots="AUC")
dev.off()

saveRDS(scenicOptions, file="SCENIC/N_scenicOptions.Rds") # To save status


