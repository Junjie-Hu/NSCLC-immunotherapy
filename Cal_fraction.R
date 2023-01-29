library(tidyr)
library(reshape)
library(ggthemes)
library(rcompanion)
library(GGally)
library(ggrepel)
library(qdapTools)
library(REdaS)

### calculate cellular fraction
cal_fraction <- function(SeuratObj,cluster){
meta.temp <- SeuratObj@meta.data[,c(cluster, "Group")]
# Create list to store frequency tables 
prop.table.error <- list()
for(i in 1:length(unique(SeuratObj@meta.data$Group))){
vec.temp <- meta.temp[SeuratObj@meta.data$Group==unique(SeuratObj@meta.data$Group)[i],cluster]
# Convert to counts and calculate 95% CI 
# Store in list 
table.temp <- freqCI(vec.temp, level = c(.95))
prop.table.error[[i]] <- print(table.temp, percent = TRUE, digits = 3)
# 
}
# Name list 
names(prop.table.error) <- unique(meta.temp$Group)
# Convert to data frame 
tab.1 <- as.data.frame.array(do.call(rbind, prop.table.error))
# Add Group column 
b <- c()
a <- c()
for(i in names(prop.table.error)){
  a <- rep(i,nrow(prop.table.error[[1]]))
  b <- c(b,a)
}
tab.1$Group <- b
# Add common cell names 
tab.1$cell <- gsub("\\.[0-9]","",row.names(tab.1))
tab.1$cell <- gsub("\\."," ",tab.1$cell)
# Resort factor Group 
tab.1$Group <- factor(tab.1$Group, levels = c("TN", "MPR", "NMPR"))
# Rename percentile columns 
colnames(tab.1)[1] <- "lower"
colnames(tab.1)[3] <- "upper"
return(tab.1)
}