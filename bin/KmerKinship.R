#!/usr/bin/R
suppressMessages(library(gplots))
suppressMessages(library(proxy))
message("\nLibraries gplots and proxy loaded successfully")
args <- commandArgs(trailingOnly = TRUE)

#Import data
message("\nBeginning data import, this may take a few minutes")
dat <- read.table(args[1],header=F,sep='\t')
#Remove the marker information, not needed for clustering
dat <- dat[c(-1,-2)]
#Import header
header <- readLines(args[2])
#Combine
colnames(dat) <- header
message("\nData import complete, calculating pairwise distance matrix")

#Heatmap calculation
mat_data <- data.matrix(dat)
dMat <- dist(mat_data, diag=TRUE, upper=TRUE, pairwise=TRUE, by_rows=FALSE)
ToPlot <- data.matrix(dMat)
row.hc <- hclust(dMat)
row.dd <- as.dendrogram(row.hc)

#Plotting
png(args[3], width=1960, height=1960)
heatmap.2(ToPlot, notecol="black", density.info="none", trace="none", dendrogram="row", key.xlab ="EuclideanDistance", col=bluered(100), keysize = 0.75, cexCol = 1.5, cexRow = 1.5, margins=c(8,8))
null <- dev.off()
message("Kmer kinship plotting complete.")
noquote(paste("Please check ",args[3]," and modify your pedigree file as needed\n"))
message("Note labels may not be well resolved. The top 20 clostest individuals are:\n")

#Obtain 20 clostest siblings
xy <- t(combn(colnames(ToPlot), 2))
xy <- data.frame(xy, dist=ToPlot[xy])
head(xy[order(xy$dist),],n=20)

write.table(xy[order(xy$dist),], file= args[4], quote=FALSE, sep='\t', col.names=FALSE) 
noquote(paste("The full list is printed in",args[4]))

