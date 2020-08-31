library(ggplot2)
library(ggrepel)
args <- commandArgs(trailingOnly = TRUE)

dat <- read.table(args[1], sep='\t', header=FALSE)

png(args[2],width=600,height=600)
ggplot(dat, aes(x=V3,y=V2,label=V1))+
geom_point()+
geom_text_repel()+
theme_minimal()+
labs(x="K-mer coverage", y="Marker Count")
dev.off()
