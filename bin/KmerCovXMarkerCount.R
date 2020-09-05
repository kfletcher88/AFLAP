library(ggplot2)
library(ggrepel)
args <- commandArgs(trailingOnly = TRUE)

dat <- read.table(args[1], sep='\t', header=FALSE)
dat$label <- dat$V1
dat$label[dat$V3 > 3] <- NA
png(args[2],width=600,height=600)
ggplot(dat, aes(x=V3,y=V2,label=label))+
geom_point()+
geom_text_repel(na.rm=TRUE)+
theme_minimal()+
labs(x="K-mer coverage", y="Marker Count")
null <- dev.off()

