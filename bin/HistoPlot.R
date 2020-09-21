#!/usr/bin/R
null <- library(ggplot2)
null <- library(dplyr)
args <- commandArgs(trailingOnly = TRUE)

His <- read.table(args[1], sep=' ', header=FALSE)

#Set ymax
His2 <- His %>% filter(between(V1, as.numeric(args[2]), as.numeric(args[3])))
YM <- signif(max(His2$V2)*1.1, digits=2)

#Set xmax
XM <- round(as.numeric(args[3])*2.5)

jpeg(args[4])
ggplot(His, aes(x=V1,y=V2))+geom_line(colour="black")+coord_cartesian(xlim=c(0,XM), ylim=c(0,YM))+theme_minimal()+geom_vline(xintercept = c(as.numeric(args[3]), as.numeric(args[2])), linetype="dashed")+labs(x="Coverage", y="Frequency")
null <- dev.off()
