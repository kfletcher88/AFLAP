#!/usr/bin/R
library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)

His <- read.table(args[1], sep=' ', header=FALSE)
jpeg(args[4])
ggplot(His, aes(x=V1,y=V2))+geom_line(colour="black")+coord_cartesian(xlim=c(0,250), ylim=c(0,4e6))+theme_minimal()+geom_vline(xintercept = c(as.numeric(args[3]), as.numeric(args[2])), linetype="dashed")+labs(x="Coverage", y="Frequency")
null <- dev.off()
