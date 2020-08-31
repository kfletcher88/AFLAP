#!/usr/bin/R
library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)

M61 <- read.table(args[1], sep=' ', header=FALSE)
M62 <- read.table(args[2], sep=' ', header=FALSE)
Mall <- read.table(args[3], sep=' ', header=FALSE)
M61sum <- sum(M61$V2)
M61$den <- M61$V2/M61sum
M62sum <- sum(M62$V2)
M62$den <- M62$V2/M62sum
Mallsum <- sum(Mall$V2)
Mall$den <- Mall$V2/Mallsum


png(args[4],width=350,height=250)
ggplot()+geom_point(data=M61, aes(x=V1,y=den, color="k+k-1"))+
geom_point(data=M62, aes(x=V1,y=den, color=">k+k-1"))+
geom_line(data=Mall, aes(x=V1, y=den), color="black")+
theme_minimal()+
labs(x="Marker Presence", y="Marker Density", Color="Marker Class")
null <- dev.off()
