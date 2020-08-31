#!/bin/R
packages <- c("ggplot2")
if (length(setdiff(packages, rownames(installed.packages())))>0){write("One or more packages missing", stderr())}
