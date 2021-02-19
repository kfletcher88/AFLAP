#!/bin/R
packages <- c("ggplot2","ggrepel", "dplyr")
if (length(setdiff(packages, rownames(installed.packages())))>0){
write("One or more packages missing, please check ggplot2, ggrepel and dplyr are all installed", stderr())
} else {
write("All necessary packages found. Good to run AFLAP!", stderr())
}
