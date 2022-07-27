#!usr/bin/env Rscript

cat("Loading libraries...")
library(ggplot2)
library(gridExtra)

args = commandArgs(trailingOnly = T)

cat("Processing inputs...")

if(length(args)<){
  stop("Need input file")
}

inputfile = args[1]

df = read.csv(inputfile)

p1.1 = ggplot(data = df)+
  geom_line(aes(x = rho, y = mh, col = t))+
  facet_wrap(~nseed)
  
  filename = paste("test-",ifelse(halo>0,yes = "repel-",""),".csv",sep = "")
  res.factor = 3
  png(filename, height = 1200*res.factor, width = 1200*res.factor, res = 72*res.factor)
  p1.1
  dev.off()
