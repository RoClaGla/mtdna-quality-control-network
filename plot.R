#!usr/bin/env Rscript

cat("Loading libraries...")
library(ggplot2)
library(gridExtra)

args = commandArgs(trailingOnly = T)

cat("Processing inputs...")

if(length(args)<2){
  stop("Need input file and repulsive halo (nonrepulsive = 0) !")
}

inputfile = args[1]
halo = as.numeric(args[2])

df = read.csv(inputfile, header = T)

plot.df = df[df$halo == halo,]

p1.1 = ggplot(data = plot.df)+
  geom_line(aes(x = rho, y = mh, col = t))+
  facet_wrap(~nseed)
  
  filename = paste("test-",ifelse(halo>0,yes = "repel-",""),".csv",sep = "")
  res.factor = 3
  png(filename, height = 1200*res.factor, width = 1200*res.factor, res = 72*res.factor)
  p1.1
  dev.off()
