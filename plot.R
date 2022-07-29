#!usr/bin/env Rscript

cat("Loading libraries...\n")
library(ggplot2)
library(gridExtra)

args = commandArgs(trailingOnly = T)

cat("Processing inputs...\n")

if(length(args)<2){
  stop("Need input file and repulsive halo (nonrepulsive = 0)!\n")
}

inputfile = args[1]
halo = as.numeric(args[2])

df = read.csv(inputfile, header = T)

p1.1 = ggplot(data = plot.df)+
  geom_line(aes(x = p, y = q, col = as.factor(t)))+
  facet_wrap(rho~nseed)
  
filename = paste("mh-",ifelse(halo>0,yes = "repel-",""),".png",sep = "")
res.factor = 3
png(filename, height = 1200*res.factor, width = 1200*res.factor, res = 72*res.factor)
grid.arrange(p1.1,p1.2,p1.3,nrow = 3)
dev.off()
