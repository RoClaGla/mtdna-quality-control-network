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

plot.df.1 = df[df$halo == halo & df$p == 0 & df$q == 0, ]
plot.df.2 = df[df$halo == halo & df$p == 1 & df$q == 1, ]
plot.df.3 = df[df$halo == halo & df$p == 1 & df$q == 0, ]

p1.1 = ggplot(data = plot.df.1)+
  geom_point(aes(x = rho, y = mh, col = as.factor(t)))+
  facet_wrap(~nseed)
p1.2 = ggplot(data = plot.df.2)+
  geom_point(aes(x = rho, y = mh, col = as.factor(t)))+
  facet_wrap(~nseed)
p1.3 = ggplot(data = plot.df.3)+
  geom_point(aes(x = rho, y = mh, col = as.factor(t)))+
  facet_wrap(~nseed)
  
  filename = paste("test-",ifelse(halo>0,yes = "repel-",""),".png",sep = "")
  res.factor = 3
  png(filename, height = 1200*res.factor, width = 1200*res.factor, res = 72*res.factor)
  grid.arrange(p1.1,p1.2,p1.3,nrow = 3)
  dev.off()
