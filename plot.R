#!usr/bin/env Rscript

cat("Loading libraries...\n")
library(ggplot2)
library(gridExtra)

args = commandArgs(trailingOnly = T)

cat("Processing inputs...\n")

if(length(args)<4){
  stop("Need input file, population size, initial heteroplasmy, and repulsive halo (nonrepulsive = 0)!\n")
}

inputfile = args[1]
n = as.numeric(args[2])
h = as.numeric(args[3])
halo = as.numeric(args[4])

df = read.csv(inputfile, header = T)
MUT_RATE = unique(df$mut_rate)
TO_RATE = unique(df$to_rate)

plot.df.1 = df[df$rho == 0.05 & df$mut_rate == MUT_RATE[1] & df$to_rate == TO_RATE[1],]
plot.df.2 = df[df$rho == 0.05 & df$mut_rate == MUT_RATE[1] & df$to_rate == TO_RATE[2],]

vhmax = max(max(plot.df.1$vh),max(plot.df.2$vh),max(plot.df.3$vh))
bl = h*(1-h)/n
fn = scale_color_gradientn(colors = c("black","blue","white","red","black"), values = c(0,bl/(2*vhmax),bl/vhmax,3*bl/(2*vhmax),1), limits = c(0,vhmax))

p1.1 = ggplot(data = plot.df.1)+#fn+
  geom_tile(aes(x = p, y = q, fill = vh))+
  facet_wrap(~nseed)
p1.2 = ggplot(data = plot.df.2)+#fn+
  geom_tile(aes(x = p, y = q, fill = vh))+
  facet_wrap(~nseed)
  
filename = paste("vh",ifelse(halo>0,yes = "-repel",""),".png",sep = "")
res.factor = 3
png(filename, height = 1200*res.factor, width = 1200*res.factor, res = 72*res.factor)
grid.arrange(p1.1,p1.2,nrow = 3)
dev.off()
