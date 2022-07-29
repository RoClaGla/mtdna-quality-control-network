#!usr/bin/env Rscript

cat("Loading libraries...\n")
library(ggplot2)
library(gridExtra)

args = commandArgs(trailingOnly = T)

cat("Processing inputs...\n")

if(length(args)<3){
  stop("Need input file, population size, and repulsive halo (nonrepulsive = 0)!\n")
}

inputfile = args[1]
n = as.numeric(args[2])
halo = as.numeric(args[3])

df = read.csv(inputfile, header = T)
MUT_RATE = unique(df$mut_rate)
TO_RATE = unique(df$to_rate)

plot.df.1 = df[df$rho == 0.16 & df$mut_rate == MUT_RATE[1] & df$to_rate == TO_RATE[1],]
plot.df.2 = df[df$rho == 0.16 & df$mut_rate == MUT_RATE[2] & df$to_rate == TO_RATE[2],]
plot.df.3 = df[df$rho == 0.16 & df$mut_rate == MUT_RATE[3] & df$to_rate == TO_RATE[3],]

vhmax = max(max(plot.df.1$vh),max(plot.df.2$vh),max(plot.df.3$vh))
bl = 1/n
fn = scale_color_gradientn(colors = c("black","blue","white","red","black"), values = c(0,0.1,0.2,0.5,1), limits = c(0,vhmax))

p1.1 = ggplot(data = plot.df.1)+fn+
  geom_tile(aes(x = p, y = q, fill = vh/(mh*(1-mh))))+
  facet_wrap(to_rate~nseed)
p1.2 = ggplot(data = plot.df.2)+fn+
  geom_tile(aes(x = p, y = q, fill = vh/(mh*(1-mh))))+
  facet_wrap(to_rate~nseed)
p1.3 = ggplot(data = plot.df.3)+fn+
  geom_tile(aes(x = p, y = q, fill = vh/(mh*(1-mh))))+
  facet_wrap(to_rate~nseed)
  
filename = paste("vh",ifelse(halo>0,yes = "-repel",""),".png",sep = "")
res.factor = 3
png(filename, height = 1200*res.factor, width = 1200*res.factor, res = 72*res.factor)
grid.arrange(p1.1,p1.2,p1.3,nrow = 3)
dev.off()
