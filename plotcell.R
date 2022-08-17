#!/usr/bin/env Rscript
# code to plot individual cell snapshots that are output from simulation

cat("Loading libraries...")

library(ggplot2)
library(gridExtra)

# helper function to draw a circle using ggplot
# this is joran's solution from https://stackoverflow.com/questions/6862742/draw-a-circle-with-ggplot2
circleFun <- function(center = c(0,0),diameter = 2, npoints = 100){
    r = diameter / 2
    tt <- seq(0,2*pi,length.out = npoints)
    xx <- center[1] + r * cos(tt)
    yy <- center[2] + r * sin(tt)
    return(data.frame(x = xx, y = yy))
}
circledat <- circleFun()

# function to produce a visualisation given a particular label for simulation output 
makeplot = function(string, rho, titlestr="", tag = "") {
  # read network and mtDNA position data
  net.df = read.csv(paste(c("network-", string, ".csv"), collapse=""))
  dna.df = read.csv(paste(c("mtdna-", string, ".csv"), collapse=""))
  
  circledat2 = circleFun(center = c(dna.df$x[1],dna.df$y[1]), diameter = 2*rho)
    
  # return the plot
  return(ggplot() +
    geom_path(data = circledat, aes(x,y)) +                                         # cell boundary
    geom_path(data = circledat2, aes(x,y)) +                                         # cell boundary
    geom_segment(data = net.df, aes(x=xs,y=ys,xend=xe,yend=ye), color="#888888") +  # network
    geom_point(data = dna.df, aes(x=x, y=y, color=factor(type))) +                  # mtDNAs
    theme_void() +
    theme(legend.position = "none") +
    ggtitle(titlestr)+
    labs(tag = tag)+
    theme(plot.title = element_text(size = rel(1.25))))
}

# strings in output filenames follow the format
# h, nseed, p, q, lambda, halo, perturb
plot.1 = makeplot("100-4-50-0.50-0.50-0.10-0.25", 0.25, "Small s,\n p=q=0.5,\n repel,\n rho=0.25", "A")
plot.2 = makeplot("100-4-50-1.00-0.00-0.00-0.25", 0.25, "Small s,\n p=1,q=0,\n no repel\n, rho=0.25", "B")
plot.3 = makeplot("100-16-50-0.50-0.50-0.10-0.25", 0.25, "Medium s,\n p=q=0.5,\n repel,\n rho=0.25", "C")
plot.4 = makeplot("100-16-50-1.00-1.00-0.00-0.25", 0.25, "Medium s,\n p=q=0.5,\n no repel,\n rho=0.25", "B")
plot.5 = makeplot("100-64-50-0.00-0.00-0.00-0.00", 0   , "Large s,\n p=0,q=0,\n no repel,\n rho=0", "C")
plot.6 = makeplot("100-64-50-0.50-0.50-0.10-0.25", 0.25, "Large s,\n p=q=0.5,\n repel,\n rho=0.25", "D")
#plot.7 = makeplot("0.50-64-0.50-0.50-0.00-0.00-2", "Large s,\n p=q=0.5", "E")
#plot.8 = makeplot("0.50-16-1.00-1.00-0.00-0.00-2", "Medium s,\n p=q=1", "F")
#plot.9 = makeplot("0.50-16-1.00-1.00-0.04-0.00-2", "Medium s,\n p=q=1,\n Medium λ", "G")
#plot.10 = makeplot("0.50-16-1.00-1.00-0.10-0.00-2", "Medium s,\n p=q=1,\n Large λ","H")
#plot.11 = makeplot("0.50-16-1.00-1.00-0.00-0.10-2", "Medium s,\n p=q=1,\n Repel","I")

# bump to output file
res.factor = 3
png("plotcell.png", width=1300*res.factor, height=500*res.factor, res=72*res.factor)
grid.arrange(plot.1, plot.2, plot.3, plot.4, plot.5, plot.6, nrow = 2)#, plot.7, plot.8, plot.9, nrow=2)
dev.off()
