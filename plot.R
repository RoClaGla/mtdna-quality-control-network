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

# Function to predict V(h) from E(h^2)-E(h)^2 with statistical simulations
vhest = function(p, q, h, n, pc, alpha, beta){
  wn0 = round(p*(1-h)*n)
  mn0 = round(q*h*n)
  wc0 = round((1-p)*(1-h)*n)
  mc0 = round((1-q)*h*n)
  
  Wn = 0:wn0
  Mn = 0:mn0
  Wc = 0:wc0
  Mc = 0:mc0
  
  INET = length(Wn)
  JNET = length(Mn)
  ICYT = length(Wc)
  JCYT = length(Mc)
  
  Eh  = 0
  Eh2 = 0
  En  = 0
  En2 = 0
  
  # Precalculate this, and grab element icyt,jcyt as you go
  dist.cyt =  dbinom(x = Wc, size = wc0, prob = pc)%*%
    t(dbinom(x = Mc, size = mc0, prob = pc))
  # For summing over u
	u = seq(from = 0, to = 1, length.out = 100)
  dist.u = dbeta(u,alpha,beta)
  deltau = u[2]-u[1]
  for(inet in 1:INET){
    dist.wn = dbinom(Wn[inet],wn0,u)
    for(jnet in 1:JNET){
      dist.mn = dbinom(Mn[jnet],mn0,u)
      for(icyt in 1:ICYT){
        for(jcyt in 1:JCYT){
          if(inet + jnet + icyt + jcyt == 4){
            # do nothing here!
          }else{
            dist.net = sum(deltau*dist.wn*dist.mn*dist.u)
            prefactor = (Mn[jnet]+Mc[jcyt])/(Wn[inet]+Wc[icyt]+Mn[jnet]+Mc[jcyt])
            
            Eh  = Eh  + prefactor*dist.cyt[icyt,jcyt]*dist.net
            Eh2 = Eh2 + prefactor^2*dist.cyt[icyt,jcyt]*dist.net
            
            En2 = En2 + (Wn[inet]+Wc[icyt]+Mn[jnet]+Mc[jcyt])^2*dist.net*dist.cyt[icyt,jcyt]
            En  = En  + (Wn[inet]+Wc[icyt]+Mn[jnet]+Mc[jcyt])*dist.net*dist.cyt[icyt,jcyt]
          }
        }
      }
    }
  }
  vn      = En2-En^2
  vhprime = (Eh2-Eh^2)/(h*(1-h))
  return(c(alpha, beta, vhprime, vn))
}

df = read.csv(inputfile, header = T)
MUT_RATE = unique(df$mut_rate)
TO_RATE = unique(df$to_rate)

plot.df.1 = df[df$t == 0 & df$rho == 0.05 & df$mut_rate == MUT_RATE[1] & df$to_rate == TO_RATE[1] & df$K == 0,]
plot.df.2 = df[df$t == 0 & df$rho == 0.05 & df$mut_rate == MUT_RATE[1] & df$to_rate == TO_RATE[1] & df$K == 10,]
plot.df.3 = df[df$t == 0 & df$rho == 0.05 & df$mut_rate == MUT_RATE[1] & df$to_rate == TO_RATE[1] & df$K == 20,]

plot.df.1$vhprime = plot.df.1$vh/(plot.df.1$mh*(1-plot.df.1$mh))
plot.df.2$vhprime = plot.df.2$vh/(plot.df.2$mh*(1-plot.df.2$mh))
plot.df.3$vhprime = plot.df.3$vh/(plot.df.3$mh*(1-plot.df.3$mh))

vhmax = max(max(plot.df.1$vhprime),max(plot.df.2$vhprime),max(plot.df.3$vhprime))
colfn = scale_fill_gradientn(colors = c("black","blue","white","red","black"),values = c(0,vhl/2,vhl,2*vhl,1),  limits = c(0,vhmax))

#vnl = nullret[4]/vnmax
#colfn2 = scale_fill_gradientn(colors = c("black","blue","white","red","black"),values = c(0,vnl/2,vnl,2*vnl,1), limits = c(0,vnmax))


fn = scale_fill_gradient2(low = "blue", high = "red", midpoint = 0.01)

p1.1 = ggplot(data = plot.df.1)+
  geom_tile(aes(x = p, y = q, fill = vhprime))+fn+
  facet_wrap(~nseed)
#p1.2 = ggplot(data = plot.df.1)+fn2+
#  geom_tile(aes(x = p, y = q, fill = vn))+
#  facet_wrap(~nseed)
p2.1 = ggplot(data = plot.df.2)+
  geom_tile(aes(x = p, y = q, fill = vhprime))+fn+
  facet_wrap(~nseed)
#p2.2 = ggplot(data = plot.df.2)+fn2+
#  geom_tile(aes(x = p, y = q, fill = vn))+
#  facet_wrap(~nseed)
p3.1 = ggplot(data = plot.df.3)+
  geom_tile(aes(x = p, y = q, fill = vhprime))+fn+
  facet_wrap(~nseed)

filename = paste("vh",ifelse(halo>0,yes = "-repel",""),".png",sep = "")
res.factor = 3
png(filename, height = 1200*res.factor, width = 1200*res.factor, res = 72*res.factor)
grid.arrange(p1.1,p2.1,p3.1,nrow = 3)
dev.off()
