library(EFDR)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)

miss_val = -999
ymin = 36.469
xmin = -14.124
xco <- 100
yco <- 0
dx = 0.352
dy = 0.234
yx <- as.data.frame(cbind(expand.grid(y = seq(ymin, ymin + 127*dy,length=128),
                                      x = seq(xmin, xmin + 127*dx,length=128))))

yx$breaks_x <- cut(yx$x,seq(xmin-0.001, xmin + 127*dx + 0.001,length=64),labels=FALSE)
yx$breaks_y <- cut(yx$y,seq(ymin-0.001, ymin + 127*dy + 0.001,length=64),labels=FALSE)

## LOAD MAPS
Europe <- map_data("world") %>%
    mutate(id = group, x=long,y=lat)

yx <- yx %>%
    mutate(z = read.table("../data/CH4_emissions_scaled_natural",skip=10)[,1], t = 1)

Z <- sqrt(abs(EFDR::df.to.mat(yx[c("x","y","z")])))

alpha =     0.05        # 5% significant level
wf =        "la8"       # wavelet name
J =         3           # 3 resolutions
b =         11          # 11 neighbours for EFDR tests
n.hyp =     c(5,10,20,30,40,60,80,100,400,800,1200,2400,4096)  # find optimal number of tests in 
# EFDR from this list
iteration = 100         # number of iterations in MC run when finding optimal n in n.hyp
# distance interpolation
parallel =  parallel::detectCores()/2 # use half the number of available cores
set.seed(20)     

Z1 <- test.bonferroni(Z)
Z2 <- test.fdr(Z)
Z3 <- test.efdr(Z,wf="la8",J=J, alpha = alpha, n.hyp = n.hyp, 
                b=b,iteration=iteration,parallel = parallel)

yx <- mutate(yx,
             Z1 = c(t(Z1$Z)),
             Z2 = c(t(Z2$Z)),
             Z3 = c(t(Z3$Z)))

g0 <- ggplot() + geom_tile(data=yx,aes(x=x,y=y,fill=sqrt(abs(z)))) +
    geom_path(data=Europe,aes(x=long,y=lat,group=group)) +
    scale_fill_gradient2(low="blue",high="red", mid="#FFFFCC", guide = guide_legend(title="sqrt(g/s)")) +
    coord_fixed(xlim=c(xmin,30.58),ylim=c(36.469,66.187)) +
    ggtitle("Emissions")


g1 <- ggplot() + geom_tile(data=yx,aes(x=x,y=y,fill=Z1)) +
    geom_path(data=Europe,aes(x=long,y=lat,group=group)) +
    scale_fill_gradient2(low="white",high="red", mid="#FFFFCC", guide = guide_legend(title="sqrt(g/s)")) +
    coord_fixed(xlim=c(xmin,30.58),ylim=c(36.469,66.187)) +
    ggtitle("Bonferroni")

g2 <- ggplot() + geom_tile(data=yx,aes(x=x,y=y,fill=Z2)) +
    geom_path(data=Europe,aes(x=long,y=lat,group=group)) +
    scale_fill_gradient2(low="white",high="red", mid="#FFFFCC", guide = guide_legend(title="sqrt(g/s)")) +
    coord_fixed(xlim=c(xmin,30.58),ylim=c(36.469,66.187)) +
    ggtitle("FDR")

g3 <- ggplot() + geom_tile(data=yx,aes(x=x,y=y,fill=Z3)) +
    geom_path(data=Europe,aes(x=long,y=lat,group=group)) +
    scale_fill_gradient2(low="white",high="red", mid="#FFFFCC", guide = guide_legend(title="sqrt(g/s)")) +
    coord_fixed(xlim=c(xmin,30.58),ylim=c(36.469,66.187)) +
    ggtitle("EFDR")

g_all <- arrangeGrob(g0,g1,g2,g3,ncol=2)
