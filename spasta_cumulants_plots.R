library(atminv)
library(Matrix)
library(ggplot2)
rm(list=ls())
set.seed(20)

save_images <- 0

## Interaction function
b <- function(s,u,p) {
    absp <- max(abs(p),0.2)
    absp*sqrt(2*pi) * dnorm(u,mean = s, sd =absp) *
        ((sign(p) == sign(u-s)) | (u-s) == -1)
}

## space axis
s_axis <- seq(-10,10,length=200)
ds <- mean(diff(s_axis))

## Choose point of reference (origin)
s_o = 0

## index of origin (halfway in domain)
c_idx <- which.min(abs(s_axis - s_o))

## number of spatial locations
ns <- length(s_axis)

## construct spatial precision matrix to simulate random field for interaction function
Q_s <- GMRF_RW(n = length(s_axis),
               order = 2,
               precinc = 2000)@Q + 
    0.001*Diagonal(length(s_axis))  # spatial precision

## and construct GMRF object for interaction function
G <- GMRF(mu = matrix(rep(0,nrow(Q_s))),   # GMRF with final precision
          Q = Q_s,
          n=nrow(Q_s))

## Now sample field from GMRF
p <- sample_GMRF(G,reps = 1)

## Interaction matrix initialisation
B <- matrix(0,ns,ns)

## The interaction matrix is b evaluated at the gridcell centre * the bin width
for(i in seq_along(s_axis)) {
    B[i,] <-  b(s = s_axis[i], 
                u = s_axis,
                p = p[i])*ds
}

## Interaction matrix at origin (s = 0)
b_o <- B[c_idx,]

## Now define the covariance function for the latent field
flux_cov_fn <- function(D,scale,nu) {
    exp(-scale*(D^nu))    
}

## and find the distance matrix
D <- sp::spDists(matrix(s_axis),
                 matrix(s_axis))

## Construct the mean vector on the Gaussian scale (mean of -2)
mu_f_trans_true <- matrix(rep(-2,length(s_axis))) 

## Construct the covariance matrix on the Gaussian scale
S_f_trans_true <- flux_cov_fn(D, scale = 0.8,nu=1.7)

## Construct the mean vector on the lognormal scale
mu_f_exp <- exp(mu_f_trans_true + diag(S_f_trans_true)/2)

## Construct the covariance matrix on the lognormal scale
S_f_exp <- tcrossprod(mu_f_exp)*(exp(S_f_trans_true) - 1)

## Visualise C22
image(B %*% S_f_exp %*% t(B))

## Now custruct the skewness matrix kappa_{211}(0,u2,u3)
##
## kappa_{211}(0,u2,u3) = \int b(0,u_1) k_{111}(u_1,u_2,u_3) d u_1
##
## kappa_{211}(0,u2,u3) \approx \sum_i b_o(0,u_i)k_{111}(u_{1i},u_2,u_3)
## 
## where b_o is the row of the interaction matrix extracted at the origin (includes
## within it the integration weight, ds
##
## In the loop below we evaluate each component of the sum, that is we evaluate the summand for each u_{1i}
skewness211 <- matrix(0,ns,ns)

## for each element in sum
for(i in 1:ns) {
    ## k_{111}(u_{1i},u_2,u_3) = k_1(u_{1i})k_1(u_2)k_1(u_3) * (
    ##                           exp(C_{11}(u_{1i},u_2) + C_{11}(u_{1i},u_3) + C_{11}(u_2,u_3)) -
    ##                           exp(C_{11}(u_{1i},u_2)) - exp(C_{11}(u_{1i},u_3)) -
    ##                           exp(C_{11}(u_2,u_3)) + 2)
    ## C_{11}(u_{1i},u_2) has repeated columns
    C11_u1_u2 <- do.call("rbind",lapply(1:ns,function(x) S_f_trans_true[i,]))
    ## C_{11}(u_{1i},u_3) has repeated rows
    C11_u1_u3 <- do.call("cbind",lapply(1:ns,function(x) S_f_trans_true[i,]))
    
    ## Now just compute the summand for this iteration i
    skewness211 <- skewness211 + b_o[i]*(tcrossprod(mu_f_exp) * mu_f_exp[i]*
                                  (exp(S_f_trans_true  +  C11_u1_u2 + C11_u1_u3) -
                                       exp(S_f_trans_true) - exp(C11_u1_u2) - exp(C11_u1_u3) + 2))
}

## Now custruct the skewness matrix kappa_{222}(0,s2,s3). To do this we write out the integral as
## k_{222}(0,s_2,s_3) = \int b(0,u_1) [ J ] d u_1
## where J is a double integral which is a function of u_1, namely
## J = \int b(s_2,u_2)b(s_3,u_3)k_{111}(u_1,u_2,u_3) d u_2 d u_3
##
## After discretisation, for each u_1 this is simply B * K * t(B) where K is the 
## discretised cumulant function as above. We then expand the outer integral as a 
## sum to obtain something that is very similar in form as skewness211
skewness222 <- matrix(0,ns,ns)
for(i in 1:ns) {
    C11_u1_u2 <- do.call("rbind",lapply(1:ns,function(x) S_f_trans_true[i,]))
    C11_u1_u3 <- do.call("cbind",lapply(1:ns,function(x) S_f_trans_true[i,]))
    
    skewness222 <- skewness222 + b_o[i]*B %*% (tcrossprod(mu_f_exp) * mu_f_exp[i]*
                                             (exp(S_f_trans_true  +  C11_u1_u2 + C11_u1_u3) -
                                                  exp(S_f_trans_true) - exp(C11_u1_u2) - exp(C11_u1_u3) + 2)) %*% t(B)
}

## Plot results
library(fields)
bw<- designer.colors(6, c( "white", "black"))
df_s1 <- data.frame(s = rep(s_axis,2))
df_s1$lines <- c(as.vector(b_o),as.vector(b_o %*% S_f_exp))
df_s1$linetype <- c(rep("b(0,u)",ns),rep("C21(0,u)",ns))
obs_pt <- geom_point(data=df_s1[c_idx,],
                     aes(x=s, y=0),
                     pch=5,cex=4)
g0 <- LinePlotTheme() + 
    geom_line(data=df_s1,aes(x=s,y=lines,linetype=linetype)) + 
    ylab("") + 
    scale_linetype(guide=guide_legend(title=""),
                   labels=c(expression(b(0,u[2])), 
                            expression({kappa^2}[Y[2]][Y[1]](0,u[2])))) +     
    #scale_linetype_manual(guide=guide_legend(title="")) +
    xlab(expression(u[2])) + obs_pt
    


df_s2 <- expand.grid(s1 = s_axis,s2=s_axis)
df_s2$skew211 <- c(skewness211)
df_s2$skew222 <- c(skewness222)
df_s2$b <- as.vector(B)

g02 <-LinePlotTheme() + geom_tile(data=df_s2,aes(s1,s2,fill=b)) + 
    scale_fill_gradientn(colours=bw,guide=guide_legend(title=expression(b(s,u)))) +
    coord_fixed(xlim=c(-10,10),ylim=c(-10,10)) +
    xlab(expression(s)) + ylab(expression(u))
g1 <- LinePlotTheme() + geom_tile(data=df_s2,aes(s1,s2,fill=skew211)) + 
    scale_fill_gradientn(colours=bw,guide=guide_legend(title=expression({kappa^3}[Y[2]][Y[1]][Y[1]](0,u[2],u[3])))) +
    coord_fixed(xlim=c(-10,10),ylim=c(-10,10)) +
    xlab(expression(u[2])) + ylab(expression(u[3])) + obs_pt
g2<- LinePlotTheme() + geom_tile(data=df_s2,aes(s1,s2,fill=skew222)) + 
    scale_fill_gradientn(colours=bw,guide=guide_legend(title=expression({kappa^3}[Y[2]][Y[2]][Y[2]](0,s[2],s[3])))) +
    coord_fixed(xlim=c(-10,10),ylim=c(-10,10)) +
    xlab(expression(s[2])) + ylab(expression(s[3])) + obs_pt
library(gridExtra)
g_all <- grid.arrange(g02,g1,g0,g2,ncol=2)
if(save_images) ggsave(filename = "./img/Fig1_cumulants_plots.pdf",plot=g_all,width=12,height=8)

## Deprecated kappa construction
# skewness211 <- matrix(0,ns,ns)
# for(i in 1:ns) {
#     C11_s_v <- do.call("rbind",lapply(1:ns,function(x) S_f_trans_true[i,]))
#     C11_u_v <- do.call("cbind",lapply(1:ns,function(x) S_f_trans_true[i,]))
#     skewness211[i,] <- b %*% (tcrossprod(mu_f_exp) * mu_f_exp[i] *
#                                     (exp(S_f_trans_true  +  C11_s_v + C11_u_v) -
#                                          exp(S_f_trans_true) - exp(C11_s_v) - exp(C11_u_v) + 2))
# }

## Deprecated Plot 
# plot.new()
# par(mfrow=c(1,3))
# plot(s_axis,b,type="l",xlim=c(-10,10),ylab="",xlab="u")
# lines(s_o,0,type="o",pch=23,cex=2,lwd=4)
# #arrows(s_obs[10,]$s + 0.2,0,y1=200)
# par(new=TRUE)
# plot(s_axis,t(b %*% S_f_exp),,type="l",lty=5,xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(-10,10))
# axis(4)
# mtext("C_21(0,u)",side=4,line=3)
# mtext("b(0,u)",side=2,line=3)
#legend(c(-4,120000),lty=c(1,5),legend=c("b(0,u)","C_21(0,u)"))
# Note the peak is to the left of the observation location (asymmetry). There is nonnegativity to the right because the spatial field is correlated, therefore although not unobserved, flux field to the right affects the mole fraction too.
#S_f_trans <- flux_cov_fn(D, scale = 0.2,nu=1)

