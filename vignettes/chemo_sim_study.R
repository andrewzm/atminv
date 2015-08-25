## ----include=FALSE--------------------
library(knitr)
opts_knit$set(width=40)

## ----installPkgs,eval=FALSE-----------
## library(devtools)
## install_github("andrewzm/hmc")
## install_github("andrewzm/atminv")

## ----loadPkgs,message=FALSE-----------
library(ggplot2)
library(grid)
library(dplyr)
library(Matrix)
library(tidyr)
library(gstat)
library(hmc)
library(atminv)
# load_all("../../../pkgs//hmc")
#library(devtools)
#load_all("..")
#library(scales) # for format_format
#load_all("../../../../../CurrentProjects/PostDoc Bristol/R Code/pkg/MVST")

## ----ST-grid--------------------------
###------------------
### Parameters
###------------------
model = "full" ## Either sparse or full or full_big or diag
misspecification = 1
#set.seed(25) # 25 , T = 400 or 15/200
ds <- 0.2            # 0.2 spacing. 
                     # NB: If we change this we need to 
                     # change the round() command further down
smin = -10 + ds/2    # first gridcell centre
smax =  10 - ds/2    # last gridcell centre

s_axis <- round(seq(smin,smax,by=ds),1) # create s-axis
if(model %in% c("full","diag")) {
  t_axis <- 1:100                       # create t-axis
  m_obs <- 6                            # number of obs. (including val.)
} else {
  t_axis <- 1:100                       # create t-axis
  m_obs <- 1000                         # number of obs.
}

ns <- length(s_axis)                    # no. of gridcells
nt <- length(t_axis)                    # no. of time points
st_grid <- expand.grid(s=s_axis,        # ST-grid (long format)
                       t=t_axis) %>% 
            data.frame()

sigma_eps_true <- 10                    # observation error std.
if(model %in% c("full","full_big")) {
  sigma_zeta_true <- 50                 # discrepancy marginal std.
} else {
  sigma_zeta_true <- 10                 # discrepancy marginal std.
}

theta_t_true <- 0.8     # temporal correlation parameter ('a' in text)
theta_s_true <- 1       # spatial range parameter ('d' in text)

## ----bfun-----------------------------
###------------------
### Transition kernel 
###------------------
# p is a ST process and reflects the std of the truncated Gaussian.
# This problem ONLY works if b is of relatively local scope. Once we 
# have b which has a very large scope we get oscillations/instability

b <- function(s,u,p) {
  absp <- max(abs(p),0.2)
  absp*sqrt(2*pi) * dnorm(u,mean = s, sd =absp) *
    ((sign(p) == sign(u-s)) | (u-s) == 0) *
    (abs(u - s) < absp)
}

## ----Wind-----------------------------
## Sample the "wind" vector
Q_s <- GMRF_RW(n = length(s_axis),
               order = 2,
               precinc = 2000)@Q + 
          0.001*.symDiagonal(length(s_axis))  # spatial precision
Q_t <- GMRF_RW(n = length(t_axis),
               order = 1,
               precinc = 20)@Q +   
          0.1*.symDiagonal(length(t_axis))    # temporal precision

Q_full = as(kronecker(Q_t,Q_s),"dgCMatrix")   # spatio-temporal precision

G <- GMRF(mu = matrix(rep(0,nrow(Q_full))),   # GMRF with final precision
          Q = Q_full,n=nrow(Q_full))

# Load the seed we used to simulate this parameter
data(sim.Random.seed)
#load("~/Desktop/Chemometrics_results/sim.Random.seed.rda")

# Now simulate this parameter by sampling from the GMRF
st_grid$p <- sample_GMRF(G,reps = 1)

## ----upsilon-plot1,fig.keep=FALSE,eval=FALSE----
## print(LinePlotTheme() +
##         geom_tile(data=st_grid,aes(s,t,fill=p)) +
##         scale_fill_gradient2(low="blue",high="red",mid="white")
##       )

## ----upsilon-plot2,fig.height=6,fig.width=5----
print(LinePlotTheme() + 
        geom_tile(data=st_grid,aes(s,t,fill=p)) +
        scale_fill_gradient2(low="blue",high="red",mid="white")
      )

## ----SRR1,include=FALSE,eval=FALSE----
## ## Simulate the SRR at one location, s = 9.1
## s_sim <- 9.1
## Z <- matrix(0,nt,ns)
## for( i in 1 : nrow(Z)) {
##   Z[i,] <- b(s = s_sim,
##              u = s_axis,
##              p = as.numeric(subset(st_grid,s== s_sim & t==i)$p))
## }

## ----sim-Yf---------------------------
###------------------
### Lognormal flux field
###------------------
set.seed(1)
variogram_model <- vgm(range = 3.334,            # Construct spherical semi-variogram
                              nugget = 0.00533,
                              psill = 0.80429,
                              model= "Sph")
S_f_log <- variogramLine(variogram_model,        # Find covariance matrix Sigma_f
                         dist_vector = sp::spDists(matrix(s_axis),
                                                 matrix(s_axis)),
                          covariance = TRUE)
mu_f_log <- matrix(rep(5,length(s_axis)))        # Construct mu_f
Yf_sim <- exp(mu_f_log + t(chol(S_f_log)) %*%    # Simulate Y_f
                rnorm(n = length(s_axis)))
st_grid <- st_grid %>%                           # Append Y_f to data frame
          left_join(data.frame(s=s_axis,
                               Yf = Yf_sim))

if(misspecification)  {                           # If we are assuming misspecification
    S_f_log <- diag(diag(S_f_log))                # We over-write Sigma_f to be diagonal
}

## -------------------------------------
###------------------
### Mole fraction field
###------------------
## Since we cannot put the discrepancy everywhere (too large), 
## we will just put it at the observation locations
st_grid <- st_grid %>%
  group_by(t,s,p) %>% # for each space-time location
  summarise(Yf = Yf,  # find the mole-fraction by finding the SRR
            Ym = sum(b(s =s, u = s_axis,p = p) * 
                       Yf_sim * ds))

## -------------------------------------
###------------------
### Observations
###------------------
if(model %in% c("full","diag")) {
      s_obs <- data.frame(s = sample(s_axis[-c(1:10,(ns-10):ns)],
                                     size = m_obs,
                                     replace=F),
                          m = 1:m_obs  )
    new_obs <- 0.3
    s_obs[6,]$s <- new_obs
} else {
  s_obs <- data.frame(s = sample(s_axis,size = m_obs, replace=T),
                      m = 1:m_obs)
}

## -------------------------------------
s_obs <- s_obs %>%
            left_join(st_grid) %>%
            mutate(z = Ym + rnorm(n = length(Ym), 
                                  sd = sigma_eps_true)) %>%
            arrange(t,s)
Qobs <- sigma_eps_true^(-2) * .symDiagonal(nrow(s_obs))

## -------------------------------------
## Now add the discrepancy
if(model %in% c("full","diag")) {

    corr_zeta_true <- corr_zeta_fn(s_obs$s[1:m_obs],
                                   t_axis,
                                   theta_t_true, 
                                   theta_s_true)
    S_zeta_true <- sigma_zeta_true^2 * corr_zeta_true
    chol_S_zeta_true <- chol(S_zeta_true)

    s_obs <- s_obs %>%
      mutate(dis = t(chol_S_zeta_true) %*%  rnorm(n = nrow(s_obs)),
             z = z + dis)
    
    C_m <- .symDiagonal(nrow(s_obs))
} else if(model == "full_big") {

    corr_s_mat <- function(theta_s) corr_s(s = s_axis,theta_s = theta_s)
    corr_t_mat <- function(theta_t) corr_t(t = t_axis,theta_t = theta_t)
    d_corr_s_mat <- function(theta_s) d_corr_s(s = s_axis,theta_s = theta_s)
    d_corr_t_mat <- function(theta_t) d_corr_t(t = t_axis,theta_t = theta_t)

    C_idx <- st_grid %>%
      as.data.frame() %>%
      select(s,t) %>%
      mutate(n = 1:nrow(st_grid)) %>%
      left_join(s_obs,.)

    C_m <- sparseMatrix(i=1:nrow(C_idx),
                        j = C_idx$n,
                        x=1,
                        dims=c(nrow(s_obs),nrow(st_grid)))

    chol_S_zeta_true <-  sigma_zeta_true * 
                              kronecker(chol(corr_t_mat(theta_t_true)),
                                        chol(corr_s_mat(theta_s_true)))

    s_obs <- s_obs %>%
        mutate(dis = as.vector(C_m %*% 
                                 (t(chol_S_zeta_true) %*% 
                                    rnorm(n = ns*nt))),
               z = z + dis)
}

## ----deprecated-sparse,include=FALSE----
# Not relavant to vignette
if (model=="sparse") {

  Q_s_mat <- function() Q_s(s = s_axis)
  Q_t_mat <- function(theta_t) Q_t(t = t_axis,theta_t = theta_t)
  d_Q_s_mat <- function() d_Q_s(s = s_axis)
  d_Q_t_mat <- function(theta_t) d_Q_t(t = t_axis,theta_t = theta_t)

  Q_zeta_true <- sigma_zeta_true^(-2) * Q_zeta_fn(s_axis,t_axis,theta_t_true)
  Q_zeta_true_sim <- sigma_zeta_true^(-2) * Q_zeta_fn(s_axis,t_axis,theta_t_true, add = TRUE)
  chol_Q_zeta_sim <- chol(Q_zeta_true_sim)

  C_idx <- st_grid %>%
       as.data.frame() %>%
       select(s,t) %>%
       mutate(n = 1:nrow(st_grid)) %>%
       left_join(s_obs,.)

  C_m <- sparseMatrix(i=1:nrow(C_idx),j = C_idx$n,x=1,dims=c(nrow(s_obs),nrow(st_grid)))

  s_obs <- s_obs %>%
    mutate(dis = as.numeric(C_m %*% backsolve(chol_Q_zeta_sim, rnorm(n = nrow(st_grid)))),
           z = z + dis)

}

## -------------------------------------
if(model == "full") {
    X <- group_by(st_grid,s) %>%
        summarise(Yf = Yf[1],Ym_av = mean(Ym)) %>%
        gather(process,value,-s)
    g <- LinePlotTheme() + 
          geom_line(data=subset(X,!(s==new_obs)),
                    aes(x=s,y=value,linetype=as.factor(process)))+
          geom_segment(data=s_obs[1:5,],
                       aes(x=s, xend=s, y = 0, yend = 50),
                        arrow=arrow(length=unit(0.1,"cm"))) +
        ylab("") +
        scale_linetype_discrete(guide=guide_legend(title="process"),
                                labels=c("Yf (g/s/degree)","Ym (ppb)")) +
        xlab("s (degrees)")
    ggsave(g,filename = "../../Sim_plot.png",width=10,height=4)
}

## ----YfYm-plot,fig.height=5,fig.width=10----
print(g)

## ----deprecated-old-plot,include=FALSE,eval=FALSE----
## ## Old plot, deprecated
## if(model == "full") {
##     source("../../R/ggplot_dual_axis.R")
##     g1 <- LinePlotTheme() + geom_line(data=subset(X,!(s==new_obs) & process=="Yf"),aes(x=s,y=value)) + ylab("Yf (g/s/degree)") + theme(text = element_text(size=15),axis.title.y=element_text(vjust=0.6)) + xlab("s (degrees)")
##     g2 <- LinePlotTheme() + geom_line(data=subset(X,!(s==new_obs) & process=="Ym_av"),aes(x=s,y=value),linetype="dashed") + ylim(0,900) + theme(text = element_text(size=15),axis.title.y=element_text(vjust=0.6)) +
##         geom_segment(data=s_obs[1:5,],aes(x=s, xend=s, y = 0, yend = 50),
##                      arrow=arrow(length=unit(0.1,"cm"))) + ylab("Ym (ppb)")
##     G <- ggplot_dual_axis(g1,g2)
## }

## ----deprecated-emulation,eval=FALSE,include=FALSE----
## ###----------------------------------
## ### Emulation (deprecated!)
## ###----------------------------------
## if(F) {
## 
##   # Assume we have 10 design points we want to emulate from, equally spaced in the domain
##   m_sim = 20
##   model_sim_obs <- expand.grid(u = s_axis,
##                                #s = s_axis[seq(5,95,by=5)]) %>%
##                                s = s_axis[round(seq(10,90,length=m_sim))]) %>%
##     data.frame() %>%
##     left_join(data.frame(st_grid)) %>%
##     mutate(s0 =s) %>%
##     select(-s) %>%
##     mutate(h =round(u - s0,1)) %>%
##     group_by(t,s0,u) %>%
##     summarise(b = b(s = s0, u = u,p = p), h=h)
##   ggplot(subset(model_sim_obs,b>0)) + geom_tile(aes(x=u,y=t,fill=as.factor(s0),alpha=pmax(b,0.7)))
## 
##   #b(s = s0[1], u = u[1],p = p[1]))
##   ###---------------------------------------
##   ### Atmospheric model prediction locations
##   ###---------------------------------------
##   model_pred <- expand.grid(up = s_axis,
##                             s0p = unique(s_obs$s),
##                             tp = t_axis) %>%
##     data.frame() %>%
##     mutate(hp = round(up - s0p,1))
## 
##   ###----------------------------------
##   ### IDW function (group by t and h)
##   ###----------------------------------
##   idw <- function(s0p,hp,model_obs_df) {
##     df <- filter(model_obs_df,h == hp)
##     B <- matrix(df$b,ncol = length(t_axis))
##     if(nrow(B) == 0) {
##       rep(0,length(t_axis))
##     } else {
##       D <- fields::rdist(df$s0[1:nrow(B)],s0p)
##       lambda <- as.numeric(1/(D)^2)
##       colSums(diag(lambda,ncol=nrow(B)) %*% B)/sum(lambda)
##     }
##   }
## 
## 
##   model_pred_em <- subset(model_pred) %>%
##     group_by(s0p,hp) %>%
##     mutate(zp = idw(s0p = s0p[1],
##                     hp = hp[1],
##                     model_obs_df = model_sim_obs)) %>%
##     mutate(zp = zp*(zp>0.5)) %>%
##     as.data.frame() %>%
##     mutate(t = tp, s = s0p, s_fine = round(s0p + hp,1))
## 
##   warning("Zeroing out very small values")
##   model_pred_em$zp <- model_pred_em$zp*(model_pred_em$zp > 0.5)
## 
##   B <- plyr::ddply(model_pred_em,c("tp","s0p"),function (X){
##     t(X$zp)
##   }) %>% select(-tp,-s0p) %>%
##     as.matrix() * ds
## 
## 
## 
##   ###----------------------------------
##   ### Plot Emulation Results
##   ###----------------------------------
##   ggplot(subset(model_pred_em,zp>0)) +
##     geom_tile(aes(x=round(s0p + hp,1),y=tp,fill=as.factor(s0p),alpha=pmax(zp,0.7))) +
##     scale_alpha_continuous(range = c(0,1))
## }

## ----eval=TRUE,include=TRUE,message=FALSE----
if(model %in% c("full","diag")) {
  df_for_B <- s_obs
} else {
  df_for_B <- st_grid
}

# TRUE B
B_true <- plyr::ddply(df_for_B,c("t","s"),function(df) {
  b = b(s=df$s[1],u=s_axis,p=df$p[1])}) %>%
  select(-s,-t) %>%
  as.matrix()*ds

B <- B_true
if(model == "sparse") {
  B <- as(B,"dgCMatrix")
}

B_true_df <- plyr::ddply(df_for_B,c("t","s"),function(df) {
  b = b(s=df$s[1],u=s_axis,p=df$p[1])}) %>%
  gather(s_grid,b,-t,-s) %>%
  separate(s_grid, into = c("V","s_fine"),sep="V") %>%
  select(-V) %>%
  mutate(s_fine = round(s_axis[as.numeric(s_fine)],1))# %>%
#left_join(model_pred_em,by = c("t","s","s_fine"))


if(model == "full") {
  ## Plot the source-receptor relationship at each observation

  ## The following code uses colours and shows the SRR on one plot  
  #   B_plot <- LinePlotTheme() + 
  #             geom_tile(data=subset(B_true_df ,b>0),
  #                       aes(x=s_fine,y=t,fill=as.factor(s),alpha=b)) +
  #     scale_alpha_continuous(range=c(0,1)) + xlab("u") +
  #     scale_fill_discrete(guide=guide_legend(title="s")) +
  #     coord_fixed(xlim=c(-10,10),ratio = 0.2)
  
  ## The following code shows one SRR per plot  
  B_plot <- LinePlotTheme() + 
    geom_tile(data=subset(B_true_df,b>0 & !(s==new_obs)),
              aes(x=s_fine,y=t,alpha=b),fill="black") +
    scale_alpha_continuous(guide=guide_legend(title="s/ng")) +
    scale_y_reverse()+
    scale_fill_discrete(guide=guide_legend(title="s")) +
    coord_fixed(xlim=c(-10,10),ratio = 0.5)  +
    facet_grid(~s) +
    theme(panel.margin = unit(1.5, "lines")) + 
    xlab("u (degrees)") + 
    ylab("t (2 h steps)")
  ggsave(filename = "../../B_plot.png",width=12)
}

## ----B-plot,fig.height=6,fig.width=10----
print(B_plot)

## ----deprecated-other-figure,include=FALSE,eval=FALSE----
## if(model == "full") {
##   BT1 <- data.frame(s = s_axis,t(B[141:145,])) %>% gather(obs,b,-s)
##   base_plot <- LinePlotTheme() + geom_line(data = subset(st_grid,t==1),aes(x=s,y=Yf)) +
##     geom_point(data=subset(s_obs,t==1),aes(x=s,y=0),colour="red",size=3) +
##     geom_ribbon(data=subset(BT1,b>0), aes(x = s,ymax = b*100,ymin=0,group=obs,fill=obs),alpha=0.5,stat="identity")
## }

## ----deprecated-Bayes-Linear,include=FALSE,eval=FALSE----
## ###----------------------------------------------------
## ### Bayes Linear
## ###----------------------------------------------------
## if(F) {
##   C_f <- Zeromat(nrow(s_obs),ns)
##   C_full <- cBind(C_f,C_m)
## 
##   corr_zeta <- crossprod(kronecker(chol_corr_zeta_t_true,chol_corr_zeta_s_true))
##   #corr_zeta <- diag(rep(1,nrow(s_obs)))
## 
## 
##   BLupdate <- function(theta) {
## 
##     S_zeta <- theta*corr_zeta
##     S_reg <- diag(rep(sigma_reg^2,nrow(B)))
##     S_m <- B %*% S_f %*% t(B) + S_zeta + S_reg ## Added S_reg for regularization
##     S <- cbind(rbind(S_f,S_mf),rbind(S_fm,S_m))
##     Cov_YZ <- S %*% t(C_full)
##     Cov_Z <- as.matrix(crossprod(chol(S) %*% t(C_full)) + chol2inv(chol(Qobs)))
##     mean <- mu + Cov_YZ %*% chol2inv(chol(Cov_Z)) %*% (s_obs$z - C_full %*% mu)
##     cov <- S - tcrossprod(Cov_YZ %*% solve(chol(Cov_Z)))
##     return(list(mean = mean,cov = cov))
##   }
## 
##   sigma2_zeta_update <- function(Ystat) {
##     Yf <- Ystat$mean[1:ns]
##     Ym <- Ystat$mean[-(1:ns)]
## 
##     cov <- as.matrix(Ystat$cov)
##     S_ff <- cov[1:ns,1:ns]
##     S_mm <- cov[-(1:ns),-(1:ns)]
##     S_fm <- cov[1:ns,-(1:ns)]
##     S_mf <- t(S_fm)
## 
##     Psi <-find_Psi(Ym,Yf,S_mm,S_ff,S_fm)
## 
##     sum(diag(solve(corr_zeta, Psi))) / nrow(B)
## 
##   }
## 
##   sigma2_zeta <- 1
##   for(i in 1:10) {
##     X <- BLupdate(sigma2_zeta)
##     sigma2_zeta <- sigma2_zeta_update(X)
##     print(sigma2_zeta)
##   }
## }

## ----mu-f-----------------------------
mu_f <- matrix(exp(mu_f_log + 0.5*diag(S_f_log)))

## ----EM-settings----------------------
###-----------------------------------------------------
### Laplace method -- use with caution because of mode close to zero
###-----------------------------------------------------
n_EM <- 100

if(model == "full") {
  s_mol = s_obs$s[1:m_obs]         # prediction locs for mol fraction
  Y_init = c(mu_f,s_obs$z)         # initial expectation value for Y
  theta_init = c(1000,0.2,0.2)     # initial parameter vector
                                   # where theta = [sigma_zeta^2, theta_t, theat_s]
  
  # Keep station out for validation
  rm_idx <- seq(6,nrow(C_m),by=6) # index to be removed
  s_obs_old <- s_obs              # save old observation location data frame
  s_obs <- s_obs[-rm_idx,]        # new observation data frame
  C_m <- C_m[-rm_idx,]            # new incidence matrix
  Qobs <- Qobs[-rm_idx,-rm_idx]   # new (observation) precision matrix

} else if (model == "full_big") {
    s_mol = s_axis                # prediction locs for mol fraction
    Y_init = c(mu_f,st_grid$Ym)   # initial expectation value for Y
    theta_init = c(1000,0.2,0.2)  # initial parameter vector
                                  # where theta = [sigma_zeta^2, theta_t, theat_s]
}

## ----EM-settings2,include=FALSE-------
if (model == "diag") {
  s_mol = s_obs$s[1:m_obs]
  Y_init = c(mu_f,s_obs$z)
  theta_init = c(10)
} else if (model == "sparse") {
  s_mol = s_axis
  Y_init = c(mu_f,st_grid$Ym)
  theta_init = c(0.1,0.2)
}

## ----EM-setup-------------------------
EM_alg <- EM(s_obs = s_obs,                  # observation data frame
        C_m = C_m,                           # incidence matrix
        Qobs = Qobs,                         # observation precision matrix
        B = B,                               # SRR matrix
        t_mol = t_axis,                      # time prediction locs
        s_mol = matrix(s_mol),               # spatial prediction locs
        S_f_log = S_f_log,                   # cov. matrix of log(Yf)
        mu_f_log = mu_f_log,                 # expectation of log(Yf)
        Yf_thresh = 1e-4,                    # drop nodes with Yf < 1e-4
        Y_init =  Y_init,                    # initialise Yf
        theta_init = theta_init,             # initialise theta
        ind = which(!(colSums(B) == 0)),     # only consider indices with SRR>0
        n_EM = n_EM,                         # number of EM iterations
        model = model)                       # model we are using

## ----EM,cache=TRUE,results="hide"-----
for(i in 1:(n_EM-2)) {
  X <- EM_alg(max_E_it = 1e6,
              max_M_it = 50,
              fine_tune_E = (i==0))
}

## ----cov-compute----------------------
if(model == "full") {
  corr_zeta <- corr_zeta_fn(c(s_obs$s[1:(m_obs-1)],new_obs),t_axis,X$theta[2,n_EM],X$theta[3,n_EM])
  S_zeta <- X$theta[1,n_EM] * corr_zeta
  Q_zeta <- chol2inv(chol(S_zeta))
}

## ----deprecated-cov-compute,include=FALSE,eval=FALSE----
## if (model == "diag") {
##   S_zeta <- X$theta[1,n_EM] * .symDiagonal(nrow(s_obs))
##   Q_zeta <- chol2inv(chol(S_zeta))
## } else if (model == "sparse") {
##   stop("Didn't sample for sparse yet")
## }

## -------------------------------------
lap2 <- Yf_marg_approx_fns(s = s_obs,            # obs. locations
                           C_m = C_m,            # incidence matrix
                           Qobs = Qobs,          # obs. precision matrix
                           B = B,                # SRR matrix
                           S_zeta = S_zeta,      # Sigma_zeta
                           mu_f_log = mu_f_log,  # expectation of log(Y_f)
                           S_f_log = S_f_log,    # covariance of log(Y_f)
                           ind=X$ind)            # only consider indices with SRR>0

## ----HMC-setup------------------------
M <- diag(1/(X$lap_approx$Yf^2)) # scaling matrix: puts variables on same scale

if(model == "full"){             # Function for generating Delta
  eps_gen <- function() 
    runif(n=1,
          min = 0.066, 
          max = 0.068)
} else if (model =="diag") {
  eps_gen <- function() 
    runif(n=1,
          min = 0.0065, 
          max = 0.0067)
}
L <- 10L                   # step-size
N <- 10000                 # number of samples
q <- matrix(0,nrow(M),N)   # matrix for storing samples
qsamp <- X$lap_approx$Yf   # first sample
dither <- i <- count <- 1  # no dithering, initialise counts

## ----HMC-init-------------------------
sampler <- hmc_sampler(U = lap2$logf,                 # log-likelihood
                       dUdq = lap2$gr_logf,           # gr. of log-likelihood 
                       M = M,                         # scaling matrix
                       eps_gen = eps_gen,             # step-size generator
                       L = L,                         # number of steps
                       lower = rep(0,length(X$ind)))  # lower limits

## ----HMC-run,results="hide",cache=TRUE----
while(i < N) {
  qsamp <- sampler(q = qsamp)
  
  if(count  == dither) {
    q[,i] <- qsamp
    count = 1
    i <- i + 1
    print(paste0("Sample: ",i," Acceptance rate: ",(nrow(unique(t(q)))-1)/i))
  } else {
    count <- count + 1
  }
}

## ----HMC-acceptance-------------------
print(paste0("Sample: ",i," Acceptance rate: ",(nrow(unique(t(q)))-1)/i))

## ----flux-plot,include=FALSE,eval=FALSE----
##  Q <- as.data.frame(t(q[1:length(X$ind),-c(1:1000,i)])) %>%
##   tidyr::gather(t,z) %>%
##   separate(t, into=c("V","s"),sep="V") %>%
##   select(-V) %>%
##   mutate(s = s_axis[as.numeric(s)])
## 
## Q2 <- Q %>%
##       mutate(s = as.factor(s))
## g <- LinePlotTheme() +
##   geom_boxplot(data=Q2,aes(x=s,y=z))  +
##   scale_x_discrete(breaks = s_axis[seq(1,length(s_axis),length=7)]) +
##   geom_point(data=subset(st_grid, s < 5),
##              aes(x=as.factor(s),y = Yf),
##              colour='black',size=3,shape=4) +
##   geom_segment(data=s_obs,aes(x=as.factor(s_obs$s),
##                               xend=as.factor(s_obs$s),
##                               y = 0, yend = 50),
##                  arrow=arrow(length=unit(0.1,"cm"))) +
##   geom_line(data=data.frame(s=rep(3.9,2),
##                             z=c(-200,1300)),
##             aes(as.factor(s),z),linetype="dashed") +
##   geom_line(data=data.frame(s=rep(-5.3,2),
##                             z=c(-200,1300)),
##             aes(as.factor(s),z),linetype="dashed") +
##   coord_fixed(ratio = 0.03, ylim=c(0,1050)) +
##   ylab("Yf (g/s/degree)") + xlab("s (degrees)")
## if(!misspecification)
##   ggsave("../../Sim1_samples.png",plot = g,width=10)

## ----HMC-plot,fig.height=6,fig.width=5----
print(g)

## ----include=FALSE,eval=FALSE---------
## comp_density <- function(j,xu,xl=0,yu=0.01) {
##   x <- seq(xl,xu,by=1)
##   hist(q[j,1:(N-1)],xlab=c("flux (g/s/degree)"),
##        ylab="[Yf | Zm]",
##        main="",
##        freq=F,
##        xlim=c(xl,xu),
##        ylim=c(0,yu))
##   lines(x,dnorm(x,mean=X$lap_approx$Yf[j],
##                 sd=sqrt(X$lap_approx$S_ff[j,j])),
##         lty=2)
## }
## 
## if(!misspecification) {
##   png("../../density_sim1.png", width=4, height=4, units="in", res=300)
##   comp_density(70,800,-200,yu=0.005); dev.off()
##   png("../../density_sim10.png", width=4, height=4, units="in", res=300)
##   comp_density(24,400,50,yu=0.014); dev.off()
## }

## ----fig.height=4,fig.width=7---------
par(mfrow=c(1,2))
comp_density(70,800,-200,yu=0.005);
comp_density(24,400,50,yu=0.014);

## ----mf-samps-------------------------
## Now get the mole fraction samples
mf_samp <- matrix(0,nrow(B),N)                # initialise mole-fraction samples array
mu_post <- matrix(0,nrow(B),N)                # initialise mole-fraction samples array
Q_post <- t(C_m) %*% Qobs %*% C_m + Q_zeta    # conditional precision of Ym
S_post <- chol2inv(chol(Q_post))              # conditional covariance of Ym
L_S_post <- t(chol(S_post))                   # Cholesky decomposition of above
SQB <- S_post %*% Q_zeta %*% B[,X$ind]        # weight given to flux sample
StCQoz <- S_post %*% t(C_m) %*% Qobs %*% s_obs$z # influence of observation
for (i in 1:N){
 mu_post[,i] <- as.vector(SQB %*% q[,i] +     # conditional mean
                            StCQoz)           
 mf_samp[,i] <- as.vector(mu_post[,i] +       # generate sample
                            L_S_post %*% rnorm(nrow(B)))
}

## ----mf-validation--------------------
stat1 <-  seq(6,600,by=6)
mu <- apply(mf_samp[stat1,-c(1:1000,i)],1,median)
uq <- apply(mf_samp[stat1,-c(1:1000,i)],1,quantile,0.75)
lq <- apply(mf_samp[stat1,-c(1:1000,i)],1,quantile,0.25)
uuq <- apply(mf_samp[stat1,-c(1:1000,i)],1,quantile,0.95)
llq <- apply(mf_samp[stat1,-c(1:1000,i)],1,quantile,0.05)
df <- data.frame(mean = mu, uq=uq,lq=lq,uuq=uuq,llq=llq,t=1:length(mu))

g <- LinePlotTheme() + 
  geom_ribbon(data=df,aes(x=t,ymax=uq,ymin=lq),alpha=0.6) +
  geom_ribbon(data=df,aes(x=t,ymax=uuq,ymin=llq),alpha=0.3) +
  geom_point(data=subset(s_obs_old,s==new_obs),
             aes(x=t,y = z),colour='black',size=3,shape=4)+
    ylab("Ym (ppb)") + xlab("t (2 h steps)")

if(!misspecification) 
  ggsave("../../MF_samples.png",plot = g,width=10,height=3)

## ----fig.height=4,fig.width=7---------
print(g)

## ----stats----------------------------
# Flux errors
residual <-(st_grid$Yf[1:75] - apply(q,1,mean))
post_unc <- apply(q,1,var)
print(paste0("S1f = ", sqrt(mean(residual^2))))
print(paste0("S2f = ", sqrt(mean(residual^2) /mean(post_unc))))

# MF errors
residual <- (apply(mf_samp[rm_idx,],1,mean) - subset(st_grid,s == new_obs)$Ym)
post_unc <- apply(mf_samp[rm_idx,],1,var)
print(paste0("S1m = ", sqrt(mean(residual^2))))

