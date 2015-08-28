## This problem ONLY works if b is of relatively local scope. Once we have b which has a very large scope we get oscillations/instability

library(devtools)
library(ggplot2)
library(grid)
library(dplyr)
library(Matrix)
library(tidyr)
#library(scales) # for format_format
#load_all("../../../../../CurrentProjects/PostDoc Bristol/R Code/pkg/MVST")
load_all("../../../pkgs//hmc")
load_all("../atminv")

rm(list=ls())

###------------------
### Parameters
###------------------
model = "full" ## Either sparse or full or full_big or diag
misspecification = 0
set.seed(25) # 25 , T = 400 or 15/200
ds <- 0.2 ## If we change this we need to change the round() command further down
smin = -10 + ds/2
smax =  10 - ds/2
s_axis <- round(seq(smin,smax,by=ds),1)
if(model %in% c("full","diag")) {
  t_axis <- 1:100 #372
  m_obs <- 6
} else {
  t_axis <- 1:100
  m_obs <- 1000
}
ns <- length(s_axis)
nt <- length(t_axis)
st_grid <- expand.grid(s=s_axis,t=t_axis) %>% data.frame()

# Observation parameters
sigma_eps_true <- 10 # 0.1


# 'full' model parameters
if(model %in% c("full","full_big")) {
  sigma_zeta_true <- 50 #1
} else {
  sigma_zeta_true <- 10
}

theta_t_true <- 0.8
theta_s_true <- 1

# 'sparse' model parameters
theta_t_true <- 0.8
theta_s_true <- 1

###------------------
### Transition kernel: p is a ST process and reflects the std of the truncated Gaussian.
###------------------
b <- function(s,u,p) {
  absp <- max(abs(p),0.2)
  absp*sqrt(2*pi) * dnorm(u,mean = s, sd =absp) *
    ((sign(p) == sign(u-s)) | (u-s) == 0) *
    (abs(u - s) < absp)
}

## Sample the "wind" vector
Q_s <- GMRF_RW(n = length(s_axis),order = 2,precinc = 2000)@Q + 0.001*.symDiagonal(length(s_axis))
Q_t <- GMRF_RW(n = length(t_axis),order = 1,precinc = 20)@Q + 0.1*.symDiagonal(length(t_axis))
Q_full = as(kronecker(Q_t,Q_s),"dgCMatrix")
G <- GMRF(mu = matrix(rep(0,nrow(Q_full))),Q = Q_full,n=nrow(Q_full))

load("~/Desktop/Chemometrics_results/sim.Random.seed.rda")
st_grid$p <- sample_GMRF(G,reps = 1)
#fields::image.plot(matrix(st_grid$p,ns,nt))

## Simulate the atmosphere
s_sim <- 9.1
Z <- matrix(0,nt,ns)
for( i in 1 : nrow(Z)) {
  Z[i,] <- b(s = s_sim, u = s_axis,p = as.numeric(subset(st_grid,s== s_sim & t==i)$p))
}
###------------------
### Lognormal flux field
###------------------
set.seed(1) # 1 for sparse obs
#S_f_log <- Matern(r = as.matrix(dist(s_axis)),nu = 1.5,var = 1,kappa = kappa_from_l(3,nu = 1.5))
variogram_model <- gstat::vgm(range = 3.334,nugget = 0.00533,psill = 0.80429,model= "Sph")
S_f_log <- gstat::variogramLine(variogram_model,
                           dist_vector = sp::spDists(matrix(s_axis),
                                                 matrix(s_axis)),
                           covariance = TRUE)
mu_f_log <- matrix(rep(5,length(s_axis)))
Yf_sim <- exp(mu_f_log + t(chol(S_f_log)) %*% rnorm(n = length(s_axis)))
st_grid <- st_grid %>%
left_join(data.frame(s=s_axis,Yf = Yf_sim))

if(misspecification)  {
    S_f_log <- diag(diag(S_f_log))
}

###------------------
### Mole fraction field
###------------------
### Since we cannot put the discrepancy everywhere (too large), we will just put it at the observation locations
st_grid <- st_grid %>%
  group_by(t,s,p) %>%
  summarise(Yf = Yf,
            Ym = sum(b(s =s, u = s_axis,p = p) * Yf_sim * ds))

###------------------
### Observations
###------------------
if(model %in% c("full","diag")) {
      s_obs <- data.frame(s = sample(s_axis[-c(1:10,(ns-10):ns)],size = m_obs,
                                     replace=F),m = 1:m_obs  )
    new_obs <- 0.3
    s_obs[6,]$s <- new_obs
#     s_obs <- data.frame(s = sample(c(-7.1, -3.1, 0.1, 2.5, 7.1, 8.5),size = m_obs,
#                                    replace=F),m = 1:m_obs)
} else {
  s_obs <- data.frame(s = sample(s_axis,size = m_obs, replace=T),m = 1:m_obs)
}
  #s_obs <- data.frame(s = round(c(-8.1,-4.1,0.1,4.1,8.1),1),m = 1:m_obs) %>%
  s_obs <- s_obs %>%
            left_join(st_grid) %>%
            mutate(z = Ym + rnorm(n = length(Ym), sd = sigma_eps_true)) %>%
            arrange(t,s)

## Now add the discrepancy (stupid to add it after but too large to do before)
if(model %in% c("full","diag")) {

#     corr_s_mat <- function(theta_s) corr_s(s = s_obs$s[1:m_obs],theta_s = theta_s)
#     corr_t_mat <- function(theta_t) corr_t(t = t_axis,theta_t = theta_t)
#     d_corr_s_mat <- function(theta_s) d_corr_s(s = s_obs$s[1:m_obs],theta_s = theta_s)
#     d_corr_t_mat <- function(theta_t) d_corr_t(t = t_axis,theta_t = theta_t)

    corr_zeta_true <- corr_zeta_fn(s_obs$s[1:m_obs],t_axis,theta_t_true, theta_s_true)
    S_zeta_true <- sigma_zeta_true^2 * corr_zeta_true
    chol_S_zeta_true <- chol(S_zeta_true)

    C_m <- .symDiagonal(nrow(s_obs))

    s_obs <- s_obs %>%
      mutate(dis = t(chol_S_zeta_true) %*% rnorm(n = nrow(s_obs)),
             z = z + dis)

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

    C_m <- sparseMatrix(i=1:nrow(C_idx),j = C_idx$n,x=1,dims=c(nrow(s_obs),nrow(st_grid)))



    chol_S_zeta_true <-  sigma_zeta_true * kronecker(chol(corr_t_mat(theta_t_true)),
                                                     chol(corr_s_mat(theta_s_true)))

    s_obs <- s_obs %>%
        mutate(dis = as.vector(C_m %*% (t(chol_S_zeta_true) %*% rnorm(n = ns*nt))),
               z = z + dis)

} else if (model=="sparse") {

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

C_f <- Zeromat(nrow(s_obs),ns)
C_full <- cBind(C_f,C_m)
Qobs <- sigma_eps_true^(-2) * .symDiagonal(nrow(s_obs))

if(model == "full") {
    X <- group_by(st_grid,s) %>%
        summarise(Yf = Yf[1],Ym_av = mean(Ym)) %>%
        gather(process,value,-s)
    g <- LinePlotTheme() + geom_line(data=subset(X,!(s==new_obs)),aes(x=s,y=value,linetype=as.factor(process)))+
        geom_segment(data=s_obs[1:5,],aes(x=s, xend=s, y = 0, yend = 50),
                     arrow=arrow(length=unit(0.1,"cm"))) +
        ylab("") +
        scale_linetype_discrete(guide=guide_legend(title="process"),labels=c("Yf (g/s/degree)","Ym (ppb)")) +
        xlab("s (degrees)")
    ggsave(g,filename = "../Sim_plot.png",width=10,height=4)

    source("ggplot_dual_axis.R")
    g1 <- LinePlotTheme() + geom_line(data=subset(X,!(s==new_obs) & process=="Yf"),aes(x=s,y=value)) + ylab("Yf (g/s/degree)") + theme(text = element_text(size=15),axis.title.y=element_text(vjust=0.6)) + xlab("s (degrees)")
    g2 <- LinePlotTheme() + geom_line(data=subset(X,!(s==new_obs) & process=="Ym_av"),aes(x=s,y=value),linetype="dashed") + ylim(0,900) + theme(text = element_text(size=15),axis.title.y=element_text(vjust=0.6)) +
        geom_segment(data=s_obs[1:5,],aes(x=s, xend=s, y = 0, yend = 50),
                     arrow=arrow(length=unit(0.1,"cm"))) + ylab("Ym (ppb)")
    G <- ggplot_dual_axis(g1,g2)
}


###----------------------------------
### Emulation
###----------------------------------
if(F) {

  # Assume we have 10 design points we want to emulate from, equally spaced in the domain
  m_sim = 20
  model_sim_obs <- expand.grid(u = s_axis,
                               #s = s_axis[seq(5,95,by=5)]) %>%
                               s = s_axis[round(seq(10,90,length=m_sim))]) %>%
    data.frame() %>%
    left_join(data.frame(st_grid)) %>%
    mutate(s0 =s) %>%
    select(-s) %>%
    mutate(h =round(u - s0,1)) %>%
    group_by(t,s0,u) %>%
    summarise(b = b(s = s0, u = u,p = p), h=h)
  ggplot(subset(model_sim_obs,b>0)) + geom_tile(aes(x=u,y=t,fill=as.factor(s0),alpha=pmax(b,0.7)))

  #b(s = s0[1], u = u[1],p = p[1]))
  ###---------------------------------------
  ### Atmospheric model prediction locations
  ###---------------------------------------
  model_pred <- expand.grid(up = s_axis,
                            s0p = unique(s_obs$s),
                            tp = t_axis) %>%
    data.frame() %>%
    mutate(hp = round(up - s0p,1))

  ###----------------------------------
  ### IDW function (group by t and h)
  ###----------------------------------
  idw <- function(s0p,hp,model_obs_df) {
    df <- filter(model_obs_df,h == hp)
    B <- matrix(df$b,ncol = length(t_axis))
    if(nrow(B) == 0) {
      rep(0,length(t_axis))
    } else {
      D <- fields::rdist(df$s0[1:nrow(B)],s0p)
      lambda <- as.numeric(1/(D)^2)
      colSums(diag(lambda,ncol=nrow(B)) %*% B)/sum(lambda)
    }
  }


  model_pred_em <- subset(model_pred) %>%
    group_by(s0p,hp) %>%
    mutate(zp = idw(s0p = s0p[1],
                    hp = hp[1],
                    model_obs_df = model_sim_obs)) %>%
    mutate(zp = zp*(zp>0.5)) %>%
    as.data.frame() %>%
    mutate(t = tp, s = s0p, s_fine = round(s0p + hp,1))

  warning("Zeroing out very small values")
  model_pred_em$zp <- model_pred_em$zp*(model_pred_em$zp > 0.5)

  B <- plyr::ddply(model_pred_em,c("tp","s0p"),function (X){
    t(X$zp)
  }) %>% select(-tp,-s0p) %>%
    as.matrix() * ds



  ###----------------------------------
  ### Plot Emulation Results
  ###----------------------------------
  ggplot(subset(model_pred_em,zp>0)) +
    geom_tile(aes(x=round(s0p + hp,1),y=tp,fill=as.factor(s0p),alpha=pmax(zp,0.7))) +
    scale_alpha_continuous(range = c(0,1))
}

###----------------------------------------------------
### Joint covariance matrix -- for now ignore emulator
###----------------------------------------------------

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
#   B_plot <- LinePlotTheme() + geom_tile(data=subset(B_true_df ,b>0),
#                                         aes(x=s_fine,y=t,fill=as.factor(s),alpha=b)) +
#     scale_alpha_continuous(range=c(0,1)) + xlab("u") +
#     scale_fill_discrete(guide=guide_legend(title="s")) +
#     coord_fixed(xlim=c(-10,10),ratio = 0.2)
  B_plot <- LinePlotTheme() + geom_tile(data=subset(B_true_df ,b>0 & !(s==new_obs)),
                                        aes(x=s_fine,y=t,alpha=b),fill="black") +
    scale_alpha_continuous(guide=guide_legend(title="s/ng")) +
      scale_y_reverse()+
    scale_fill_discrete(guide=guide_legend(title="s")) +
    coord_fixed(xlim=c(-10,10),ratio = 0.5)  +
        facet_grid(~s) +
    theme(panel.margin = unit(1.5, "lines")) + xlab("u (degrees)") + ylab("t (2 h steps)")
  ggsave(filename = "../B_plot.png",width=12)
}

mu_f <- matrix(exp(0.5*diag(S_f_log)))
mu_m <- as.matrix(B %*% mu_f)
mu <- rbind(mu_f,mu_m)


S_f <- exp(S_f_log[1,1]) *(exp(S_f_log) - 1)
S_fm <- S_f %*% t(B)
S_mf <- B %*% S_f

if(model == "full") {
  BT1 <- data.frame(s = s_axis,t(B[141:145,])) %>% gather(obs,b,-s)
  base_plot <- LinePlotTheme() + geom_line(data = subset(st_grid,t==1),aes(x=s,y=Yf)) +
    geom_point(data=subset(s_obs,t==1),aes(x=s,y=0),colour="red",size=3) +
    geom_ribbon(data=subset(BT1,b>0), aes(x = s,ymax = b*100,ymin=0,group=obs,fill=obs),alpha=0.5,stat="identity")
}


###----------------------------------------------------
### Bayes Linear
###----------------------------------------------------
if(F) {
  corr_zeta <- crossprod(kronecker(chol_corr_zeta_t_true,chol_corr_zeta_s_true))
  #corr_zeta <- diag(rep(1,nrow(s_obs)))


  BLupdate <- function(theta) {

    S_zeta <- theta*corr_zeta
    S_reg <- diag(rep(sigma_reg^2,nrow(B)))
    S_m <- B %*% S_f %*% t(B) + S_zeta + S_reg ## Added S_reg for regularization
    S <- cbind(rbind(S_f,S_mf),rbind(S_fm,S_m))
    Cov_YZ <- S %*% t(C_full)
    Cov_Z <- as.matrix(crossprod(chol(S) %*% t(C_full)) + chol2inv(chol(Qobs)))
    mean <- mu + Cov_YZ %*% chol2inv(chol(Cov_Z)) %*% (s_obs$z - C_full %*% mu)
    cov <- S - tcrossprod(Cov_YZ %*% solve(chol(Cov_Z)))
    return(list(mean = mean,cov = cov))
  }

  sigma2_zeta_update <- function(Ystat) {
    Yf <- Ystat$mean[1:ns]
    Ym <- Ystat$mean[-(1:ns)]

    cov <- as.matrix(Ystat$cov)
    S_ff <- cov[1:ns,1:ns]
    S_mm <- cov[-(1:ns),-(1:ns)]
    S_fm <- cov[1:ns,-(1:ns)]
    S_mf <- t(S_fm)

    Psi <-find_Psi(Ym,Yf,S_mm,S_ff,S_fm)

    sum(diag(solve(corr_zeta, Psi))) / nrow(B)

  }

  sigma2_zeta <- 1
  for(i in 1:10) {
    X <- BLupdate(sigma2_zeta)
    sigma2_zeta <- sigma2_zeta_update(X)
    print(sigma2_zeta)
  }
}

###-----------------------------------------------------
### Laplace method -- use with caution because of mode close to zero
###-----------------------------------------------------
if(model == "full") {
  n_EM <- 100
  s_mol = s_obs$s[1:m_obs]
#   Zmat <- matrix(s_obs$z,nrow=5)
#   Zmat <- rbind(Zmat,Zmat[5,])
#   Y_init = c(mu_f,c(Zmat))
  Y_init = c(mu_f,s_obs$z)
  theta_init = c(1000,0.2,0.2)

  # Assume a station has failed
  rm_idx <- seq(6,nrow(C_m),by=6)
  s_obs_old <- s_obs
  s_obs <- s_obs[-rm_idx,]
  C_m <- C_m[-rm_idx,]
  Qobs <- Qobs[-rm_idx,-rm_idx]

} else if (model == "diag") {
  n_EM <- 100
  s_mol = s_obs$s[1:m_obs]
  Y_init = c(mu_f,s_obs$z)
  theta_init = c(10)
} else if (model == "sparse") {
  n_EM <- 100
  s_mol = s_axis
  Y_init = c(mu_f,st_grid$Ym)
  theta_init = c(0.1,0.2)
} else if (model == "full_big") {
    n_EM <- 100
    s_mol = s_axis
    Y_init = c(mu_f,st_grid$Ym)
    theta_init = c(1000,0.2,0.2)
}


EM_alg <- EM(s_obs = s_obs,
        C_m = C_m,
        Qobs = Qobs,
        B = B,
        t_mol = t_axis,
        s_mol = matrix(s_mol),
        S_f_log = S_f_log,
        mu_f_log = mu_f_log,
        Yf_thresh = 1e-4,
        Y_init =  Y_init,
        theta_init = theta_init,
        ind = which(!(colSums(B) == 0)),
        n_EM = n_EM,
        model = model)

breakme
readline("Press enter to run EM algorithm")
for(i in 1:(n_EM-2)) {
  X <- EM_alg(max_E_it = 1e6,
              max_M_it = 50,
              fine_tune_E = (i==0))
}
#theta = [sigma_zeta^2, theta_t, theat_s]

### Now use HMC sampler with fixed paramaters


### MARGINAL STATE
if(model == "full") {
  corr_zeta <- corr_zeta_fn(c(s_obs$s[1:(m_obs-1)],new_obs),t_axis,X$theta[2,n_EM],X$theta[3,n_EM])
  S_zeta <- X$theta[1,n_EM] * corr_zeta
  Q_zeta <- chol2inv(chol(S_zeta))

} else if (model == "diag") {
  S_zeta <- X$theta[1,n_EM] * .symDiagonal(nrow(s_obs))

} else if (model == "sparse") {
  stop("Didn't sample for sparse yet")
}

lap2 <- Yf_marg_approx_fns(s = s_obs,
                           C_m = C_m,
                           Qobs = Qobs,
                           B = B,
                           S_zeta = S_zeta, #sigma_zeta_true^2 * diag(rep(1,nrow(s_obs)))
                           mu_f_log = mu_f_log,
                           S_f_log = S_f_log,
                           ind=X$ind)


###-----------------------------------------------------
### HMC Sampler
###-----------------------------------------------------

readline("Press enter to run HMC sampler")
if(1) {
  # M <- diag(rep(0.3, length(X$ind)))
  # M <- chol2inv(chol(gr2[1:length(ind),1:length(ind)])) ## WARNING: Laplace is horrible where there are no observations, observations too tight
  #M <- solve(diag(diag(gr2[1:length(ind),1:length(ind)])))
  M <- diag(1/(X$lap_approx$Yf^2))
  if(model == "full"){
    #eps_gen <- function() runif(n=1,min = 0.035, max = 0.037)
    eps_gen <- function() runif(n=1,min = 0.066, max = 0.068)
    L <- 10L
  } else if (model =="diag") {
    eps_gen <- function() runif(n=1,min = 0.0065, max = 0.0067)
    L <- 10L
  }
  #eps_gen <- function() runif(n=1,min = 0.001, max = 0.002)

  ### GOOD SETTINGS L = 100, eps = 0.03--0.032
  sampler <- hmc_sampler(U = lap2$logf,dUdq = lap2$gr_logf,
                         M = M,
                         eps_gen = eps_gen,
                         L = L,
                         lower = rep(0,length(X$ind)))
  N <- 10000
  q <- matrix(0,nrow(M),N)
  qsamp <- X$lap_approx$Yf
  #   qsamp <- rep(1,ns)
  #   qsamp[ind] <- mode[1:length(ind)]
  dither <- 1
  i <- count <- 1

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


}

 Q <- as.data.frame(t(q[1:length(X$ind),-c(1:1000,i)])) %>%
# Q <- as.data.frame(t(q[1:length(X$ind),1:(i-1)])) %>%
# Q <- as.data.frame(t(q[,-c(1:1000,i)])) %>%
  tidyr::gather(t,z) %>%
  separate(t, into=c("V","s"),sep="V") %>%
  select(-V) %>%
  mutate(s = s_axis[as.numeric(s)])

Q2 <- Q %>% #rbind(Q,data.frame(s = s_axis[setdiff(1:ns,X$ind)], z = 0)) %>%
      mutate(s = as.factor(s))
g <- LinePlotTheme() + geom_boxplot(data=Q2,aes(x=s,y=z))  +
  #scale_x_discrete(breaks = s_axis[seq(1,length(s_axis),length=10)]) +
  scale_x_discrete(breaks = s_axis[seq(1,length(s_axis),length=7)]) +
  #geom_point(data=subset(st_grid, s < 4.4),aes(x=as.factor(s),y = Yf),colour='red',size=3) +
  geom_point(data=subset(st_grid, s < 5),aes(x=as.factor(s),y = Yf),colour='black',size=3,shape=4) +
  geom_segment(data=s_obs,aes(x=as.factor(s_obs$s), xend=as.factor(s_obs$s), y = 0, yend = 50),
                 arrow=arrow(length=unit(0.1,"cm"))) +
  geom_line(data=data.frame(s=rep(3.9,2),z=c(-200,1300)),aes(as.factor(s),z),linetype="dashed") +
  geom_line(data=data.frame(s=rep(-5.3,2),z=c(-200,1300)),aes(as.factor(s),z),linetype="dashed") +
#  geom_point(data=subset(s_obs,t==1),aes(x=as.factor(s),y = 0),colour='blue',size=3) +
  coord_fixed(ratio = 0.03, ylim=c(0,1050)) + ylab("Yf (g/s/degree)") + xlab("s (degrees)")
g
ggsave("../Sim1_samples.png",plot = g,width=10)


comp_density <- function(j,xu,xl=0,yu=0.01) {
  x <- seq(xl,xu,by=1)
  hist(q[j,1:(N-1)],xlab=c("flux (g/s/degree)"),ylab="[Yf | Zm]",main="",freq=F,xlim=c(xl,xu),ylim=c(0,yu))
  lines(x,dnorm(x,mean=X$lap_approx$Yf[j],sd=sqrt(X$lap_approx$S_ff[j,j])),lty=2)
}

png("../density_sim1.png", width=4, height=4, units="in", res=300)
comp_density(70,800,-200,yu=0.005); dev.off()
png("../density_sim10.png", width=4, height=4, units="in", res=300)
comp_density(24,400,50,yu=0.014); dev.off()


## Now get the mole fraction samples
mf_samp <- matrix(0,nrow(B),N)
mu_post <- matrix(0,nrow(B),N)
Q_post <- t(C_m) %*% Qobs %*% C_m + Q_zeta
S_post <- chol2inv(chol(Q_post))
L_S_post <- t(chol(S_post))
SQB <- S_post %*% Q_zeta %*% B[,X$ind]
StCQoz <- S_post %*% t(C_m) %*% Qobs %*% s_obs$z
for (i in 1:N){
 mu_post[,i] <- as.vector(SQB %*% q[,i] + StCQoz)
 mf_samp[,i] <- as.vector(mu_post[,i] + L_S_post %*% rnorm(nrow(B)))
}


stat1 <-  seq(6,600,by=6)
mu <- apply(mf_samp[stat1,-c(1:1000,i)],1,median)
uq <- apply(mf_samp[stat1,-c(1:1000,i)],1,quantile,0.75)
lq <- apply(mf_samp[stat1,-c(1:1000,i)],1,quantile,0.25)
uuq <- apply(mf_samp[stat1,-c(1:1000,i)],1,quantile,0.95)
llq <- apply(mf_samp[stat1,-c(1:1000,i)],1,quantile,0.05)
df <- data.frame(mean = mu, uq=uq,lq=lq,uuq=uuq,llq=llq,t=1:length(mu))

g <- LinePlotTheme() + geom_ribbon(data=df,aes(x=t,ymax=uq,ymin=lq),alpha=0.6) +
  geom_ribbon(data=df,aes(x=t,ymax=uuq,ymin=llq),alpha=0.3) +
  geom_point(data=subset(s_obs_old,s==new_obs),aes(x=t,y = z),colour='black',size=3,shape=4)+
    ylab("Ym (ppb)") + xlab("t (2 h steps)")
ggsave("../MF_samples.png",plot = g,width=10,height=3)

# MF <- as.data.frame(t(mf_samp[rm_idx,-c(1:1000,i)])) %>%
#     # Q <- as.data.frame(t(q[1:length(X$ind),1:(i-1)])) %>%
#     # Q <- as.data.frame(t(q[,-c(1:1000,i)])) %>%
#     tidyr::gather(s,z) %>%
#     separate(s, into=c("V","t"),sep="V") %>%
#     select(-V) %>%
#     mutate(t = t_axis[as.numeric(t)])
#
# MF2 <- MF %>% #rbind(Q,data.frame(s = s_axis[setdiff(1:ns,X$ind)], z = 0)) %>%
#     mutate(t = as.factor(t))
# g <- LinePlotTheme() + geom_boxplot(data=MF2,aes(x=t,y=z))  +
#     scale_x_discrete(breaks = t_axis[seq(1,length(t_axis),length=7)]) +
#     #geom_point(data=subset(st_grid, s < 4.4),aes(x=as.factor(s),y = Yf),colour='red',size=3) +
#     geom_point(data=subset(st_grid, s == new_obs),aes(x=as.factor(t),y = Ym),colour='black',size=3,shape=4)+
#     coord_fixed(ratio = 0.025, ylim=c(-200,1600)) + ylab("Yf")
# g
# ggsave("../../MVCF/Chemometrics/MF_samples.png",plot = g,width=10)

# Flux errors
residual <-(st_grid$Yf[1:75] - apply(q,1,mean))
post_unc <- apply(q,1,var)
print(paste0("S1f = ", sqrt(mean(residual^2))))
print(paste0("S2f = ", sqrt(mean(residual^2) /mean(post_unc))))
print(paste0("C1 = ", mean(residual) / (1/length(post_unc) * sqrt(sum(post_unc)))))

# MF errors
residual <- (apply(mf_samp[rm_idx,],1,mean) - subset(st_grid,s == new_obs)$Ym)
post_unc <- apply(mf_samp[rm_idx,],1,var)
print(paste0("S1m = ", sqrt(mean(residual^2))))
print(paste0("S2m = ", sqrt(mean(residual^2) /mean(post_unc))))
