## This problem ONLY works if b is of relatively local scope. Once we have b which has a very large scope we get oscillations/instability

library(devtools)
library(ggplot2)
library(dplyr)
library(actuar)
library(Matrix)
library(tidyr)
load_all("../../../../CurrentProjects/PostDoc Bristol/R Code/pkg/MVST")
load_all("../..//pkgs//hmc")
load_all("../../pkgs//atminv")

rm(list=ls())

###------------------
### Parameters
###------------------
model = "sparse" ## Either sparse or full
set.seed(25) # 25 , T = 400 or 15/200
ds <- 0.2 ## If we change this we need to change the round() command further down
smin = -10 + ds/2
smax =  10 - ds/2
s_axis <- round(seq(smin,smax,by=ds),1)
if(model=="full") {
  t_axis <- 1:100 #372
  m_obs <- 5
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
if(model=="full") {
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
### Transition kernel
###------------------

b <- function(s,u,p) {
  absp <- max(abs(p),0.2)
  absp*sqrt(2*pi) * dnorm(u,mean = s, sd =absp) * 
    ((sign(p) == sign(u-s)) | (u-s) == 0) *
    (abs(u - s) < absp)
}

## Sample the "wind" vector
Q_s <- GMRF_RW(n = length(s_axis),order = 2,precinc = 2000,name="space")@Q + 0.001*Imat(length(s_axis))
Q_t <- GMRF_RW(n = length(t_axis),order = 1,precinc = 20,name="time")@Q + 0.1*Imat(length(t_axis))
Q_full = as(kronecker(Q_t,Q_s),"dgCMatrix")
G <- GMRF(mu = matrix(rep(0,nrow(Q_full))),Q = Q_full)
st_grid$p <- sample_GMRF(G,reps = 1)
fields::image.plot(matrix(st_grid$p,ns,nt))

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
S_f_log <- Matern(r = as.matrix(dist(s_axis)),nu = 1.5,var = 1,kappa = kappa_from_l(3,nu = 1.5))
mu_f_log <- matrix(rep(5,length(s_axis)))
warning("Setting var = 1 and not var = 5")
Yf_sim <- exp(mu_f_log + t(chol(S_f_log)) %*% rnorm(n = length(s_axis)))
st_grid <- st_grid %>%
  left_join(data.frame(s=s_axis,Yf = Yf_sim))


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
if(model == "full") {
  s_obs <- data.frame(s = sample(s_axis[-c(1:10,(ns-10):ns)],size = m_obs,
                                 replace=F),m = 1:m_obs)   
} else {
  s_obs <- data.frame(s = sample(s_axis,size = m_obs, replace=T),m = 1:m_obs)
}
  #s_obs <- data.frame(s = round(c(-8.1,-4.1,0.1,4.1,8.1),1),m = 1:m_obs) %>%
  s_obs <- s_obs %>%
            left_join(st_grid) %>%
            mutate(z = Ym + rnorm(n = length(Ym), sd = sigma_eps_true)) %>%
            arrange(t,s)

## Now add the discrepancy (stupid to add it after but too large to do before)
if(model == "full") {

    corr_s_mat <- function(theta_s) corr_s(s = s_obs$s[1:m_obs],theta_s = theta_s)
    corr_t_mat <- function(theta_t) corr_t(t = t_axis,theta_t = theta_t)
    d_corr_s_mat <- function(theta_s) d_corr_s(s = s_obs$s[1:m_obs],theta_s = theta_s)
    d_corr_t_mat <- function(theta_t) d_corr_t(t = t_axis,theta_t = theta_t)
    
    corr_zeta_true <- corr_zeta_fn(s_obs$s[1:m_obs],t_axis,theta_t_true, theta_s_true)
    S_zeta_true <- sigma_zeta_true^2 * corr_zeta_true
    chol_S_zeta_true <- chol(S_zeta_true)
    
    C_m <- Imat(nrow(s_obs))
    
    s_obs <- s_obs %>%
      mutate(dis = t(chol_S_zeta_true) %*% rnorm(n = nrow(s_obs)),
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
Qobs <- sigma_eps_true^(-2) * Imat(nrow(s_obs))

X <- group_by(st_grid,s) %>% 
  summarise(Yf = Yf[1],Ym_av = mean(Ym)) %>%
  gather(process,value,-s)
g <- LinePlotTheme() + geom_line(data=X,aes(x=s,y=value,colour=as.factor(process)))+ 
  geom_point(data=s_obs,aes(x=s_obs$s,y=0),colour="red",size=3) +
  ylab("") +
  scale_colour_discrete(guide=guide_legend(title="process"))
if(model == "full") ggsave(g,filename = "../../MVCF/Chemometrics/Sim_plot.png")

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

if(model == "full") {
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

B_true_df <- plyr::ddply(df_for_B,c("t","s"),function(df) { 
  b = b(s=df$s[1],u=s_axis,p=df$p[1])}) %>%
  gather(s_grid,b,-t,-s) %>%
  separate(s_grid, into = c("V","s_fine"),sep="V") %>%
  select(-V) %>%
  mutate(s_fine = round(s_axis[as.numeric(s_fine)],1))# %>%
#left_join(model_pred_em,by = c("t","s","s_fine"))

if(model == "full") { 
  ## Plot the source-receptor relationship at each observation
  B_plot <- LinePlotTheme() + geom_tile(data=subset(B_true_df ,b>0),
                                        aes(x=s_fine,y=t,fill=as.factor(s),alpha=b)) +
    scale_alpha_continuous(range=c(0,1)) + xlab("u") +
    scale_fill_discrete(guide=guide_legend(title="s")) + 
    coord_fixed(xlim=c(-10,10),ratio = 0.2)
  ggsave(filename = "../../MVCF/Chemometrics/B_plot.png")
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
  s_mol = s_obs$s[1:m_obs]
  Y_init = c(mu_f,s_obs$z)
  theta_init = c(0.1,0.2,0.2)
} else {
  s_mol = s_axis
  Y_init = c(mu_f,st_grid$Ym)
  theta_init = c(0.1,0.2)
}
n_EM <- 100
EM_alg <- EM(s_obs = s_obs,
        C_m = C_m,
        Qobs = Qobs,
        B = B,
        t_mol = t_axis,
        s_mol = s_mol,
        S_f_log = S_f_log,
        mu_f_log = mu_f_log,
        Yf_thresh = 1e-4,
        Y_init =  Y_init,
        theta_init = theta_init,
        ind = which(!(colSums(B) == 0)),
        model = model)

for(i in 1:n_EM) {
  X <- EM_alg(max_E_it = 1e6,
              max_M_it = 50,
              fine_tune_E = (i==1))
}


### Now use HMC sampler with fixed paramaters
### FULL STATE
if(0) {
  M <- diag(rep(1, length(mode)))
  eps_gen <- function() runif(n=1,min = 0.05, max = 0.15)
  L <- 10L
  sampler <- hmc_sampler(U = lap$logf,dUdq = lap$gr_logf, M = chol2inv(chol(gr2)),eps_gen = eps_gen,L = L)
  N <- 1000
  dither <- 10
  q <- matrix(0,nrow(M),N)
  qsamp <- mode
  i <- count <- 1
  
  
  while(i < N) {
    qsamp <- sampler(q = qsamp)
    
    if(count  == dither) {
      q[,i] <- qsamp
      count = 1
      i <- i + 1
      print(i)
    } else {
      count <- count + 1
    }
  }
  
  Q <- as.data.frame(t(q[1:100,1:(i-1)])) %>%
    tidyr::gather(t,z) %>%
    separate(t, into=c("V","s"),sep="V") %>%
    select(-V) %>% 
    mutate(s = s_axis[as.numeric(s)])
  
  ggplot(Q) + geom_boxplot(aes(x=as.factor(s),y=z)) + 
    geom_point(data=st_grid,aes(x=factor(s),y = Yf),colour='red',size=3)
  
}

### MARGINAL STATE
corr_zeta <- corr_zeta_fn(s_obs$s[1:m_obs],t_axis,X$theta[2,n_EM],X$theta[3,n_EM]) 
S_zeta <- X$theta[1,n_EM] * corr_zeta

lap2 <- Yf_marg_approx_fns(s = s_obs, 
                           C_m = C_m,
                           Qobs = Qobs,
                           B = B, 
                           S_zeta = S_zeta, #sigma_zeta_true^2 * diag(rep(1,nrow(s_obs)))
                           mu_f_log = mu_f_log,
                           S_f_log = S_f_log,
                           ind=X$ind)
  
if(0) {
  # M <- diag(rep(0.3, length(X$ind)))
  # M <- chol2inv(chol(gr2[1:length(ind),1:length(ind)])) ## WARNING: Laplace is horrible where there are no observations, observations too tight
  #M <- solve(diag(diag(gr2[1:length(ind),1:length(ind)])))
  M <- 1/diag(mode[1:length(X$ind)]^2)
  eps_gen <- function() runif(n=1,min = 0.0001, max = 0.0002)
  #eps_gen <- function() runif(n=1,min = 0.001, max = 0.002)
  L <- 100L
  sampler <- hmc_sampler(U = lap2$logf,dUdq = lap2$gr_logf, 
                         M = M,
                         eps_gen = eps_gen,
                         L = L,
                         lower = rep(0,length(X$ind)))
  N <- 1000
  q <- matrix(0,nrow(M),N)
  qsamp <- mode[1:length(X$ind)]
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
  
  Q <- as.data.frame(t(q[1:length(ind),1:(i-1)])) %>%
    tidyr::gather(t,z) %>%
    separate(t, into=c("V","s"),sep="V") %>%
    select(-V) %>% 
    mutate(s = s_axis[ind[as.numeric(s)]])
  
  ggplot(Q) + geom_boxplot(aes(x=as.factor(s),y=z)) + 
    geom_point(data=subset(st_grid,s %in%s_axis[ind]),aes(x=as.factor(s),y = Yf),colour='red',size=3) 
}

slice_sampler <- slice::oaatslice(d = length(X$ind), scale=X$lap_approx$Yf)
#  slice_sampler <- slice::slice(100,scale = 50)
# slice_sampler <- slice::blockslice(list(ind, setdiff(1:ns,ind)),scale = c(4,4))
  
log_f <- function(x) -lap2$logf(x)
  
N <- 10000
q <- matrix(0,length(X$ind),N)
qsamp <- pmax(X$lap_approx$Yf,1)
dither <- 1
i <- count <- 1

while(i < N) {
  qsamp <- slice_sampler(qsamp,log_f,learn = i <= 200)
  
  
  if(count  == dither) {
    q[,i] <- qsamp
    count = 1
    i <- i + 1
    print(i)
  } else { 
    count <- count + 1
  }
}
  
Q <- as.data.frame(t(q[1:length(X$ind),1:(i-1)])) %>%
  tidyr::gather(t,z) %>%
  separate(t, into=c("V","s"),sep="V") %>%
  select(-V) %>% 
  mutate(s = s_axis[X$ind[as.numeric(s)]])

Q2 <- rbind(Q,data.frame(s = s_axis[setdiff(1:ns,X$ind)], z = 0))
g <- LinePlotTheme() + geom_boxplot(data=Q2,aes(x=s,y=z))  +
  scale_x_discrete(breaks = s_axis[seq(1,length(s_axis),length=10)]) +
  geom_point(data=subset(st_grid),aes(x=as.factor(s),y = Yf),colour='red',size=3) +
  geom_point(data=subset(s_obs,t==1),aes(x=as.factor(s),y = 0),colour='blue',size=3) +
  coord_fixed(ratio = 0.05) + ylab("Yf")
ggsave("../../MVCF/Chemometrics/Sim1_samples.png",plot = g,width=10)


  geom_point(data=data.frame(s=s_axis[X$ind],y=X$lap_approx$Yf),aes(x=factor(s),y = y),colour='green',size=3)   

  geom_point(data=data.frame(s=s_axis[X$ind],y=X$lap_approx$Yf,sd=sqrt(diag(X$lap_approx$S_ff))),aes(x=factor(s),y = y + sd),colour='blue',size=3) + 
  geom_point(data=data.frame(s=s_axis[X$ind],y=mode[1:length(X$ind)],sd=sqrt(diag(gr2)[1:length(X$ind)])),aes(x=factor(s),y = y-sd),colour='blue',size=3)
 