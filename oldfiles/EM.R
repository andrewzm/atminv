## Laplace functions for the ATCP problem

### EM method needs input
# s_obs -- observation df
# C_m -- C matrix relating Z to MF
# Qobs -- observation precision
# B -- atmospheric model matrix
# corr_t -- temporal correlation function with one free parameter
# corr_s -- spatial correlation function with one free parameter
# S_f_log -- covariance function of log(Yf)
# mu_f_log -- exepectation of log(Yf)
# Yf_thresh -- threshold at which to start dropping Yf nodes
# Y_init -- initial Y values
# theta_init -- initial theta values
# n_EM -- number of EM iterations

EM <- function(s_obs,
               C_m,
               Qobs,
               B,
               S_f_log,
               mu_f_log,
               Yf_thresh,
               Y_init,
               theta_init,
               ind = 1:ncol(B),
               s_mol = NULL,
               t_mol = NULL,
               n_EM = 100,
               model="full") {
  
    
    
  if(model=="full") {
    browser()
    corr_s_mat <- function(theta_s) corr_s(s = s_mol,theta_s = theta_s)
    corr_t_mat <- function(theta_t) corr_t(t = t_mol,theta_t = theta_t)
    d_corr_s_mat <- function(theta_s) d_corr_s(s = s_mol,theta_s = theta_s)
    d_corr_t_mat <- function(theta_t) d_corr_t(t = t_mol,theta_t = theta_t)
    
    S_zeta_t <- corr_t_mat(theta_init[2])
    S_zeta_s <- corr_s_mat(theta_init[3])
    Q_zeta_t <- chol2inv(chol(S_zeta_t))
    Q_zeta_s <- chol2inv(chol(S_zeta_s))

    corr_zeta <- kronecker(S_zeta_t,S_zeta_s)
    S_zeta <- theta_init[1] * corr_zeta
    Q_zeta <- 1/theta_init[1] * kronecker(Q_zeta_t,Q_zeta_s)
    Q_zeta <- chol2inv(chol(S_zeta))
    
  } else if (model == "diag") {
    corr_zeta <- kronecker(Imat(nrow(s_mol)),Imat(length(t_mol)))
    S_zeta <- theta_init[1] * corr_zeta
    Q_zeta <- 1/theta_init[1] * corr_zeta
  } else if(model == "sparse") {
    Q_s_mat <- function() Q_s(s = s_mol)
    Q_t_mat <- function(theta_t) Q_t(t = t_mol,theta_t = theta_t)
    d_Q_s_mat <- function() d_Q_s(s = s_mol)
    d_Q_t_mat <- function(theta_t) d_Q_t(t = t_mol,theta_t = theta_t)
    
    Q_zeta <- theta_init[1]^(-1) * kronecker(Q_t_mat(theta_init[2]),Q_s_mat())
  }
  
  remove <- setdiff(1:ncol(B),ind)
  mode <- Y_init
  if(length(remove) > 0)  mode <- mode[-remove]
  theta <- matrix(theta_init,length(theta_init),n_EM)
  i <- 2
  
  function(max_E_it = 1e6,
           max_M_it = 10,
           fine_tune_E = F) {
    
    
    OK <- FALSE
    while(!OK) {
#       lap <- laplace_approx_fns(s_obs = s_obs, 
#                                 C_m = C_m, 
#                                 Qobs = Qobs,
#                                 B = B,
#                                 Q_zeta = Q_zeta,
#                                 mu_f_log = mu_f_log,
#                                 S_f_log = S_f_log,
#                                 ind = ind)
        
        lap <- laplace_approx_fns(s_obs = s_obs, 
                                  C_m = C_m, 
                                  Qobs = Qobs,
                                  B = B,
                                  prec_zeta = 1/theta_init[1],
                                  Q_zeta_s = Q_zeta_s,
                                  Q_zeta_t = Q_zeta_t,
                                  mu_f_log = mu_f_log,
                                  S_f_log = S_f_log,
                                  ind = ind)
      mode <<- optim(mode,
                    fn = lap$logf,
                    gr = lap$gr_logf,
                    method = "BFGS",
                    control=list(trace=4,reltol=1e-8,maxit=max_E_it))$par
      
      ## Fine tune at peak if needed
      if(fine_tune_E) mode <<- optim(mode,fn = lap$logf,
                             gr = lap$gr_logf,
                             method = "BFGS",
                             control=list(trace=4,reltol=1e-12,maxit=1e6,parscale = abs(mode)))$par
      
      if (any(mode[1:length(ind)] < Yf_thresh)) {
        remove <- which(mode[1:length(ind)] < Yf_thresh)
        ind <<- ind[-remove]
        mode <<- mode[-remove]
      
      } else {
        OK <- TRUE
      }
    }  
    
    
    if(model %in% c("diag","sparse")) {
      inf <- lap$gr2_logf(mode) 
      XX <- cholPermute(inf)
      gr2 <- Takahashi_Davis(Q = inf,cholQp = XX$Qpermchol, P = XX$P)
      rm(XX)
      rm(inf)
    } else {
      gr2 <- chol2inv(chol(as.matrix(lap$gr2_logf(mode)))) 
    }
    lap_approx <- summ_stats(mode,gr2,ind=ind)
    rm(gr2)
    #sigma2_zeta[i+1] <- sigma2_zeta_update_only(Estep = lap_approx,
    #                                            corr_zeta = corr_zeta_true)
    #print(sigma2_zeta[i+1])
    if(model=="full") {
      
      M <- M_step_fns(Estep = lap_approx,
                      corr_t=corr_t_mat,
                      corr_s=corr_s_mat,
                      d_corr_t = d_corr_t_mat, 
                      d_corr_s = d_corr_s_mat, 
                      B=B, 
                      ind = ind)
    } else if (model == "sparse") {
      M <- M_step_fns_sparse(Estep = lap_approx,
                             Q_t=Q_t_mat,
                             Q_s=Q_s_mat,
                             d_Q_t = d_Q_t_mat, 
                             d_Q_s = d_Q_s_mat, 
                             B=B, 
                             ind = ind)
    }
    
    if(!(model == "diag")) {
      theta[,(i+1)] <<- optim(theta[,i],f = M$logf,gr = M$gr_logf,method = "BFGS",
                             control=list(trace=4,reltol=1e-8,maxit=max_M_it))$par
    } else {
      theta[,(i+1)] <<- sigma2_zeta_update_only(Estep = lap_approx, B=B,corr_zeta = corr_zeta,ind=ind)
    }
    
    if(model=="full") {
      corr_zeta <<- kronecker(corr_t_mat(theta[2,i+1]),corr_s_mat(theta[3,i+1]))
      S_zeta <<- theta[1,i+1] * corr_zeta
      Q_zeta <<- chol2inv(chol(S_zeta))
    } else if (model == "sparse") {
      Q_zeta <<-  theta[1,i+1]^(-1) * kronecker(Q_t_mat(theta[2,i+1]),Q_s_mat())
    } else if (model == "diag") {
      S_zeta <<- theta[1,i+1] * corr_zeta
      Q_zeta <<- 1/theta[1,i+1] * corr_zeta
      
    }
    i <<- i + 1
    
    return(list(lap_approx=lap_approx,
                theta=theta,
                ind=ind))
  }
 
}

## Returns the negative function, negative gradient and negative Hessian on current parameters
## Uses Kronecker form only (not that important since we have a problem with the Hessian)
laplace_approx_fns <- function(s_obs,
                               C_m,
                               Qobs,
                               B,
                               Q_zeta_t,
                               Q_zeta_s,
                               prec_zeta,
                               mu_f_log,
                               S_f_log,
                               ind)
{

  mu_f_log <- matrix(mu_f_log[ind])
  B <- B[,ind]
  Q_f_log <- chol2inv(chol(S_f_log[ind,ind]))
  
  cholQo <- chol(Qobs)
  cholQf_log <- chol(Q_f_log)
  
  Q_zeta_x <- function(X) {
      prec_zeta * c(Q_zeta_s %*% matrix(as.vector(X),ns,nt) %*% Q_zeta_t)
  }
  

  ns <- nrow(Q_zeta_s)
  nt <- nrow(Q_zeta_t)
  tC_Qo <- t(C_m) %*% Qobs
  tC_Qo_C <- tC_Qo %*% C_m

  nf <- length(ind)

  browser()
  logf <- function(x) {
    
    x <- as.matrix(x)
    Ym = x[-(1:nf)]
    Yf = x[1:nf]
    

    if(any(Yf <= 0)) {
      return(Inf)
    } else {
      -(
        -0.5 * crossprod(cholQo %*% (s_obs$z - C_m %*% Ym)) -
          0.5 * t(Ym - B %*% Yf) %*% Q_zeta_x(Ym - B %*% Yf)-
          #0.5 * t(Ym - B %*% Yf) %*% Q_m %*% (Ym - B %*% Yf) -
          0.5 * crossprod(cholQf_log %*% (log(Yf) - mu_f_log)) - 
          sum(log(Yf))
      ) %>% 
        as.numeric() %>%
        return()
    }
  }
  
  gr_logf <- function(x) {
    x <- as.matrix(x)
    Ym = x[-(1:nf)]
    Yf = x[1:nf]
    
    grYf <- as.matrix(t(B) %*% Q_zeta_x(Ym) - 
      t(B) %*% Q_zeta_x(B %*% Yf)) - 
      diag(1/Yf) %*% Q_f_log %*% log(Yf) + 
      diag(1/Yf) %*% Q_f_log %*% mu_f_log -
      1/Yf 
    
    grYm <- as.matrix(-tC_Qo_C %*% Ym + tC_Qo %*% s_obs$z -
                          Q_zeta_x(Ym) + Q_zeta_x(B %*% Yf))
  
    - rbind(grYf, grYm) %>% as.numeric()
  }

  
  gr2_logf <- function(x) {
    Ym = x[-(1:nf)]
    Yf = x[1:nf]
        
    Qm <- prec_zeta*kronecker(Q_zeta_t,Q_zeta_s)
    
    grYmYm <- -tC_Qo_C - prec_zeta * krone
    grYmYf <- Qm %*% B
    grYfYm <- t(B) %*% Qm
    grYfYf <- -t(B) %*% Qm %*% B  -
      diag(as.numeric(Q_f_log %*% mu_f_log)) %*% diag(1/Yf^2) -
      diag(1/Yf) %*% Q_f_log %*% diag(1/Yf) + 
      diag(as.numeric(Q_f_log %*% log(Yf))) %*% diag(1/Yf^2) +
      diag(1/Yf^2)
    
    -rBind(cBind(grYfYf, grYfYm),cBind(grYmYf,grYmYm))
  }
  
  list(logf = logf,
       gr_logf = gr_logf,
       gr2_logf = gr2_logf)
  
}  

# ## Returns the negative function, negative gradient and negative Hessian on current parameters
# laplace_approx_fns <- function(s_obs,
#                                C_m,
#                                Qobs,
#                                B,
#                                Q_zeta,
#                                mu_f_log,
#                                S_f_log,
#                                ind)
# {
#   mu_f_log <- matrix(mu_f_log[ind])
#   B <- B[,ind]
#   Q_m <- Q_zeta
#   Q_f_log <- chol2inv(chol(S_f_log[ind,ind]))
#   
#   cholQo <- chol(Qobs)
#   cholQf_log <- chol(Q_f_log)
#   #cholQm <- chol(Q_m)
#   
#   tB_Qm <- t(B) %*% Q_m
#   tB_Qm_B <- tB_Qm %*% B
#   tC_Qo <- t(C_m) %*% Qobs
#   tC_Qo_C <- tC_Qo %*% C_m
#   
#   Qm_B <- Q_m %*% B
#   
#   nf <- length(ind)
#   
#   logf <- function(x) {
#     
#     x <- as.matrix(x)
#     Ym = x[-(1:nf)]
#     Yf = x[1:nf]
#     
#     if(any(Yf <= 0)) {
#       return(Inf)
#     } else {
#       -(
#         -0.5 * crossprod(cholQo %*% (s_obs$z - C_m %*% Ym)) -
#           #0.5 * crossprod(cholQm %*% (Ym - B %*% Yf)) - 
#           0.5 * t(Ym - B %*% Yf) %*% Q_m %*% (Ym - B %*% Yf) -
#           0.5 * crossprod(cholQf_log %*% (log(Yf) - mu_f_log)) - 
#           sum(log(Yf))
#       ) %>% 
#         as.numeric() %>%
#         return()
#     }
#   }
#   
#   gr_logf <- function(x) {
#     x <- as.matrix(x)
#     Ym = x[-(1:nf)]
#     Yf = x[1:nf]
#     
#     grYf <- as.matrix(tB_Qm %*% Ym - 
#       tB_Qm_B %*% Yf) - 
#       diag(1/Yf) %*% Q_f_log %*% log(Yf) + 
#       diag(1/Yf) %*% Q_f_log %*% mu_f_log -
#       1/Yf 
#     
#     grYm <- as.matrix(-tC_Qo_C %*% Ym + tC_Qo %*% s_obs$z -
#                         Q_m %*% Ym + Qm_B %*% Yf)
#   
#     - rbind(grYf, grYm) %>% as.numeric()
#   }
# 
#   
#   gr2_logf <- function(x) {
#     Ym = x[-(1:nf)]
#     Yf = x[1:nf]
#         
#     grYmYm <- -tC_Qo_C - Q_m
#     grYmYf <- Qm_B
#     grYfYm <- tB_Qm
#     grYfYf <- -tB_Qm_B -
#       diag(as.numeric(Q_f_log %*% mu_f_log)) %*% diag(1/Yf^2) -
#       diag(1/Yf) %*% Q_f_log %*% diag(1/Yf) + 
#       diag(as.numeric(Q_f_log %*% log(Yf))) %*% diag(1/Yf^2) +
#       diag(1/Yf^2)
#     
#     -rBind(cBind(grYfYf, grYfYm),cBind(grYmYf,grYmYm))
#   }
#   
#   list(logf = logf,
#        gr_logf = gr_logf,
#        gr2_logf = gr2_logf)
#   
# }  



summ_stats <- function(mode,cov,ind=ns) {
  nf <-length(ind)
  list(Yf = mode[1:nf],
       Ym = mode[-(1:nf)],
       S_ff = cov[1:nf,1:nf],
       S_mm = cov[-(1:nf),-(1:nf)],
       S_fm = cov[1:nf,-(1:nf)],
       S_mf = cov[-(1:nf),1:nf])
}

tr <- function(A) { sum(diag(A))}

sigma2_zeta_update_only <- function(Estep,B,corr_zeta,ind) {
  B <- B[,ind]
  Psi <-find_Psi(Estep$Ym,Estep$Yf,Estep$S_mm,Estep$S_ff,Estep$S_fm,B = B)
  tr(solve(corr_zeta, Psi)) / nrow(corr_zeta)
}

## For the EM algorithm
find_Psi <- function(Ym,Yf,S_mm,S_ff,S_fm, B) {
  S_mm + outer(Ym,Ym) + B %*% (S_ff + outer(Yf,Yf)) %*% t(B) - 
    2 * B %*% (S_fm + outer(Yf,Ym))
}

# find_Psi2 <- function(Ym,Yf,S_mm,S_ff,S_fm, B) {
#   S_mm_T <- as(S_mm,"dgTMatrix")
#   Ym_outer <- S_mm
#   Ym_outer@x <- Ym[S_mm_T@i+1] * Ym[S_mm_T@j+1]
#   
#   X <- S_mm + Ym_outer + B %*% (S_ff + outer(Yf,Yf)) %*% t(B) - 
#     2 * B %*% (S_fm + outer(Yf,Ym))
# }

M_step_fns <- function(Estep,corr_t,corr_s,d_corr_t,d_corr_s, B, ind=ind) {
  ### s_axis has to be s_obs$s[1:m_obs]
  Psi <-find_Psi(Estep$Ym,Estep$Yf,Estep$S_mm,Estep$S_ff,Estep$S_fm,B = B[,ind])
  ns <- nrow(corr_s(0.2))
  nt <- nrow(corr_t(0.2))
  n_tot <- ns*nt
  
  f_YinvX <- function(cholY,X) {
    X %>%
      forwardsolve(t(cholY),.) %>%
      backsolve(cholY,.) 
  }
  
  
  logf <- function(theta) {
    
    sigma2_zeta <- theta[1]
    prec_zeta <- 1/theta[1]
    theta_t <- theta[2]
    theta_s <- theta[3]
    
    
    if( (abs(theta_t) > 1) | (theta_s <= 0) | (sigma2_zeta <= 0)) {
      return(Inf)
    } else {
      
      corr_zeta_t <- corr_t(theta_t = theta_t)
      corr_zeta_s <- corr_s(theta_s = theta_s)
      chol_corr_zeta_t <- chol(corr_zeta_t)
      chol_corr_zeta_s <- chol(corr_zeta_s)
      
      chol_corr_zeta <- kronecker(chol_corr_zeta_t,chol_corr_zeta_s)
      
      logdet_corr_zeta <- ns*logdet(chol_corr_zeta_t) + nt*logdet(chol_corr_zeta_s)
      corr_zeta_inv_Psi <- f_YinvX(chol_corr_zeta,Psi)
      
      -(
        -0.5 * n_tot * log(sigma2_zeta) -
          0.5 * logdet_corr_zeta -
          0.5 * prec_zeta * tr(corr_zeta_inv_Psi)
      )
    }
  }
  
  gr_logf <- function(theta) {
    
    sigma2_zeta <- theta[1]
    prec_zeta <- 1/theta[1]
    theta_t <- theta[2]
    theta_s <- theta[3]
    
    corr_zeta_t <- corr_t(theta_t = theta_t)
    corr_zeta_s <- corr_s(theta_s = theta_s)
    d_corr_zeta_t <- d_corr_t(theta_t = theta_t)
    d_corr_zeta_s <- d_corr_s(theta_s = theta_s)
    chol_corr_zeta_t <- chol(corr_zeta_t)
    chol_corr_zeta_s <- chol(corr_zeta_s)
    
    chol_corr_zeta <- kronecker(chol_corr_zeta_t,chol_corr_zeta_s)
    
    logdet_corr_zeta <- ns*logdet(chol_corr_zeta_t) + nt*logdet(chol_corr_zeta_s)
    
    corr_zeta_inv_Psi <- f_YinvX(chol_corr_zeta,Psi)
    #corr_zeta_inv_dt_s <- f_YinvX(chol_corr_zeta,kronecker(d_corr_zeta_t, corr_zeta_s))
    #corr_zeta_inv_t_ds <- f_YinvX(chol_corr_zeta,kronecker(corr_zeta_t, d_corr_zeta_s))
    corr_zeta_inv_dt_s <- kronecker(f_YinvX(chol_corr_zeta_t,d_corr_zeta_t),diag(nrow(corr_zeta_s)))
    corr_zeta_inv_t_ds <- kronecker(diag(nrow(corr_zeta_t)),f_YinvX(chol_corr_zeta_s,d_corr_zeta_s))
    
    
    gr_theta1 <- -0.5 * n_tot * prec_zeta +
      0.5 * prec_zeta^2 * tr(corr_zeta_inv_Psi)
    
    gr_theta2 <-  -0.5 * tr(corr_zeta_inv_dt_s) +
      0.5 * prec_zeta * tr(corr_zeta_inv_dt_s %*% corr_zeta_inv_Psi)
    
    gr_theta3 <-  -0.5 * tr(corr_zeta_inv_t_ds) +
      0.5 * prec_zeta * tr(corr_zeta_inv_t_ds %*% corr_zeta_inv_Psi)
    
    -(
      c(gr_theta1,gr_theta2,gr_theta3)
    )
  }
  list(logf = logf,
       gr_logf = gr_logf)
}



M_step_fns_sparse <- function(Estep,Q_t,Q_s,d_Q_t,d_Q_s, B, ind=ind) {
  ### s_axis has to be s_obs$s[1:m_obs]
  Psi <-find_Psi(Estep$Ym,Estep$Yf,Estep$S_mm,Estep$S_ff,Estep$S_fm,B = B[,ind])
  
  S_mm_symb <- Estep$S_mm
  S_mm_symb@x <- 1
  warning("Sparsing Psi using partial matrix inversion -- still need to prove I can do this")
  Psi <- S_mm_symb*Psi
  
  ns <- nrow(Q_s())
  nt <- nrow(Q_t(0.2))
  n_tot <- ns*nt
  
  logf <- function(theta) {
    
    sigma2_zeta <- theta[1]
    prec_zeta <- 1/theta[1]
    theta_t <- theta[2]
    
    if( (abs(theta_t) > 1) | (sigma2_zeta <= 0)) {
      return(Inf)
    } else {
      
      Q_zeta_t <- Q_t(theta_t = theta_t)
      Q_zeta_s <- Q_s()
      Q_zeta_ts <- kronecker(Q_zeta_t,Q_zeta_s)
      chol_Q_zeta_t <- chol(Q_zeta_t)
      
      logdet_Q_zeta <- ns*logdet(chol_Q_zeta_t) + 0
      Q_zeta_Psi <- Q_zeta_ts %*% Psi
      
      -(
        0.5 * n_tot * log(prec_zeta) +
          0.5 * logdet_Q_zeta -
          0.5 * prec_zeta * tr(Q_zeta_Psi)
      )
    }
  }
  
  gr_logf <- function(theta) {
    
    sigma2_zeta <- theta[1]
    prec_zeta <- 1/theta[1]
    theta_t <- theta[2]
    
    Q_zeta_t <- Q_t(theta_t = theta_t)
    Q_zeta_s <- Q_s()
    Q_zeta_ts <- kronecker(Q_zeta_t,Q_zeta_s)
    
    d_Q_zeta_t <- d_Q_t(theta_t = theta_t)
    d_Q_zeta_s <- d_Q_s()
    
    chol_Q_zeta_t <- chol(Q_zeta_t)
    
    logdet_Q_zeta <- ns*logdet(chol_Q_zeta_t) + 0
    
    Q_zeta_Psi <- Q_zeta_ts %*% Psi
    Q_zeta_dt_s <- kronecker(d_Q_zeta_t,Q_zeta_s)
    
    gr_theta1 <- -0.5 * n_tot * prec_zeta +
      0.5 * prec_zeta^2 * tr(Q_zeta_Psi)
    
    gr_theta2 <-  0.5 * tr(kronecker(solve(Q_zeta_t,d_Q_zeta_t),Imat(nrow(Q_zeta_s)))) -
      0.5 * prec_zeta * tr(Q_zeta_dt_s %*% Psi)
    
    -(
      c(gr_theta1,gr_theta2)
    )
  }
  list(logf = logf,
       gr_logf = gr_logf)
}

Yf_marg_approx_fns <- function(s_obs,
                               C_m,
                               Qobs,
                               B,
                               S_zeta,
                               mu_f_log,
                               S_f_log,
                               ind = 1:ns)
{
  B <- B[,ind]
  mu_f_log <- matrix(mu_f_log[ind])
  #S_obs_zeta <- solve(Qobs) + S_zeta
  S_obs_zeta <- solve(Qobs) + C_m %*% S_zeta %*% t(C_m)
  Q_obs_zeta <- chol2inv(chol(S_obs_zeta))
  cholQo_zeta <- chol(Q_obs_zeta)
  Q_f_log <- chol2inv(chol(S_f_log[ind,ind]))
  cholQf_log <- chol(Q_f_log)
  tB_tC_Qozeta_z <- t(B) %*% t(C_m) %*% Q_obs_zeta %*% s_obs$z
  tB_tC_Qozeta_C_B <- t(B) %*% t(C_m) %*% Q_obs_zeta %*% C_m %*% B
  
  logf <- function(x) {
    Yf <- as.matrix(x)
    
    if(any(Yf <= 0)) {
      return(Inf)
    } else {
      -(
        -0.5 * crossprod(cholQo_zeta %*% (s_obs$z - C_m %*% B %*% Yf)) -
          0.5 * crossprod(cholQf_log %*% (log(Yf) - mu_f_log)) - 
          sum(log(Yf))
      ) %>% 
        as.numeric() %>%
        return()
    }
  }
  
  gr_logf <- function(x) {
    Yf <- as.numeric(x)
    
    grYf <- tB_tC_Qozeta_z - 
        tB_tC_Qozeta_C_B %*% Yf - 
      diag(1/Yf) %*% Q_f_log %*% log(Yf) +
      diag(1/Yf) %*% Q_f_log %*% mu_f_log - 
      1/Yf
    
    -as.numeric(grYf[,1])
  }
  
  list(logf = logf,
       gr_logf = gr_logf)
}

Q_AR1 <- function(n=10,a=0.1) {
  i <- c(1,1,rep(2:(n-1),each=3),n,n)
  j <- c(1,2,c(outer(1:3,0:(n-3),FUN = "+")),(n-1),n)
  x <- c(1, -a, rep(c(-a,1+a^2,-a),(n-2)),-a,1)
  Q <- sparseMatrix(i=i,j=j,x=x)
}

dQda_AR1 <- function(n=10,a=0.1) {
  i <- c(1,1,rep(2:(n-1),each=3),n,n)
  j <- c(1,2,c(outer(1:3,0:(n-3),FUN = "+")),(n-1),n)
  x <- c(0, -1, rep(c(-1,2*a,-1),(n-2)),-1,0)
  Q <- sparseMatrix(i=i,j=j,x=x)
}

corr_t <- function(t,theta_t) {theta_t^as.matrix(dist(matrix(t)))}
d_corr_t <- function(t,theta_t) {as.matrix(dist(matrix(t))) * corr_t(t,theta_t) * theta_t^(-1) }
corr_s <- function(s,theta_s) {exp(-as.matrix(dist(s)/theta_s))}
d_corr_s <- function(s,theta_s) {theta_s^(-2)*as.matrix(dist(s)) * corr_s(s,theta_s)}

corr_zeta_fn <- function(s,t,theta_t,theta_s) {
  corr_zeta_t <- corr_t(t,theta_t)
  corr_zeta_s <- corr_s(s,theta_s)
  chol_corr_zeta_t <- chol(corr_zeta_t)
  chol_corr_zeta_s <- chol(corr_zeta_s)
  corr_zeta <- crossprod(kronecker(chol_corr_zeta_t,chol_corr_zeta_s))
}

Q_t <- function(t,theta_t) {Q_AR1(n=length(t),a=theta_t)}
d_Q_t <- function(t,theta_t) {dQda_AR1(n=length(t),a=theta_t)}

Q_s <- function(s,add=F) {Q_AR1(n=length(s), a=1) + add}
d_Q_s <- function(s,add=F) {Q_AR1(n=length(s), a=1) + add}

Q_zeta_fn <- function(s,t,theta_t,add=F) {
  Q_zeta_t <- Q_t(t,theta_t)
  Q_zeta_s <- Q_s(s,add=add)
  Q_zeta <- kronecker(Q_zeta_t,Q_zeta_s)
}