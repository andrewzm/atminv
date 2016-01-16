#' @title Box Cox transformations and related
#' @description These functions are used to transform from the original space to the the transformed space. In general, if \eqn{Z} are the non-Gaussian observations, then \eqn{g(Z)} is assumed to be Gaussian
#' @param Z the observations (these have to be positive)
#' @param lambda the Box-Cox parameter
#' @details The function \code{g} is the transformation, \code{g_dash} is the derivative of \code{g}, \code{g_2dash} is the second derivative of \code{g}, and \code{ginv} is the inverse of \code{g}. \code{J_lambda} is the Jacobian and is a diagonal matrix with the derivative of \code{g} evaluated at the observations as its elements.
#' @examples 
#' Z <- 1:10 
#' g(Z,0.4)
#' g_dash(Z,0.4)
#' g_2dash(Z,0.4)
#' J_lambda(Z,0.4)
#' ginv(g(Z,0.4),0.4)
#'@export
g <- function(Z,lambda) {
    if(lambda == 0) {
        log(Z)
    } else {
        (Z^lambda -1)/lambda
    }
}

#' @rdname g
#' @export
g_dash <- function(Z,lambda) {
    if(lambda == 0) {
        1/Z
    } else {
        Z^(lambda-1)
    }
}

#' @rdname g
#' @export
g_2dash <- function(Z,lambda) {
    if(lambda == 0) {
        -1/(Z^2)
    } else {
        (lambda-1)*Z^(lambda-2)
    }
}

#' @rdname g
#' @export
ginv <- function(Y,lambda) {
    if(lambda == 0) {
        exp(Y)
    } else {
        (lambda * Y + 1) ^(1/lambda)
    }
}

#' @rdname g
#' @export
J_lambda <- function(Z,lambda) {
    Diagonal(x=as.numeric(g_dash(Z,lambda)))
}


#' @title Conditional for Yf
#' @description Returns the conditional distribution and its derivative of the flux field given the flux field parameters and the mole-fraction observations
#' @param Z matrix of mole fraction observations (of size m x 1)
#' @param C_m observation incidence matrix (of size m x n_m)
#' @param Qobs observation precision matrix (of size m x m)
#' @param B source-receptor relationship matrix (of size n_m x n_f)
#' @param X covariates matrix (of size n_f x p)
#' @param Q_zeta precision matrix of discrepancy field (of size n_m x n_m). Should be sparse 
#' @param S_f_trans covariance matrix of flux field in transformed (Gaussian) space (of size n_f x n_f)
#' @param lambda Box-Cox transformation parameter
#' @param ind vector of indices on which to consider the flux (some unidentifiable grid cells may be omitted)
#' @export 
#' @examples 
#' n_f <- 20
#' n_m <- m_obs <- 100
#' cond_fns <- Yf_BC_conditionals(Z = 200 + matrix(rnorm(n_m,sd=100),n_m,1),
#'                                C_m = .symDiagonal(n_m),
#'                                Qobs = Diagonal(n_m),
#'                                B = matrix(rnorm(n_m*n_f),n_m,n_f),
#'                                X = matrix(rep(1,n_f)),
#'                                Q_zeta = n_m*Diagonal(n_m),
#'                                S_f_trans = Diagonal(n_f),
#'                                lambda = 0.1
#'                              )
#' Yf <- matrix(rpois(n_f,20),n_f,1)                              
#' cond_fns$logf(Yf)
#' cond_fns$gr_logf(Yf)
Yf_BC_conditionals <- function(Z,
                               Zinv,
                               C_m,
                               Qobs,
                               B,
                               X,
                               Q_zeta,
                               S_f_trans,
                               lambda=0,
                               lambda_fix=NA,
                               Pericchi=TRUE,
                               ind = 1:nrow(S_f_trans))               {
                                   
                                   # CtQoC_Q_zeta <- crossprod(chol(Qobs) %*% C_m) + Q_zeta
                                   # chol_CtQoC_Q_zeta <- chol(CtQoC_Q_zeta)
                                   # RmT_Q_zeta <- backsolve(chol_CtQoC_Q_zeta,
                                   # Q_zeta %*% B,
                                   # transpose = TRUE)  # X = solve(t(R)) %*% Q_zeta     
                                   # Mat12 <- t(B) %*% Q_zeta %*% B  -  crossprod(RmT_Q_zeta)
                                   # Mat2 <- t(B) %*% Q_zeta %*% cholsolve(chol_CtQoC_Q_zeta,t(C_m) %*% Qobs %*% Z)
                                   
                                   
                                   stopifnot(is.matrix(Z) | is(Z,"Matrix"))
                                   stopifnot(is.matrix(Zinv))
                                   stopifnot(is.matrix(C_m) | is(C_m,"Matrix"))
                                   stopifnot(is.matrix(Qobs) | is(Qobs,"Matrix"))
                                   stopifnot(is.matrix(B) | is(B,"Matrix"))
                                   stopifnot(is.matrix(X) | is(X,"Matrix"))
                                   stopifnot(is.matrix(Q_zeta) | is(Q_zeta,"Matrix"))
                                   stopifnot(is.matrix(S_f_trans) | is(S_f_trans,"Matrix"))
                                   stopifnot(is.numeric(lambda))
                                   stopifnot(is.numeric(ind))
                                   stopifnot(is.na(lambda_fix) | is.numeric(lambda_fix))
                                   
                                   stopifnot(nrow(C_m) == nrow(Z))
                                   stopifnot(nrow(Qobs) == nrow(Z))
                                   stopifnot(nrow(B) == ncol(C_m))
                                   stopifnot(nrow(X) == nrow(S_f_trans))
                                   stopifnot(nrow(Q_zeta) == ncol(C_m))
                                   
                                   ## Mole-fraction observation contribution
                                   CtQoC_Q_zeta <- as(crossprod(chol(Qobs) %*% C_m) + Q_zeta,"dgCMatrix")
                                   chol_permuted <- cholPermute(CtQoC_Q_zeta)
                                   Mat1 <- t(B) %*% Q_zeta %*% (B  - linalg:::cholsolve(Q=CtQoC_Q_zeta,
                                                                                        y=Q_zeta %*% B,
                                                                                        cholQp = chol_permuted$Qpermchol,
                                                                                        P = chol_permuted$P))
                                   Mat2 <- t(B) %*% Q_zeta %*% linalg:::cholsolve(Q = CtQoC_Q_zeta,
                                                                                  y = t(C_m) %*% Qobs %*% Z,
                                                                                  cholQp = chol_permuted$Qpermchol,
                                                                                  P = chol_permuted$P)
                                   
                                   ## Flux prior contribution
#                                    Q_f_trans <- chol2inv(chol(S_f_trans[ind,ind]))
#                                    cholQf_trans <- chol(Q_f_trans)
#                                    XtQX <- t(X) %*% Q_f_trans %*% X
#                                    PSI <- Q_f_trans - Q_f_trans %*% X %*% solve(XtQX) %*% t(X) %*% Q_f_trans
                                   
                                   n_inv <- ncol(Zinv)                                # Number of flux-field maps
                                   Q_f_trans <- chol2inv(chol(S_f_trans[ind,ind]))
                                   Q_f_trans <- do.call("bdiag",
                                                rep(list(Q_f_trans),n_inv+1))        # Repeat prec matrix n_obs times
                                   cholQf_trans <- chol(Q_f_trans)
                                   X2 <- do.call("rbind",
                                                 rep(list(X),n_inv+1))             # Repeat covariances n_obs times
                                   XtQX <- t(X2) %*% Q_f_trans %*% X2
                                   PSI <- Q_f_trans - Q_f_trans %*% X2 %*% solve(XtQX) %*% t(X2) %*% Q_f_trans
                                   
                                   n <- nrow(X)
                                   p <- ncol(X)
                                   n2 <- nrow(X2)                             # Number of flux data points (augmented)
                                   p2 <- ncol(X2)                             # Number of regressors
                                   
                                   # We have a data-dependent prior if we are not considering the Pericchi 
                                   # alternative AND lambda is a free parameter
                                   lmbda_prior <- !Pericchi & is.na(lambda_fix) 
                                   
                                   logf <- function(x) {
                                       Yf <- as.matrix(x)
                                       Ztot <- rbind(matrix(Yf),
                                                     matrix(as.vector(Zinv))) # Put observations into one vector
                                       gZ <- g(Ztot,lambda)
                                       
                                       if(any(Yf <= 0)) {
                                           return(Inf)
                                       } else {
                                           
#                                            betahat <- solve(XtQX) %*% t(X) %*% Q_f_trans %*% g(Yf,lambda)
#                                            qtilde <- crossprod(cholQf_trans %*% (g(Yf,lambda) - X %*% betahat))
#                                            -(
#                                                #-0.5 * crossprod(cholQo_zeta %*% (s_obs$z - C_m %*% B %*% Yf)) -
#                                                - 0.5 * t(Yf) %*% Mat1 %*% Yf +
#                                                    t(Yf) %*% Mat2 -
#                                                    0.5*(n-p)*log(qtilde) + 
#                                                    sum(log(g_dash(Yf,lambda)))    
#                                            ) %>% 
#                                                as.numeric() %>%
#                                                return()
                                           betahat <- solve(XtQX) %*% t(X2) %*% Q_f_trans %*% gZ
                                           qtilde <- crossprod(cholQf_trans %*% (gZ - X2 %*% betahat))
                                           -(
                                               #-0.5 * crossprod(cholQo_zeta %*% (s_obs$z - C_m %*% B %*% Yf)) -
                                               - 0.5 * t(Yf) %*% Mat1 %*% Yf +
                                                   t(Yf) %*% Mat2 -
                                                   0.5*(n2-(!Pericchi)*p2)*log(qtilde) +    # Pericchi prior => +p/2
                                                   (1 - lmbda_prior*p2/n2)*sum(log(g_dash(Yf,lambda)))    
                                           ) %>% 
                                               as.numeric() %>%
                                               return()
                                           ## the  - is.na(lambda_fix) is because if we have fixed lambda we do not have a data-dependent prior distribution anymore
                                           
                                       }
                                   }
                                   
                                   gr_logf <- function(x) {
#                                        Yf <- as.numeric(x)
#                                        betahat <- solve(XtQX) %*% t(X) %*% Q_f_trans %*% g(Yf,lambda)
#                                        qtilde <- crossprod(cholQf_trans %*% (g(Yf,lambda) - X %*% betahat))
#                                        
#                                        grYf <- 
#                                            #tB_tC_Qozeta_z - 
#                                            #tB_tC_Qozeta_C_B %*% Yf - 
#                                            Mat2 - Mat1 %*% Yf -
#                                            as.numeric((n-p)/qtilde) * J_lambda(Yf,lambda) %*% PSI %*% g(Yf,lambda) +
#                                            g_2dash(Yf,lambda) / g_dash(Yf,lambda)
#                                        -as.numeric(grYf[,1])
                                       
                                       Yf <- as.matrix(x)
                                       Ztot <- rbind(matrix(Yf),
                                                     matrix(as.vector(Zinv))) # Put observations into one vector
                                       gZ <- g(Ztot,lambda)
                                       betahat <- solve(XtQX) %*% t(X2) %*% Q_f_trans %*% gZ
                                       qtilde <- crossprod(cholQf_trans %*% (gZ - X2 %*% betahat))
                                       
                                       grYf <- 
                                           #tB_tC_Qozeta_z - 
                                           #tB_tC_Qozeta_C_B %*% Yf - 
                                           Mat2 - Mat1 %*% Yf -
                                           as.numeric((n2-(!Pericchi)*p2)/qtilde) * cBind(J_lambda(Yf,lambda),J_lambda(Yf,lambda)*0) %*% PSI %*% gZ +  # Pericchi prior => +p/2 in exponent
                                           (1 - lmbda_prior*p2/n2)*g_2dash(Yf,lambda) / g_dash(Yf,lambda)
                                       -as.numeric(grYf[,1])
                                       
                                       
                                   }
                                   
                                   list(logf = logf,
                                        gr_logf = gr_logf)
}



#' @title Conditional for mole-fraction field discrepancy parameters
#' @description Returns the conditional distribution of the mole-fraction field discrepancy parameters given the mole-fraction observations and the flux field 
#' @param theta mole-fraction field parameters -- see details
#' @param s observation spatial locations
#' @param t observation temporal indices
#' @param C_m observation incidence matrix (of size m x n_m)
#' @param Qobs observation precision matrix (of size m x m)
#' @param Z matrix of mole fraction observations (of size m x 1)
#' @param B source-receptor relationship matrix (of size n_m x n_f)
#' @param Q_zeta_fn function which returns the inverse of the mole-fraction discrepancy correlation matrix. This function should be of the form \code{ Q_zeta_fn(s,t,d_t,d_s)} where \code{s,t} are as above, \code{d_t} is the temporal length scale parameter and \code{d_s} is the spatial length scale parameter
#' @details To be completed
#' @export 
#' @examples
#' n_m <- m_obs <- 100
#' n_f <- 20
#' Yf <- matrix(rpois(n_f,20),n_f,1)    
#' log_f_theta_m( theta = c(4,0.8,1),
#'                Yf=Yf,
#'                s = 1:10,
#'                t = 1:10,
#'                C_m = .symDiagonal(n_m),
#'                Qobs = Diagonal(n_m),
#'                Z = 200 + matrix(rnorm(n_m,sd=100),n_m,1),
#'                B = matrix(rnorm(n_m*n_f),n_m,n_f),
#'                Q_zeta_fn = function(s,t,d_t,d_s) Diagonal(length(t) * length(s)))
log_f_theta_m <- function(theta,
                          Yf,
                          s,
                          t,
                          C_m,
                          Qobs,
                          Z,
                          B,
                          Q_zeta_fn) {
    s2  <- exp(theta[1])
    d_t <- theta[2]
    d_s <- exp(theta[3])
    
    ## We use priors theta[1] in (-2,20), theta[2] in (-1,1) and thata[3] in (0.1,5)
    
    if(theta[1] < -2 | theta[1] > 20 | abs(theta[2]) > 1 | theta[3] < log(0.1) | theta[3] > log(5)) {
        return(-Inf)
    } else {
        
        Q_zeta <- s2^(-1) * Q_zeta_fn(s,
                                      t,
                                      d_t, 
                                      d_s)
        
        chol_Q_zeta <- chol(Q_zeta)
        ybar <- t(C_m) %*% Qobs %*% Z + Q_zeta %*% B %*% Yf
        
        #  CtQoC_Q_zeta <- crossprod(chol(Qobs) %*% C_m) + Q_zeta
        #  chol_CtQoC_Q_zeta <- chol(CtQoC_Q_zeta)
        #  RmT_X <- backsolve(chol_CtQoC_Q_zeta,
        #                            ybar,
        #                            transpose = TRUE)  # X = solve(t(R)) %*% X
        
        CtQoC_Q_zeta <- as(crossprod(chol(Qobs) %*% C_m) + Q_zeta,"dgCMatrix")
        chol_permuted <- cholPermute(Q=CtQoC_Q_zeta)
        CtQoC_Q_zeta_inv_ybar <- linalg:::cholsolve(Q=CtQoC_Q_zeta,y=ybar,
                                                    cholQp=chol_permuted$Qpermchol,
                                                    P=chol_permuted$P)
        
        logf <- 0.5 * determinant(Q_zeta)$modulus -
            #0.5 * atminv:::logdet(chol_CtQoC_Q_zeta) -
            0.5 * atminv:::logdet(chol_permuted$Qpermchol) -
            0.5 * crossprod(chol_Q_zeta %*% B %*% Yf) +
            #0.5 * crossprod(RmT_X)
            0.5 * t(ybar) %*% CtQoC_Q_zeta_inv_ybar
        
        
        
        return(as.numeric(logf))
    }
}

#' @title Conditional for flux field parameters
#' @description Returns the conditional distribution of the flux field parameters given the flux field and, possibly, other flux maps. Uniform priors are assumed over the scale and smoothness, the bounds of which can be specified. 
#' @param theta flux field parameters -- see details
#' @param Z matrix of flux field data (of size m x n_obs, where n_obs is generally 2; an inventory and the sample derived from the mole fraction observations)
#' @param D distance matrix of flux field spatial locations
#' @param flux_cov_fn covariance function used to generate the transformed correlation matrix. This should be of the form flux_cov_fn(D,scale,nu), where scale and nu are scale and smoothness parameters. 
#' @param X covariates matrix (of size n_f x p)
#' @param scale_ul upper limit of uniform distribution over scale parameter. Default of 2
#' @param nu_ul upper limit of uniform distribution over smoothness parameter. Default of 2
#' @param lambda_l absolute limit of uniform distribution over Box-Cox parameter. Default of 3
#' @param lambda_fix use a fixed value for lambda (do not sample from it). If this is set then \code{theta} should only contain two numbers
#' @details To be completed
#' @export 
#' @examples 
#' n_f <- 20
#' s <- 1:n_f
#' cov_fn <- function(D,scale,nu) {
#'      exp(-scale*(D^nu))    
#' }
#' print(log_f_theta_f(theta = c(1,1,0),
#'               Z = matrix(rpois(n_f,100)),
#'               D = as.matrix(dist(s)),
#'               flux_cov_fn = cov_fn,
#'               X = matrix(rnorm(n_f),n_f,1)))
log_f_theta_f <- function(theta,
                          Z,
                          D,
                          flux_cov_fn,
                          X,
                          scale_ul = 2,
                          nu_ul = 2,
                          lambda_ll = -3,
                          lambda_ul = 3,
                          lambda_fix = NA,
                          uncorrelated = FALSE,
                          Pericchi=TRUE) {
    
      if(uncorrelated & !is.na(lambda_fix)) 
        stop("lambda fixed and correlation matrix is diagonal -- Nothing to sample!")
      
      if(!uncorrelated) {
        scale <- theta[1]
        nu <- theta[2]
        if(is.na(lambda_fix)) {
          lambda <- theta[3]
        } else {
          lambda <- lambda_fix
        }
      } else { # if uncorrelated == T
        lambda <- theta[1]
        scale <-scale_ul/2 
        nu <- nu_ul / 2     # Just dummy variables to pass check later on -- ignored
      }
     
    stopifnot(scale_ul >0 & nu_ul > 0)
    
    if(scale < 0 | scale > scale_ul | nu < 0 | nu > nu_ul | lambda <lambda_ll | lambda > lambda_ul) {
        XX <- -Inf
    } else {
        if(!uncorrelated) {
          Sigma <- flux_cov_fn(D,scale,nu)  # Covariance in transformed space
        } else {
          Sigma <- diag(nrow(D))
        }
        R <- chol(Sigma)                  # Cholesky
        Q <- chol2inv(R)                  # Precision matric
        
        n_obs <- ncol(Z)                                # Number of flux-field maps
        Q <- do.call("bdiag",rep(list(Q),n_obs))        # Repeat precision matrix n_obs times
        X2 <- do.call("rbind",rep(list(X),n_obs))       # Repeat covariances n_obs times
        Ztot <- matrix(as.vector(Z))                    # Put observations into one long vector
        n2 <- nrow(X2)                                  # Number of flux data points (augmented)
        p2 <- ncol(X2)                                  # Number of regressors
        
        XtQX <- t(X2) %*% Q %*% X2
        gZ <- g(Ztot,lambda)
        
        beta_hat <- solve(XtQX) %*% t(X2) %*% Q %*% gZ
        qtilde <- t((gZ - X2 %*% beta_hat)) %*% Q %*% (gZ - X2 %*% beta_hat)
        log_J <- sum(log(g_dash(Ztot,lambda))) ## a constant if lambda_fix = TRUE
        
        lmbda_prior <- !Pericchi & is.na(lambda_fix) 
        
        XX <- as.numeric(
            #0.5*determinant(Q)$modulus -
            - n_obs*0.5*determinant(Sigma)$modulus - #recall Sigma is half the size of Q so this = -0.5|Q|
                0.5*determinant(XtQX)$modulus -
                (n2-(!Pericchi)*p2)/2*log(qtilde) +  # Pericchi prior => +p/2 in exponent
                (1 - lmbda_prior*p2/n2) * log_J
        )
    }
    XX
}











## DIRECT OBSERVATIONS OF FLUX FIELD
Yf_BC_conditionals_DO <- function(Z,
                               Qobs,
                               X,
                               S_f_trans,
                               lambda=0,
                               ind = 1:nrow(S_f_trans))               {
    
    stopifnot(is.matrix(Z) | is(Z,"Matrix"))
    stopifnot(is.matrix(Qobs) | is(Qobs,"Matrix"))
    stopifnot(is.matrix(X) | is(X,"Matrix"))
    stopifnot(is.matrix(S_f_trans) | is(S_f_trans,"Matrix"))
    stopifnot(is.numeric(lambda))
    stopifnot(is.numeric(ind))
    
    stopifnot(nrow(Qobs) == nrow(Z))
    stopifnot(nrow(X) == nrow(S_f_trans))
    
    Q_f_trans <- chol2inv(chol(S_f_trans[ind,ind]))
    cholQf_trans <- chol(Q_f_trans)
    XtQX <- t(X) %*% Q_f_trans %*% X
    PSI <- Q_f_trans - Q_f_trans %*% X %*% solve(XtQX) %*% t(X) %*% Q_f_trans
    
    n <- nrow(X)
    p <- ncol(X)
    
    logf <- function(x) {
        Yf <- as.matrix(x)
        Ztot <- Yf
        gZ <- g(Ztot,lambda)
        
        if(any(Yf <= 0)) {
            return(Inf)
        } else {
            betahat <- solve(XtQX) %*% t(X) %*% Q_f_trans %*% gZ
            qtilde <- crossprod(cholQf_trans %*% (gZ - X2 %*% betahat))
            -(
                - 0.5 * t(Yf) %*% Qobs %*% Yf +
                    t(Yf) %*% Qobs %*% Z -
                    0.5*(n2-p2)*log(qtilde) +        
                    sum(log(g_dash(Yf,lambda)))    
            ) %>% 
                as.numeric() %>%
                return()
            
        }
    }
    
    gr_logf <- function(x) {
        Yf <- as.matrix(x)
        Ztot <- Yf
        gZ <- g(Ztot,lambda)
        betahat <- solve(XtQX) %*% t(X) %*% Q_f_trans %*% gZ
        qtilde <- crossprod(cholQf_trans %*% (gZ - X %*% betahat))
        
        grYf <- 
            Qobs %*% Z - Qobs %*% Yf -
            as.numeric((n-p)/qtilde) * cBind(J_lambda(Yf,lambda),J_lambda(Yf,lambda)*0) %*% PSI %*% gZ +
            g_2dash(Yf,lambda) / g_dash(Yf,lambda)
        -as.numeric(grYf[,1])
        
        
    }
    
    list(logf = logf,
         gr_logf = gr_logf)
}
