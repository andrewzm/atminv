#### CONSTRUCTORS ###########

#' @title GMRF
#' 
#' @description This function initialises a GMRF with some mean \code{mu} and precision matrix \code{Q}. The returned object is of class \code{GMRF}
#'
#' @param mu mean, of class \code{matrix}
#' @param Q sparse precision matrix (of class \code{Matrix})
#' @param intrinsic set intrinsic level if Q is singular
#' @param n number of nodes. Note that this can be different from \code{nrow(mu)} if the system is multi-variate
#' @param rep data frame of length \code{N} with more details (for example axis, covariate information)
#' @param t_axis this is the time horizon for the system under consideration. If you are considering a spatial problem set this to zero.
#' @return Object of class GMRF
#' @keywords GMRF
#' @export
#' @examples
#'
#' require(Matrix)
#' # Create a GMRF
#' Q <- sparseMatrix(i=c(1,2,1,2),j=c(1,1,2,2),x=c(1,0.1,0.1,1))
#' mu <- matrix(c(0,0))
#' my_GMRF <- GMRF(mu=mu, Q=Q)
#' print(getPrecision(my_GMRF))
#' print(getMean(my_GMRF))
#' print(getDf(my_GMRF))
GMRF <- function(mu=NA,Q= sparseMatrix(i=c(1,2),j=c(1,2),x=4),intrinsic=0,n=NULL,t_axis=0,
                 rep=data.frame()) {
                     
                     if(any(is.na(mu))) {
                         mu <- matrix(0,nrow(Q),1)
                     }
                     stopifnot(nrow(mu) == nrow(Q))
                     return(new("GMRF",mu=mu,Q=Q,intrinsic=intrinsic,n=n,t_axis=t_axis,rep=rep)) 
                 }

#' @title Random walk GMRF
#' 
#' @description This function initialises a random walk and represents it as a Gaussian Markov Random Field with mean \code{mu} and precision matrix \code{Q}. Only random walks along the real line, first-order and second-order variants are implemented for now. Also a data frame can be specified with more details on the GMRF.
#'
#' @param n number of vertices
#' @param order 1 or 2, depending on the order of the random walk
#' @param precinc precision constant (multiples the template precision matrix)
#' @param df data frame of length \code{n} with more details (for example axis, covariate information)
#' @return Object of class GMRF with zero mean
#' @keywords GMRF, random walk
#' @export
#' @examples
#'
#' require(Matrix)
#' # Create a first-order random walk GMRF
#' my_RW <- GMRF_RW(n=10, order=1, precinc =2)
#' print(getPrecision(my_RW))
#' print(getMean(my_RW))
#' print(getDf(my_RW))
GMRF_RW <- function(n = 10,order=1,precinc = 1,df=data.frame()) {
    
    stopifnot(order %in% c(1,2))
    
    if(is.null(n)) n<- nrow(mu)
    
    mu <- matrix(0,nrow=n,ncol=1)
    i = c(1:(n-1),1:(n-1))
    j = c(1:(n-1),2:n)
    x <- numeric(length=((n-1)*2))
    x[1:(n-1)] = -1
    x[n:((2*n)-2)] = 1
    Dmat = sparseMatrix(i,j,x=x)
    R = t(Dmat)%*%Dmat
    if (order == 1) {
        Q = precinc*R
        intrinsic = 1
    }
    if (order == 2) {
        R <- R %*% R
        R[1,(1:3)] = c(1,-2,1)
        R[2,(1:4)] = c(-2,5,-4,1)
        R[(n-1),(n-3):n] = c(1,-4,5,-2)
        R[(n),(n-2):n] = c(1,-2,1)
        Q <- precinc*R
        intrinsic = 2
    }
    return(new("GMRF",mu=mu,Q=Q,intrinsic=1,n=n,t_axis=0:(n-1),rep=df))
}


#### Sample GMRF ###########
#' @title Sample from a GMRF
#' @description Takes a GMRF object and, possibly, associated permutation matrix and Cholesky factor of permuted precision matrix, to generate
#' samples.
#' @param G object of class \code{GMRF}
#' @param L Cholesky factor of the precision matrix, if \code{P} is \code{NULL} then this is treated as the factor of an unpermuted matrix
#' @param reps number of samples to generate
#' @param P permutation matrix
#' @return a matrix of size \code{n} by \code{reps}
#' @export
#' @examples
#' G <- GMRF_RW()
#' G <- setPrecision(G, getPrecision(G) + 0.01*Imat(nrow(G)))
#' X <- sample_GMRF(G,reps=10)
setGeneric("sample_GMRF", function(G,L=NULL,reps=1,P=NULL) standardGeneric("sample_GMRF"))

setMethod("sample_GMRF",signature="GMRF",function(G,L=NULL,reps=1,P=NULL) {
    n = G@n
    z <- matrix(rnorm(n*reps),n,reps)
    mu <- c(G@mu)
    if(is.null(L)) {
        L <-t(chol(G@Q))
        P <-NULL
    }
    
    # Algorithm 2.4, Rue and Held
    if (is.null(P)) {
        v <- solve(t(L),z)
        x <- mu + v
    } else {
        v <-  P %*% solve(t(L),z)
        x <- mu + v
    }
    ifelse(reps==1, return(as.vector(x)), return(x))
})