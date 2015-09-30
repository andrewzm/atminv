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
                     if(is.null(n)) n <- nrow(Q)
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


#' @title Autoregressive GMRF
#' 
#' @description This function initialises an autoregressive block and represents it as a Gaussian Markov Random Field with mean \code{mu} and precision matrix \code{Q}. Only AR1 processes along the real line  are implemented for now. Also a data frame can be specified with more details on the GMRF.
#'
#' @param n number of vertices
#' @param a autoregression coefficient
#' @param precinc precision constant (multiplies the template precision matrix)
#' @param order order of AR process (only AR1 implemented for now)
#' @param df data frame of length \code{n} with more details (for example axis, covariate information)
#' @return Object of class GMRF with zero mean
#' @keywords GMRF, random walk
#' @export
#' @examples
#'
#' require(Matrix)
#' # Create a first-order random walk GMRF
#' my_RW <- GMRF_AR1(n=10, a = 0.8, precinc =2)
#' print(getPrecision(my_RW))
#' print(getMean(my_RW))
#' print(getDf(my_RW))
GMRF_RW <- function(n = 10,a=0.8,order=1,precinc = 1,df=data.frame()) {
    
    stopifnot(order %in% c(1))
    stopifnot(abs(a) < 1)
    if(is.null(n)) n <- nrow(mu)
    
    mu <- matrix(0,nrow=n,ncol=1)
    i = c(1:(n-1),1:(n-1))
    j = c(1:(n-1),2:n)
    x <- numeric(length=((n-1)*2))
    x[1:(n-1)] = -a
    x[n:((2*n)-2)] = 1
    Dmat = sparseMatrix(i,j,x=x)
    R = t(Dmat)%*%Dmat
    if (order == 1) {
        Q = precinc*R
        intrinsic = 0
        Q[1,1] <- 1
    }
    return(new("GMRF",mu=mu,Q=Q,intrinsic=1,n=n,t_axis=0:(n-1),rep=df))
}


#' @title Observation block
#' 
#' @description This function initialises an object of class \code{Obs} which defines a an observation data set. By default, this is for observations with negligible spatial footprint. For larger supports, use \code{Obs_poly}. 
#'
#' @param df a data frame which should contain at least 5 entries, \code{x,y,t,z} and \code{std} which denote the horizontal, vertical and temporal indices of the observations, the value and error respectively. Alternatively this could be a path name.
#' @param name the name of the observation process
#' @param remove_cross_ins removes data outside a circle centred at zero with specified radius. Convenient when working with satellite data in polar stereographic projection when some cross-ins are detected.
#' @param ... other arguments passed on to \code{preprocess_obs}
#' @return Object of class \code{Obs} (which inherits from class \code{block} and is thus also a block)
#' @keywords Observations, change of support, block
#' @export
#' @examples
#' O <- Obs(df=data.frame(x=runif(5),y=runif(5),t=c(1,1,1,2,2),z=runif(5),std=runif(5)))
#' print(O)
#' plot(subset(O,t==1),"z",pt_size=4)
Obs <- function(df,name="Obs",remove_cross_ins=0,...) {
    return(new("Obs",df=df,name=name,remove_cross_ins=remove_cross_ins,...))
}


# ... is passed on to preprocess_obs
setMethod("initialize",signature="Obs",function(.Object,df,name=NA,remove_cross_ins=0,pol=NA,alpha0=NA,av_dist=NA,...) { 
    
    args<-list(...)
    args <- c(args,df=df,name=name,remove_cross_ins=remove_cross_ins,pol=pol,alpha0=alpha0,av_dist=av_dist)
    .Object@args <- args
    
    stopifnot((is.character(df)) | is.data.frame(df))
    
    if(is.character(df)) {
        cat(paste("Loading from",df),sep="\n")
        data_df <- read.table(df,header=T)    
    } else {
        data_df <- df 
    }
    
    .Object@df <- data_df
    .Object <- preprocess_obs(.Object,...)
    if (is.null(data_df$obs_name))
        .Object["obs_name"] <- as.factor(name)
    if(remove_cross_ins > 0) {
        .Object@df <- subset(.Object@df,sqrt(x^2 + y^2) > remove_cross_ins)
    }
    
    if(!is.na(pol[1])) {
        poly_points <- pol
        if (!("id" %in% names(.Object@df))) stop("Need to merge by id field which is not supplied")
        .Object@df <- merge(poly_points,.Object@df,by=c("id","t"))
        .Object@df <- arrange(.Object@df,id,t)
        .Object@df2 <- .expand_poly(.Object@df)
    }
    .Object@df$n <- 1:nrow(.Object@df)
    .Object@n <- nrow(.Object@df)
    
    if("cmweq" %in% names(.Object@df)) {
        .Object@df$z <-  as.numeric(.Object@df$cmweq)*0.01*.Object@df$area2    #Convert to Mt
        .Object@df2$z <-  as.numeric(.Object@df2$cmweq)*0.01*.Object@df2$area2    #Convert to Mt
        .Object@df$std <-  as.numeric(.Object@df$std)*0.01*.Object@df$area2    #Convert to Mt
        .Object@df2$std <-  as.numeric(.Object@df2$std)*0.01*.Object@df2$area2    #Convert to Mt 
    }
    
    
    if(!is.na(alpha0)) {
        if(is.na(alpha0)) stop("Cannot specify alpha0 without averaging distance")
        .Object@args$P <-  Find_Smooth_mat(subset(.Object@df,t==0),alpha0,av_dist)
    }
    
    callNextMethod(.Object)})

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

#' @title Pre-processes observations
#' @description This function simply takes a data frame and performs standard pre-processing functions such as removing obvious outliers and averaging data
#' in computationally efficient way over a regular grid
#'
#' @param Obs.obj The observation object
#' @param std_lim Values with error higher than \code{std_lim} are deleted
#' @param abs_lim Values which are bigger than \code{abs_lim} are deleted
#' @param avr_method If \code{mean} then the mean value and mean error on a sub-grid are used to sub-sample the observations. If \code{median} the median value and MAD of the z-values are used
#' for subsampling. If \code{mean_no_std} or \code{median_no_std}, the error is not computed ignored.
#' @param box_size The grid width over which the observations are averaged
#' @param min_pts If a grid box contains less than \code{min_pts} it is marked as empty
#' @param ... Other parameters which are ignored
#' @export
preprocess_obs <- function(Obs.obj,std_lim=NA,abs_lim=NA,avr_method=NA,box_size=10,min_pts=4,...) {
    
    if(!(is.na(std_lim))) {
        stopifnot(is.numeric(std_lim))
        Obs.obj@df <- subset(Obs.obj@df,std < std_lim)
    }
    
    if (!is.na(abs_lim)) {
        stopifnot(is.numeric(abs_lim))
        Obs.obj@df <- subset(Obs.obj@df,abs(z) < abs_lim)
    }
    
    if(any(Obs.obj@df$std == 0)) {
        warning("Data points with zero error detected. These are automatically deleted")
        Obs.obj@df <- subset(Obs.obj@df,std > 0)
    }
    
    
    if(!is.na(avr_method)) {
        stopifnot(avr_method %in% c("mean","median","mean_no_std","median_no_std"))
        
        breaksx <- seq(min(Obs.obj@df$x)-1,max(Obs.obj@df$x)+1,by=box_size)
        breaksy <- seq(min(Obs.obj@df$y)-1,max(Obs.obj@df$y)+1,by=box_size)
        Obs.obj@df$box_x <- cut(Obs.obj@df$x,breaksx,labels=F)
        Obs.obj@df$box_y <- cut(Obs.obj@df$y,breaksy,labels=F) 
        
        
        averaging_fun <- function(z) {
            x <- z[1]  
            if(length(z) > min_pts) {
                if(avr_method %in% c("median","median_no_std")) {
                    x <- median(z)
                } else if(avr_method  %in% c("mean","mean_no_std") ) {
                    x <- mean(z)
                }
            }
            return(x)
        }
        
        std_fun <- function(z,std) {
            x <- std[1]  
            if(length(z) > 4) {
                if(avr_method=="median") {
                    x <- mad(z)
                } else if(avr_method=="mean") {
                    x <- sd(z)
                }
            }
            return(x)
        }
        
        boxes <- group_by(Obs.obj@df, box_x,box_y,t)
        if(avr_method %in% c("median","mean")) { 
            Obs.obj@df <- data.frame(summarise(boxes,x=round(mean(x)),y=round(mean(y)),z2=averaging_fun(z),std=std_fun(z,std)))
        } else if(avr_method %in% c("median_no_std","mean_no_std")) {
            Obs.obj@df <- data.frame(summarise(boxes,x=round(mean(x)),y=round(mean(y)),z2=averaging_fun(z)))
        }
        Obs.obj@df$z <- Obs.obj@df$z2
        Obs.obj@df$z2 <- NULL
        
        Obs.obj@df <- arrange(Obs.obj@df,t)
        Obs.obj@n <- 1:nrow(Obs.obj@df)
    }  
    return(Obs.obj)
}

setMethod("[",signature = "GMRF",function(x,i,j) { return(x@rep[i][,])})

#' @aliases [<-,GMRF,ANY,ANY-method
#' @export
setMethod("[<-",signature = "GMRF",function(x,i,j,value) {
    x@rep[i] <- value
    return(x)
})
setMethod("[",signature = "Obs",function(x,i,j) { return(x@df[i][,])})

#' @aliases [<-,Obs,ANY,ANY-method
#' @export
setMethod("[<-",signature = "Obs",function(x,i,j,value) {
    x@df[i] <- value
    return(x)
})
