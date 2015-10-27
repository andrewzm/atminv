#### CLASS DEFINITIONS ###########

#'  @docType class
#'  @title Parent (virtual) class for all blocks (latent fields and observations)
#'
#' @description A block is the most basic component of a multi-variate spatio-temporal model. It can either consist of a latent field or an observations. Each block 
#' has a unique identifier (uid) and some cardinality n, the meaning of which depends on the specific class
#' 
#' @keywords Block
#' @rdname blockclass 
setClass("block", representation(uid="numeric",n="numeric","VIRTUAL"))

#'  @docType class
#'  @title Process block (virtual)
#'
#' @description A process block inherits from class \code{block}, however can only be used in certain ways in a graph. For example it can be linked to an observation
#' data block but not to another process block. It is the parent class of all latent processes.
#' 
#' @keywords Block, process
#' @rdname processclass
setClass("process", contains="block",representation("VIRTUAL"))

#'  @docType class
#'  @title GMRF
#'
#' @description A basic GMRF object which can also take an intrinsic value for future definitions (intrinsic GMRFs not implemented yet). This
#' class inherits from \code{process} which in turn is a block.
#' @rdname GMRFclass
setClass("GMRF",contains="process",
         representation(mu="matrix", Q="Matrix",intrinsic="numeric",rep="data.frame",
                        t_axis="numeric"),
         prototype(mu=matrix(0,2,1),
                   Q = sparseMatrix(i=c(1,2),j=c(1,2),x=1),
                   intrinsic=0,
                   n=2,rep=data.frame(),t_axis=0))

#'  @docType class
#'  @title GMRF_RW
#'
#' @description A random walk  which inherits from class \code{GMRF}. The primary difference is that this class is constructed using
#' a first-order auto-regressice structure. All random walks are intrinsic GMRFs.
#' @rdname GMRF_RWclass
setClass("GMRF_RW",contains="GMRF")


#'  @docType class
#'  @title GMRF_AR
#'
#' @description An autoregressive proces which inherits from class \code{GMRF}. The primary difference is that this class is constructed using
#' a first-order auto-regressive structure. Unlike the random walk, this is not an intrinsic GMRF.
#' @rdname GMRF_ARclass
setClass("GMRF_AR",contains="GMRF")

#'  @docType class
#'  @title observation block
#'
#' @description An observation block inherits from class \code{block}, however can only be used in certain ways in a graph. For example it can be linked to an observation
#' process block but not to another observation block. To initialise use \code{initObs}.
#' @rdname Obsclass
setClass("Obs",contains=c("block"),representation(df="data.frame",args="list"))