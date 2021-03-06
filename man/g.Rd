% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BoxCox.R
\name{g}
\alias{J_lambda}
\alias{g}
\alias{g_2dash}
\alias{g_dash}
\alias{ginv}
\title{Box Cox transformations and related}
\usage{
g(Z, lambda)

g_dash(Z, lambda)

g_2dash(Z, lambda)

ginv(Y, lambda)

J_lambda(Z, lambda)
}
\arguments{
\item{Z}{the observations (these have to be positive)}

\item{lambda}{the Box-Cox parameter}
}
\description{
These functions are used to transform from the original space to the the transformed space. In general, if \eqn{Z} are the non-Gaussian observations, then \eqn{g(Z)} is assumed to be Gaussian
}
\details{
The function \code{g} is the transformation, \code{g_dash} is the derivative of \code{g}, \code{g_2dash} is the second derivative of \code{g}, and \code{ginv} is the inverse of \code{g}. \code{J_lambda} is the Jacobian and is a diagonal matrix with the derivative of \code{g} evaluated at the observations as its elements.
}
\examples{
Z <- 1:10 
g(Z,0.4)
g_dash(Z,0.4)
g_2dash(Z,0.4)
J_lambda(Z,0.4)
ginv(g(Z,0.4),0.4)
}

