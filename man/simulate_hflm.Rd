\name{simulate_hflm}
\alias{simulate_hflm}
\title{Simulate data for HFLM}
\description{
  Function used to simulate data sets to illustrate the use of hflm. 
  \code{RcppArmadillo.package.skeleton()} 
}
\usage{
  simulate_hflm(nobs, n.tau, delta, varx, eSNR, plot = TRUE)
}
\arguments{
  \item{x}{a numeric vector}
}
\value{
  \code{ytau} A matrix of \code{nobs} functional responses.
  
  \code{ytau} A matrix of \code{nobs} functional predictors.

  
  \code{theta.s.tau} Regression surface..
  
  \code{vary} Noise s
  
}
\details{
  These are example functions which should be largely
  self-explanatory. Their main benefit is to demonstrate how to write a
  function using the Armadillo C++ classes, and to have to such a
  function accessible from R.
}
\references{
  See the documentation for Armadillo, and RcppArmadillo, for more details.
}
\examples{
  x <- sqrt(1:4)
  rcpparma_innerproduct(x)
  rcpparma_outerproduct(x)
}
\author{Dirk Eddelbuettel}

