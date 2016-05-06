\name{norm2-package}

\alias{norm2-package}
\alias{norm2}

\docType{package}

\title{
Analysis of incomplete multivariate data under a normal model
}

\description{
Functions for estimation and multiple imputation from
incomplete multivariate data under a normal model
}

\details{
The \code{norm2} package provides functions for analyzing incomplete
multivariate 
data using techniques and algorithms described by Schafer (1997).  The
name of this package derives from the assumed model for the complete
data, which is a multivariate normal model. The 
major functions are:
\preformatted{
   emNorm        EM algorithm estimating model parameters
   mcmcNorm      MCMC algorithm for simulating parameters and missing values
   impNorm       Simulate or predict missing values
   loglikNorm    Loglikelihood function
   logpostNorm   Log-posterior density function
   miInference   Combine results from analyses after multiple imputation
}
The package also includes three datasets:
\preformatted{
   cholesterol   Cholesterol levels for heart-attack patients
   flas          Foreign Language Attitude Scale
   marijuana     Changes in heart rate after marijuana use
}
}

\note{
Fortran source code written by the author for a much earlier version
called \code{norm} was ported to an R package by Alvaro A. Novo and
distributed through the Comprehensive R Archive
Network (CRAN). The old package \code{norm} is still available on
CRAN, but it has some 
major limitations (e.g., it does not work reliably when the number of
variables exceeds 30) and the author does not recommend its use.
}

\references{
Schafer, J.L. (1997) \emph{Analysis of Incomplete Multivariate
Data}. London: Chapman & Hall/CRC Press.
\cr

For more information about functions in
\code{norm2}, see  \emph{User's Guide for \code{norm2}}
in the library subdirectory \code{doc}.
}

\author{
Joseph L. Schafer <Joseph.L.Schafer@census.gov>

Maintainer: Joseph L. Schafer <Joseph.L.Schafer@census.gov>
}
\keyword{ package }
