\name{miInference}

\docType{methods}

\alias{miInference}
\alias{print.miInference}

\title{ Combine results from analyses after multiple imputation}

\description{
This function combines the results from data analyses performed after
multiple imputation using methods described by Rubin (1987) and others.
}


\usage{

miInference( est.list, std.err.list, method = "scalar",
   df.complete = NULL )

\method{print}{miInference}( x, \ldots)

}


\arguments{

  \item{est.list}{a list of estimates to be combined.  Each component of
  this list should be a scalar or vector containing point estimates 
  from the analysis of an imputed dataset.  This list should have
  \emph{M} components, where \emph{M} is the number of imputations,
  and all components should have the same length.}

  \item{std.err.list}{a list of standard errors to be combined.  Each
  component of this list should be a scalar or vector of
  standard errors associated with the estimates in \code{est.list}. }

  \item{method}{how are the estimates to be combined? At present,
  the only type allowed is \code{"scalar"}, which means that 
  estimands are treated as one-dimensional entities.  If
  \code{est.list} contains 
  vectors, inference for each element of the vector is carried out
  separately; covariances among them are not considered.}

  \item{df.complete}{degrees of freedom assumed for the complete-data
  inference.  This should be a scalar
  or a vector of the same length as the components of \code{est.list}
  and \code{std.err.list}.}

  \item{x}{ a result from \code{miInference}. }

  \item{\dots}{values to be passed to the methods.}
}

\details{

  If \code{df.complete = NULL} or \code{Inf}, the degrees of
  freedom are computed by 
  the method of Rubin (1987, Chap.3),  which assumes that if there were no
  missing data, the usual normal approximation for large samples
  would be appropriate, i.e. that 
  a 95\% interval would be computed as the estimate
  plus or minus 1.96 standard errors.  Otherwise, the degrees of
  freedom are computed by the method of Barnard and Rubin (1999),
  which assumes that an approximate 95\% interval without missing data
  would be the estimate plus or
  minus \code{qt(.975, df.complete)} standard errors.

  The result from this function is a list whose \code{class}
  attribute has been set to \code{"miInference"}.  If this list is displayed
  or printed via the generic function \code{print}, it will be
  formatted into a table resembling the output from a regression
  analysis with columns for the estimates, standard errors,
  t-ratios (estimates divided by their standard errors) and p-values
  for testing the null hypothesis that each estimate is zero.
}

\value{ a list with the following components:

\item{names}{character-string labels for the estimands. This is
derived from the \code{names} attribute, if any, of the 
components of \code{est.list}.}

\item{est}{combined estimate(s).}

\item{std.err}{standard error(s) for \code{est}.}

\item{df}{degrees of freedom for Student-t approximation. For
example, 95\% intervals can be computed as \code{est} plus or minus
\code{qt(.975,df)*std.err}.}

\item{p}{p-value(s) for testing the null hypothesis that each estimand
is zero against a two-tailed alternative.}

\item{rel.incr}{ estimated relative increase(s) in variance due to
nonresponse. }

\item{mis.inf}{ estimated rate(s) of missing information.}

}

\note{

Rubin (1987) defined the rate of missing information
as \code{rel.incr + 2/(df+3)} divided by \code{(rel.incr+1)}, which
estimates the information lost due to missing values and due to the
fact that the number of multiple imputations is finite.   We define it
as \code{rel.incr} divided by \code{(rel.incr+1)}, the information
lost due 
to missing values, which is consistent with the formulas of
Barnard and Rubin (1999).

}

\references{

Barnard, J. and Rubin, D.B. (1999) Small-sample degrees of freedom
with multiple imputation. \emph{Biometrika}, 86, 948-955.
\cr

Rubin, D.B. (1987) \emph{Multiple Imputation for
Nonresponse in Surveys}. New York: Wiley.
\cr

For more information about functions in
the NORM package, see the manual \emph{NORM Package for R, Version 2}
in the subdirectory \code{doc} where the library has been installed.

}

\author{Joe Schafer \email{Joseph.L.Schafer@census.gov} }


\examples{

## create 25 multiple imputations for the cholesterol data
## using functions from the NORM library
data(cholesterol)
emResult <- emNorm(cholesterol)
set.seed(234)
mcmcResult <- mcmcNorm(emResult, iter=2500, impute.every=100)

## get estimates and standard errors for the mean of Y1
## minus the mean of Y3
est.list <- as.list(NULL)
std.err.list <- as.list( NULL )
for( m in 1:25 ){
   yimp <- mcmcResult$imp.list[[m]]  # one imputed dataset
   diff <- yimp[,"Y1"] - yimp[,"Y3"]
   est.list[[m]] <- mean(diff)
   std.err.list[[m]] <- sqrt( var(diff) / length(diff) ) }

## combine the results by rules of Barnard and Rubin (1999)
## with df.complete = 27, because a paired t-test in a dataset with
## N=28 cases has 27 degrees of freedom
miResult <- miInference(est.list, std.err.list, df.complete=27)
print(miResult)

}

\keyword{ NA }
