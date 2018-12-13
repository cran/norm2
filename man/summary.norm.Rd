\name{summary.norm}

\alias{summary.norm}
\alias{print.summary.norm}

\title{Summarize information from EM or MCMC algorithms}

\description{
Method for summarizing the results from a call to the functions
\code{\link{emNorm}} or \code{\link{mcmcNorm}.}
}

\usage{
\method{summary}{norm}(object, show.variables = (object$method == "EM"),
   show.patterns = (object$method == "EM"),
   show.params = (object$method =="EM"), \ldots) 

\method{print}{summary.norm}(x, \ldots)
}


\arguments{
  \item{object}{an object of class \code{"norm"} which is produced as the
  result of a call to \code{emNorm} or \code{\link{mcmcNorm}}; see
  DETAILS below.} 

  \item{show.variables}{ if \code{TRUE}, then tables summarizing the
  response variables and covariates used in the model will be printed.}

  \item{show.patterns}{ if \code{TRUE}, then a table summarizing the
  patterns of missingness will be printed.}

  \item{show.params}{ if \code{TRUE}, then the final values of the
  parameters (estimates after the final iteration of EM or simulated
  values after the final iteration of MCMC) will be printed.}

  \item{x}{ a result from \code{summary.norm}. }

  \item{\dots}{values to be passed to the methods.}

}


\details{
  The result from a call to \code{emNorm} or \code{mcmcNorm}
  is an object of 
  class \code{"norm"}, which is a list containing results from the EM
  or MCMC run.  The function \code{summary.norm}, which is
  invoked through the generic method \code{summary}, summarizes the
  information contained in this object. 

  The result from \code{summary.norm} is an object of
  class \code{"summary.norm"} which can be displayed or printed via
  the generic method \code{print}.
}


\value{
  A list that includes all
  the original components of \code{obj} plus some additional summaries
  that are printed via a call to the generic method \code{print}.
  These include:

  \item{x.table}{a summary of all variables appearing in the model as
  predictors or covariates.}

  \item{y.table}{a summary of all variables appearing in the model as
  responses or outcomes.}

  \item{em.summary}{a summary of the results from the EM run,
  including: the number of iterations; whether EM converged; and an
  empirical estimate of the rate of convergence which estimates the
  worst fraction of missing information.}

  \item{mcmc.summary}{a summary of the results from the MCMC run,
  including: the number of iterations; whether imputations were created
  and, if so, how many; and whether parameter series were saved.}

}

\references{
For more information about this function and other functions in
the \code{norm2} package, see  \emph{User's Guide for \code{norm2}}
in the library subdirectory \code{doc}.
}

\author{Joe Schafer \email{Joseph.L.Schafer@census.gov} }

\seealso{
\code{\link{emNorm}}, \code{\link{mcmcNorm}}
}

\examples{
## run EM for cholesterol data and summarize
data(cholesterol)
emResult <- emNorm(cholesterol)
summary(emResult)

## run MCMC starting from the ML estimates and summarize
mcmcResult <- mcmcNorm(emResult)
summary(mcmcResult)
}

\keyword{ multivariate }
\keyword{ NA }
