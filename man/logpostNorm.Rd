\name{logpostNorm}

\docType{methods}

\alias{logpostNorm}
\alias{logpostNorm.default}
\alias{logpostNorm.formula}
\alias{logpostNorm.norm}

\title{ Observed-data log-posterior density for incomplete
multivariate normal data} 

\description{
Computes the observed-data log-posterior density function 
at given parameter values for an incomplete dataset
under a normal model.
}


\usage{
% the generic function
logpostNorm(obj, \dots)

% the default method
\method{logpostNorm}{default}(obj, x = NULL, intercept = TRUE, param, 
   prior = "uniform", prior.df = NULL, prior.sscp = NULL, \ldots)

% method for class formula
\method{logpostNorm}{formula}(formula, data, param, 
   prior = "uniform", prior.df = NULL, prior.sscp = NULL, \ldots)

% method for class norm
\method{logpostNorm}{norm}(obj, param = obj$param, prior = obj$prior,
   prior.df = obj$prior.df, prior.sscp = obj$prior.sscp, \ldots)

}


\arguments{

  \item{obj}{an object used to select a method.  It may be \code{y},
  a numeric matrix, vector or data frame of responses 
  to be  modeled as normal. Missing values (\code{NA}s) are allowed. If
  \code{y} is a data frame, any factors or ordered factors will be
  replaced by their internal codes, and a warning will be given.
  Alternatively, this first argument may be a \code{formula} as described
  below, or an object of class \code{"norm"} 
  resulting from a call to \code{emNorm} or \code{\link{mcmcNorm}};
  see DETAILS.} 

  \item{x}{a numeric matrix, vector or data frame of covariates to be
  used as predictors for \code{y}. Missing values (\code{NA}'s) are
  not allowed. If \code{x} is a matrix, it must have the same number
  of rows as \code{y}.  If \code{x} is a data frame, any factors or
  ordered factors are   replaced by their internal codes, and a
  warning is given. If \code{NULL}, it defaults to \code{x =
  rep(1,nrow(y))}, an intercept-only model.}

  \item{intercept}{if \code{TRUE}, then a column of \code{1}'s is
  appended to \code{x}.  Ignored if \code{x = NULL}.}

  \item{formula}{an object of class \code{"\link{formula}"} (or one
  that can be coerced to that class): a symbolic description of the
  model which is provided in lieu of \code{y} and \code{x}. The
  details of model specification are given 
  under DETAILS.}

  \item{data}{an optional data frame, list or environment (or object
  coercible by \code{\link{as.data.frame}} to a data frame) containing
  the variables in the model. If not found in \code{data}, the variables are
  taken from \code{environment(formula)}, typically the environment
  from which \code{logpostNorm} is called.} 

  \item{param}{assumed values for the model
  parameters. This must be a list with two named components,
  \code{beta} and \code{sigma}, which are numeric matrices with correct
  dimensions.  In most circumstances, the parameter
  values will be obtained from a run of \code{\link{emNorm}} or
  \code{mcmcNorm}; see DETAILS.}

  \item{prior}{should be \code{"uniform"}, \code{"jeffreys"},
  \code{"ridge"} or \code{"invwish"}. If \code{"ridge"} then
  \code{prior.df}  must be supplied.  If \code{"invwish"} then
  \code{prior.df} and \code{prior.sscp} must be 
  supplied.}

  \item{prior.df}{prior degrees of freedom for a ridge
  (\code{prior="ridge"}) or inverted Wishart (\code{prior="invwish"})
  prior.}

  \item{prior.sscp}{prior sums of squares and cross-products (SSCP)
  matrix for an inverted Wishart prior (\code{prior="invwish"}).}

  \item{\dots}{values to be passed to the methods.}
}


\details{ 
  The simplest way to call
  \code{logpostNorm} is to provide an object of class \code{"norm"} as its
  sole argument, where that object is the result of a call to
  \code{\link{emNorm}} or \code{\link{mcmcNorm}}. The parameter values
  stored in that object will then be passed to \code{logpostNorm} automatically.

  Alternatively, one may call \code{logpostNorm} by providing as the first
  argument \code{y}, a vector or matrix of data to be modeled as
  normal, and an optional vector or matrix of predictors \code{x}.
  Missing values \code{NA}
  are allowed in \code{y} but not in \code{x}. 

  A third way to call \code{logpostNorm} is to provide
  \code{formula},  a formula for a (typically
  multivariate) linear regression model in the manner expected by
  \code{\link{lm}}. A formula is given as \code{y ~ model}, where
  \code{y} is either a single numeric variable or a matrix of numeric
  variables bound together with the function \code{\link{cbind}}. The
  right-hand   side of the formula (everything to the right of \code{~}) is a
  linear predictor, a series of terms separated by operators \code{+},
  \code{:} or \code{*} to specify main effects and
  interactions. Factors are allowed on the right-hand side and will
  enter the model as contrasts among the \code{\link{levels}}. The
  intercept term \code{1} is included by default; to remove the
  intercept, use \code{-1}.
}

\value{ a numeric value reporting the observed-data log-posterior
density
}

\references{
Schafer, J.L. (1997) \emph{Analysis of Incomplete Multivariate
Data}. London: Chapman & Hall/CRC Press.
\cr

For more information about this function and other functions in
the \code{norm2} package, see \emph{User's Guide for \code{norm2}}
in the library subdirectory \code{doc}.
}

\author{Joe Schafer \email{Joseph.L.Schafer@census.gov} }

\seealso{ \code{\link{emNorm}}, 
\code{\link{mcmcNorm}},
\code{\link{loglikNorm}}
}

\examples{

## run EM for marijuana data with ridge prior and print the
## last value of the log-posterior density
data(marijuana)
emResult <- emNorm(marijuana, prior="ridge", prior.df=0.5)
print( emResult$logpost[ emResult$iter ] )

## compute the log-posterior density at the final estimate
## and compare it to the last value reported by emNorm
logpost.max <- logpostNorm(emResult)
print( logpost.max - emResult$logpost[ emResult$iter ] )

## The result from logpostNorm is slightly higher,
## because the last value reported by emNorm is the
## log-posterior at the BEGINNING of the last iteration

}

\keyword{ multivariate }
\keyword{ NA }
