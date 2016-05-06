\name{impNorm}

\docType{methods}

\alias{impNorm}
\alias{impNorm.default}
\alias{impNorm.formula}
\alias{impNorm.norm}

\title{ Imputation and prediction for incomplete multivariate normal data}

\description{
Simulates or predicts missing values from their predictive
distribution given the observed data 
under a normal model with fixed parameters.
}


\usage{
% the generic function
impNorm(obj, \dots)

% the default method
\method{impNorm}{default}(obj, x = NULL, intercept = TRUE, param,
   seeds = NULL, method = "random", \ldots)

% method for class formula
\method{impNorm}{formula}(formula, data, param,
   seeds = NULL, method = "random", \ldots)

% method for class norm
\method{impNorm}{norm}(obj, param = obj$param, seeds = NULL, 
   method = "random", \ldots)

}


\arguments{

  \item{obj}{an object used to select a method.  It may be \code{y},
  a numeric matrix, vector or data frame of responses 
  to be  modeled as normal. Missing values (\code{NA}s) are allowed. If 
  \code{y} is a data frame, any factors or ordered factors will be
  replaced by their internal codes, and a warning will be given.
  Alternatively, this first argument may be an object of class \code{"norm"} 
  resulting from a call to \code{emNorm} or
  \code{\link{mcmcNorm}}; see DETAILS.} 

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
  from which \code{impNorm} is called.} 

  \item{param}{assumed values for the model
  parameters. This must be a list with two named components,
  \code{beta} and \code{sigma}, which are numeric matrices with correct
  dimensions.  In most circumstances, the parameter
  values will be obtained from a run of \code{\link{emNorm}} or
  \code{mcmcNorm}; see DETAILS.}

  \item{seeds}{two integers to initialize the random number
  generator; see DETAILS.}

  \item{method}{if \code{"random"}, the missing values in each row of \code{y}
  will be simulated from their joint predictive distribution given
  \code{x} and the observed values in \code{y}. If \code{"predict"},
  missing values will be 
  replaced by regression predictions given the observed values. See
  DETAILS.} 

  \item{\dots}{values to be passed to the methods.}
}


\details{ 
  This function is used primarily in conjunction with
  \code{\link{mcmcNorm}} to draw multiple imputations by the
  multiple-chain method.  In those instances, the simplest way to call
  \code{impNorm} is to provide an object of class \code{"norm"} as its
  first argument, where that object is the result of a call to
  \code{mcmcNorm}. The parameter values stored in that object will
  then be passed to \code{impNorm} automatically.

  Alternatively, one may call \code{impNorm} by providing as the first
  argument \code{y}, a vector or matrix of data to be modeled as
  normal, and an optional vector or matrix of predictors \code{x}.
  Missing values \code{NA}
  are allowed in \code{y} but not in \code{x}. 

  A third way to call \code{impNorm} is to provide
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

  \code{norm2} functions use their own internal random number generator which
  is seeded by two integers, for example, \code{seeds=c(123,456)},
  which allows results to be reproduced in the future. If
  \code{seeds=NULL} then 
  the function will seed itself with two random
  integers from R.  Therefore, results  can also be  made reproducible by
  calling \code{\link{set.seed}} beforehand and taking \code{seeds=NULL}.

}

\value{ a data matrix resembling the original data \code{y}, but
  with \code{NA}'s replaced by simulated values or predictions.
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

\seealso{\code{\link{emNorm}}, \code{\link{mcmcNorm}}
}

\examples{

## run EM for marijuana data with ridge prior
data(marijuana)
emResult <- emNorm(marijuana, prior="ridge", prior.df=0.5)

## generate 25 multiple imputations by running 25 chains
## of 100 iterations each, starting each chain at the 
## posterior mode
set.seed(456)
imp.list <- as.list(NULL)
for(m in 1:25){
   mcmcResult <- mcmcNorm(emResult, iter=100)
   imp.list[[m]] <- impNorm(mcmcResult)}

}

\keyword{ multivariate }
\keyword{ NA }
