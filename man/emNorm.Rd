\name{emNorm}

\alias{emNorm}
\alias{emNorm.default}
\alias{emNorm.formula}
\alias{emNorm.norm}

\title{ EM algorithm for incomplete multivariate normal data}

\description{
Computes maximum likelihood estimates and posterior modes from
incomplete multivariate data under a normal model.
}


\usage{
% the generic function
emNorm(obj, \dots)

% the default method
\method{emNorm}{default}(obj, x = NULL, intercept = TRUE,
   iter.max = 1000, criterion = NULL, estimate.worst = TRUE, 
   prior = "uniform", prior.df = NULL, prior.sscp = NULL,
   starting.values = NULL, \ldots)

% method for class formula
\method{emNorm}{formula}(formula, data, iter.max = 1000, 
   criterion = NULL, estimate.worst = TRUE, prior = "uniform",
   prior.df = NULL, prior.sscp = NULL, starting.values = NULL, \ldots)

% method for class norm
\method{emNorm}{norm}(obj, iter.max = 1000, 
   criterion = obj$criterion, estimate.worst = obj$estimate.worst,
   prior = obj$prior, prior.df = obj$prior.df,
   prior.sscp = obj$prior.sscp, starting.values = obj$param, \ldots)

}


\arguments{

  \item{obj}{an object used to select a method.  It may be \code{y},
  a numeric matrix, vector or data frame containing response variables 
  to be  modeled as multivariate normal. Missing values (\code{NA}s)
  are allowed. If \code{y}  
  is a data frame, any factors or ordered factors will be
  replaced by their internal codes, and a warning will be given.
  Alternatively, this first argument may be a \code{formula} as described
  below, or an object of class \code{"norm"} 
  resulting from a call to \code{emNorm} or \code{\link{mcmcNorm}};
  see DETAILS.} 

  \item{x}{a numeric matrix, vector or data frame of covariates to be
  used as model predictors. Missing values (\code{NA}'s) are
  not allowed. If \code{x} is a matrix, it must have the same number
  of rows as \code{y}.  If \code{x} is a data frame, any factors or
  ordered factors will be replaced by their internal codes, and a
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
  from which \code{emNorm} is called.} 

  \item{iter.max}{maximum number of iterations to be performed.  Each
  iteration consists of an Expectation or E-step followed by a
  Maximization or M-step. The procedure halts if it has not converged
  by this many iterations.}

  \item{criterion}{convergence criterion.  The procedure halts if the
  maximum relative change in all parameters from one iteration to the
  next falls below this value. If \code{NULL}, then the default
  criterion of \code{1e-05} is used.}

  \item{estimate.worst}{if \code{TRUE}, then upon convergence of the EM
  algorithm, a procedure is attempted to numerically estimate the worst
  fraction of missing information and the worst linear function of the
  parameters; see DETAILS.} 

  \item{prior}{should be \code{"uniform"}, \code{"jeffreys"},
  \code{"ridge"} or \code{"invwish"}. If \code{"ridge"} then
  \code{prior.df}  must be supplied.  If \code{"invwish"} then
  \code{prior.df} and \code{prior.sscp} must be 
  supplied. For more information, see DETAILS.}

  \item{prior.df}{prior degrees of freedom for a ridge
  (\code{prior="ridge"}) or inverted Wishart (\code{prior="invwish"})
  prior.}

  \item{prior.sscp}{prior sums of squares and cross-products (SSCP)
  matrix for an inverted Wishart prior (\code{prior="invwish"}).}

  \item{starting.values}{optional starting values for the model
  parameters. This must be a list with two named components,
  \code{beta} and \code{sigma}, which are numeric matrices with correct
  dimensions; see DETAILS.}

  \item{\dots}{values to be passed to the methods.}
}


\details{ 

  There are three different ways to specify the data and model when
  calling \code{emNorm}: 
  \itemize{
     \item by directly supplying as the initial argument a matrix of
     numeric response variables \code{y}, along with an optional 
     matrix of predictor variables \code{x};
     \item by supplying a model specification
     \code{formula}, along with an optional data frame \code{data}; or
     \item by supplying an object of class
     \code{"norm"}, which was produced by an earlier call to
     \code{emNorm} or \code{\link{mcmcNorm}}.
   }

  In the first case, the matrix \code{y} is assumed to have a
  multivariate normal 
  linear regression on \code{x} with coefficients \code{beta} and
  covariance matrix \code{sigma}, where
  \code{dim(beta)=c(ncol(x),ncol(y))} and
  \code{dim(sigma)=c(ncol(y),ncol(y))}. Missing values \code{NA}
  are allowed in \code{y} but not in \code{x}.

  In the second case, \code{formula} is a formula for a (typically
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

  In the third case, the initial argument to \code{emNorm} is an
  object of class 
  \code{"norm"} returned by a previous call to \code{emNorm}
  or \code{\link{mcmcNorm}}. The value of the parameters
  carried in this object (the estimates from the last iteration of
  EM or the simulated values from the last iteration of MCMC) will be
  used as the starting values.  This
  can be useful, for example, if a previous run of EM did not converge
  by \code{max.iter} iterations.  Supplying the result from that
  EM run as the sole argument to \code{emNorm} allows the algorithm to
  continue from where it was halted.

  If \code{prior="uniform"} (the default), the EM algorithm computes a
  maximum-likelihood estimate for the parameters \code{beta} and
  \code{sigma}; otherwise it computes a posterior mode.

  If \code{prior="invwish"} then an inverted Wishart prior distribution
  is applied to \code{sigma} with hyperparameters \code{prior.df} (a
  scalar) and \code{prior.sscp} (a symmetric, positive definite matrix
  of the same dimension as \code{sigma}).  Using the device of imaginary
  results, we can interpret \code{prior.sscp/prior.df} as a prior guess
  for \code{sigma}, and \code{prior.df} as the prior degrees of
  freedom on which this guess is based.

  The usual noninformative prior for the normal regression model
  (\code{prior="jeffreys"}) is equivalent to the inverted 
  Wishart density with \code{prior.df} equal to 0 and
  \code{prior.sscp} equal to a matrix of 0's.

  The ridge prior (\code{prior="ridge"}) is a special case of the
  inverted Wishart (Schafer, 1997).  The prior
  guess for \code{sigma} is a diagonal matrix with diagonal elements
  estimated by regressing the observed values in each column of
  \code{y} on the corresponding rows of \code{x}. When
  \code{prior="ridge"}, the user must supply a value for
  \code{prior.df}, which 
  determines how strongly the estimated correlations are smoothed
  toward zero.

  If \code{estimate.worst}, then upon convergence, a procedure is run
  to numerically estimate the
  worst fraction of missing information and the worst linear function
  of the parameters. The worst fraction of missing information is
  closely related to EM's convergence rate. Values near one
  correspond to slow convergence, and values near zero indicate fast
  convergence. If there are no missing values in the response variables,
  the worst fraction of missing information is exactly zero, and EM
  converges after one step from any starting values. The worst linear
  function is a linear combination of the parameters (elements of
  \code{beta} and \code{sigma}) for which the rate of missing information
  is highest.

  For details of the EM algorithm, see the manual distributed
  with the \code{norm2} package in the library subdirectory \code{doc}.
}

\value{
  a list whose
  \code{class} attribute has been set to \code{"norm"}.
  This object may be 
  passed as the first argument in subsequent calls to \code{emNorm}, 
  \code{\link{mcmcNorm}}, \code{\link{impNorm}},
  \code{\link{loglikNorm}} or \code{\link{logpostNorm}}. 
  The object also carries the original data and specifies the prior
  distribution, so that these do not need to be provided again.
  \cr

  To see a summary of
  this object, use the generic function \code{summary},
  which passes the object to \code{\link{summary.norm}}.
  \cr

  Components of the list may also be directly accessed
  and examined by the user.  Components  which may be of interest 
  include:

  \item{iter}{number of EM iterations actually performed.}

  \item{rel.diff}{maximum relative difference between the parameters
  the last two iterations.}

  \item{converged}{logical value indicating whether the algorithm
  converged by \code{iter} iterations. Will be \code{TRUE} if
  \code{rel.diff<=criterion}.}

  \item{loglik}{a numeric vector of length \code{iter} reporting the
  logarithm of the observed-data likelihood function at the start of
  each iteration.  If \code{prior="uniform"} then the loglikelihood
  values will be non-decreasing.}

  \item{logpost}{a numeric vector of length \code{iter} reporting the
  logarithm of the observed-data posterior density function at the start of
  each iteration.  The log-posterior density values will be
  non-decreasing. If \code{prior="uniform"} then the log-posterior
  density and loglikelihood will be identical.} 

  \item{param}{a list with elements \code{beta} and \code{sigma}
  containing the estimated parameters after the final iteration of
  EM. This may be supplied as starting values to \code{emNorm} or
  \code{\link{mcmcNorm}}, or as an argument to \code{\link{impNorm}},
  \code{\link{loglikNorm}} or \code{\link{logpostNorm}}
  .}

  \item{miss.patt}{logical matrix with \code{ncol(y)} columns
  reporting the missingness patterns seen in \code{y}.  Each row of
  \code{miss.patt} corresponds to a distinct missingness pattern, with
  \code{TRUE} indicating that the \code{y}-variable is missing and
  \code{FALSE} indicating that the \code{y}-variable is observed.} 

  \item{miss.patt.freq}{integer vector of length
  \code{nrow(miss.patt)} indicating, for each missingness pattern, the
  number of cases or rows of \code{y} having that pattern.}

  \item{which.patt}{integer vector of length \code{nrow(y)} indicating
  the missingness pattern for each
  row of \code{y}. Thus \code{is.na( y[i,] )} is the same thing as
  \code{miss.patt[ which.patt[i], ]}.

 }

}

\references{
Schafer, J.L. (1997) \emph{Analysis of Incomplete Multivariate
Data}. London: Chapman & Hall/CRC Press.
\cr

For more information about this function and other functions in
\code{\link{norm2}}, see the manual \emph{NORM Package for
R, Version 2} in the library subdirectory \code{doc}.
}

\author{Joe Schafer \email{Joseph.L.Schafer@census.gov} }

\seealso{
\code{\link{mcmcNorm}}, 
\code{\link{summary.norm}},
\code{\link{impNorm}},
\code{\link{loglikNorm}},
\code{\link{logpostNorm}}
}


\examples{

## run EM for marijuana data with strict convergence criterion
data(marijuana)
result <- emNorm(marijuana, criterion=1e-06)

## re-run with ridge prior and examine results
result <- emNorm(marijuana, prior="ridge", prior.df=0.5)
summary(result)
}

\keyword{ multivariate }
\keyword{ NA }
