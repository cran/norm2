\name{mcmcNorm}

\docType{methods}

\alias{mcmcNorm}
\alias{mcmcNorm.default}
\alias{mcmcNorm.formula}
\alias{mcmcNorm.norm}

\title{ MCMC algorithm for incomplete multivariate normal data}

\description{
Simulates parameters and missing values from a joint posterior
distribution under a normal model using Markov chain Monte Carlo.
}


\usage{
% the generic function
mcmcNorm(obj, \dots)

% the default method
\method{mcmcNorm}{default}(obj, x = NULL, intercept = TRUE,
   starting.values, iter = 1000, multicycle = NULL,
   seeds = NULL, prior = "uniform",
   prior.df = NULL, prior.sscp = NULL, save.all.series = TRUE,
   save.worst.series = FALSE, worst.linear.coef = NULL,
   impute.every = NULL, \ldots)

% method for class formula
\method{mcmcNorm}{formula}(formula, data, starting.values, 
   iter = 1000, multicycle = NULL, seeds = NULL, prior = "uniform", 
   prior.df = NULL, prior.sscp = NULL, save.all.series = TRUE, 
   save.worst.series = FALSE, worst.linear.coef = NULL,
   impute.every=NULL, \ldots)

% method for class norm
\method{mcmcNorm}{norm}(obj, starting.values = obj$param,
   iter = 1000, multicycle = obj$multicycle, 
   seeds = NULL, prior = obj$prior, prior.df = obj$prior.df, 
   prior.sscp = obj$prior.sscp,
   save.all.series = !(obj$method=="MCMC" & is.null( obj$series.beta )), 
   save.worst.series = !is.null( obj$worst.linear.coef ),
   worst.linear.coef = obj$worst.linear.coef,
   impute.every = obj$impute.every, \ldots)
}


\arguments{

  \item{obj}{an object used to select a method.  It may be \code{y},
  a numeric matrix, vector or data frame containing response variables
  to be  modeled as normal. Missing values (\code{NA}s)
  are allowed. If \code{y}  
  is a data frame, any factors or ordered factors will be
  replaced by their internal codes, and a warning will be given.
  Alternatively, this first argument may be a \code{formula} as described
  below, or an object of class \code{"norm"} 
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
  from which \code{mcmcNorm} is called.} 

  \item{starting.values}{starting values for the model
  parameters. This must be a list with two named components,
  \code{beta} and \code{sigma}, which are numeric matrices with correct
  dimensions.  In most circumstances, the starting
  values will be obtained from a prior run of \code{\link{emNorm}} or
  \code{mcmcNorm}; see DETAILS.}

  \item{iter}{number of iterations to be performed.  By default, each
  iteration consists of one Imputation or I-step followed by 
  one Posterior or P-step, but this can be changed by
  \code{multicycle}.} 

  \item{multicycle}{number of cycles per iteration, with
  \code{NULL} equivalent to \code{multicycle=1}. 
  Specifying   
  \code{multicycle=}\emph{k} for some \emph{k}>1 instructs
  \code{mcmcNorm} to perform the I-step and P-step cycle \code{k}
  times within each iteration; see DETAILS.}

  \item{seeds}{two integers to initialize the random number
  generator; see DETAILS.}

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

  \item{save.all.series}{if \code{TRUE}, then the simulated values of all
  parameters at all iterations will be saved.}

  \item{save.worst.series}{if \code{TRUE}, then the simulated values
  of the worst linear function of the parameters will be saved. Under
  ordinary circumstances, this function will have been estimated by 
  \code{\link{emNorm}} after the EM algorithm converged.}

  \item{worst.linear.coef}{vector or coefficients that define the worst
  linear function of the parameters. Under ordinary circumstances,
  these are provided automatically in the result from \code{\link{emNorm}}.}

  \item{impute.every}{how many iterations to perform between
  imputations? If \code{impute.every=}\emph{k}, then the simulated
  values for the missing data after every \emph{k} iterations will be
  saved, resulting in \code{floor(iter/impute.every)} multiple
  imputations. If \code{NULL}, then no imputations will be saved.}

  \item{\dots}{values to be passed to the methods.}
}


\details{ 

  There are three different ways to specify the data and model when
  calling \code{mcmcNorm}: 
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

  In the third case, the initial argument to \code{mcmcNorm} is an
  object of class 
  \code{"norm"} returned by a previous call to \code{emNorm}
  or \code{\link{mcmcNorm}}. The value of the parameters
  carried in this object (the estimates from the last iteration of
  EM or the simulated values from the last iteration of MCMC) will be
  used as the starting values.  

  The matrix \code{y} is assumed to have a multivariate normal
  linear regression on \code{x} with coefficients \code{beta} and
  covariance matrix \code{sigma}, where
  \code{dim(beta)=c(ncol(x),ncol(y))} and
  \code{dim(sigma)=c(ncol(y),ncol(y))}.  

  Starting values for the parameters must be provided. In most cases
  these will be the result of a previous call to \code{emNorm} or
  \code{mcmcNorm}. If the starting
  values are close to the mode (i.e., if they are the result of an EM
  run that converged) then the worst linear function of the
  parameters will  be saved at each iteration. If the starting
  values are the result of a previous run of MCMC, then the new
  run will be a continuation of the same Markov chain.

  If \code{multicycle=}\emph{k} for some \emph{k}>1,
  then the length of the saved parameter
  series will be reduced by a factor of \emph{k}, and the serial
  correlation in the series will also be reduced. This option is
  useful in large problems with many parameters  and in slowly
  converging problems for which many iterations are needed.

  \code{norm2} functions use their own internal random number generator which
  is seeded by two integers, for example, \code{seeds=c(123,456)},
  which allows results to be reproduced in the future. If
  \code{seeds=NULL} then 
  the function will seed itself with two random
  integers from R.  Therefore, results  can also be  made reproducible by
  calling \code{\link{set.seed}} beforehand and taking \code{seeds=NULL}.

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

  If the first argument to \code{mcmcNorm} is an object of class
  \code{"norm"}, then the parameter values stored in that object will
  automatically be used as starting values.

  For details of the MCMC algorithm, see the manual distributed
  with the NORM package in the subdirectory \code{doc}.
}

\value{

  a list whose
  \code{class} attribute has been set to \code{"norm"}.
  This object may be 
  passed as the first argument in subsequent calls to \code{emNorm}, 
  \code{\link{mcmcNorm}}, \code{\link{impNorm}},
  \code{\link{loglikNorm}} or \code{\link{logpostNorm}}. The 
  object also carries the original data and specifies the prior
  distribution, so that these do not need to be provided again.
  \cr

  To see a summary of
  this object, use the generic function \code{summary},
  which passes the object to \code{\link{summary.norm}}.
  \cr

  Components of the list may also be directly accessed
  and examined by the user.  Components  which may be of interest 
  include:

  \item{iter}{number of MCMC iterations performed.}

  \item{param}{a list with elements \code{beta} and \code{sigma}
  containing the estimated parameters after the final iteration of
  MCMC. This may be supplied as starting values to \code{emNorm} or
  \code{\link{mcmcNorm}}, or as an argument to \code{\link{impNorm}},
  \code{\link{loglikNorm}} or \code{\link{logpostNorm}}.}

  \item{loglik}{a numeric vector of length \code{iter} reporting the
  logarithm of the observed-data likelihood function at the start of
  each iteration.}

  \item{logpost}{a numeric vector of length \code{iter} reporting the
  logarithm of the observed-data posterior density function at the start of
  each iteration.}

  \item{series.worst}{a time-series object (class \code{"ts"}) which
  contains the simulated values of the worst linear function of the
  parameters from all iterations.  This will be present if the
  starting values provided to \code{mcmcNorm} were close enough to the
  mode to provide a reliable estimate of the worst linear function.
  The dependence in this series tends to be higher than for
  other parameters, so examining the dependence by plotting the series
  with \code{\link{plot}} or its autocorrelation function with
  \code{\link{acf}} may help the user to judge how quickly the Markov
  chain achieves stationarity. For the definition of the worst linear
  function, see the manual accompanying the NORM package
  in the subdirectory \code{doc}.}

  \item{series.beta}{a multivariate time-series object (class
  \code{"ts"}) which 
  contains the simulated values of the coefficients \code{beta} 
  from all iterations.  This will present if \code{save.all.series=TRUE}.}

  \item{series.sigma}{a multivariate time-series object (class
  \code{"ts"}) which 
  contains the simulated values of the variances and
  covariances (elements of the lower triangle of \code{sigma})
  from all iterations.  This will be present if \code{save.all.series=TRUE}.}

  \item{imp.list}{a list containing the multiple imputations. Each
  component of this list is a data matrix resembling \code{y}, but
  with \code{NA}'s replaced by imputed values. The
  length of the list depends on the values of \code{iter} and
  \code{impute.every}.}

  \item{miss.patt}{logical matrix with \code{ncol(y)} columns
  reporting the missingness patterns seen in \code{y}.  Each row of
  \code{miss.patt} corresponds to a distinct missingness pattern, with
  \code{TRUE} indicating that the variable is missing and
  \code{FALSE} indicating that the variable is observed.} 

  \item{miss.patt.freq}{integer vector of length
  \code{nrow(miss.patt)} indicating, for each missingness pattern, the
  number of cases or rows of \code{y} having that pattern.}

  \item{which.patt}{integer vector of length \code{nrow(y)} indicating
  the missingness pattern for each
  row of \code{y}. Thus \code{is.na( y[i,] )} is the same thing as
  \code{miss.patt[ which.patt[i], ]}.}

}

\references{
Schafer, J.L. (1997) \emph{Analysis of Incomplete Multivariate
Data}. London: Chapman & Hall/CRC Press.
\cr

For more information about this function and other functions in
the NORM package, see \emph{User's Guide for \code{norm2}}
in the library subdirectory \code{doc}.
}

\author{Joe Schafer \email{Joseph.L.Schafer@census.gov} }

\seealso{\code{\link{emNorm}}, \code{\link{summary.norm}},
\code{\link{impNorm}},
\code{\link{loglikNorm}},
\code{\link{logpostNorm}}
}


\examples{

## run EM for marijuana data with ridge prior
data(marijuana)
emResult <- emNorm(marijuana, prior="ridge", prior.df=0.5)

## run MCMC for 5,000 iterations starting from the 
## posterior mode using the same prior
mcmcResult <- mcmcNorm(emResult, iter=5000)

## summarize and plot worst linear function
summary(mcmcResult)
plot(mcmcResult$series.worst)
acf(mcmcResult$series.worst, lag.max=50)

## generate 25 multiple imputations, taking 
## 100 steps between imputations, and look st
## the first imputed dataset
mcmcResult <- mcmcNorm(emResult, iter=2500, impute.every=100)
mcmcResult$imp.list[[1]]

}

\keyword{ multivariate }
\keyword{ NA }
