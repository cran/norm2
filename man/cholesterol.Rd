\name{cholesterol}

\docType{data}

\alias{cholesterol}

\title{Cholesterol levels for heart-attack patients}

\description{
This dataset reports cholesterol levels for 28 patients treated at a
Pennsylvania medical center. 
All patients have cholesterol recorded on day 2 and day 4
after attack, but
some have missing values on day 14.
These data were analyzed by Schafer (1997, Chap. 5).
}

\usage{data(cholesterol)}

\format{a data frame with 28 rows and 3 variables:

\describe{    

    \item{\code{Y1}}{cholesterol 2 days after heart attack.}

    \item{\code{Y2}}{cholesterol 4 days after heart attack.}

    \item{\code{Y3}}{cholesterol 14 days after heart attack.}

}
}

\source{
Ryan, B.F. and Joiner, B.L. (1994) \emph{Minitab Handbook}
(Third edition). Belmont, CA: Wadsworth.
}

\references{
Schafer, J.L. (1997) \emph{Analysis of Incomplete Multivariate
Data}. London: Chapman & Hall/CRC Press.
\cr

For example analyses of this dataset using functions in
the \code{norm2} package, see the manual \emph{NORM Package for R, Version 2}
in the library subdirectory \code{doc}.
}

\keyword{ datasets }
