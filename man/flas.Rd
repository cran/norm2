\name{flas}

\docType{data}

\alias{flas}

\title{Foreign Language Attitude Scale}

\description{
This dataset comes from a study to evaluate the
reliability of the Foreign Lanuage Attitude Scale, an instrument for
predicting success in the study or foreign languages (Raymond and
Roberts, 1983). The
questionnaire was given to 279 students enrolled in four different
language courses (French, German, Spanish, Russian) at 
Penn State University.  This dataset includes FLAS score,
final grade in the course, and other variables for predicting
achievement.  These data were analyzed by Schafer (1997,
Chap. 6).
}

\usage{data(flas)}

\format{a data frame with 279 rows and 14 variables:

\describe{    

    \item{\code{LAN2}}{1=Spanish, 0=other.}

    \item{\code{LAN3}}{1=German, 0=other.}

    \item{\code{LAN4}}{1=Russian, 0=other.}

    \item{\code{AGE}}{age group (1=less than 20, 2=20+).}

    \item{\code{PRI}}{number of prior foreign language courses
    (1=none, 2=1-2, 3=3+).}

    \item{\code{SEX}}{0=male, 1=female}

    \item{\code{MLAT}}{Modern Language Aptitude Test}

    \item{\code{FLAS}}{Foreign Language Attitude Scale}

    \item{\code{SATV}}{Scholastic Aptitude Test, verbal score}

    \item{\code{SATM}}{Scholastic Aptitude Test, math score}

    \item{\code{ENG}}{score on Penn State English placement exam}

    \item{\code{HGPA}}{high school grade point average}

    \item{\code{CGPA}}{current college grade point average}

    \item{\code{GRD}}{final grade in foreign language course (1=B or
    lower, 2=A)}
}
}

\source{
Schafer, J.L. (1997) \emph{Analysis of Incomplete Multivariate
Data}. London: Chapman & Hall/CRC Press.
}

\references{

Raymond, M.R. and Roberts, D.M. (1983) Development and
validation of a foreign language attitude scale. \emph{Educational and
Psychological Measurement}, 43, 1239-1246.
\cr

For example analyses of this dataset using functions in
the \code{norm2} package, see the manual \emph{NORM Package for R, Version 2}
in the library subdirectory \code{doc}.
}

\keyword{ datasets }
