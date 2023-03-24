\name{griffing56}
\alias{griffing56}
\docType{data}
\title{
Griffing's dataset for diallel analysis
}
\description{
Data for a diallel in maize, with no selfs and no selfed parents. Data are the means of several replicates.
}
\usage{
data("griffing56")
}
\format{
  A data.frame with 36 observations on the following 5 variables
  \describe{
    \item{\code{Par1}}{male parent, a factor with 8 levels}
    \item{\code{Par2}}{female parent, a factor with 8 levels}
    \item{\code{Yield}}{Maize Yield}
    \item{\code{Cob}}{Cob weight}
    \item{\code{Shelled}}{Shelled corn weight}
  }
}
\author{
Andrea Onofri, Niccolo' Terzaroli, Luigi Russi}
\source{
\cite{Griffing, B., 1956. Concept of general and specific combining ability in relation to diallel crossing systems. Australian Journal of Biological Science 9, 463--493.}
}
\references{
\cite{Onofri, A., Terzaroli, N. & Russi, L. Linear models for diallel crosses: a review with R functions. Theor Appl Genet (2020). https://doi.org/10.1007/s00122-020-03716-8}
}
\examples{
data("griffing56")
dMod2 <- lm.diallel(Yield ~ Par1 + Par2,
                    data = griffing56, fct = "GRIFFING4")
anova(dMod2, MSE = 21.05, dfr = 2558)
summary(dMod2, MSE = 21.05, dfr = 2558)
}
\keyword{datasets}
