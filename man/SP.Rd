\name{SP}
\alias{SP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Selfed Parents effect
}
\description{
SP effect to fit GE3 model with \code{lm} function
}
\usage{
SP(P1, P2, type = "fix", data)
}
\arguments{
  \item{P1}{\code{a variable for the first parent}}
  \item{P2}{\code{a variable for the second parent}}
  \item{type}{\code{a variable for model selection. May be "fix" (fixed model) or "random" (random model).}}
  \item{data}{\code{a 'data.frame' where to look for explanatory variables}}
  }
\details{
a design matrix of n-1 parentals
}
\references{
\cite{Onofri, A., Terzaroli, N. & Russi, L. Linear models for diallel crosses: a review with R functions. Theor Appl Genet (2020). https://doi.org/10.1007/s00122-020-03716-8}
}
\author{
Andrea Onofri \email{(andrea.onofri@unipg.it)}, Niccolo' Terzaroli \email{(n.terzaroli@gmail.com)}, Luigi Russi \email{(luigi.russi@unipg.it)}
}
\examples{
data("hayman54")
SP(Par1,Par2, data=hayman54)
}

\keyword{ ~diallel }
\keyword{ ~genetic effects }
