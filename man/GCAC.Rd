\name{GCAC}
\alias{GCAC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
General Combining Ability without considering the selfed parents
}
\description{
Design matrix for GCAC, useful to fit Gardner & Eberhart model 3 (GE3) with \code{lm} function
}
\usage{
GCAC(P1,P2,data)
}
\arguments{
 \item{P1}{\code{a variable for the first parent}}
  \item{P2}{\code{a variable for the second parent}}
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
GCAC(Par1,Par2, data=hayman54)
}
\keyword{ ~diallel }
\keyword{ ~genetic effects }
