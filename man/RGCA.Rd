\name{RGCA}
\alias{RGCA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reciprocal General Combining Ability
}
\description{
RGCA effect to fit Hayman1 & 2 models with \code{lm} function
}
\usage{
RGCA(P1,P2,data)
}
%- maybe also 'usage' for other objects documented here.
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
dMod<-lm(Ftime~ Block + GCA(Par1,Par2)
         + tSCA(Par1, Par2) + RGCA(Par1, Par2)
         + RSCA(Par1,Par2), data = hayman54)
anova(dMod)
}
\keyword{ ~diallel }
\keyword{ ~genetic effects }
