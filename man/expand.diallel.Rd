\name{expand.diallel}
\alias{expand.diallel}
\title{
Create a Data Frame from All Combinations of Parentals
}
\description{
This is a modification of the 'expand.grid()' function working specifically with diallel experiments. It creates a data frame from all combinations of the supplied vector of parents, depending on the mating scheme.
}
\usage{
expand.diallel(pars, mating = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{pars}{\code{pars} a vector of parentals}
\item{mating}{The type of mating scheme. 1: full diallel experiment; 2: no reciprocals; 3: no selfs; 4: no reciprocals and no selfs}}

\value{
returns a data.frame object}
\references{
\cite{Onofri, A., Terzaroli, N. & Russi, L. Linear models for diallel crosses: a review with R functions. Theor Appl Genet (2020). https://doi.org/10.1007/s00122-020-03716-8}
}
\author{
Andrea Onofri \email{(andrea.onofri@unipg.it)}, Niccolo' Terzaroli \email{(n.terzaroli@gmail.com)}, Luigi Russi \email{(luigi.russi@unipg.it)}
}
\examples{
pars <- LETTERS[1:4]
expand.diallel(pars, mating = 3)
}
