\name{int.matrix}
\alias{int.matrix}
\alias{checkScheme}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Utilities for fitting diallel models.
}
\description{
 These functions are used internally by the package.
}
\usage{
int.matrix(Xa, Xb)
checkScheme(P1, P2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Xa}{Incidence matrix of genetic effects}
  \item{Xb}{Incidence matrix for an external factor}
  \item{P1}{A vector with parentals}
  \item{P2}{A vector of parentals}
}

\details{
No more detail, so far.
}
\value{
The function 'int.matrix()' returns an incidence matrix. The function 'checkScheme()' returns a list, containing the main traits of the mating scheme.
}
\references{
\cite{Onofri, A., Terzaroli, N. & Russi, L. Linear models for diallel crosses: a review with R functions. Theor Appl Genet (2020). https://doi.org/10.1007/s00122-020-03716-8}
}
\author{
Andrea Onofri
}
\note{
No further notes
}


\examples{
data(griffing56)
head(griffing56)
checkScheme(griffing56$Par1, griffing56$Par2)
}
