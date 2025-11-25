\name{int.matrix}
\alias{int.matrix}
\alias{checkScheme}
\alias{emm.diallel}

\title{
Utilities for fitting diallel models.
}
\description{
 These functions are used internally by the package, but they can also called from the outside, to fit specific needs
}
\usage{
int.matrix(Xa, Xb)
checkScheme(P1, P2)
emm.diallel(obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Xa}{Incidence matrix of genetic effects}
  \item{Xb}{Incidence matrix for an external factor}
  \item{P1}{A vector with parentals}
  \item{P2}{A vector with parentals}
  \item{obj}{A glht object}
}

\details{
The function 'int.matrix()' produces the incidence matrix for the interaction between two main effects; 'Xa' and 'Xb' are two incidence matrices for two main effects.
The function 'checkScheme()' takes two vectors containing the codings for parentals (P1 and P2), retrieves the mating scheme and detects whether there are missing crosses.
The function 'emm.diallel()' is used with multi-environment diallel experiments to obtain the expected marginal means for genetic effects across environments.
}
\value{
The function 'int.matrix()' returns an incidence matrix. The function 'checkScheme()' returns a list, containing the main traits of the mating scheme. The function 'emm.diallel()' retrns a data.frame with the marginal means, standard errors and t-test statistics.
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
