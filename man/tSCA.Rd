\name{tSCA}
\alias{tSCA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Total Specific Combining Ability
}
\description{
Total SCA to fit Hayman1, Griffing1 and Griffing2 models with \code{lm} function
}
\usage{
tSCA(P1,P2,data)
}
\arguments{
 \item{P1}{\code{a variable for the first parent}}
  \item{P2}{\code{a variable for the second parent}}
  \item{data}{\code{a 'data.frame' where to look for explanatory variables}}
}
\details{
a design matrix of all possible combinations between parentals with selfs but no reciprocals
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
