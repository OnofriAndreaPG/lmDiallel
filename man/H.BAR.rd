\name{H.BAR}
\alias{H.BAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Average heterosis effect
}
\description{
H.BAR effect to fit GE2 and GE3 models with \code{lm} function
}
\usage{
H.BAR(P1, P2, data)
}
\arguments{
  \item{P1}{\code{a variable for the first parent}}
  \item{P2}{\code{a variable for the second parent}}
  \item{data}{\code{a 'data.frame' where to look for explanatory variables}}
}
\references{
\cite{Onofri, A., Terzaroli, N. & Russi, L. Linear models for diallel crosses: a review with R functions. Theor Appl Genet (2020). https://doi.org/10.1007/s00122-020-03716-8}
}
\author{
Andrea Onofri \email{(andrea.onofri@unipg.it)}, Niccolo' Terzaroli \email{(n.terzaroli@gmail.com)}, Luigi Russi \email{(luigi.russi@unipg.it)}
}
\examples{
data("zhang05")
dMod <- lm(Yield ~ Env/Block + H.BAR(Par1, Par2) + VEi(Par1, Par2) +
                   Hi(Par1, Par2) + SCA(Par1, Par2) +
                   H.BAR(Par1, Par2):Env + VEi(Par1, Par2):Env +
                   Hi(Par1, Par2):Env + SCA(Par1, Par2):Env, data = zhang05)
anova(dMod)
}
\keyword{ ~diallel }
\keyword{ ~genetic effects }
