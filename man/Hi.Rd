\name{Hi}
\alias{Hi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Average heterosis effect
}
\description{
H.i effect to fit GE2 Model with \code{lm} function
}
\usage{
Hi(P1, P2, type = "fix", data)
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
Andrea Onofri, Niccolo' Terzaroli , Luigi Russi}
\value{
A design matrix for the Hi effect
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
