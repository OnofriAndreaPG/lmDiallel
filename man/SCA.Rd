\name{SCA}
\alias{SCA}
\alias{SCA.G3}
\alias{SCA.GE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Specific Combining Ability
}
\description{
SCA effect to fit Hayman2, Griffing3 (SCA.G3), GE2 and GE3 (SCA.GE) models with \code{lm} function
}
\usage{
SCA(P1, P2, type = "fix", data)
SCA.G3(P1, P2,  type = "fix", data)
}
\arguments{
\item{P1}{\code{a variable for the first parent}}
  \item{P2}{\code{a variable for the second parent}}
  \item{type}{\code{a variable for model selection. May be "fix" (fixed model) or "random" (random model).}}
  \item{data}{\code{a 'data.frame' where to look for explanatory variables}}
}
\details{
a design matrix of all possible combinations between parentals with no selfs and no reciprocals
}
\references{
\cite{Onofri, A., Terzaroli, N. & Russi, L. Linear models for diallel crosses: a review with R functions. Theor Appl Genet (2020). https://doi.org/10.1007/s00122-020-03716-8}
}
\author{
Andrea Onofri, Niccolo' Terzaroli, Luigi Russi}
\value{
A design matrix for the SCA effect
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
