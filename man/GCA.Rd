\name{GCA}
\alias{GCA}
\alias{GCAmis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
General Combining Ability effect
}
\description{
GCA effect to fit Hayman1 & 2 and Griffing 1 & 2 models with \code{lm} function
}
\usage{
GCA(P1,P2,type = "fix", data)
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
\examples{
data("hayman54")
contrasts(hayman54$Block) <- "contr.sum"
dMod<-lm(Ftime ~ Block + GCA(Par1,Par2)
         + tSCA(Par1, Par2) + RGCA(Par1, Par2)
         + RSCA(Par1,Par2), data = hayman54)
anova(dMod)
}
\value{
A design matrix for the GCA effect
}
\keyword{ ~diallel }
\keyword{ ~genetic effects }
