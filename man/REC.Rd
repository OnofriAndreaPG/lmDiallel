\name{REC}
\alias{REC}
\alias{REC.G3}
\title{
Reciprocal Effect not parted into RGCA and RSCA
}
\description{
Build incidence matrix to fit reciprocal effects in Griffing's model 1, 2, 4 (REC) and 3 (REC.G3) with \code{lm} function
}
\usage{
REC(P1, P2, type = "fix", data)
}
\arguments{
\item{P1}{\code{a variable for the first parent}}
  \item{P2}{\code{a variable for the second parent}}
  \item{type}{\code{a variable for model selection. May be "fix" (fixed model) or "random" (random model).}}
  \item{data}{\code{a 'data.frame' where to look for explanatory variables}}}

\details{
a design matrix of n-1 parentals
}
\references{
\cite{Onofri, A., Terzaroli, N. & Russi, L. Linear models for diallel crosses: a review with R functions. Theor Appl Genet (2020). https://doi.org/10.1007/s00122-020-03716-8}
}
\author{
Andrea Onofri, Niccolo' Terzaroli, Luigi Russi}
\value{
A design matrix for the reciprocal (REC) effect
}
\examples{
data("hayman54")
dMod<-lm(Ftime~ Block + GCA(Par1,Par2)
         +tSCA(Par1, Par2)+REC(Par1, Par2)
         , data = hayman54)
anova(dMod)
}
\keyword{ ~diallel }
\keyword{ ~genetic effects }
