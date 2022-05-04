\name{MDD}
\alias{MDD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mean Dominance Deviation effect
}
\description{
It relates to the difference between the average yield of selfed parents and the average yield of crosses. DD effect to fit Hayman2 model with \code{lm} function
}
\usage{
MDD(P1, P2, type = "fix", data)
}
\arguments{
  \item{P1}{\code{a variable for the first parent}}
  \item{P2}{\code{a variable for the second parent}}
  \item{type}{\code{a variable for the model}}
  \item{data}{\code{a 'data.frame' where to look for explanatory variables}}
}
\references{
\cite{Onofri, A., Terzaroli, N. & Russi, L. Linear models for diallel crosses: a review with R functions. Theor Appl Genet (2020). https://doi.org/10.1007/s00122-020-03716-8}
}
\author{
Andrea Onofri, Niccolo' Terzaroli, Luigi Russi}
\value{
A design matrix for the MDD effect
}
\examples{
data("hayman54")
MDD(Par1, Par2, data = hayman54)
}
\keyword{ ~diallel }
\keyword{ ~genetic effects }
