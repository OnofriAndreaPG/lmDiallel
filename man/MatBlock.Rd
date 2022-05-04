\name{matBlock}
\alias{matBlock}
\title{
Design matrix for blocks
}
\description{
It creates a disign matrix for block effects (with sum-to-zero constraint). It is used internally
}
\usage{
matBlock(formula)
}
\arguments{
  \item{formula}{\code{a formula containing the block variables}}
}
\references{
\cite{Onofri, A., Terzaroli, N. & Russi, L. Linear models for diallel crosses: a review with R functions. Theor Appl Genet (2020). https://doi.org/10.1007/s00122-020-03716-8}
}
\author{
Andrea Onofri, Niccolo' Terzaroli, Luigi Russi}
\value{
A design matrix for the block effect
}
\examples{
data("hayman54")
matBlock(~hayman54$Block)
}
\keyword{ ~diallel }
\keyword{ ~genetic effects }
