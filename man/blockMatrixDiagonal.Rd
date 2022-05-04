\name{blockMatrixDiagonal}
\alias{blockMatrixDiagonal}
\title{
Creates block diagonal matrix. It is used internally.
}
\description{
This function takes a list of matrices and creates a block diagonal matrix. It is used to fit multi-environment diallel models
}
\usage{
blockMatrixDiagonal(matList)
}
\arguments{
  \item{matList}{
It is a list of matrices to be combined
}
}

\value{
Returns a matrix object
}

\references{
\cite{Onofri, A., Terzaroli, N. & Russi, L. Linear models for diallel crosses: a review with R functions. Theoretical Applied Genetics (2020). https://doi.org/10.1007/s00122-020-03716-8}
}
\author{
Andrea Onofri, Niccolo' Terzaroli, Luigi Russi}

\examples{
a <- matrix(1:16, 8, 2)
b <- matrix(1:9, 3, 3)
c <- list(a, b)
blockMatrixDiagonal(c)
}

