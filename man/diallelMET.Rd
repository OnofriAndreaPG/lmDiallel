\name{diallelMET}
\alias{diallelMET}
\docType{data}
\title{
Factitious dataset for Diallel analysis
}
\description{
Multi-environment half-diallel dataset with six parentals, in five blocks and ten environments; the dataset is factitious and was obtained by Monte Carlo simulation.
}
\usage{data("diallelMET")}
\format{
  A data.frame with 1050 observations on the following 5 variables.
  \describe{
    \item{\code{Env}}{environment, a factor with 10 levels}
    \item{\code{Block}}{block, a factor with 5 levels}
    \item{\code{Par1}}{male parent, a factor with 6 levels}
    \item{\code{Par2}}{female parent, a factor with 6 levels}
    \item{\code{Yield}}{yield, a numeric vector}
  }
}
\author{
Andrea Onofri \email{(andrea.onofri@unipg.it)}, Niccolo' Terzaroli \email{(n.terzaroli@gmail.com)}, Luigi Russi \email{(luigi.russi@unipg.it)}
}
\source{
Onofri, A., Terzaroli, N. & Russi, L. Linear models for diallel crosses: a review with R functions. Theor Appl Genet (2020). https://doi.org/10.1007/s00122-020-03716-8
}
\references{
Onofri, A., Terzaroli, N. & Russi, L. Linear models for diallel crosses: a review with R functions. Theor Appl Genet (2020). https://doi.org/10.1007/s00122-020-03716-8
}
\examples{
data("diallelMET")
}
\keyword{datasets}
