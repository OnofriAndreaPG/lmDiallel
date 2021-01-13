\name{diallelMod-methods}
\alias{glht.diallelMod}
\alias{diallel.eff}
\alias{G1.eff}
\alias{G2.eff}
\alias{G3.eff}
\alias{G4.eff}
\alias{hayman1.eff}
\alias{hayman2.eff}
\title{Functions to retreive full list of genetical effects}
\description{
  Diallel model parameters are estimated under a set of restrictions and, therefore, the methods 'coef' and 'summary' do not return the full list of genetical parameters. Therefore, the 'glht.diallelMod' method can be used, which works by way of a series of helper functions, providing the necessary contrast matrices.
}
\usage{
\method{glht}{diallelMod}(model, linfct, ...)
}
\arguments{
  \item{model}{a model object (OPTIONAL)}
  \item{linfct}{ a diellel.eff() function}
  \item{...}{ Other optional arguments}
}
\details{
...
}
\value{
  \code{summary} Returns the full list of genetical parameters
}
\references{
\cite{Onofri, A., Terzaroli, N. & Russi, L. Linear models for diallel crosses: a review with R functions. Theor Appl Genet (2020). https://doi.org/10.1007/s00122-020-03716-8}
}
\examples{
library(lmDiallel)

# Data with replicates
data("hayman54")

fit <- lm.diallel(Ftime ~ Par1 + Par2, data = hayman54,
                  fct = "HAYMAN1")
summary(fit)
anova(fit)
gh <- glht(linfct = diallel.eff(fit), adjust = "none")
}

\keyword{diallel}
