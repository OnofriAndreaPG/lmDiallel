\name{lmDiallel-methods}
\alias{vcov.diallel}
\alias{summary.diallel}
\alias{anova.diallel}
\alias{predict.diallel}
\alias{model.matrix.diallel}
\alias{glht.diallel}
\title{ Methods for diallel model fitting }
\description{
  The object returned by the 'lm.diallel()' function is of classes 'lm' and 'diallel'. Specific methods were devised to explore the 'diallel' object.
}
\usage{
\method{summary}{diallel}(object, MSE, dfr, ...)
\method{vcov}{diallel}(object, MSE, ...)
\method{anova}{diallel}(object, MSE, dfr, ...)
\method{model.matrix}{diallel}(object, ...)
}
\arguments{
  \item{object}{ an object of class diallel.}
  \item{MSE}{ Mean Square Error, when it cannot be derived from model fit}
  \item{dfr}{ Residual degrees of freedom, when they cannot be derived from model fit}
  \item{...}{ Other optional arguments}
}
\details{
To be defined
}
\value{
  \code{summary} computes (adjusted) p values for general linear hypotheses,
}
\references{
\cite{Onofri, A., Terzaroli, N. & Russi, L. Linear models for diallel crosses: a review with R functions. Theor Appl Genet (2020). https://doi.org/10.1007/s00122-020-03716-8}
}
\examples{
library(lmDiallel)
library(multcomp)
# Data with replicates
data("hayman54")

fit <- lm.diallel(Ftime ~ Par1 + Par2, data = hayman54,
                  fct = "HAYMAN1")
summary(fit)
anova(fit)
gh <- glht(linfct = diallel.eff(fit), adjust = "none")
}

\keyword{diallel}
