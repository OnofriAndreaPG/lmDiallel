\name{lmDiallel-methods}
\alias{vcov.diallel}
\alias{summary.diallel}
\alias{anova.diallel}
\alias{predict.diallel}
\alias{glht.diallel}
\title{ Methods for diallel model fitting }
\description{
  The object returned by the 'lm.diallel()' function is of classes 'lm' and 'diallel'. Specific methods were devised to explore the 'diallel' object.
}
\usage{
\method{summary}{diallel}(object, MSA, dfr)
\method{vcov}{diallel}(object, MSA, dfr)
\method{anova}{diallel}(object, MSA, dfr)
}
\arguments{
  \item{object}{ an object of class \code{\link{diallel}}.}
  \item{MSA}{ Mean Square Error, when it cannot be derived from model fit}
  \item{dfr}{ Residual degrees of freedom, when they cannot be derived from model fit}
}
\details{

}
\value{

  \code{summary} computes (adjusted) p values for general linear hypotheses,
  \code{\link{confint}} computes (adjusted) confidence intervals.
  \code{\link{coef}} returns estimates of the linear function \eqn{K \theta}
  and \code{\link{vcov}} its covariance.

}
\references{


}
\examples{

  }

\keyword{diallel}
