\name{lonnquist61}
\alias{lonnquist61}
\docType{data}
\title{
Half diallel of maize dataset
}
\description{
Diallel experiment with six maize varieties and no reciprocals. The data here are means adjusted for block effects.
}
\usage{data("lonnquist61")}
\format{
  A data.frame with 21 observations on the following 3 variables.
  \describe{
    \item{\code{Par1}}{male parent, a factor with 6 levels}
    \item{\code{Par2}}{female parent, a factor with 6 levels}
    \item{\code{Yield}}{mean across blocks, a numeric vector}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\author{
Andrea Onofri \email{(andrea.onofri@unipg.it)}, Niccolo' Terzaroli \email{(n.terzaroli@gmail.com)}, Luigi Russi \email{(luigi.russi@unipg.it)}
}
\source{
\cite{J. H. Lonnquist, C. O. Gardner. (1961) Heterosis in Intervarietal Crosses in Maize and Its Implication in Breeding Procedures. Crop Science, 1, 179-183. Table 1.}
}
\references{
\cite{Onofri, A., Terzaroli, N. & Russi, L. Linear models for diallel crosses: a review with R functions. Theor Appl Genet (2020). https://doi.org/10.1007/s00122-020-03716-8}

\cite{Mohring, Melchinger, Piepho. (2011). REML-Based Diallel Analysis. Crop Science, 51, 470-478. http://doi.org/10.2135/cropsci2010.05.0272}

\cite{C. O. Gardner and S. A. Eberhart. 1966. Analysis and Interpretation of the Variety Cross Diallel and Related Populations. Biometrics, 22, 439-452. http://doi.org/10.2307/2528181}
}
\examples{
data(lonnquist61)
dMod <- lm(Yield ~ H.BAR(Par1, Par2) + VEi(Par1, Par2) +
             Hi(Par1, Par2) + SCA(Par1, Par2),
                  data = lonnquist61)
summary.diallel(dMod, MSE = 7.10, dfr = 60)
anova.diallel(dMod, MSE = 7.10, dfr = 60)
}
\keyword{datasets}
