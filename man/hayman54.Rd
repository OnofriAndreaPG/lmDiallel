\name{hayman54}
\alias{hayman54}
\docType{data}
\title{
Hayman dataset for diallel analysis
}
\description{
Data for a diallel in tobacco with 2 reps
}
\usage{data(hayman54)}
\format{
  A data.frame with 128 observations on the following 4 variables
  \describe{

    \item{\code{Block}}{block, a factor with 2 levels}
    \item{\code{Par1}}{male parent, a factor with 8 levels}
    \item{\code{Par2}}{female parent, a factor with 8 levels}
    \item{\code{Ftime}}{mean flowering time (days), a numeric vector}
  }
}
\details{

}
\author{
Andrea Onofri \email{(andrea.onofri@unipg.it)}, Niccolo' Terzaroli \email{(n.terzaroli@gmail.com)}, Luigi Russi \email{(luigi.russi@unipg.it)}
}

\source{
\cite{B. I. Hayman (1954a). The Analysis of Variance of Diallel Tables. Biometrics, 10, 235-244. Table 5, page 241. http://doi.org/10.2307/3001877}
}
\references{
\cite{Onofri, A., Terzaroli, N. & Russi, L. Linear models for diallel crosses: a review with R functions. Theor Appl Genet (2020). https://doi.org/10.1007/s00122-020-03716-8}
}
\examples{
data(hayman54)

contrasts(hayman54$Block) <- c("contr.sum")
dMod <- lm(Ftime ~ Block  + GCA(Par1, Par2)
             + tSCA(Par1, Par2) + RGCA(Par1, Par2)
             + RSCA(Par1, Par2), data = hayman54)
anova(dMod)
#or
dMod2 <- lm.diallel(Ftime ~ Par1 + Par2, Block = Block,
                  data = hayman54,
                  fct = "HAYMAN1")
anova(dMod2)
}

\keyword{datasets}
