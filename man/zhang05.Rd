\name{zhang05}
\alias{zhang05}
\docType{data}
\title{
Data for diallel analysis from Zhang (2005)
}
\description{
Data collected in XXX with 5 parents, 2 reps and 2 environments
}
\usage{data("zhang05")}
\format{
  A data.frame with 60 observations on the following 6 variables.
  \describe{
    \item{\code{Par1}}{male parent, a factor with 5 levels}
    \item{\code{Par2}}{female parent, a factor with 5 levels}
    \item{\code{Block}}{block, a factor with 2 levels}
    \item{\code{Combination}}{combination between environment and block,  an integer vector}
    \item{\code{Env}}{environment, a factor with 2 levels}
    \item{\code{Yield}}{yield, a numeric vector}
  }
}
\author{
Andrea Onofri \email{(andrea.onofri@unipg.it)}, Niccolo' Terzaroli \email{(n.terzaroli@gmail.com)}, Luigi Russi \email{(luigi.russi@unipg.it)}
}
\source{
\cite{Zhang, Y., Kang, M.S. and Lamkey, K.R. (2005), DIALLEL-SAS05: A Comprehensive Program for Griffing's and Gardner&Eberhart Analyses. Agron. J., 97: 1097-1106. https://doi.org/10.2134/agronj2004.0260}
}
\references{
\cite{Onofri, A., Terzaroli, N. & Russi, L. Linear models for diallel crosses: a review with R functions. Theor Appl Genet (2020). https://doi.org/10.1007/s00122-020-03716-8}
}
\examples{
data("zhang05")
contrasts(zhang05$Block) <- c("contr.sum")
contrasts(zhang05$Env) <- c("contr.sum")
dMod <- lm(Yield ~ Env/Block + H.BAR(Par1, Par2) + VEi(Par1, Par2) +
                   Hi(Par1, Par2) + SCA(Par1, Par2) +
                   H.BAR(Par1, Par2):Env + GCA(Par1, Par2):Env +
                   Hi(Par1, Par2):Env + SCA(Par1, Par2):Env, data = zhang05)
anova(dMod)
#or
dMod2 <- lm(Yield ~ Env/Block + H.BAR(Par1, Par2) + VEi(Par1, Par2) +
                   Hi(Par1, Par2) + SCA(Par1, Par2), data = zhang05)
summary(dMod2)$coefficients
}
\keyword{datasets}
