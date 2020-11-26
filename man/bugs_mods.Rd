\name{bugs_mods}
\alias{bugs_mods}
\docType{data}
\title{
Definitions for a JAGS model
}
\description{
The definition of a JAGS model may require some experience; to ease the task, we provide several template files as a list in the 'bugs_mods' dataset
}
\usage{data("bugs_mods")}
\format{
  A list with the following 3 objects.
  \describe{
    \item{\code{mod_diallel.fix}}{a character}
    \item{\code{mod_diallel.ranB}}{a character}
    \item{\code{mod_GE2.ranEnvB}}{a character}
    }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\author{
Andrea Onofri \email{(andrea.onofri@unipg.it)}, Niccolo' Terzaroli \email{(n.terzaroli@gmail.com)}, Luigi Russi \email{(luigi.russi@unipg.it)}
}
\source{
\cite{Onofri, A., Terzaroli, N. & Russi, L. Linear models for diallel crosses: a review with R functions. Theor Appl Genet (2020). https://doi.org/10.1007/s00122-020-03716-8}
}
\references{
%%
}
\examples{
data(bugs_mods)
modelDef <- bugs_mods$mod_GE2.ranEnvB
writeLines(modelDef, "modelDef.txt")
}
\keyword{datasets}
