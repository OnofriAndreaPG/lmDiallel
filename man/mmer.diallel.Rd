\name{mmer.diallel}
\alias{mmer.diallel}
\title{
Fitting random diallel linear models
}
\description{
Wrapper function for the function 'mmer()' in the 'sommer' package. It can be used to fit random diallel models and retreive variance components for main effects.
}
\usage{
mmer.diallel(formula, Block, Env, fct, data, type = "all")
}
\arguments{
  \item{formula}{
     an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.'formula' uses the regular R syntax to specify the response variable and the two variables for parentals \code{(e.g., Yield ~ Par1 + Par2)}
}
\item{Block}{used to specify an optional variable coding for blocks}
\item{Env}{used to specify an optional variable coding for environments}
\item{data}{a 'data.frame' where to look for explanatory variables}
\item{fct}{a string variable coding for the selected model. 8 main diallel models: Hayman's model 1 (="HAYMAN1"), Hayman's model 2 (="HAYMAN2"), Griffing's model 1 (="GRIFFING1"), Griffing's model 2 (="GRIFFING2"), Griffing's model 3 (="GRIFFING3"), Griffing's model 4 (="GRIFFING4"), Gardner-Eberhart model 2 (="GE2")  and Gardner-Eberhart model 3 (="GE3"). The strings "GE2r" and "GE3r" can be used to specify the 'enhanced' GE2 and GE3 models, including the effect of reciprocals (REC).}
\item{type}{a string variable coding for the selected model. It is only used for multi-environment experiments and it is equal to "all" when both the environment and genetical effects are random or "environment" when the environment is random and genetical effects are fixed.}
}

\value{
Returns a data frame of variance components with standard errors
}
\references{
\cite{Covarrubias-Pazaran, G., 2016. Genome-Assisted Prediction of Quantitative Traits Using the R Package sommer. PLOS ONE 11, e0156744. https://doi.org/10.1371/journal.pone.0156744}
}
\author{
  Andrea Onofri \email{(andrea.onofri@unipg.it)}, Niccolo' Terzaroli \email{(n.terzaroli@gmail.com)}, Luigi Russi \email{(luigi.russi@unipg.it)}
}
\examples{
data("hayman54")
rMod <- mmer.diallel(Ftime ~ Par1 + Par2, Block = Block,
                  data = hayman54,
                  fct = "HAYMAN1")
rMod
}

\keyword{ diallel }
\keyword{ genetic effects }
