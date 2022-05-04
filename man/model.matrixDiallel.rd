\name{model.matrixDiallel}
\alias{model.matrixDiallel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Design matrix for Diallel model parametrisation
}
\description{
\code{model.matrixDiallel} is useful to build design matrices, according to the user-defined (or default) parameterisation for \code{lm} function. It shares the same syntax of the \code{lm.diallel} function.
}
\author{
Andrea Onofri, Niccolo' Terzaroli, Luigi Russi}
\arguments{
  \item{formula}{
     an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.'formula' uses the regular R syntax to specify the response variable and the two variables for parentals \code{(e.g., Yield ~ Par1 + Par2)}
}
\item{Block}{used to specify an optional variable coding for blocks}
\item{Env}{used to specify an optional variable coding for environments}
\item{data}{a 'data.frame' where to look for explanatory variables}
\item{fct}{a string variable coding for the selected model. 6 main diallel models: Hayman's model 1 (="HAYMAN1"), Hayman's model 2 (="HAYMAN2"), Griffing's model 1 (="GRIFFING1"), Griffing's model 2 (="GRIFFING2"), Gardner-Eberhart model 2 (="GE2")  and Gardner-Eberhart model 3 (="GE3"). The strings "GE2r" and "GE3r" can be used to specify the 'enhanced' GE2 and GE3 models, including the effect of reciprocals (REC).}
}
\details{
model.matrixDiallel creates a design matrix for a diallel model, as specified in the 'fct' argument.
}
\value{
The design matrix for a diallel model as specified in the 'fct' argument.
}
\references{
\cite{Onofri, A., Terzaroli, N. & Russi, L. Linear models for diallel crosses: a review with R functions. Theor Appl Genet (2020). https://doi.org/10.1007/s00122-020-03716-8}
}
\examples{
data("diallelMET")
ModMat <- model.matrixDiallel(Yield ~ Par1 + Par2,
                   Env, Block, fct= "GE3",
                   data = diallelMET)
}
\keyword{ ~diallel }
\keyword{ ~genetic effects }
