\name{model.matrixDiallel}
\alias{model.matrixDiallel}
\alias{model.matrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Design matrix for Diallel model parametrisation
}
\description{
\code{model.matrixDiallel} is useful to build design matrices, according to the user-defined (or default) parameterisation for \code{lm} function. It shares the same syntax of the \code{lm.diallel} function.
}
\usage{
model.matrixDiallel(formula, Block, Env, data, fct)
}
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
RIGUARDARE
model.matrixDiallel creates a design matrix from the description given in terms(object), using the data in data which must supply variables with the same names as would be created by a call to model.frame(object) or, more precisely, by evaluating attr(terms(object), "variables"). If data is a data frame, there may be other columns and the order of columns is not important. Any character variables are coerced to factors. After coercion, all the variables used on the right-hand side of the formula must be logical, integer, numeric or factor.

If contrasts.arg is specified for a factor it overrides the default factor coding for that variable and any "contrasts" attribute set by C or contrasts. Whereas invalid contrasts.args have been ignored always, they are warned about since R version 3.6.0.

In an interaction term, the variable whose levels vary fastest is the first one to appear in the formula (and not in the term), so in ~ a + b + b:a the interaction will have a varying fastest.

By convention, if the response variable also appears on the right-hand side of the formula it is dropped (with a warning), although interactions involving the term are retained.
}
\value{
RIGUARDARE
The design matrix for a regression-like model with the specified formula and data.

There is an attribute "assign", an integer vector with an entry for each column in the matrix giving the term in the formula which gave rise to the column. Value 0 corresponds to the intercept (if any), and positive values to terms in the order given by the term.labels attribute of the terms structure corresponding to object.

If there are any factors in terms in the model, there is an attribute "contrasts", a named list with an entry for each factor. This specifies the contrasts that would be used in terms in which the factor is coded by contrasts (in some terms dummy coding may be used), either as a character vector naming a function or as a numeric matrix.
}
\references{
\cite{Onofri, A., Terzaroli, N. & Russi, L. Linear models for diallel crosses: a review with R functions. Theor Appl Genet (2020). https://doi.org/10.1007/s00122-020-03716-8}
}
\author{
Andrea Onofri \email{(andrea.onofri@unipg.it)}, Niccolo' Terzaroli \email{(n.terzaroli@gmail.com)}, Luigi Russi \email{(luigi.russi@unipg.it)}
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
data("diallelMET")
ModMat <- model.matrixDiallel(Yield ~ Par1 + Par2,
                   Env, Block, fct= "GE3",
                   data = diallelMET)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
