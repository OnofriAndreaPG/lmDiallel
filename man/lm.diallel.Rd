\name{lm.diallel}
\alias{lm.diallel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fitting diallel linear models
}
\description{
Wrapper function for lm.fit and diallel models. It can be used to carry out several powerful methods for linear models, such as 'summary()', anova() or 'glht()' in the 'multcomp' package.
}
\usage{
lm.diallel(formula, Block, Env, data, fct)
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
Notations for the 8 models
\tabular{lc}{
Model name in 'lm.diallel()' \tab 	Model notation in 'lm()' \cr
HAYMAN1 	\tab Y ~ GCA(Par1, Par2) + tSCA(Par1, Par2) + RGCA(Par1, Par2) + RSCA(Par1, Par2) \cr
GRIFFING1	\tab Y ~ GCA(Par1, Par2) + tSCA(Par1, Par2) + REC(Par1, Par2) \cr
GRIFFING2	\tab Y ~ GCA(Par1, Par2) + tSCA(Par1, Par2) \cr
HAYMAN2	\tab Y ~ GCA(Par1, Par2) + MDD(Par1, Par2) + DD(Par1, Par2) + SCA(Par1, Par2) + RGCA(Par1, Par2) + RSCA(Par1, Par2) \cr
GE2	\tab Y ~ H.BAR(Par1, Par2) + VE.i(Par1, Par2) + H.i(Par1, Par2) + SCA(Par1, Par2) \cr
GE3	\tab Y ~ H.BAR(Par1, Par2) + SP(Par1, Par2) + GCAC(Par1, Par2) + SCA(Par1, Par2) \cr
GE2r \tab Y ~ H.BAR(Par1, Par2) + VE.i(Par1, Par2) + H.i(Par1, Par2) + SCA(Par1, Par2) + RGCA(Par1, Par2) + RSCA(Par1, Par2) \cr
GE3r	\tab Y ~ H.BAR(Par1, Par2) + SP(Par1, Par2) + GCAC(Par1, Par2) + SCA(Par1, Par2) + RGCA(Par1, Par2) + RSCA(Par1, Par2)}
}
\value{
RIGUARDARE TUTTO
lm.diallel returns an object of class "lm" or for multiple responses of class c("mlm", "lm").

The functions summary and anova are used to obtain and print a summary and analysis of variance table of the results. The generic accessor functions coefficients, effects, fitted.values and residuals extract various useful features of the value returned by lm.

An object of class "lm" is a list containing at least the following components:
\describe{
    \item{\code{coefficients}}{a named vector of coefficients}
    \item{\code{residuals}}{the residuals, that is response minus fitted values}
    \item{\code{fitted.values}}{the fitted mean values}
    \item{\code{rank}}{the numeric rank of the fitted linear models}
    \item{\code{weights}}{(only for weighted fits) the specified weights}
    \item{\code{df.residual}}{the residual degrees of freedom}
    \item{\code{call}}{the matched call}
    \item{\code{terms}}{the terms object used}
    \item{\code{contrasts	}}{(only where relevant) the contrasts used}
    \item{\code{xlevels}}{(only where relevant) a record of the levels of the factors used in fitting}
    \item{\code{call}}{the matched call}
    \item{\code{offset}}{the offset used (missing if none were used)}
    \item{\code{y}}{if requested, the response used}
    \item{\code{x}}{if requested, the model matrix used}
    \item{\code{model}}{if requested (the default), the model frame used}
    \item{\code{na.action}}{(where relevant) information returned by model.frame on the special handling of NAs}
    }
In addition, non-null fits will have components assign, effects and (unless not requested) qr relating to the linear fit, for use by extractor functions such as summary and effects.
}
\references{
 \cite{Onofri, A., Terzaroli, N. & Russi, L. Linear models for diallel crosses: a review with R functions. Theor Appl Genet (2020). https://doi.org/10.1007/s00122-020-03716-8}
}
\author{
  Andrea Onofri \email{(andrea.onofri@unipg.it)}, Niccolo' Terzaroli \email{(n.terzaroli@gmail.com)}, Luigi Russi \email{(luigi.russi@unipg.it)}
}
\note{
  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\describe{
\item{\code{\code{\link{model.matrixDiallel}}}}
\item{\code{\code{\link{RSCA}}}}
\item{\code{\code{\link{tSCA}}}}
\item{\code{\code{\link{SCA}}}}
\item{\code{\code{\link{GCA}}}}
\item{\code{\code{\link{GCAC}}}}
\item{\code{\code{\link{RGCA}}}}
\item{\code{\code{\link{REC}}}}
\item{\code{\code{\link{SP}}}}
\item{\code{\code{\link{DD}}}}
\item{\code{\code{\link{MDD}}}}
\item{\code{\code{\link{Hi}}}}
\item{\code{\code{\link{H.BAR}}}}
\item{\code{\code{\link{Vei}}}}
\item{\code{\code{\link{summary.diallel}}}}
\item{\code{\code{\link{anova.diallel}}}}
\item{\code{\code{\link{matBlock}}}}
}}}
}
\examples{
data(hayman54)
dMod <- lm.diallel(Ftime ~ Par1 + Par2, Block = Block,
                  data = hayman54,
                  fct = "HAYMAN1")
anova(dMod)

#or
data(zhang05)
contrasts(zhang05$Block) <- c("contr.sum")
dMod2 <- lm.diallel(Yield ~ Par1 + Par2, Env = Env, Block = Block,
            data = zhang05, fct = "GE2")

#or
dMod3 <- lm.diallel(Yield ~ Par1 + Par2,
                  data = lonnquist61,
                  fct = "GE2")
summary(dMod3, MSE = 7.10, dfr = 60)
anova(dMod3, MSE = 7.10, dfr = 60)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ diallel }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
