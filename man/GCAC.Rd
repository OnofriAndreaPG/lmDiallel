\name{GCAC}
\alias{GCAC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
General Combining Ability without considering the selfed parents
}
\description{
Design matrix for GCAC, useful to fit Gardner & Eberhart model 3 (GE3) with \code{lm} function
}
\usage{
GCAC(x,y,data)
}

\arguments{
 \item{x}{\code{a variable for the first parent}}
  \item{y}{\code{a variable for the second parent}}
  \item{data}{\code{a 'data.frame' where to look for explanatory variables}}
}
\details{
a design matrix of n-1 parentals
}
\value{

%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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

\seealso{
\describe{
\item{\code{\code{\link{model.matrixDiallel}}}}
\item{\code{\code{\link{lm.diallel}}}}
\item{\code{\code{\link{SCA}}}}
\item{\code{\code{\link{tSCA}}}}
\item{\code{\code{\link{RSCA}}}}
\item{\code{\code{\link{GCA}}}}
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
}}
\examples{
data("hayman54")
GCAC(Par1,Par2, data=hayman54)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
