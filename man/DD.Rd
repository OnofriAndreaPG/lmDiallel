\name{DD}
\alias{DD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Dominant Deviation effect
}
\description{
DD effect to fit Hayman2 model with \code{lm} function
}
\usage{
DD(x, y, data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{\code{a variable for the first parent}}
  \item{y}{\code{a variable for the second parent}}
  \item{data}{\code{a 'data.frame' where to look for explanatory variables}}
}
\details{
a design matrix of n-1 parentals
}
\value{
%%  ~Describe the value returned
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

\examples{
data("hayman54")
contrasts(hayman54$Block) <- "contr.sum"
dMod <- lm(Ftime ~ Block + GCA(Par1, Par2) + MDD(Par1, Par2) +
             DD(Par1, Par2) + SCA(Par1, Par2) +
             RGCA(Par1, Par2) + RSCA(Par1, Par2), data = hayman54)
summary(dMod)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~diallel }% use one of  RShowDoc("KEYWORDS")
