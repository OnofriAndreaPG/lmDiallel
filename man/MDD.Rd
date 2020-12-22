\name{MDD}
\alias{MDD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mean Dominance Deviation effect
}
\description{
It relates to the difference between the average yield of selfed parents and the average yield of crosses. DD effect to fit Hayman2 model with \code{lm} function
}
\usage{
MDD(x, y, data)
}
\arguments{
  \item{x}{\code{a variable for the first parent}}
  \item{y}{\code{a variable for the second parent}}
  \item{data}{\code{a 'data.frame' where to look for explanatory variables}}
}
\details{
??? Only a column???
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
MDD(Par1, Par2, data = hayman54)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~diallel }% use one of  RShowDoc("KEYWORDS")
