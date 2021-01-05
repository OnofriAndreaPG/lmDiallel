\name{diallelMET}
\alias{diallelMET}
\docType{data}
\title{
Factitious dataset for Diallel analysis
}
\description{
Multi-environment half-diallel dataset with six parentals, in five blocks and ten environments; the dataset is factitious and was obtained by Monte Carlo simulation.
}
\usage{data("diallelMET")}
\format{
  A data.frame with 1050 observations on the following 5 variables.
  \describe{
    \item{\code{Env}}{environment, a factor with 10 levels}
    \item{\code{Block}}{block, a factor with 5 levels}
    \item{\code{Par1}}{male parent, a factor with 6 levels}
    \item{\code{Par2}}{female parent, a factor with 6 levels}
    \item{\code{Yield}}{yield, a numeric vector}
  }
}
\author{
Andrea Onofri \email{(andrea.onofri@unipg.it)}, Niccolo' Terzaroli \email{(n.terzaroli@gmail.com)}, Luigi Russi \email{(luigi.russi@unipg.it)}
}
\source{
Onofri, A., Terzaroli, N. & Russi, L. Linear models for diallel crosses: a review with R functions. Theor Appl Genet (2020). https://doi.org/10.1007/s00122-020-03716-8
}
\references{
Onofri, A., Terzaroli, N. & Russi, L. Linear models for diallel crosses: a review with R functions. Theor Appl Genet (2020). https://doi.org/10.1007/s00122-020-03716-8
}
\examples{
data("diallelMET")
contrasts(diallelMET$Block) <- c("contr.sum")
contrasts(diallelMET$Env) <- c("contr.sum")

data("bugs_mods")
modelDef <- bugs_mods$mod_GE2.ranEnvB
writeLines(modelDef, "modelDef.txt")
# Define design matrices
modMat <- model.matrix(~ H.BAR(Par1, Par2) + VEi(Par1, Par2) +
                   Hi(Par1, Par2) + SCA(Par1, Par2) + Env + Env:Block +
                   H.BAR(Par1, Par2):Env + VEi(Par1, Par2):Env +
                   Hi(Par1, Par2):Env + SCA(Par1, Par2):Env,
                   data = diallelMET)
X <- modMat[,attr(modMat, "assign") <= 4]
Z.1 <- modMat[,attr(modMat, "assign") == 5]
Z.2 <- modMat[,attr(modMat, "assign") == 6]
Z.3 <- modMat[,attr(modMat, "assign") == 7]
Z.4 <- modMat[,attr(modMat, "assign") == 8]
Z.5 <- modMat[,attr(modMat, "assign") == 9]
Z.6 <- modMat[,attr(modMat, "assign") == 10]


# Get starting values for fixed effects
start <- lm(Yield ~ H.BAR(Par1, Par2) + VEi(Par1, Par2) +
                   Hi(Par1, Par2) + SCA(Par1, Par2), data = diallelMET)
betaInit  <- coef(start)

# Create input lists for WinBugs
dataMod <- list(Y = diallelMET$Yield, X = X, Z.1 = Z.1, Z.2 = Z.2,
                Z.3 = Z.3, Z.4 = Z.4,
                Z.5 = Z.5, Z.6 = Z.6)
init <- list(beta = betaInit, sigma = 0.1, sigma.1 = 0.1,
                 sigma.2 = 0.1, sigma.3 = 0.1,
                 sigma.4 = 0.1, sigma.5 = 0.1,
                 sigma.6 = 0.1)

# Start sampler (you need JAGS installed on your computer)
library(rjags)
mcmc <- jags.model("modelDef.txt",
                   data = dataMod, inits = init,
                   n.chains = 4, n.adapt = 1000)
# Specify the parameters to be sampled from posterior
params <- c("beta", "sigma2", "sigma2.1", "sigma2.2", "sigma2.3",
                "sigma2.4", "sigma2.5", "sigma2.6")
res <- coda.samples(mcmc, params, n.iter = 10000)
out <- summary(window(res, start = 1000))
out$quantiles[,c(1,3,5)]
}
\keyword{datasets}
