rm(list=ls())

# BUGS model definition 1 (Diallel fixed) ##############
mod_diallel.fix <- "
data {
n <- length(Y)
np <- dim(X)
}

model {

# Priors
beta[1] ~ dunif(0, 1000000)
for (i in 2:np[2]){
beta[i] ~ dnorm(0, 0.0001)
}
sigma ~ dunif(0, 100)

# Likelihood
for (i in 1:n) {
   expected[i] <- inprod(X[i,], beta)
   Y[i] ~ dnorm(expected[i], tau)
}

# Derived quantities
tau <- 1 / ( sigma * sigma)
}
"

# BUGS model definition 2 (Diallel fixed, with random blocks) ##############
mod_diallel.ranB <- "
data{
n <- length(Y)
nf <- dim(X)
nr <- dim(Z)
}

model {

# Priors
beta[1] ~ dunif(0, 1000000)
for (i in 2:nf[2]){
beta[i] ~ dnorm(0, 0.001)
}
for (i in 1:nr[2]){
b[i] ~ dnorm(0, tau.r)
}

sigma ~ dunif(0, 100)
sigma.r ~ dunif(0, 100)


# Likelihood
for (i in 1:n) {
   expected[i] <- inprod(X[i,], beta) + inprod(Z[i,], b)
   Y[i] ~ dnorm(expected[i], tau)
}

# Derived quantities
sigma2 <- sigma * sigma
sigma2.r <- sigma.r * sigma.r
tau <- 1 / sigma2
tau.r <- 1 / sigma2.r
}
"

# BUGS model definition 3 (Diallel fixed, with random blocks and Env) ##############
mod_GE2.ranEnvB <- "
data {
n <- length(Y)
nf <- dim(X)
nb1 <- dim(Z.1)
nb2 <- dim(Z.2)
nb3 <- dim(Z.3)
nb4 <- dim(Z.4)
nb5 <- dim(Z.5)
nb6 <- dim(Z.6)
}

model {

# Definition of priors
beta[1] ~ dunif(0, 1000000)
for (i in 2:nf[2]){
beta[i] ~ dnorm(0, 0.0001)
}
for (i in 1:nb1[2]){
b1[i] ~ dnorm(0, tau.1)
}
for (i in 1:nb2[2]){
b2[i] ~ dnorm(0, tau.2)
}
for (i in 1:nb3[2]){
b3[i] ~ dnorm(0, tau.3)
}
for (i in 1:nb4[2]){
b4[i] ~ dnorm(0, tau.4)
}
for (i in 1:nb5[2]){
b5[i] ~ dnorm(0, tau.5)
}
for (i in 1:nb6[2]){
b6[i] ~ dnorm(0, tau.6)
}

sigma ~ dunif(0, 500)
sigma.1 ~ dunif(0, 500)
sigma.2 ~ dunif(0, 500)
sigma.3 ~ dunif(0, 500)
sigma.4 ~ dunif(0, 500)
sigma.5 ~ dunif(0, 500)
sigma.6 ~ dunif(0, 500)

# Likelihood
for (i in 1:n) {
   expected[i] <- inprod(X[i,], beta) + inprod(Z.1[i,], b1) + inprod(Z.2[i,], b2) + inprod(Z.3[i,], b3) + inprod(Z.4[i,], b4) + inprod(Z.5[i,], b5) + inprod(Z.6[i,], b6)
   Y[i] ~ dnorm(expected[i], tau)
}

# Derived quantities
sigma2 <- sigma * sigma
sigma2.1 <- sigma.1 * sigma.1
sigma2.2 <- sigma.2 * sigma.2
sigma2.3 <- sigma.3 * sigma.3
sigma2.4 <- sigma.4 * sigma.4
sigma2.5 <- sigma.5 * sigma.5
sigma2.6 <- sigma.6 * sigma.6
tau <- 1 / sigma2
tau.1 <- 1 / sigma2.1
tau.2 <- 1 / sigma2.2
tau.3 <- 1 / sigma2.3
tau.4 <- 1 / sigma2.4
tau.5 <- 1 / sigma2.5
tau.6 <- 1 / sigma2.6
}
"
#modelDef <- mod_GE2.ranEnvB
bugs_mods <- list(mod_diallel.fix=mod_diallel.fix,
                  mod_diallel.ranB = mod_diallel.ranB,
                  mod_GE2.ranEnvB = mod_GE2.ranEnvB)
save(bugs_mods, file = "data/bugs_mods.rda")
