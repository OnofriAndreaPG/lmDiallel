
MLfit1 <- function(X, Y, z, Par1, Par2){
  # Maximum likelihood fit for HAYMAN1
  # Building new incidence matrix
  #z <- fit
  #mat <- fit$modMatrix
  #Par1 <- df$Par1; Par2 <- df$Par2
  mat <- X
  P1 <- factor(Par1)
  P2 <- factor(Par2)
  p <- length(levels(P1))
  contrasts(P1) <- c("contr.sum")
  contrasts(P2) <- c("contr.sum")
  Z1 <- model.matrix(~P1)[,-1]
  Z2 <- model.matrix(~P2)[,-1]
  colnames(Z1) <- paste("k_", c(1:(p-1)), sep="")
  colnames(Z2) <- colnames(Z1)
  X <- cbind(mat[,1:p], Z1, Z2, mat[,(2*p):length(mat[1,])])
  groups <- 0:3
  attr(X, "namEff") <- c("Intercept", "GCA", "RGCA", "SCA", "RSCA")
  asgn <- rep(groups, c(1, p - 1,
                        2 * (p - 1),
                        (p^2 - p) / 2))
  attr(X, "assign") <- asgn
  sigma <- sqrt(sum(z$residuals^2)/z$df.residual)
  paramInit <- c(z$coefficients, "sigma" = sigma)

  logLik <- function(param, x, y, p){
    # x è la matrice di incidenza corretta
    # param deve contenere il sigma in fondo
    sigmaI <- param[length(param)]
    param <- param[-length(param)]
    mu <- Hayman1.fun(x, param, p)
    res <- dnorm(y, mu, sigmaI, log=T)
    .value <- -sum(res)
    return(.value)
  }

  # Maximum likelihood fit
  opt <- optim(par = paramInit,
               fn = logLik,
               method = "BFGS",
               x = X,
               y = Y, p = p,
               hessian = T)
  estpar <- opt$par
  coefs <- estpar
  coefs <- coefs[c(1:p, (p+1):(p + p - 1), (p+1):(p + p -1), (p+p):(length(coefs)-1))]
  coefs[(2*p):(3*p-2)] <- -coefs[(2*p):(3*p-2)]
  sesML <- sqrt(diag(solve(opt$hessian)))
  tab <- cbind(estpar, sesML)
  #row.names(tab) <- c(names(z$coefficients), "sigma")
  #print(names(z$coefficients))
  #print(row.names(tab))
  exp <- Hayman1.fun(X, estpar[-length(estpar)], p)
  res <- df$Ftime - exp
  RSS <- sum(res^2)
  dfreed <- length(X[,1]) - length(opt$par) + 1
  sigma.mod <- sqrt(RSS/dfreed)
  sigma.ML <- tab[length(tab[,1]), 1]
  tabLS <- data.frame(Estimate = estpar[-length(estpar)],
                      SE = sesML[-length(sesML)] * sigma.mod/sigma.ML)
  #print(row.names(tabLS))
  returnList<- list(MLcoefficients = tab, fitted = exp,
                    residuals = res,
              deviance = RSS, df.residual = dfreed,
              sigma = sigma.mod, LScoefficients = tabLS,
              X = X, hessian = opt$hessian, assign = asgn,
              coefficients = coefs)

  return(returnList)
}

MLfit1b <- function(X, Y, z, Par1, Par2){

  Hayman1b.fun <- function(x, param, p, assign1, assign2){
  # Param contiene il sigma!!!!!!!
  #print(assign1); print(assign2)
  lcoef <- split(param, assign1)
  modMat <- lapply( split(as.data.frame(t(x)), assign2), t)
  separ <- length(modMat$`4`[1,])/2
  .expr0 <- modMat$`0` %*% lcoef$`0`
  .expr1 <- modMat$`1` %*% lcoef$`1`
  .expr2 <- modMat$`2` %*% lcoef$`2`
  .expr3 <- modMat$`3` %*% lcoef$`3`
  .expr4 <- modMat$`4`[,1:separ] %*% lcoef$`4`
  .expr5 <- modMat$`4`[,(separ+1):(separ * 2)] %*% (- lcoef$`4`)
  .expr6 <-modMat$`5` %*% lcoef$`5`
  exp <- .expr0 + .expr1 + .expr2 +.expr3 +.expr4 +.expr5 +.expr6
  exp
}

  logLik <- function(param, x, y, p, assign1, assign2){
    # x è la matrice di incidenza corretta
    # param deve contenere il sigma in fondo,
    # che viene rimosso prima del calcolo
    # della risposta attesa

    sigmaI <- param[length(param)]
    lcoef <- split(param, assign1)
    modMat <- lapply( split(as.data.frame(t(x)), assign2), t)
    separ <- length(modMat$`4`[1,])/2
    mu <- Hayman1b.fun(x, param, p, assign1, assign2)
    res <- dnorm(y, mu, sigmaI, log=T)
    .value <- -sum(res)
    return(.value)
  }

  # Maximum likelihood fit for HAYMAN1
  p <- length(levels(factor(Par1)))
  sigma <- sqrt(sum(z$residuals^2)/z$df.residual)
  paramInit <- c(z$coefficients, "sigma" = sigma)
  paramInit[is.na(paramInit)] <- 0
  assign1 <- c(z$assign, max(z$assign) + 1)
  assign2 <- attr(X, "assign")

  # Maximum likelihood fit
  opt <- optim(par = paramInit,
               fn = logLik,
               method = "BFGS",
               x = X,
               y = Y, p = p, assign1 = assign1,
               assign2 = assign2,
               hessian = T, control = list(trace = T))
  estpar <- opt$par
  coefs <- estpar
  sesML <- sqrt(diag(solve(opt$hessian)))
  # lcoefs <- split(estpar, assign1)
  # lcoefs$`4` <- c(lcoefs$`4`, -lcoefs$`4`)
  tab <- cbind(estpar, sesML)
  tab
  exp <- Hayman1b.fun(X, estpar, p, assign1, assign2)
  res <- df$Ftime - exp
  RSS <- sum(res^2)
  dfreed <- length(X[,1]) - length(opt$par) + 1
  sigma.mod <- sqrt(RSS/dfreed)
  sigma.ML <- tab[length(tab[,1]), 1]
  tabLS <- data.frame(Estimate = estpar[-length(estpar)],
                      SE = sesML[-length(sesML)] * sigma.mod/sigma.ML)
  #print(row.names(tabLS))
  returnList<- list(MLcoefficients = tab, fitted = exp,
                    residuals = res,
              deviance = RSS, df.residual = dfreed,
              sigma = sigma.mod, LScoefficients = tabLS,
              X = X, hessian = opt$hessian, assign = assign2,
              coefficients = coefs)

  return(returnList)
}

MLfit2 <- function(X, Y, z, Par1, Par2){
  # Maximum likelihood fit for HAYMAN2

  # Building new incidence matrix
  # z <- fit; y <- df$Ftime
  # mat <- fit$modMatrix
  # Par1 <- df$Par1; Par2 <- df$Par2
  mat <- X
  P1 <- factor(Par1)
  P2 <- factor(Par2)
  p <- length(levels(P1))
  contrasts(P1) <- c("contr.sum")
  contrasts(P2) <- c("contr.sum")
  Z1 <- model.matrix(~P1)[,-1]
  Z2 <- model.matrix(~P2)[,-1]
  colnames(Z1) <- paste("k_", c(1:(p-1)), sep="")
  colnames(Z2) <- colnames(Z1)
  Z <- mat[,1:(length(mat[1,]) - (p-1))]
  X <- cbind(Z, Z1, Z2)

  attr(X, "namEff") <- attr(mat, "namEff")
  asgn <- attr(mat, "assign")
  asgn <- c(asgn, rep(asgn[length(asgn)], p - 1))
  attr(X, "assign") <- asgn
  sigma <- sqrt(sum(z$residuals^2)/z$df.residual)
  paramInit <- c(z$coefficients, "sigma" = sigma)

  logLik2 <- function(param, x, y, p){
    # x è la matrice di incidenza corretta
    # param deve contenere il sigma in fondo
    sigmaI <- param[length(param)]
    param <- param[-length(param)]
    mu <- Hayman2.fun(x, param, p)
    res <- dnorm(y, mu, sigmaI, log=T)
    .value <- -sum(res)
    return(.value)
  }

  # Maximum likelihood fit
  opt <- optim(par = paramInit,
               fn = logLik2,
               method = "BFGS",
               x = X,
               y = Y, p = p,
               hessian = T)

  estpar <- opt$par
  coefs <- estpar
  add <- coefs[(length(coefs)-p + 1):(length(coefs) - 1)]
  coefs <- c(coefs[-length(coefs)], -add)
  sesML <- sqrt(diag(solve(opt$hessian)))
  tab <- cbind(estpar, sesML)
  #row.names(tab) <- c(names(z$coefficients), "sigma")
  #print(names(z$coefficients))
  #print(row.names(tab))
  exp <- Hayman2.fun(X, estpar[-length(estpar)], p)
  res <- df$Ftime - exp
  RSS <- sum(res^2)
  dfreed <- length(X[,1]) - length(opt$par) + 1
  sigma.mod <- sqrt(RSS/dfreed)
  sigma.ML <- tab[length(tab[,1]), 1]
  tabLS <- data.frame(Estimate = estpar[-length(estpar)],
                      SE = sesML[-length(sesML)] * sigma.mod/sigma.ML)
  #print(row.names(tabLS))
  returnList<- list(MLcoefficients = tab, fitted = exp,
                    residuals = res,
                    deviance = RSS, df.residual = dfreed,
                    sigma = sigma.mod, LScoefficients = tabLS,
                    X = X, hessian = opt$hessian, assign = asgn,
                    coefficients = coefs)

  return(returnList)
}

Hayman1.fun <- function(x, param, p){
  # Param non contiene il sigma!!!!!!!
  Z <- x[,2:p]
  Z1 <- x[,(p+1):(2*p-1)]; Z2 <- x[,(2*p):(3*p - 2)];
  SCA <- x[,(3*p - 1):length(x[1,])]
  int <- param[1]; c1 <- param[2:p]; c2 <- param[(p+1):(2*p-1)]
  c3 <- param[(2*p):length(param)]
  int + Z %*% c1 + Z1 %*% c2 + Z2 %*% -c2 + SCA %*% c3
}


Hayman2.fun <- function(x, param, p){
  limit <- (2 + 2*(p - 1) + 1/2*p*(p - 3))
  Z <- x[,1:limit]
  Z1 <- x[,(limit+1):(limit+p-1)]
  Z2 <- x[,(limit+p):(limit+p+6)];
  param1 <- param[1:limit]
  param2 <- param[(limit+1):(limit+p-1)]
  Z %*% param1 + Z1 %*% param2 + Z2 %*% - param2
}

