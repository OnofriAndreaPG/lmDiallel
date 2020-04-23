lm.diallel <- function(formula, Block = NULL, Env = NULL,
                       fct = "GRIFFING2", data,
                       ML = FALSE){
  #formula Ftime + Block + Par1 + Par2, data = df
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)

  m <- match(c("formula", "Block", "Env", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)

  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  Y <- model.response(mf, "numeric")

  mlm <- is.matrix(Y)
  ny <- if (mlm)
        nrow(Y)
  else length(Y)

  bName <- deparse(substitute(Block))  # storing name of Blocks
  Block <- model.extract(mf, "Block")

  eName <- deparse(substitute(Env))  # storing name of Env
  Env <- model.extract(mf, "Env")

  pars <- attr(mt, "term.labels")
  Par1 <- data[[pars[1]]]
  Par2 <- data[[pars[2]]]
  X <- model.matrixDiallel(~Par1 + Par2, Block, Env = Env, fct = fct)
  #print(head(X))
  z <- lm.fit(X, Y)
  if(ML == T){
    if(fct == "HAYMAN1") {
      if(is.null(Block)) {
        res <- MLfit1(X, Y, z, Par1, Par2)
      } else {
        X2 <- model.matrixDiallel(~Par1 + Par2, Block = Block,
                                   fct = fct, ML = T)
        res <- MLfit1b(X2, Y, z, Par1, Par2)
      }
    } else if (fct == "HAYMAN2"){
      res <- MLfit2(X, Y, z, Par1, Par2)
    }
    z$coefficients <- res$coefficients
    z$residuals <- res$residuals
    z$fitted.values <- res$fitted
    z$assign <- res$assign
    z$MLcoefficients <- res$MLcoefficients
    z$LScoefficients <- res$LScoefficients
    X <- res$X
    }
  class(z) <- c(if (mlm) "mlm", "lm")
  z$response <- Y

  z$fct <- fct
  z$Env <- ifelse(is.null(Env), F, T)
  z$na.action <- attr(mf, "na.action")
  z$offset <- NULL
    #z$contrasts <- attr(x, "contrasts")
  z$xlevels <- .getXlevels(mt, mf)
  z$call <- cl
  z$terms <- mt
  z$model <- mf
  z$namEff <- attr(X, "namEff")
  z$modMatrix <- X
  z$ML <- ML

  # if (ret.x)
  #       z$x <- x
  # if (ret.y)
  #       z$y <- y
  # if (!qr)
  #       z$qr <- NULL
  class(z) <- c("diallel", "lm")
  return(z)
}


summary.diallel <- function (object, correlation = FALSE, symbolic.cor = FALSE,
    ...)
{ #print(is.null(object$MLcoefficients))
  if(object$ML == T){
  tab <- object$LScoefficients
  tab$"t value" <- tab[,1]/tab[,2]
  tab$"Pr(>|t|)" <- 2 * pt(abs(tab$"t value"), object$df.residual, lower.tail = F)
  return(tab)
  }else{
  z <- object
  class(z) <- "lm"
  summary(z)
  }
}

vcov.diallel <- function(object, ...)
{
    so <- summary(object)
    so$sigma^2 * so$cov.unscaled
}

anova.diallel <- function(object, MSE = NULL, dfr = NULL, ...)
{
    if(object$ML == F & is.null(MSE) & object$Env == F
       & object$fct != "HAYMAN1" & object$fct != "HAYMAN2"){
      ## Analisi con blocco: uso dell'errore residuo ####
      ## Do not copy this: anova.lmlist is not an exported object.
      ## See anova.glm for further comments.
    if(length(list(object, ...)) > 1L) return(anova.lmlist(object, ...))
    # object <- fit
    #if(!inherits(object, "lm"))
	  #warning("calling anova.lm(<fake-lm-object>) ...")
    w <- object$weights
    ssr <- sum(if(is.null(w)) object$residuals^2 else w*object$residuals^2)
    mss <- sum(if(is.null(w)) object$fitted.values^2 else w*object$fitted.values^2)
    if(ssr < 1e-10*mss)
        warning("ANOVA F-tests on an essentially perfect fit are unreliable")
    dfr <- df.residual(object)
    p <- object$rank
    if(p > 0L) {
        p1 <- 1L:p
        comp <- object$effects[p1]
        asgn <- object$assign[object$qr$pivot][p1]
        nmeffects <- c("(Intercept)", attr(object$terms, "term.labels"))
        tlabels <- object$namEff
        ss <- c(unlist(lapply(split(comp^2,asgn), sum)), ssr)
        df <- c(lengths(split(asgn,  asgn)), dfr)
    } else {
        ss <- ssr
        df <- dfr
        tlabels <- character()
    }
    ms <- ss/df
    f <- ms/(ssr/dfr)
    P <- pf(f, df, dfr, lower.tail = FALSE)
    table <- data.frame(df, ss, ms, f, P)
    table[length(P), 4:5] <- NA

    dimnames(table) <- list(c(tlabels),
                            c("Df","Sum Sq", "Mean Sq", "F value", "Pr(>F)"))
    if(attr(object$terms,"intercept")) table <- table[-1, ]
    structure(table, heading = c("Analysis of Variance Table\n",
		     paste("Response:", deparse(formula(object)[[2L]]))),
	      class = c("anova", "data.frame"))# was "tabular"

    }else if(object$Env == F & (object$fct == "GE2" | object$fct == "GE2r")) {
       ## Analisi senza blocco, per GE2 e GE2r ####
       ## Deve ricalcolare in modo diverso
       ssr <- sum(object$residuals^2)
       if(!is.null(MSE)) dfr1 <- dfr
    #print(dfr1)
    dfr <- df.residual(object)
    #if(!is.null(MSE)) ssr1 <- MSE * dfr1
    mss <- sum(object$fitted.values^2)
    p <- object$rank
    p1 <- 1L:p
    comp <- object$effects[p1]
    asgn <- object$assign[object$qr$pivot][p1]
    nmeffects <- c("(Intercept)", attr(object$terms, "term.labels"))
    tlabels <- object$namEff
    ss <- c(unlist(lapply(split(comp^2,asgn), sum)), ssr)
    df <- c(lengths(split(asgn,  asgn)), dfr)
    ms <- ss/df

    table <- data.frame(df, ss, ms)#, f, P)
    row.names(table) <- fit$namEff
    if(!is.null(MSE)){
        f <- ms/MSE
        P <- pf(f, df, dfr1, lower.tail = FALSE)
        table <- data.frame(df, ss, ms, dfr1, f, P)
        #table[length(P), 4:5] <- NA
        tlabels <- object$namEff
        dimnames(table) <- list(c(tlabels),
                                c("Df","Sum Sq", "Mean Sq", "Den df", "F value", "Pr(>F)"))
        if(attr(object$terms,"intercept")) table <- table[-1, ]
        structure(table, heading = c("Analysis of Variance Table\n",
                                     paste("Response:", deparse(formula(object)[[2L]]))),
                  class = c("anova", "data.frame"))# was "tabular"
      }
      table
    }else if(object$Env == T) {
      # Analisi con anno ############
      X <- object$modMatrix
      Y <- object$response
      namEff <- object$namEff
      numEff <- length(namEff) - 2
      namEff <- namEff[-1]
      asgn <- attr(X, "assign")
      fct <- object$fct
      dataset <- object$model
      names(dataset)[4:5] <- c("Block", "Env")
      dataset$Block <- factor(dataset$Block)
      dataset$Env <- factor(dataset$Env)
      matsOr <- model.matrix.diallel(~dataset[,2]+dataset[,3],
                           dataset$Block,
                           fct = fct)
      asgn2 <- attr(matsOr, "assign")
      ss <- c()
      dfr <- c(); labTab <- c()
      rss <- sum(object$residuals^2)
      ss[1] <- rss
      resdf <- object$df.residual
      dfr[1] <- resdf
      labTab[1] <- "Residuals"
      cont <- 2

      for(i in numEff:1){
      # Model with common effect (G:E)
      sel <- asgn2 == i
      sel2 <- (asgn >=( i + numEff) & asgn <= 2*numEff) | (asgn >= i & asgn <= numEff)
      sel2 <- ifelse(sel2==T, F, T)
      df <- matsOr[,sel]
      X2 <- as.matrix( cbind(X[, sel2], df) )
      if(i != 1){
        reg2 <- lm.fit(X2, Y)
        ssGE <- sum(reg2$residuals^2)
        dfGE <- reg2$df.residual
        ss[cont] <- ssGE; dfr[cont] <- dfGE
        labTab[cont] <- paste(namEff[i], "Env", sep = ":")
        cont <- cont + 1
        }

      # Model with no effects
      X3 <- as.matrix( X[, sel2] )
      reg3 <- lm.fit(X3, Y)
      ssG <- sum(reg3$residuals^2)
      dfG <- reg3$df.residual
      ss[cont] <- ssG; dfr[cont] <- dfG; labTab[cont] <- namEff[i]
      cont <- cont + 1
      }
      reg.null <- lm(Y ~ 1)
      totss <- deviance(reg.null)
      totdf <- reg.null$df.residual
      ss[cont] <- totss; dfr[cont] <- totdf
      ss <- diff(ss); dfr <- diff(dfr)
      ss <- c(rev(ss), rss); dfr <- c(rev(dfr), resdf)
      labTab <- c("Environment", rev(labTab))
      ms <- ss/dfr
      MSE <- ms[length(ms)]
      dfr1 <- dfr[length(dfr)]
      f <- ms/MSE
      P <- pf(f, dfr, dfr1, lower.tail = FALSE)
      table <- data.frame(dfr, ss, ms, f, P)
      colnames(table) <- c("Df","Sum Sq", "Mean Sq", "F value", "Pr(>F)")
      row.names(table) <- labTab
      structure(table, heading = c("Analysis of Variance Table\n",
                                     paste("Response:", deparse(formula(object)[[2L]]))),
                  class = c("anova", "data.frame"))# was "tabular"
      table
      } else {
      # Hyman #############################
      rss <- c()
      fit <- object
      asgn <- fit$assign
      X <- fit$modMatrix
      coefs <- fit$coefficients
      y <- fit$response
      ngroup <- length(levels(factor(fit$assign)))
      rss[1] <- deviance(lm(y ~ 1))
      for(i in 1:ngroup){
        #i <- 2
        val <- fit$assign <= i
        exp <- X[,val] %*% as.matrix(coefs[val])
        res <- y - exp
        rss[i+1] <- sum(res^2)
      }
      rss[i + 1] <- 0
      ss <- rss[1:ngroup] - rss[2:(ngroup+1)]
      ss <- c(rss[1], ss)
      df <- c(lengths(split(asgn,  asgn)), fit$df.residual)
      if(object$fct == "HAYMAN1" & object$ML == T) {
        df[3] <- df[3]/2
      } else if(object$fct == "HAYMAN2" & object$ML == T){
        df[6] <- df[6]/2
      }
      ms <- ss/df
      table <- data.frame(df, ss, ms)#, f, P)
      row.names(table) <- fit$namEff
      if(!is.null(MSE)){
        f <- ms/MSE
        P <- pf(f, df, dfr, lower.tail = FALSE)
        table <- data.frame(df, ss, ms, dfr, f, P)
        #table[length(P), 4:5] <- NA
        tlabels <- object$namEff
        dimnames(table) <- list(c(tlabels),
                                c("Df","Sum Sq", "Mean Sq", "Den df", "F value", "Pr(>F)"))
        if(attr(object$terms,"intercept")) table <- table[-1, ]
        structure(table, heading = c("Analysis of Variance Table\n",
                                     paste("Response:", deparse(formula(object)[[2L]]))),
                  class = c("anova", "data.frame"))# was "tabular"

      }
      table
      }
    }

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




