bugs.diallel <- function(formula, Block = NULL, Env = NULL,
                       fct = "GRIFFING2", data,
                       random = "fixed", n.iter = 10000,
                       burn.in = 1000){
  
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
  if(missing(data) == T){
    Par1 <- mf[,2]
    Par2 <- mf[,3]
  } else {
    Par1 <- data[[pars[1]]]
    Par2 <- data[[pars[2]]]
  }
  
  # Create model matrix
  X <- model.matrixDiallel(~Par1 + Par2, Block=Block, Env = Env, fct = fct)
  
  if(random == "fixed"){
    #Assegnazione valori iniziali
    start <- lm.diallel(Y ~ Par1 + Par2, Block = Block,
                  fct = fct)

    sigmaInit <- summary(start)$sigma
    betaInit <- coef(start)
    np <- length(betaInit); n <- length(Y)

    #Creazione delle liste per WINBUGS
    dataMod <- list(Y = Y, X = X)
    init <- list(beta = betaInit, sigma = sigmaInit)

    #Lancio del campionatore
    mcmc <- jags.model("mod_diallel.fix.txt", 
                   data = dataMod, inits = init, 
                   n.chains = 4, n.adapt = 100)

    params <- c("beta", "sigma") 
  
    } else if(random == "design" & is.null(Env)){
    if(is.null(Block)){print("ERROR: Block variable is needed"); stop()}  
    
    # Get the fixed effects matrix  
    X <- model.matrixDiallel(~Par1 + Par2,
                         fct = fct)
    B <- factor(Block)
    contrasts(B) <- "contr.sum"
    Z <- matrix(model.matrix(~B)[,-1]) 

    # Get starting values (method of moments)
    start <- lm.diallel(Y ~ Par1 + Par2, Block=Block,
                         fct = "HAYMAN1")
    aovMat <- anova(start)
    MSe <- summary(start)$sigma^2
    MSb <- aovMat[1,3]

    n <- mean( tapply(Y, Block, length) ) 
    sigma2.b <- max(c( (MSb - MSe)/n, 0.0001))
    sigmaInit <- summary(start)$sigma
    sigmab <- sqrt(sigma2.b)
    betaInit <- coef(start)[-2]
    
    # Create input for WinBugs
    dataMod <- list(Y = Y, X = X, Z = matrix(Z))
    init <- list(beta = betaInit, sigma = sigmaInit, 
             sigma.r = sigmab)

    # Start sampler
    mcmc <- jags.model("mod_diallel.ran.txt", 
                   data = dataMod, inits = init, 
                   n.chains = 4, n.adapt = 100)

    params <- c("beta", "sigma", "sigma.r", "sigma2", "sigma2.r")
  }
  

  res <- coda.samples(mcmc, params, n.iter = n.iter)
  burn.in <- burn.in
  out <- summary(window(res, start=burn.in))
  out
}
