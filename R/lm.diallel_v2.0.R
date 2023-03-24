# Wrapper for lm and diallel models
# Last edited: 20/03/2023
# Added check of scheme and balance
lm.diallel <- function(formula, Block = NULL, Env = NULL,
                       fct = "GRIFFING2", data){
  # formula Ftime + Block + Par1 + Par2, data = df
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
    Par1 <- mf[,2]
    Par2 <- mf[,3]
  }

  # Before starting, the mating scheme and balance is checked
  # Updated on 20/03/2023
  Par1 <- factor(Par1)
  Par2 <- factor(Par2)
  chk_design <- checkScheme(Par1, Par2)
  if(!is.null(chk_design$missingCrosses) & (fct != "GRIFFING4" | chk_design$matingScheme != 4)){
    # Se ci sono missing crosses, lavora solo con GRIFFING4
    stop("Missing cells (crosses/selfs) were detected. In this condition, only 'GRIFFING4' model (with mating scheme 4) is allowed. Enhancements are on their way.")
  }
  if(!is.null(chk_design$missingCrosses) & fct == "GRIFFING4" & length(chk_design$parNoMis) == 0){
    # Nessun parent non ha missing crosses
    stop("Missing crosses can be handled with GRIFFING4 only when at least one parent has no missing crosses!")
  }

  if(chk_design$matingScheme == 1 & (fct == "GRIFFING2" | fct == "GRIFFING3" | fct == "GRIFFING3") ) {
    wrnmsg <- paste("Mating scheme", chk_design$matingScheme, "has been detected. The model you selected may not be appropriate for this scheme",
                    sep = " ")
    warning(wrnmsg, call. = FALSE)
  }
  if((fct == "HAYMAN1" | fct == "HAYMAN2") & chk_design$matingScheme != 1){
    msg <- paste("Mating scheme ", chk_design$matingScheme, " detected. ", sep = "")
    stop(paste(msg, "The HAYMAN1 and HAYMAN2 models can only be fitted to full diallel experiments (mating scheme 1)", sep = ""))
  }
  if((fct == "GE2" | fct == "GE3") & chk_design$matingScheme != 2){
    msg <- paste("Mating scheme ", chk_design$matingScheme, " detected. ", sep = "")
    stop(paste(msg, "The GE2 and GE3 models can only be fitted to diallel experiments with mating scheme 2", sep = ""))
  }
  if((fct == "GE2r" | fct == "GE3r") & chk_design$matingScheme != 1){
      msg <- paste("Mating scheme ", chk_design$matingScheme, " detected. ", sep = "")
      stop("The GE2r and GE3r models can only be fitted to diallel experiments with mating scheme 1")
  }


  X <- model.matrixDiallel(~ Par1 + Par2, Block=Block, Env = Env,
                           fct = fct, type = "nested")
  z <- lm.fit(X, Y)
  class(z) <- c(if (mlm) "mlm", "lm")
  z$response <- Y
  z$fct <- fct
  z$Env <- ifelse(is.null(Env), F, T)
  z$Block <- ifelse(is.null(Block), F, T)
  z$na.action <- attr(mf, "na.action")
  z$offset <- NULL

  z$xlevels <- .getXlevels(mt, mf)
  z$call <- cl
  z$terms <- mt
  z$model <- mf
  z$namEff <- attr(X, "namEff")
  z$modMatrix <- X
  z$matingScheme <- chk_design$matingScheme
  z$chk_design <- chk_design
  class(z) <- c("diallel", "lm")
  return(z)
}


summary.diallel <- function (object, # correlation = FALSE, symbolic.cor = FALSE,
                             MSE = NULL, dfr = NULL, ...)
{ #print(is.null(object$MLcoefficients))
  if(!is.null(MSE) & !is.null(dfr)){
  sigma <- sqrt(MSE)
  if(any(class(object) == "diallel") == T) {X <- object$modMatrix
  } else { X <- model.matrix(object)}
  ses <- sqrt(diag(solve( t(X) %*% X ))) * sigma
  tab <- data.frame("Estimate" = object$coef, "SE" = ses)
  tab$"t value" <- tab[,1]/tab[,2]
  tab$"Pr(>|t|)" <- 2 * pt(abs(tab$"t value"), dfr, lower.tail = F)
  return(tab)
  }else{
  z <- object
  class(z) <- "lm"
  summary(z)
  }
}

vcov.diallel <- function(object, MSE = NULL, ...)
{
    so <- summary(object)
    if(is.na(so$sigma) & is.null(MSE)){
      cat("No residual variance estimate is available")
      stop() #retVal <- NA
    } else if(is.na(so$sigma) == T & is.null(MSE) == F){
      retVal <- MSE * so$cov.unscaled
    } else if(is.na(so$sigma) == F & is.null(MSE) == F){
      retVal <- MSE * so$cov.unscaled
    } else if(is.na(so$sigma) == F & is.null(MSE) == T){
      retVal <- so$sigma^2 * so$cov.unscaled
    }
    retVal
}

predict.diallel <- function(object, ...){
  fitted(object)
}



