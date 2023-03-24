anova.diallel <- function(object, MSE = NULL, dfr = NULL, type = 1, ...){
  # The problem is the parameterisation of the model with environments,
  # that is not fully correct, because it is nested. Therefore, we need
  # a preliminary step.
  if(is.null(object$Env)) {object$Env <- FALSE }
  if(object$Env == T){
    # for METs, I have to recover the elements from model call
    # Define a crossed incidence matrix and refit the model
    # This part only works for diallel objects, where we have an
    # Env argument
    cl <- getCall(object)
    tmp <- as.data.frame(eval(cl$data))
    Block <- tmp[, colnames(tmp) %in% as.character(cl$Block)]
    Env <- tmp[, colnames(tmp) %in% as.character(cl$Env)]
    Par1 <- tmp[,all.vars(cl$formula)[2]]
    Par2 <- tmp[,all.vars(cl$formula)[3]]
    Y <- tmp[,all.vars(cl$formula)[1]]
    X <- model.matrixDiallel.MET(Par1, Par2, Block, Env, fct = object$fct,
                                 type = "crossed")
    z <- lm.fit(X, Y)
    mlm <- is.matrix(Y)
    class(z) <- c(if (mlm) "mlm", "lm")
    z$response <- Y
    z$fct <- object$fct
    z$Env <- ifelse(is.null(Env), F, T)
    z$Block <- ifelse(is.null(Block), F, T)

    # z$na.action <- attr(mf, "na.action")
    z$offset <- NULL
    z$xlevels <- object$xlevels
    z$call <- cl
    z$terms <- object$terms
    z$model <- object$terms
    # attr(object$terms, "term.labels") Da manipolare!!!!!!!!
    z$namEff <- attr(X, "namEff")# [-c(1, length(attr(X, "namEff")))]
    z$modMatrix <- X
    z$matingScheme <- object$matingScheme
    z$chk_design <- object$chk_design
    class(z) <- c("diallel", "lm")
    object <- z
    rm(z)
  } else {
    X <- model.matrix(object)
  }

  ## Se non trova MSE esplicito: uso del residuo come errore
  # if(length(list(object, ...)) > 1L) return(anova.lmlist(object, ...))
  # object <- fit
  #if(!inherits(object, "lm"))
  #warning("calling anova.lm(<fake-lm-object>) ...")
  w <- object$weights
  ssr <- sum(if(is.null(w)) object$residuals^2 else w*object$residuals^2)
  mss <- sum(if(is.null(w)) object$fitted.values^2 else w*object$fitted.values^2)
  if(ssr < 1e-10*mss & is.null(MSE))
      warning("ANOVA F-tests on an essentially perfect fit are unreliable")
  if( is.null(dfr) ) dfr <- df.residual(object)
  p <- object$rank
  p1 <- 1L:p
  comp <- object$effects[p1]
  asgn <- object$assign[object$qr$pivot][p1]
  nmeffects <- c("(Intercept)", attr(object$terms, "term.labels"))
  if(any(class(object) == "diallel") == T) {
    tlabels <- object$namEff
  } else {
    tlabels <- nmeffects[1 + unique(asgn)]
    # nmeffects
  }
  # print(tlabels)
  if(type == 1 | type == "I"){
     # Type I Sum of Squares

    ss <- c(unlist(lapply(split(comp^2,asgn), sum)), ssr)
    df <- c(lengths(split(asgn,  asgn)), dfr)
    # } else {
    #     ss <- ssr
    #     df <- dfr
    #     tlabels <- character()
    # }
    if(is.null(MSE)){
      ms <- ss/df
      f <- ms/(ssr/dfr)
      P <- pf(f, df, dfr, lower.tail = FALSE)
      # P[length(P)] <- NULL
      # f[length(f)] <- NULL
      table <- data.frame(df, ss, ms, f, P)
      table[length(P), 4:5] <- NA
    } else {
      ms <- ss/df
      f <- ms/MSE
      P <- pf(f, df, dfr, lower.tail = FALSE)
      # P[length(P)] <- NULL
      # f[length(f)] <- NULL
      table <- data.frame(df, ss, ms, f, P)
      table[length(P), 2:5] <- NA
      table[length(P), 3] <- MSE
    }
  } else if (type == 3 | type == "III"){
      # Type III Sum of Squares
      tlabels <- tlabels[-1]
      ss <- c()
      df <- c()

      if(!any(inherits(object, "diallel"))){
        num.eff <- length(tlabels)
        # X <- model.matrix(object)
        cl <- getCall(object)
        tmp <- as.data.frame(eval.parent(cl$data, n = 3))
        Y <- tmp[,all.vars(cl$formula)[1]]
      } else {
        num.eff <- length(object$namEff) - 2
        Y <- object$response
      }
      mod <- lm(Y ~ X - 1)
      for(i in 0:num.eff){
        Xred <- X[,attr(X, "assign") != i]
        modRed <- lm(Y ~ Xred - 1)
        tab <- anova(modRed, mod)
        ss[i] <- tab$`Sum of Sq`[2]
        df[i] <- tab$Df[2]
      }
      ss <- c(ss, ssr)
      df <- c(df, dfr)
      if(is.null(MSE)){
        ms <- ss/df
        f <- ms/(ssr/dfr)
        P <- pf(f, df, dfr, lower.tail = FALSE)
        table <- data.frame(df, ss, ms, f, P)
        table[length(P), 4:5] <- NA
      } else {
        ms <- ss/df
        f <- ms/MSE
        P <- pf(f, df, dfr, lower.tail = FALSE)
        # P[length(P)] <- NULL
        # f[length(f)] <- NULL
        table <- data.frame(df, ss, ms, f, P)
        table[length(P), 2:5] <- NA
        table[length(P), 3] <- MSE
    }
  }
  if(all(class(object) != "diallel")) tlabels <- c(tlabels, "Residuals")
  dimnames(table) <- list(c(tlabels),
                          c("Df","Sum Sq", "Mean Sq", "F value", "Pr(>F)"))
  if(attr(object$terms,"intercept") & (type == 1 | type == "I")) table <- table[-1, ]
  # structure(table, heading = c(paste("Analysis of Variance Table\n", "Mating Scheme: ",
  #                                    object$matingScheme, "\n",
  #                                        "Diallel model:", object$fct, sep = ""),
  structure(table, heading = c("Analysis of Variance Table\n",
       paste("Response:", deparse(formula(object)[[2L]]))),
     class = c("anova", "data.frame"))# was "tabular"
  class(table) <- c("anova", "data.frame")
  print(table)
  invisible(table)
  }


  # } else if(object$Env == F &
  #           (object$fct == "GE2" | object$fct == "GE2r")) {
  # ## Analisi senza blocco, per GE2 e GE2r
  # ## Deve ricalcolare in modo diverso
  # ssr <- sum(object$residuals^2)
  # if(!is.null(MSE)) dfr1 <- dfr
  # dfr <- df.residual(object)
  # mss <- sum(object$fitted.values^2)
  # p <- object$rank
  # p1 <- 1L:p
  # comp <- object$effects[p1]
  # asgn <- object$assign[object$qr$pivot][p1]
  # nmeffects <- c("(Intercept)", attr(object$terms, "term.labels"))
  # tlabels <- object$namEff
  # ss <- c(unlist(lapply(split(comp^2,asgn), sum)), ssr)
  # df <- c(lengths(split(asgn,  asgn)), dfr)
  # ms <- ss/df
  # table <- data.frame(df, ss, ms)#, f, P)
  # row.names(table) <- fit$namEff
  # if(!is.null(MSE)){
  #     f <- ms/MSE
  #     P <- pf(f, df, dfr1, lower.tail = FALSE)
  #     table <- data.frame(df, ss, ms, dfr1, f, P)
  #     tlabels <- object$namEff
  #     dimnames(table) <- list(c(tlabels),
  #                             c("Df","Sum Sq", "Mean Sq", "Den df", "F value", "Pr(>F)"))
  #     if(attr(object$terms,"intercept")) table <- table[-1, ]
  #     structure(table, heading = c("Analysis of Variance Table\n",
  #                                  paste("Response:", deparse(formula(object)[[2L]]))),
  #               class = c("anova", "data.frame"))# was "tabular"
  #   }
  #   #table
  #
  # } else if(object$Env == T) {
  #   # Analisi poliennale: deve rifittare il modello
  #   # con parametrizzazione crossed (not nested)
  #   MSEor <- MSE; dfrOr <- dfr
  #   X <- object$modMatrix
  #   Y <- object$response
  #   namEff <- object$namEff
  #   numEff <- length(namEff) - 2
  #   namEff <- namEff[-1]
  #   asgn <- attr(X, "assign")
  #   fct <- object$fct
  #   dataset <- object$model
  #
  #   # Crea una matrice nella quale sono confusi gli ambienti
  #   if(object$Block == T){
  #     names(dataset)[4:5] <- c("Block", "Env")
  #     dataset$Block <- factor(dataset$Block)
  #     dataset$Env <- factor(dataset$Env)
  #     matsOr <- model.matrixDiallel(~dataset[,2] + dataset[,3],
  #                        dataset$Block,
  #                        fct = fct)
  #   } else {
  #     names(dataset)[4] <- c("Env")
  #     dataset$Env <- factor(dataset$Env)
  #     matsOr <- model.matrixDiallel(~dataset[,2] + dataset[,3],
  #                                   fct = fct)
  #   }
  #
  #   asgn2 <- attr(matsOr, "assign")
  #   resdf <- object$df.residual
  #   ss <- c()
  #   dfr <- c(); labTab <- c()
  #   rss <- sum(object$residuals^2)
  #   ss[1] <- rss
  #   dfr[1] <- resdf
  #   labTab[1] <- "Residuals"
  #   cont <- 2
  #
  #   for(i in numEff:1){
  #   # Model with common effect (G:E)
  #   sel <- asgn2 == i
  #   sel2 <- (asgn >=( i + numEff) & asgn <= 2*numEff) | (asgn >= i & asgn <= numEff)
  #   sel2 <- ifelse(sel2==T, F, T)
  #   df <- matsOr[,sel]
  #   X2 <- as.matrix( cbind(X[, sel2], df) )
  #   if(i != 1){
  #     reg2 <- lm.fit(X2, Y)
  #     ssGE <- sum(reg2$residuals^2)
  #     dfGE <- reg2$df.residual
  #     ss[cont] <- ssGE; dfr[cont] <- dfGE
  #     labTab[cont] <- paste(namEff[i], "Env", sep = ":")
  #     cont <- cont + 1
  #   } else {
  #       if(object$Block == F) {
  #         reg2 <- lm.fit(X2, Y)
  #         ssGE <- sum(reg2$residuals^2)
  #         dfGE <- reg2$df.residual
  #         ss[cont] <- ssGE; dfr[cont] <- dfGE
  #         labTab[cont] <- paste(namEff[i], "Env", sep = ":")
  #         cont <- cont + 1
  #       }
  #     }
  #
  #   # Model with no effects
  #   X3 <- as.matrix( X[, sel2] )
  #   reg3 <- lm.fit(X3, Y)
  #   ssG <- sum(reg3$residuals^2)
  #   dfG <- reg3$df.residual
  #   ss[cont] <- ssG; dfr[cont] <- dfG; labTab[cont] <- namEff[i]
  #   cont <- cont + 1
  #   }
  #   reg.null <- lm(Y ~ 1)
  #   totss <- deviance(reg.null)
  #   totdf <- reg.null$df.residual
  #   ss[cont] <- totss; dfr[cont] <- totdf
  #   ss <- diff(ss); dfr <- diff(dfr)
  #   ss <- c(rev(ss), rss); dfr <- c(rev(dfr), resdf)
  #   labTab <- c("Environment", rev(labTab))
  #   labTab[labTab=="Block"] <- "Env:Block"
  #   ms <- ss/dfr
  #   if(is.null(MSEor)){
  #     MSE <- ms[length(ms)]
  #   } else {
  #     MSE <- MSEor
  #     ms[length(ms)] <- MSE
  #   }
  #   if(!is.null(dfrOr)){
  #     dfr1 <- dfrOr
  #     dfr[length(dfr)] <- dfrOr
  #     ss[length(ss)] <- NA
  #   } else {
  #     dfr1 <- dfr[length(dfr)]
  #   }
  #
  #   f <- ms/MSE
  #   P <- pf(f, dfr, dfr1, lower.tail = FALSE)
  #   table <- data.frame(dfr, ss, ms, f, P)
  #   table[length(P), 4:5] <- NA
  #   colnames(table) <- c("Df","Sum Sq", "Mean Sq", "F value", "Pr(>F)")
  #   row.names(table) <- labTab
  #   structure(table, heading = c("Analysis of Variance Table\n",
  #                                  paste("Response:", deparse(formula(object)[[2L]]))),
  #               class = c("anova", "data.frame"))# was "tabular"
  #
  #   } else {
  #   # Non uso il residuo come errore,
  #   ## ma quello fornito
  #   # print("OK")
  #   rss <- c()
  #   fit <- object
  #   asgn <- fit$assign
  #   X <- fit$modMatrix
  #   coefs <- fit$coefficients
  #   y <- fit$response
  #   ngroup <- length(levels(factor(fit$assign)))
  #   rss[1] <- deviance(lm(y ~ 1))
  #   for(i in 1:ngroup){
  #     #i <- 2
  #     val <- fit$assign <= i
  #     exp <- X[,val] %*% as.matrix(coefs[val])
  #     res <- y - exp
  #     rss[i+1] <- sum(res^2)
  #   }
  #   rss[i + 1] <- 0
  #   ss <- rss[1:ngroup] - rss[2:(ngroup+1)]
  #   ss <- c(rss[1], ss)
  #   df <- c(lengths(split(asgn,  asgn)), fit$df.residual)
  #   tlabels <- object$namEff
  #
  #   if(fit$df.residual == 0){
  #      ss <- ss[-length(ss)]
  #      df <- df[-length(df)]
  #      tlabels <- object$namEff[-length(object$namEff)]
  #   }
  #   #   df[3] <- df[3]/2
  #   # } else if(object$fct == "HAYMAN2" & object$ML == T){
  #   #   df[6] <- df[6]/2
  #   # }
  #   ms <- ss/df
  #   table <- data.frame(df, ss, ms)#, f, P)
  #   #row.names(table) <- fit$namEff
  #   #if(!is.null(MSE)){
  #   f <- ms/MSE
  #   P <- pf(f, df, dfr, lower.tail = FALSE)
  #   table <- data.frame(df, ss, ms, dfr, f, P)
  #
  #   dimnames(table) <- list(c(tlabels),
  #                             c("Df","Sum Sq", "Mean Sq", "Den df", "F value", "Pr(>F)"))
  #     if(attr(object$terms,"intercept")) table <- table[-1, ]
  #     structure(table, heading = c("Analysis of Variance Table\n",
  #                                  paste("Response:", deparse(formula(object)[[2L]]))),
  #               class = c("anova", "data.frame"))# was "tabular"
  #
  #   #}
      # table
      # }

# old.anova.diallel <- function(object, MSE = NULL, dfr = NULL, ...){
#   # Different sections
#   if(is.null(object$Env)) {object$Env <- FALSE }
#   if(object$Env == F){
#     ## Se non trova MSE esplicito: uso del residuo come errore
#     # if(length(list(object, ...)) > 1L) return(anova.lmlist(object, ...))
#     # object <- fit
#     #if(!inherits(object, "lm"))
# 	  #warning("calling anova.lm(<fake-lm-object>) ...")
#     w <- object$weights
#     ssr <- sum(if(is.null(w)) object$residuals^2 else w*object$residuals^2)
#     mss <- sum(if(is.null(w)) object$fitted.values^2 else w*object$fitted.values^2)
#     if(ssr < 1e-10*mss & is.null(MSE))
#         warning("ANOVA F-tests on an essentially perfect fit are unreliable")
#     if( is.null(dfr) ) dfr <- df.residual(object)
#     p <- object$rank
#     if(p > 0L) {
#       # Type I Sum of Squares
#         p1 <- 1L:p
#         comp <- object$effects[p1]
#         asgn <- object$assign[object$qr$pivot][p1]
#         nmeffects <- c("(Intercept)", attr(object$terms, "term.labels"))
#         if(any(class(object) == "diallel") == T) {tlabels <- object$namEff
#         } else {tlabels <- nmeffects[1 + unique(asgn)]
#         nmeffects}
#         ss <- c(unlist(lapply(split(comp^2,asgn), sum)), ssr)
#         df <- c(lengths(split(asgn,  asgn)), dfr)
#     } else {
#         ss <- ssr
#         df <- dfr
#         tlabels <- character()
#     }
#     if(is.null(MSE)){
#       ms <- ss/df
#       f <- ms/(ssr/dfr)
#       P <- pf(f, df, dfr, lower.tail = FALSE)
#       # P[length(P)] <- NULL
#       # f[length(f)] <- NULL
#       table <- data.frame(df, ss, ms, f, P)
#       table[length(P), 3:5] <- NA
#     } else {
#       ms <- ss/df
#       f <- ms/MSE
#       P <- pf(f, df, dfr, lower.tail = FALSE)
#       # P[length(P)] <- NULL
#       # f[length(f)] <- NULL
#       table <- data.frame(df, ss, ms, f, P)
#       table[length(P), 2:5] <- NA
#       table[length(P), 3] <- MSE
#     }
#     if(all(class(object) != "diallel")) tlabels <- c(tlabels, "Residuals")
#     dimnames(table) <- list(c(tlabels),
#                             c("Df","Sum Sq", "Mean Sq", "F value", "Pr(>F)"))
#     if(attr(object$terms,"intercept")) table <- table[-1, ]
#     # structure(table, heading = c(paste("Analysis of Variance Table\n", "Mating Scheme: ",
#     #                                    object$matingScheme, "\n",
#     #                                        "Diallel model:", object$fct, sep = ""),
# 		structure(table, heading = c("Analysis of Variance Table\n",
#          paste("Response:", deparse(formula(object)[[2L]]))),
# 	     class = c("anova", "data.frame"))# was "tabular"
#
#     } else if(object$Env == F &
#               (object$fct == "GE2" | object$fct == "GE2r")) {
#     ## Analisi senza blocco, per GE2 e GE2r
#     ## Deve ricalcolare in modo diverso
#     ssr <- sum(object$residuals^2)
#     if(!is.null(MSE)) dfr1 <- dfr
#     dfr <- df.residual(object)
#     mss <- sum(object$fitted.values^2)
#     p <- object$rank
#     p1 <- 1L:p
#     comp <- object$effects[p1]
#     asgn <- object$assign[object$qr$pivot][p1]
#     nmeffects <- c("(Intercept)", attr(object$terms, "term.labels"))
#     tlabels <- object$namEff
#     ss <- c(unlist(lapply(split(comp^2,asgn), sum)), ssr)
#     df <- c(lengths(split(asgn,  asgn)), dfr)
#     ms <- ss/df
#     table <- data.frame(df, ss, ms)#, f, P)
#     row.names(table) <- fit$namEff
#     if(!is.null(MSE)){
#         f <- ms/MSE
#         P <- pf(f, df, dfr1, lower.tail = FALSE)
#         table <- data.frame(df, ss, ms, dfr1, f, P)
#         tlabels <- object$namEff
#         dimnames(table) <- list(c(tlabels),
#                                 c("Df","Sum Sq", "Mean Sq", "Den df", "F value", "Pr(>F)"))
#         if(attr(object$terms,"intercept")) table <- table[-1, ]
#         structure(table, heading = c("Analysis of Variance Table\n",
#                                      paste("Response:", deparse(formula(object)[[2L]]))),
#                   class = c("anova", "data.frame"))# was "tabular"
#       }
#       #table
#
#     } else if(object$Env == T) {
#       # Analisi poliennale: deve rifittare il modello
#       # con parametrizzazione crossed (not nested)
#       MSEor <- MSE; dfrOr <- dfr
#       X <- object$modMatrix
#       Y <- object$response
#       namEff <- object$namEff
#       numEff <- length(namEff) - 2
#       namEff <- namEff[-1]
#       asgn <- attr(X, "assign")
#       fct <- object$fct
#       dataset <- object$model
#
#       # Crea una matrice nella quale sono confusi gli ambienti
#       if(object$Block == T){
#         names(dataset)[4:5] <- c("Block", "Env")
#         dataset$Block <- factor(dataset$Block)
#         dataset$Env <- factor(dataset$Env)
#         matsOr <- model.matrixDiallel(~dataset[,2] + dataset[,3],
#                            dataset$Block,
#                            fct = fct)
#       } else {
#         names(dataset)[4] <- c("Env")
#         dataset$Env <- factor(dataset$Env)
#         matsOr <- model.matrixDiallel(~dataset[,2] + dataset[,3],
#                                       fct = fct)
#       }
#
#       asgn2 <- attr(matsOr, "assign")
#       resdf <- object$df.residual
#       ss <- c()
#       dfr <- c(); labTab <- c()
#       rss <- sum(object$residuals^2)
#       ss[1] <- rss
#       dfr[1] <- resdf
#       labTab[1] <- "Residuals"
#       cont <- 2
#
#       for(i in numEff:1){
#       # Model with common effect (G:E)
#       sel <- asgn2 == i
#       sel2 <- (asgn >=( i + numEff) & asgn <= 2*numEff) | (asgn >= i & asgn <= numEff)
#       sel2 <- ifelse(sel2==T, F, T)
#       df <- matsOr[,sel]
#       X2 <- as.matrix( cbind(X[, sel2], df) )
#       if(i != 1){
#         reg2 <- lm.fit(X2, Y)
#         ssGE <- sum(reg2$residuals^2)
#         dfGE <- reg2$df.residual
#         ss[cont] <- ssGE; dfr[cont] <- dfGE
#         labTab[cont] <- paste(namEff[i], "Env", sep = ":")
#         cont <- cont + 1
#       } else {
#           if(object$Block == F) {
#             reg2 <- lm.fit(X2, Y)
#             ssGE <- sum(reg2$residuals^2)
#             dfGE <- reg2$df.residual
#             ss[cont] <- ssGE; dfr[cont] <- dfGE
#             labTab[cont] <- paste(namEff[i], "Env", sep = ":")
#             cont <- cont + 1
#           }
#         }
#
#       # Model with no effects
#       X3 <- as.matrix( X[, sel2] )
#       reg3 <- lm.fit(X3, Y)
#       ssG <- sum(reg3$residuals^2)
#       dfG <- reg3$df.residual
#       ss[cont] <- ssG; dfr[cont] <- dfG; labTab[cont] <- namEff[i]
#       cont <- cont + 1
#       }
#       reg.null <- lm(Y ~ 1)
#       totss <- deviance(reg.null)
#       totdf <- reg.null$df.residual
#       ss[cont] <- totss; dfr[cont] <- totdf
#       ss <- diff(ss); dfr <- diff(dfr)
#       ss <- c(rev(ss), rss); dfr <- c(rev(dfr), resdf)
#       labTab <- c("Environment", rev(labTab))
#       labTab[labTab=="Block"] <- "Env:Block"
#       ms <- ss/dfr
#       if(is.null(MSEor)){
#         MSE <- ms[length(ms)]
#       } else {
#         MSE <- MSEor
#         ms[length(ms)] <- MSE
#       }
#       if(!is.null(dfrOr)){
#         dfr1 <- dfrOr
#         dfr[length(dfr)] <- dfrOr
#         ss[length(ss)] <- NA
#       } else {
#         dfr1 <- dfr[length(dfr)]
#       }
#
#       f <- ms/MSE
#       P <- pf(f, dfr, dfr1, lower.tail = FALSE)
#       table <- data.frame(dfr, ss, ms, f, P)
#       table[length(P), 4:5] <- NA
#       colnames(table) <- c("Df","Sum Sq", "Mean Sq", "F value", "Pr(>F)")
#       row.names(table) <- labTab
#       structure(table, heading = c("Analysis of Variance Table\n",
#                                      paste("Response:", deparse(formula(object)[[2L]]))),
#                   class = c("anova", "data.frame"))# was "tabular"
#
#       } else {
#       # Non uso il residuo come errore,
#       ## ma quello fornito
#       # print("OK")
#       rss <- c()
#       fit <- object
#       asgn <- fit$assign
#       X <- fit$modMatrix
#       coefs <- fit$coefficients
#       y <- fit$response
#       ngroup <- length(levels(factor(fit$assign)))
#       rss[1] <- deviance(lm(y ~ 1))
#       for(i in 1:ngroup){
#         #i <- 2
#         val <- fit$assign <= i
#         exp <- X[,val] %*% as.matrix(coefs[val])
#         res <- y - exp
#         rss[i+1] <- sum(res^2)
#       }
#       rss[i + 1] <- 0
#       ss <- rss[1:ngroup] - rss[2:(ngroup+1)]
#       ss <- c(rss[1], ss)
#       df <- c(lengths(split(asgn,  asgn)), fit$df.residual)
#       tlabels <- object$namEff
#
#       if(fit$df.residual == 0){
#          ss <- ss[-length(ss)]
#          df <- df[-length(df)]
#          tlabels <- object$namEff[-length(object$namEff)]
#       }
#       #   df[3] <- df[3]/2
#       # } else if(object$fct == "HAYMAN2" & object$ML == T){
#       #   df[6] <- df[6]/2
#       # }
#       ms <- ss/df
#       table <- data.frame(df, ss, ms)#, f, P)
#       #row.names(table) <- fit$namEff
#       #if(!is.null(MSE)){
#       f <- ms/MSE
#       P <- pf(f, df, dfr, lower.tail = FALSE)
#       table <- data.frame(df, ss, ms, dfr, f, P)
#
#       dimnames(table) <- list(c(tlabels),
#                                 c("Df","Sum Sq", "Mean Sq", "Den df", "F value", "Pr(>F)"))
#         if(attr(object$terms,"intercept")) table <- table[-1, ]
#         structure(table, heading = c("Analysis of Variance Table\n",
#                                      paste("Response:", deparse(formula(object)[[2L]]))),
#                   class = c("anova", "data.frame"))# was "tabular"
#
#       #}
#       table
#       }
#     }

