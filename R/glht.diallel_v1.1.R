hayman2.eff <- function(obj){
  # Get the data
    assign <- attr(model.matrix(obj), "assign")
    P1 <- obj$model[,2]
    P2 <- obj$model[,3]
    fct <- obj$fct

    # Intercept
    temp <- matrix(0, 1, length(assign))
    X <- 1
    temp[,assign == 0] <- X
    row.names(temp) <- "Intercept"
    i <- 0
    if(obj$Block == T) {i <- 1}

    # MD
    temp.m <- matrix(0, 1, length(assign))
    X <- 1
    temp.m[,assign == i + 1] <- X
    row.names(temp.m) <- "m"

    # GCA
    P1 <- factor(as.character(P1))
    P2 <- factor(as.character(P2))
    levs <- c(levels(P1), levels(P2))
    levs <- levels(factor(levs))
    levs <- factor(levs)
    temp1 <- matrix(0, length(levs), length(assign))
    temp3 <- temp1
    contrasts(levs) <- "contr.sum"
    X <- model.matrix(~levs)[,-1]
    temp1[,assign == i + 2] <- X
    row.names(temp1) <- paste("g", levs, sep = "_")
    # rbind(temp, temp.m, temp1)

    # DD
    temp4 <- temp3
    temp3[,assign == i + 3] <- X
    row.names(temp3) <- paste("d", levs, sep = "_")
    #rbind(temp, temp.m, temp1, temp3)

    # SCA
    expl <- expand.grid(levs,levs)
    X2 <- SCA(expl[,2], expl[,1])
    temp2 <- matrix(0, length(X2[,1]), length(assign))
    temp2[,assign == i + 4] <- X2
    row.names(temp2) <- paste("s", "_", expl[,2], ":", expl[,1], sep = "")
    # rbind(temp, temp.m, temp1, temp3, temp2)

    # RGCA
    temp4[,assign == i + 5] <- X
    row.names(temp4) <- paste("j", levs, sep = "_")
    # rbind(temp, temp.m, temp1, temp3, temp2, temp4)

    # RSCA
    X <- RSCA(expl[,2], expl[,1])
    temp5 <- matrix(0, length(data.frame(X)[,1]), length(assign))
    temp5[,assign == i + 6] <- X
    row.names(temp5) <- paste("rs", "_", expl[,2], ":", expl[,1], sep = "")
    # rimuovere le righe senza elementi non-zero
    X <- rbind(temp, temp.m, temp1, temp3, temp2, temp4, temp5)
    X <- X[apply(X, 1, function(x) !all(x==0)),]
    return(X)
  }

hayman1.eff <- function(obj){
  # Get the data
    assign <- attr(model.matrix(obj), "assign")
    P1 <- obj$model[,2]
    P2 <- obj$model[,3]
    fct <- obj$fct
    # Intercept
    temp <- matrix(0, 1, length(assign))
    X <- 1
    temp[,assign == 0] <- X
    row.names(temp) <- "Intercept"
    i <- 0
    if(obj$Block == T) {i <- 1}

    # GCA
    P1 <- factor(as.character(P1))
    P2 <- factor(as.character(P2))
    levs <- c(levels(P1), levels(P2))
    levs <- levels(factor(levs))
    levs <- factor(levs)
    temp1 <- matrix(0, length(levs), length(assign))
    temp3 <- temp1
    contrasts(levs) <- "contr.sum"
    X <- model.matrix(~levs)[,-1]
    temp1[,assign == i + 1] <- X
    row.names(temp1) <- paste("g", levs, sep = "_")
    # RGCA
    temp3[,assign == i + 3] <- X
    row.names(temp3) <- paste("j", levs, sep = "_")
    # tSCA
    expl <- expand.grid(levs,levs)
    X <- tSCA(expl[,2], expl[,1])
    temp2 <- matrix(0, length(X[,1]), length(assign))
    temp2[,assign == i + 2] <- X
    row.names(temp2) <- paste("ts", "_", expl[,2], ":", expl[,1], sep = "")
    # RSCA
    X <- RSCA(expl[,2], expl[,1])
    temp4 <- matrix(0, length(data.frame(X)[,1]), length(assign))
    temp4[,assign == i + 4] <- X
    row.names(temp4) <- paste("rs", "_", expl[,2], ":", expl[,1], sep = "")
    # rimuovere le righe senza elementi non-zero
    X <- rbind(temp, temp1, temp2, temp3, temp4)
    X <- X[apply(X, 1, function(x) !all(x==0)),]
    return(X)
  }

diallel.eff <- function(obj, MSE = NULL, dfr = NULL) {
  if(all(class(obj) != "diallel")) {
    cat("This method works only with diallel objects")
    stop()
  }
    if(obj$fct == "HAYMAN1"){
       k <- hayman1.eff(obj)
    } else if(obj$fct == "HAYMAN2"){
       k <- hayman2.eff(obj)
    }
    linfct.list <- list(linfct = k, MSE = MSE, dfr = dfr, obj = obj)
    class(linfct.list) <- "diallelMod"
    return(linfct.list)
}


### multiple comparison procedures
glht.diallelMod <- function(linfct, ...) {

    ### extract factors and contrast matrices from `model'
    obj <- linfct$obj
    MSE <- linfct$MSE
    dfr <- ifelse(is.null(linfct$dfr), obj$df.residual, linfct$dfr)
    k <- linfct$linfct

    coefMod <- coef(obj)
    vcovMod <- vcov(obj, MSE = MSE)
    args <- list(coef = coefMod, vcov = vcovMod, df = 26)
    class(args) <- "parm"
    ret <- multcomp::glht(args, k)
    return(ret)

}
