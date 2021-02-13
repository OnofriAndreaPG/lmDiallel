GE3r.eff <- function(obj){
  # Get the data
  # obj <- dMod2; coef(dMod2)
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
  # h.bar
  X <- 1
  temp1 <- matrix(0, 1, length(assign))
  temp1[,assign == i + 1] <- X
  row.names(temp1) <- "h.bar"
  # SPi
  P1 <- factor(as.character(P1))
  P2 <- factor(as.character(P2))
  levs <- c(levels(P1), levels(P2))
  levs <- levels(factor(levs))
  levs <- factor(levs)
  temp2 <- matrix(0, length(levs), length(assign))
  contrasts(levs) <- "contr.sum"
  X <- model.matrix(~levs)[,-1]
  temp2[,assign == i + 2] <- X
  row.names(temp2) <- paste("sp", levs, sep = "_")
  # GC
  temp3 <- matrix(0, length(levs), length(assign))
  temp3[,assign == i + 3] <- X
  row.names(temp3) <- paste("gc", levs, sep = "_")
  
  # SCA
  expl <- expand.grid(levs,levs)
  X2 <- SCA(expl[,2], expl[,1])
  temp4 <- matrix(0, length(X2[,1]), length(assign))
  temp4[,assign == i + 4] <- X2
  row.names(temp4) <- paste("s", "_", expl[,2], ":", expl[,1], sep = "")
  # REC
  X <- REC(expl[,2], expl[,1])
  temp5 <- matrix(0, length(data.frame(X)[,1]), length(assign))
  temp5[,assign == i + 5] <- X
  row.names(temp5) <- paste("r", "_", expl[,2], ":", expl[,1], sep = "")
  
  X <- rbind(temp, temp1, temp2, temp3, temp4)
  X <- X[apply(X, 1, function(x) !all(x==0)),]
  return(X)
}

GE2r.eff <- function(obj){
  # Get the data
  # obj <- dMod2; coef(dMod2)
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
  # h.bar
  X <- 1
  temp1 <- matrix(0, 1, length(assign))
  temp1[,assign == i + 1] <- X
  row.names(temp1) <- "h.bar"
  # VEi
  P1 <- factor(as.character(P1))
  P2 <- factor(as.character(P2))
  levs <- c(levels(P1), levels(P2))
  levs <- levels(factor(levs))
  levs <- factor(levs)
  temp2 <- matrix(0, length(levs), length(assign))
  contrasts(levs) <- "contr.sum"
  X <- model.matrix(~levs)[,-1]
  temp2[,assign == i + 2] <- X
  row.names(temp2) <- paste("v", levs, sep = "_")
  # Hi
  # P1 <- factor(as.character(P1))
  # P2 <- factor(as.character(P2))
  # levs <- c(levels(P1), levels(P2))
  # levs <- levels(factor(levs))
  # levs <- factor(levs)
  temp3 <- matrix(0, length(levs), length(assign))
  #contrasts(levs) <- "contr.sum"
  #X <- model.matrix(~levs)[,-1]
  temp3[,assign == i + 3] <- X
  row.names(temp3) <- paste("h", levs, sep = "_")
  
  # SCA
  expl <- expand.grid(levs,levs)
  X2 <- SCA(expl[,2], expl[,1])
  temp4 <- matrix(0, length(X2[,1]), length(assign))
  temp4[,assign == i + 4] <- X2
  row.names(temp4) <- paste("s", "_", expl[,2], ":", expl[,1], sep = "")
  
  # REC
  X <- REC(expl[,2], expl[,1])
  temp5 <- matrix(0, length(data.frame(X)[,1]), length(assign))
  temp5[,assign == i + 5] <- X
  row.names(temp5) <- paste("r", "_", expl[,2], ":", expl[,1], sep = "")
  
  X <- rbind(temp, temp1, temp2, temp3, temp4, temp5)
  
  # tSCA - NO
  # expl <- expand.diallel(as.character(levs), 2)
  # X <- tSCA(expl[,1], expl[,2])
  # temp2 <- matrix(0, length(X[,1]), length(assign))
  # temp2[,assign == i + 2] <- X
  # row.names(temp2) <- paste("ts", "_", expl[,1], ":", expl[,2], sep = "")
  # rimuovere le righe senza elementi non-zero
  # X <- rbind(temp, temp1, temp2)
  X <- X[apply(X, 1, function(x) !all(x==0)),]
  return(X)
}

GE3.eff <- function(obj){
  # Get the data
  # obj <- dMod2; coef(dMod2)
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
  # h.bar
  X <- 1
  temp1 <- matrix(0, 1, length(assign))
  temp1[,assign == i + 1] <- X
  row.names(temp1) <- "h.bar"
  # SPi
  P1 <- factor(as.character(P1))
  P2 <- factor(as.character(P2))
  levs <- c(levels(P1), levels(P2))
  levs <- levels(factor(levs))
  levs <- factor(levs)
  temp2 <- matrix(0, length(levs), length(assign))
  contrasts(levs) <- "contr.sum"
  X <- model.matrix(~levs)[,-1]
  temp2[,assign == i + 2] <- X
  row.names(temp2) <- paste("sp", levs, sep = "_")
  # GC
  temp3 <- matrix(0, length(levs), length(assign))
  temp3[,assign == i + 3] <- X
  row.names(temp3) <- paste("gc", levs, sep = "_")
  
  # SCA
  expl <- expand.grid(levs,levs)
  X2 <- SCA(expl[,2], expl[,1])
  temp4 <- matrix(0, length(X2[,1]), length(assign))
  temp4[,assign == i + 4] <- X2
  row.names(temp4) <- paste("s", "_", expl[,2], ":", expl[,1], sep = "")
  X <- rbind(temp, temp1, temp2, temp3, temp4)
  
  # tSCA - NO
  # expl <- expand.diallel(as.character(levs), 2)
  # X <- tSCA(expl[,1], expl[,2])
  # temp2 <- matrix(0, length(X[,1]), length(assign))
  # temp2[,assign == i + 2] <- X
  # row.names(temp2) <- paste("ts", "_", expl[,1], ":", expl[,2], sep = "")
  # rimuovere le righe senza elementi non-zero
  # X <- rbind(temp, temp1, temp2)
  X <- X[apply(X, 1, function(x) !all(x==0)),]
  return(X)
}

GE2.eff <- function(obj){
  # Get the data
  # obj <- dMod2; coef(dMod2)
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
  # h.bar
  X <- 1
  temp1 <- matrix(0, 1, length(assign))
  temp1[,assign == i + 1] <- X
  row.names(temp1) <- "h.bar"
  # VEi
  P1 <- factor(as.character(P1))
  P2 <- factor(as.character(P2))
  levs <- c(levels(P1), levels(P2))
  levs <- levels(factor(levs))
  levs <- factor(levs)
  temp2 <- matrix(0, length(levs), length(assign))
  contrasts(levs) <- "contr.sum"
  X <- model.matrix(~levs)[,-1]
  temp2[,assign == i + 2] <- X
  row.names(temp2) <- paste("v", levs, sep = "_")
  # Hi
  # P1 <- factor(as.character(P1))
  # P2 <- factor(as.character(P2))
  # levs <- c(levels(P1), levels(P2))
  # levs <- levels(factor(levs))
  # levs <- factor(levs)
  temp3 <- matrix(0, length(levs), length(assign))
  #contrasts(levs) <- "contr.sum"
  #X <- model.matrix(~levs)[,-1]
  temp3[,assign == i + 3] <- X
  row.names(temp3) <- paste("h", levs, sep = "_")
  
  # SCA
  expl <- expand.grid(levs,levs)
  X2 <- SCA(expl[,2], expl[,1])
  temp4 <- matrix(0, length(X2[,1]), length(assign))
  temp4[,assign == i + 4] <- X2
  row.names(temp4) <- paste("s", "_", expl[,2], ":", expl[,1], sep = "")
  X <- rbind(temp, temp1, temp2, temp3, temp4)
  
  # tSCA - NO
  # expl <- expand.diallel(as.character(levs), 2)
  # X <- tSCA(expl[,1], expl[,2])
  # temp2 <- matrix(0, length(X[,1]), length(assign))
  # temp2[,assign == i + 2] <- X
  # row.names(temp2) <- paste("ts", "_", expl[,1], ":", expl[,2], sep = "")
  # rimuovere le righe senza elementi non-zero
  # X <- rbind(temp, temp1, temp2)
  X <- X[apply(X, 1, function(x) !all(x==0)),]
  return(X)
}

G4.eff <- function(obj){
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
    # temp3 <- temp1
    contrasts(levs) <- "contr.sum"
    X <- model.matrix(~levs)[,-1]
    temp1[,assign == i + 1] <- X
    row.names(temp1) <- paste("g", levs, sep = "_")
    # SCA
    expl <- expand.diallel(as.character(levs), 3)
    X <- SCA.G3(expl[,1], expl[,2])
    temp2 <- matrix(0, length(X[,1]), length(assign))
    temp2[,assign == i + 2] <- X
    row.names(temp2) <- paste("s", "_", expl[,1], ":", expl[,2], sep = "")
    # rimuovere le righe senza elementi non-zero
    X <- rbind(temp, temp1, temp2)
    X <- X[apply(X, 1, function(x) !all(x==0)),]
    return(X)
  }

G3.eff <- function(obj){
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
    # temp3 <- temp1
    contrasts(levs) <- "contr.sum"
    X <- model.matrix(~levs)[,-1]
    temp1[,assign == i + 1] <- X
    row.names(temp1) <- paste("g", levs, sep = "_")
    # SCA
    expl <- expand.diallel(as.character(levs), 3)
    X <- SCA.G3(expl[,1], expl[,2])
    temp2 <- matrix(0, length(X[,1]), length(assign))
    temp2[,assign == i + 2] <- X
    row.names(temp2) <- paste("s", "_", expl[,1], ":", expl[,2], sep = "")
    # REC
    X <- REC.G3(expl[,1], expl[,2])
    temp4 <- matrix(0, length(data.frame(X)[,1]), length(assign))
    temp4[,assign == i + 3] <- X
    row.names(temp4) <- paste("r", "_", expl[,1], ":", expl[,2], sep = "")
    # rimuovere le righe senza elementi non-zero
    X <- rbind(temp, temp1, temp2, temp4)
    X <- X[apply(X, 1, function(x) !all(x==0)),]
    return(X)
  }

G2.eff <- function(obj){
  # Get the data
  # obj <- dMod2
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
    # temp3 <- temp1
    contrasts(levs) <- "contr.sum"
    X <- model.matrix(~levs)[,-1]
    temp1[,assign == i + 1] <- X
    row.names(temp1) <- paste("g", levs, sep = "_")
    # tSCA
    expl <- expand.diallel(as.character(levs), 2)
    X <- tSCA(expl[,1], expl[,2])
    temp2 <- matrix(0, length(X[,1]), length(assign))
    temp2[,assign == i + 2] <- X
    row.names(temp2) <- paste("ts", "_", expl[,1], ":", expl[,2], sep = "")
    # rimuovere le righe senza elementi non-zero
    X <- rbind(temp, temp1, temp2)
    X <- X[apply(X, 1, function(x) !all(x==0)),]
    return(X)
  }

G1.eff <- function(obj){
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
    # temp3 <- temp1
    contrasts(levs) <- "contr.sum"
    X <- model.matrix(~levs)[,-1]
    temp1[,assign == i + 1] <- X
    row.names(temp1) <- paste("g", levs, sep = "_")
    # tSCA
    expl <- expand.grid(levs,levs)
    X <- tSCA(expl[,2], expl[,1])
    temp2 <- matrix(0, length(X[,1]), length(assign))
    temp2[,assign == i + 2] <- X
    row.names(temp2) <- paste("ts", "_", expl[,2], ":", expl[,1], sep = "")
    # REC
    X <- REC(expl[,2], expl[,1])
    temp4 <- matrix(0, length(data.frame(X)[,1]), length(assign))
    temp4[,assign == i + 3] <- X
    row.names(temp4) <- paste("r", "_", expl[,2], ":", expl[,1], sep = "")
    # rimuovere le righe senza elementi non-zero
    X <- rbind(temp, temp1, temp2, temp4)
    X <- X[apply(X, 1, function(x) !all(x==0)),]
    return(X)
  }

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
    } else if(obj$fct == "GRIFFING1"){
       k <- G1.eff(obj)
    } else if(obj$fct == "GRIFFING2"){
       k <- G2.eff(obj)
    } else if(obj$fct == "GRIFFING3"){
       k <- G3.eff(obj)
    } else if(obj$fct == "GRIFFING4"){
       k <- G4.eff(obj)
    } else if(obj$fct == "GE2"){
      k <- GE2.eff(obj)
    } else if(obj$fct == "GE3"){
      k <- GE3.eff(obj)
    } else if(obj$fct == "GE2r"){
      k <- GE2r.eff(obj)
    } else if(obj$fct == "GE3r"){
      k <- GE3r.eff(obj)
    }

    linfct.list <- list(linfct = k, MSE = MSE, dfr = dfr, obj = obj)
    class(linfct.list) <- "diallelMod"
    return(linfct.list)
}


### multiple comparison procedures
glht.diallelMod <- function(model, linfct, ...) {

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

expand.diallel <- function(pars, mating = 1){
  pars <- sort(pars)
  # print(pars)
  if(mating == 1){
    Par1 <- rep(pars, each = length(pars))
    Par2 <- rep(pars, length(pars))
  } else if(mating == 2){
    Par1 <- rep(pars, c(length(pars):1))
    Par2 <- pars
    for(i in 1:length(pars)){
    Par2 <- c(Par2, pars[-c(1:i)])}
  } else if(mating == 3) {
    Par1i <- rep(pars, each = length(pars))
    Par2i <- rep(pars, length(pars))
    Par1 <- Par1i[Par1i != Par2i]
    Par2 <- Par2i[Par1i != Par2i]
  } else if(mating == 4) {
    Par1i <- rep(pars, c(length(pars):1))
    Par2i <- pars
    for(i in 1:length(pars)){
    Par2i <- c(Par2i, pars[-c(1:i)])}
    Par1 <- Par1i[Par1i != Par2i]
    Par2 <- Par2i[Par1i != Par2i]
  }
df <- data.frame(Par1, Par2)
return(df)
}
