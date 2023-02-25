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
    # X <- SCA.G3(expl[,1], expl[,2])
    X <- SCA(expl[,1], expl[,2]) # Corrected on 2/7/21
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
    # X <- SCA.G3(expl[,1], expl[,2])
    X <- SCA(expl[,1], expl[,2]) # Corrected on 2/7/2021
    temp2 <- matrix(0, length(X[,1]), length(assign))
    temp2[,assign == i + 2] <- X
    row.names(temp2) <- paste("s", "_", expl[,1], ":", expl[,2], sep = "")
    # REC
    # X <- REC.G3(expl[,1], expl[,2])
    X <- REC(expl[,1], expl[,2]) # Corrected on 2/2/2021
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
    row.names(temp3) <- paste("rg", levs, sep = "_")
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
MET1.eff <- function(obj){
  Y <- obj$model[,1]
  P1 <- factor(obj$model[,2])
  P2 <- factor(obj$model[,3])
  Blk <- factor(obj$model$`(Block)`)
  Env <- factor(obj$model$`(Env)`)
  fct <- obj$fct
  if(length(Blk) == 0){
    temp <- data.frame(Y, P1, P2, Env)
    mods <- plyr::dlply(temp, c("Env"),
      function(df) lm.diallel(Y ~ P1 + P2,
                              fct = fct))
  } else {
    temp <- data.frame(Y, P1, P2, Blk, Env)
    mods <- plyr::dlply(temp, c("Env"),
      function(df) lm.diallel(Y ~ P1 + P2,
                              Block = Blk,
                              fct = fct))
  }

  k_env <- function(ll){
    res <- diallel.eff(ll)
    res <- res$linfct
    colnames(res) <- names(coef(ll))
    res
    }
  mats <- lapply(mods, k_env)
  # Rimuove l'intercetta
  mats <- lapply(mats, function(x) x[-1, -1])
  for(i in 1:length(levels(temp$Env))) colnames(mats[[i]]) <- paste(colnames(mats[[i]]), names(mats)[i], sep = ":")
  for(i in 1:length(levels(temp$Env))) rownames(mats[[i]]) <- paste(rownames(mats[[i]]), names(mats)[i], sep = ":")
  colNames <- unlist(lapply(mats, colnames))
  rowNames <- unlist(lapply(mats, rownames))

  mats <- lmDiallel::blockMatrixDiagonal(mats)
  colnames(mats) <- colNames
  rownames(mats) <- rowNames
  mats2 <- matrix(0, length(mats[,1]), length(levels(temp$Env)))
  k <- cbind(mats2, mats)
  return(k)
}

MET2.eff <- function(obj){
  Y <- obj$model[,1]
  P1 <- factor(obj$model[,2])
  P2 <- factor(obj$model[,3])
  Blk <- factor(obj$model$`(Block)`)
  Env <- factor(obj$model$`(Env)`)
  fct <- obj$fct
    if(length(Blk) == 0){
    temp <- data.frame(Y, P1, P2, Env)
    mods <- plyr::dlply(temp, c("Env"),
      function(df) lm.diallel(Y ~ P1 + P2,
                              fct = fct))
  } else {
    temp <- data.frame(Y, P1, P2, Blk, Env)
    mods <- plyr::dlply(temp, c("Env"),
      function(df) lm.diallel(Y ~ P1 + P2,
                              Block = Blk,
                              fct = fct))
  }
  k_env <- function(ll){
    res <- diallel.eff(ll)
    res <- res$linfct
    colnames(res) <- names(coef(ll))
    res
    }
  mats <- lapply(mods, k_env)

  # Rimuove l'intercetta
  mats <- lapply(mats, function(x) x[-1, -1])
  for(i in 1:length(levels(temp$Env))) colnames(mats[[i]]) <- paste(colnames(mats[[i]]), names(mats)[i], sep = ":")
  for(i in 1:length(levels(temp$Env))) rownames(mats[[i]]) <- paste(rownames(mats[[i]]), names(mats)[i], sep = ":")
  colNames <- unlist(lapply(mats, colnames))
  rowNames <- unlist(lapply(mats, rownames))
  k <- mats[[1]]
  for(i in 2:length(mats)){
    k <- cbind(k, mats[[i]])
  }
  k
  mats2 <- matrix(0, length(k[,1]), length(mats))
  k <- cbind(mats2, k)/length(mats)
  return(k)
}

MET3.eff <- function(obj){
  Y <- obj$model[,1]
  P1 <- factor(obj$model[,2])
  P2 <- factor(obj$model[,3])
  Blk <- factor(obj$model$`(Block)`)
  Env <- factor(obj$model$`(Env)`)
  fct <- obj$fct
  if(length(Blk) == 0){
    temp <- data.frame(Y, P1, P2, Env)
    mod <- lm.diallel(Y ~ P1 + P2, data = temp,
                       Block = Env,
                       fct = fct)
  } else {
    BlockEnv <- factor(paste(Blk, Env, sep = ":"))
    temp <- data.frame(Y, P1, P2, BlockEnv)
    mod <- lm.diallel(Y ~ P1 + P2, data = temp,
                       Block = BlockEnv,
                       fct = fct)
  }
  k <- diallel.eff(mod)
  return(k)
}

diallel.eff <- function(obj, MSE = NULL, dfr = NULL, type = "all") {
  if(all(class(obj) != "diallel")) {
    cat("This method works only with diallel objects")
    stop()
  }
  if(obj$Env == F){
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
  } else {
    # Multiambiente: refitta il modello per ambiente
    if(type == "all"){
      k <- MET1.eff(obj)
    } else if(type == "means"){
      k <- MET2.eff(obj)
    } else if(type == "reduced"){
      k <- MET3.eff(obj)
    }else {
      print("Argument type may be either 'all' or 'means' or 'reduced'")
      stop()
    }
  }



    linfct.list <- list(linfct = k, MSE = MSE, dfr = dfr, obj = obj)
    class(linfct.list) <- "diallelMod"
    return(linfct.list)
}


### multiple comparison procedures
# glht.diallelMod <- function(model, linfct, ...) {
#     obj <- linfct$obj
#     if(obj$Env == F){
#     ### extract factors and contrast matrices from `model'
#     # obj <- linfct$obj
#     MSE <- linfct$MSE
#     dfr <- ifelse(is.null(linfct$dfr), obj$df.residual, linfct$dfr)
#     k <- linfct$linfct
#
#     coefMod <- coef(obj)
#     vcovMod <- vcov(obj, MSE = MSE)
#     args <- list(coef = coefMod, vcov = vcovMod, df = dfr)
#     class(args) <- "parm"
#     ret <- multcomp::glht(args, k)
#     return(ret)
#     } else {
#     newData <- data.frame(Yield = obj$model[,1],
#                       Par1 = obj$model[,2],
#                       Par2 = obj$model[,3],
#                       BlockEnv = paste(obj$model$`(Block)`,
#                                        obj$model$`(Env)`, sep = ":"))
#     newFit <- lm.diallel(Yield ~ Par1 + Par2, data = newData,
#                       Block = BlockEnv,
#                       fct = obj$fct)
#     ret <- glht(linfct = diallel.eff(newFit))
#     return(ret)
#     }
#
# }


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
