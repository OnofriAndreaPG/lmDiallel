# This functions create model matrices for diallel models
# Date of last edit: 23/6/2020
model.matrix.diallel <- function(obj){
  return(obj$modMatrix)
}
model.matrixDiallel <- function(formula, Block = NULL, Env = NULL,
                                 fct = NULL, data = NULL, ML = F){
  if(is.null(fct)){
    # This is a formula based output ###############
    X1 <- model.matrix(formula, data)
    X <- X1[,-1]
    attr(X, "assign") <- attr(X1, "assign")[-1]
  } else {
  # fct based output #############################
  mf <- match.call(expand.dots = FALSE) # Riprende la chiamata, con i nomi
  m <- match(c("formula", "Block", "Env", "data"), names(mf), 0L) # Trova nella chiamata la formula. m è la posizione della formula nella chiamata
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  pars <- attr(mt, "term.labels")

  bName <- deparse(substitute(Block))  # storing name of Blocks
  Block <- model.extract(mf, "Block")
  eName <- deparse(substitute(Env))  # storing name of Env
  Env <- model.extract(mf, "Env")

  if(missing(data) == T){
    Par1 <- mf[,1]
    Par2 <- mf[,2]
  } else {
    Par1 <- data[[pars[1]]]
    Par2 <- data[[pars[2]]]
  }

if(is.null(Env) == T){
      n <- length(Par1)
      P1 <- factor(as.character(Par1))
      P2 <- factor(as.character(Par2))
      nGroups <- 0
      groups <- c(nGroups)
      X <- cbind(Intercept = rep(1, n))
      reps <- c(1)
      namEffs <- c("Intercept")

      if(is.null(Block) == F) {
        Block <- factor(Block)
        B <- matBlock(~Block)
        colnames(B) <- paste("Block", levels(Block)[-length(levels(Block))], sep = "")
        nGroups <- 1
        groups <- c(groups, nGroups)
        X <- cbind(X, B)
        reps <- c(reps, length(B[1,]))
        namEffs <- c(namEffs, "Block")
      }

    if(fct == "HAYMAN1"){
      # HAYMAN1 - con tSCA #############################
      # 6/04/2020

      # Matrix for GCA
      Z <- GCA(P1, P2)


      # Matrix tSCA
      SCA <- tSCA(P1, P2)

      #Matrix for RGCA
      RGCA <- RGCA(P1, P2)


      #Matrix for RSCA
      rec <- RSCA(P1, P2)


      # Building matrix (0:5)
      X <- cbind(X, Z, SCA, RGCA, rec)
      groups <- c(groups, seq(nGroups+1, nGroups+4, 1))
      reps <- c(reps, length(Z[1,]), length(SCA[1,])
                              ,length(RGCA[1,])
                              ,length(rec[1,]))
      levs <- rep(groups, reps)
      attr(X, "assign") <- levs
      namEffs <- c(namEffs, "GCA","tSCA", "RGCA", "RSCA",
                                 "Residuals")
      attr(X, "namEff") <- namEffs

      } else if(fct == "HAYMAN2"){
      # HAYMAN2 - SCA decomposta ###########
      # 23/03/2020

      # Matrix for crosses
      crM <- MDD(P1, P2)


      # Matrix for GCA
      Z <- GCA(P1, P2)
      # nams <- paste("gca_", levels(P1)[1:(length(levels(P1))-1)], sep="")
      # colnames(Z) <- c(nams)

      # Matrix for h.i
      H <- DD(P1, P2)

      # Matrix for sca
      SCA <- SCA(P1, P2)
      #colnames(SCA) <- paste("sca_", colnames(SCA), sep = "")

      # Matrix for RGCA
      RGCA <- RGCA(P1, P2)
      # nams <- paste("rgca_", levels(P1)[1:length(levels(P1))-1], sep="")
      # colnames(RGCA) <- c(nams)

      # Matrix for RSCA
      rec <- RSCA(P1, P2)
      # colnames(rec) <- paste("rsca_", colnames(rec), sep = "")

      # Building incidence matrix (0:7)
      X <- cbind(X, crM, Z, H, SCA, RGCA, rec)
      groups <- c(groups, seq(nGroups+1, nGroups+6, 1))
      reps <- c(reps, 1, length(Z[1,]), length(H[1,])
                              ,length(SCA[1,])
                              ,length(RGCA[1,])
                              ,length(rec[1,]))
      levs <- rep(groups, reps)
      attr(X, "assign") <- levs
      namEffs <- c(namEffs, "MDD", "GCA", "DD", "SCA",
                   "RGCA", "RSCA", "Residuals")
      attr(X, "namEff") <- namEffs

      } else if(fct == "GRIFFING1"){

      # GRIFFING 1 - Reciprocal effects #########################################
      #B <- matBlock(~Block)

      Z <- GCA(P1, P2)
      # nams <- paste("gca_", levels(P1)[1:(length(levels(P1))-1)], sep="")
      # colnames(Z) <- c(nams)

      SCA <- tSCA(P1, P2)
      # colnames(SCA) <- paste("sca_", colnames(SCA), sep = "")

      rec <- REC(P1, P2)
      # colnames(rec) <- paste("rec_", colnames(rec), sep = "")

      # Building incidence matrix (0:4)
      X <- cbind(X, Z, SCA, rec)
      groups <- c(groups, seq(nGroups+1, nGroups+3, 1))
      reps <- c(reps, length(Z[1,]),
                      length(SCA[1,])
                    , length(rec[1,]))
      levs <- rep(groups, reps)
      attr(X, "assign") <- levs
      namEffs <- c(namEffs, "GCA", "SCA", "Reciprocals",
                                 "Residuals")
      attr(X, "namEff") <- namEffs

    } else if(fct == "GRIFFING2"){
      # GRIFFING 2 - No reciprocals #########################################
      # 23/03/2020

      #B <- matBlock(~Block)

      Z <- GCA(P1, P2)
      # nams <- paste("gca_", levels(P1)[1:(length(levels(P1))-1)], sep="")
      # colnames(Z) <- c(nams)
      SCA <- tSCA(P1, P2)
      # colnames(SCA) <- paste("sca_", colnames(SCA), sep = "")

      # Building incidence matrix (0:3)
      X <- cbind(X, Z, SCA)
      groups <- c(groups, seq(nGroups+1, nGroups+2, 1))
      reps <- c(reps, length(Z[1,]),
                      length(SCA[1,]))
      levs <- rep(groups, reps)
      attr(X, "assign") <- levs
      namEffs <- c(namEffs, "GCA", "SCA",
                                 "Residuals")
      attr(X, "namEff") <- namEffs

    } else if(fct == "GE2"){
      # GE2 - Senza reciproci ####################
      # 23/03/2020
      # B <- matBlock(~Block)

      # Matrix for bar_h
      crM <- H.BAR(P1, P2)
      # colnames(crM) <- "h.bar"

      # Matrix for nu.i
      Z <- VEi(P1, P2)
      # nams <- paste("ve_", levels(P1)[1:(length(levels(P1))-1)], sep="")
      # colnames(Z) <- c(nams)

      # Matrix for h.i
      H <- Hi(P1, P2)


      # Matrix for sca
      SCA <- SCA(P1, P2)
      # colnames(SCA) <- paste("sca_", colnames(SCA), sep = "")

      # Building incidence matrix (0:5)
      X <- cbind(X, crM, Z, H, SCA)
      groups <- c(groups, seq(nGroups+1, nGroups+4, 1))
      reps <- c(reps,  1 ,length(Z[1,])
                         ,length(H[1, ])
                         ,length(SCA[1,]))
      levs <- rep(groups, reps)
      attr(X, "assign") <- levs
      namEffs <- c(namEffs, "h.bar", "Variety", "h.i", "SCA",
                                 "Residuals")
      attr(X, "namEff") <- namEffs

    } else if(fct == "GE2r"){
      # GE2r - With reciprocals ####################
      # 23/03/2020

      #B <- matBlock(~Block)

      # Matrix for bar_h
      crM <- H.BAR(P1, P2)
      # colnames(crM) <- "h.bar"

      # # Matrix for crosses
      # slM <- H.BAR(crosses)
      # colnames(slM) <- "Selfs"

      Z <- VEi(P1, P2)
      # nams <- paste("ve_", levels(P1)[1:(length(levels(P1))-1)], sep="")
      # colnames(Z) <- c(nams)

      H <- Hi(P1, P2)
      # nams <- paste("h_", levels(P1)[1:(length(levels(P1))-1)], sep="")
      # colnames(H) <- c(nams)

      # Matrix for sca
      SCA <- SCA(P1, P2)
      # colnames(SCA) <- paste("sca_", colnames(SCA), sep = "")

      rec <- REC(P1, P2)
      # colnames(rec) <- paste("rec_", colnames(rec), sep = "")

      # Building incidence matrix (0:6)
      X <- cbind(X, crM, Z, H, SCA, rec)
      groups <- c(groups, seq(nGroups+1, nGroups+5, 1))
      reps <- c(reps,  1 ,length(Z[1,])
                         ,length(H[1, ])
                         ,length(SCA[1,]),
                          length(rec[1,]))
      levs <- rep(groups, reps)
      attr(X, "assign") <- levs
      namEffs <- c(namEffs, "h.bar", "Variety", "h.i", "SCA", "Reciprocals",
                                 "Residuals")
      attr(X, "namEff") <- namEffs

    } else if(fct == "GE3"){
      # GE3 - Senza reciproci ###############################################
      # 23/03/2020

      # Matrix for crosses
      #crM <- matrix(crosses, n, 1)
      crM <- H.BAR(P1, P2)
      # colnames(crM) <- "h.bar"

      # Matrix for selfs
      # slM <- H.BAR(crosses)
      # colnames(slM) <- "Selfs"
      #
      # Matrix for GCA
      H <- SP(P1, P2)


      Z <- GCAC(P1, P2)


      # Matrix for sca
      SCA <- SCA(P1, P2)
      # colnames(SCA) <- paste("sca_", colnames(SCA), sep = "")

      # Building incidence matrix (0:5)
      X <- cbind(X, crM, H, Z, SCA)
      groups <- c(groups, seq(nGroups+1, nGroups+4, 1))
      reps <- c(reps,  1 ,length(H[1,])
                         ,length(Z[1, ])
                         ,length(SCA[1,]))
      levs <- rep(groups, reps)
      attr(X, "assign") <- levs
      namEffs <- c(namEffs, "h.bar",
                             "Selfed par.",
                             "Varieties", "SCA",
                                 "Residuals")
      attr(X, "namEff") <- namEffs

    } else if(fct == "GE3r"){
      # GE3r - Con reciproci ###############################################
      # 23/03/2020

      # Matrix for crosses
      #crM <- matrix(crosses, n, 1)
      crM <- H.BAR(P1, P2)
      # colnames(crM) <- "h.bar"

      Z <- GCAC(P1, P2)
      # nams <- paste("gcac_", levels(P1)[1:(length(levels(P1))-1)], sep="")
      # colnames(Z) <- c(nams)

      H <- SP(P1, P2)
      # nams <- paste("sp_", levels(P1)[1:(length(levels(P1))-1)], sep="")
      # colnames(H) <- c(nams)

      # Matrix for sca
      SCA <- SCA(P1, P2)
      # colnames(SCA) <- paste("sca_", colnames(SCA), sep = "")

      rec <- REC(P1, P2)
      #  colnames(rec) <- paste("rec_", colnames(rec), sep = "")

      # Building incidence matrix (0:6)
      X <- cbind(X, crM, Z, H, SCA, rec)
      groups <- c(groups, seq(nGroups+1, nGroups+5, 1))
      reps <- c(reps,  1 ,length(Z[1,])
                         ,length(H[1, ])
                         ,length(SCA[1,]),
                length(rec[1,]))
      levs <- rep(groups, reps)
      attr(X, "assign") <- levs
      namEffs <- c(namEffs, "h.bar",
                             "gcac",
                             "Selfed par.", "SCA", "Reciprocals",
                                 "Residuals")
      attr(X, "namEff") <- namEffs

    }else{
         stop("Not yet implemented")
      }
  } else {
    # GE Data ####################################
    Block <- factor(Block)
    Env <- factor(Env)
    datasetS <- data.frame(Id= 1:length(Block), Env, Block, Par1, Par2)
    datasetS <- datasetS[order(datasetS$Env, datasetS$Par1,
                   datasetS$Par2, datasetS$Block), ]

    # Creazione matrice incidenza (per ogni ambiente)
    matsOr <- plyr::dlply(datasetS, c("Env"), function(df){
              model.matrixDiallel(~ df$Par1+df$Par2, df$Block,
                  fct = fct)})
    mats <- matsOr
    mats <- lapply(mats, function(x) x[, -1])
    for(i in 1:length(levels(datasetS$Env))) colnames(mats[[i]]) <- paste(colnames(mats[[i]]), names(mats)[i], sep = ":")
    colNames <- unlist(lapply(mats, colnames))
    mats <- blockMatrixDiagonal(mats)
    colnames(mats) <- colNames
    mats2 <- model.matrix(~ Env - 1, data = datasetS)
    X <- cbind(mats2, mats)
    X <- X[order(datasetS$Id), ]

    # Creating the submatrices
    asgnList <- lapply(matsOr, function(x) attr(x, "assign"))
    asgnList <- lapply(asgnList, function(x) unlist(x)[-1])
    addVal <- max(unlist(asgnList[1]))
    asgn <- c(unlist(asgnList[1]),unlist(lapply(asgnList[-1], function(x) unlist(x) + addVal)) )
    asgn <- as.numeric(c(rep(0, length(levels(datasetS$Env))), asgn))
    attr(X, "assign") <- asgn
    attr(X, "namEff") <- as.character( unlist( lapply(matsOr, function(x) attr(x, "namEff"))[1] ) )
  }}
  return(X)
}

blockMatrixDiagonal<-function(matList){
  dimensionsRow <- sapply(matList, FUN=function(x) dim(x)[1])
  dimensionsCol <- sapply(matList, FUN=function(x) dim(x)[2])

  finalDimensionRow <- sum(dimensionsRow)
  finalDimensionCol <- sum(dimensionsCol)
  finalMatrix<-matrix(0, nrow=finalDimensionRow, ncol=finalDimensionCol)
  indexRow <- 1; indexCol <- 1
  for(k in 1:length(dimensionsRow)){
    #print(paste(k, indexRow, (indexRow + dimensionsRow[k]-1), sep="-"))

    finalMatrix[indexRow:(indexRow + dimensionsRow[k]-1),indexCol:(indexCol+dimensionsCol[k]-1)] <- matList[[k]]
    indexRow <- indexRow + dimensionsRow[k]
    indexCol <- indexCol + dimensionsCol[k]
    }
    finalMatrix
}

matBlock <- function(formula){
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  nameFac <- attr(mt, "term.labels")
  fac <- factor( mf[[1]])
  n <- length(fac)
  contrasts(fac) <- c("contr.sum")
  B <-  model.matrix(~fac)
  B <- B[,-1]
  if(is.vector(B) == T) B <- matrix(B, n, 1)
  colnames(B) <- paste(nameFac, levels(fac)[-length(levels(fac))], sep = "")
  B
  }

# GCA.old <- function(P1, P2){
#   P1 <- factor(as.character(P1))
#   P2 <- factor(as.character(P2))
#   contrasts(P1) <- c("contr.sum") #Parte da qui
#   contrasts(P2) <- c("contr.sum")
#   Z1 <- model.matrix(~P1)
#   Z2 <- model.matrix(~P2)
#   Z <- (Z1 + Z2)
#   Z <- Z[,-1]
#   Z
# }

GCA <- function(P1, P2, data = NULL){
  # This is modified to work with mating design 4
  if(!is.null(data)){
    P1Name <- deparse(substitute(P1))
    P2Name <- deparse(substitute(P2))
    P1 <- data[[P1Name]]
    P2 <- data[[P2Name]]
  }
  P1 <- factor(as.character(P1))
  P2 <- factor(as.character(P2))
  levs <- c(levels(P1), levels(P2))
  levs <- levels(factor(levs))
  Z1n <- factor(P1, levels = levs, ordered = T)
  Z2n <- factor(P2, levels = levs)
  contrasts(Z1n) <- c("contr.sum")
  contrasts(Z2n) <- c("contr.sum")
  Z1 <- model.matrix(~ Z1n)
  Z2 <- model.matrix(~ Z2n)
  # contrasts(P1) <- c("contr.sum") #Parte da qui
  # contrasts(P2) <- c("contr.sum")
  # Z1 <- model.matrix(~P1)
  # Z2 <- model.matrix(~P2)
  Z <- (Z1 + Z2)
  Z <- Z[,-1]
  #nams <- paste("g_", levels(P1)[1:(length(levels(P1))-1)], sep="")
  nams <- paste("g_", levs[1:(length(levs)-1)], sep="")
  colnames(Z) <- c(nams)
  Z
}

VEi <- function(P1, P2, data = NULL){
    if(!is.null(data)){
    P1Name <- deparse(substitute(P1))
    P2Name <- deparse(substitute(P2))
    P1 <- data[[P1Name]]
    P2 <- data[[P2Name]]
  }
  # For GE2 models
  P1 <- factor(as.character(P1))
  P2 <- factor(as.character(P2))
  contrasts(P1) <- c("contr.sum")
  contrasts(P2) <- c("contr.sum")
  Z1 <- model.matrix(~P1)
  Z2 <- model.matrix(~P2)
  Z <- (Z1 + Z2)/2
  Z <- Z[,-1]
  nams <- paste("v_", levels(P1)[1:(length(levels(P1))-1)], sep="")
  colnames(Z) <- c(nams)
  Z
}

SP <- function(P1, P2, data = NULL){
    if(!is.null(data)){
    P1Name <- deparse(substitute(P1))
    P2Name <- deparse(substitute(P2))
    P1 <- data[[P1Name]]
    P2 <- data[[P2Name]]
  }
  # For GE3 models
  P1 <- factor(as.character(P1))
  P2 <- factor(as.character(P2))
  P1c <- as.character(P1)
  P2c <- as.character(P2)
  selfs <- ifelse(P1c == P2c, 1, 0)
  contrasts(P1) <- c("contr.sum")
  contrasts(P2) <- c("contr.sum")
  Z1 <- model.matrix(~P1)
  Z2 <- model.matrix(~P2)
  Z <- (Z1 + Z2)/2 * selfs
  Z <- Z[,-1]
  nams <- paste("sp_", levels(P1)[1:(length(levels(P1))-1)], sep="")
  colnames(Z) <- c(nams)
  Z
}

RGCA <- function(P1, P2, data = NULL){
  if(!is.null(data)){
    P1Name <- deparse(substitute(P1))
    P2Name <- deparse(substitute(P2))
    P1 <- data[[P1Name]]
    P2 <- data[[P2Name]]
  }
  P1 <- factor(as.character(P1))
  P2 <- factor(as.character(P2))
  contrasts(P1) <- c("contr.sum")
  contrasts(P2) <- c("contr.sum")
  # p <- length(levels(P1))
  p <- levels(P1)[length(levels(P1))] #Correction 14/11/20. AO
  P1c <- as.character(P1)
  P2c <- as.character(P2)
  dr <- ifelse(P1c == P2c, 0, ifelse(P1c < P2c, -1, 1))
  Z3 <- model.matrix(~P1 - 1)
  Z4 <- model.matrix(~P2 - 1)
  RGCA <- (Z3 - Z4) #* -dr
  RGCA[P1==p,] <- RGCA[P1==p,] - 1
  RGCA[P2==p,] <- RGCA[P2==p,] + 1
  # RGCA <- RGCA[,-p]
  RGCA <- RGCA[,-length(levels(P1))] ##Correction 14/11/20. AO
  nams <- paste("rg_", levels(P1)[1:length(levels(P1))-1], sep="")
  colnames(RGCA) <- c(nams)
  RGCA
}

# RGCA2 <- function(P1, P2){
#   # Per ML estimation?
#   P1 <- factor(as.character(P1))
#   P2 <- factor(as.character(P2))
#   contrasts(P1) <- c("contr.sum")
#   contrasts(P2) <- c("contr.sum")
#   p <- length(levels(P1))
#   P1c <- as.character(P1)
#   P2c <- as.character(P2)
#   dr <- ifelse(P1c == P2c, 0, ifelse(P1c < P2c, -1, 1))
#   Z1 <- model.matrix(~P1)[,-1]
#   Z2 <- model.matrix(~P2)[,-1]
#   colnames(Z1) <- paste("k_", c(1:(p-1)), sep="")
#   colnames(Z2) <- colnames(Z1)
#   RGCA <- cbind(Z1, Z2)
# }

tSCA <- function(P1, P2, data = NULL){
    if(!is.null(data)){
    P1Name <- deparse(substitute(P1))
    P2Name <- deparse(substitute(P2))
    P1 <- data[[P1Name]]
    P2 <- data[[P2Name]]
  }
  # Matrix tSCA: final version: 6/5/2020
  P1 <- factor(as.character(P1))
  P2 <- factor(as.character(P2))
  P1c <- as.character(P1); P2c <- as.character(P2)
  P1n <- as.numeric(P1); P2n <- as.numeric(P2)
  combination <- factor(apply(cbind(P1n*10 + P2n, P2n * 10 + P1n), 1, min))
  combLev <- factor( ifelse(P1c < P2c, paste(P1c, P2c, sep = ":"), paste(P2c, P1c, sep = ":") ) )
  mating <- factor(P1n*10 + P2n)
  p <- length(levels(factor(c(levels(P1), levels(P2)) )))
  n <- length(combination)

  # Step 1. gets the parameters to be estimated, removing
  # the unnecessary combinations
  last <- seq(10, p*10, 10) + p
  levs <- as.numeric(levels(combination))
  idx <- c() # Identifica la posizione degli ultimi
   for(i in 1:length(last)){
       #i <- 2
       y <- which(levs == last[i])
       idx[i] <- y
  }
  levs <- as.character(levs[-idx])
  SCA <- matrix(0, nrow = n, ncol = length(levs))
  colnames(SCA) <- levs
  colNamsOrd <- levels(combLev)[-idx]
  # Step 2. Insert 1s for all levels, but the last one
  for(i in 1:length(levs)){
    cond <- (combination == colnames(SCA)[i])*1
    SCA[, i] <- cond
   }
  # Step 3. Insert the -1s for the last level. The last level of
  # Par2, within each level of Par1. The last level of Par1
  # requires another step
  for(i in 1:(length(last) - 1)){
      start <- ceiling(ifelse(i == 1, 1, last[i-1])/10) * 10 +1
      arrival <- last[i]
      sel <- seq(start, arrival, 1)
      tmp <- as.character( sel[1:i] ) # Se necessario, inverte i reciproci
      splits <- strsplit(tmp, "")
      reversed <- lapply(splits, rev)
      tmp <- as.character(lapply(reversed, paste, collapse = ""))
      sel[1:i] <- as.numeric(tmp)
      sel
      idx <- c() # Identifica la posizione di quelli da scrivere
      for(j in 1:length(sel)){
           #i <- 7
           y <- which(levs == sel[j])
           if(length(y) > 0) idx[j] <- y
      }
      idx
      SCA[,idx]
      idx1 <- last[i]
      SCA[combination == idx1, idx] <- -1
  }
  # Scrive il self dell'ultimo livello
  #SCA[combination == last[p],] <- - apply(SCA, 2, sum)
  SCA[combination == last[p],] <- 2
  for(i in 1:p) {SCA[combination == last[p], colnames(SCA) == paste(i, i, sep = "")] <- 1}
  #colnames(SCA) <- paste("ts_", colnames(SCA), sep = "")
  colnames(SCA) <- paste("ts_", colNamsOrd, sep = "")
  row.names(SCA) <- c(1:length(SCA[,1]))
  SCA
}

# tSCA <- function(P1, P2){
#   # Matrix tSCA: final version: 6/5/2020
#   # P1 <- df$Par1; P2 <- df$Par2
#   P1 <- factor(as.character(P1))
#   P2 <- factor(as.character(P2))
#   P1n <- as.numeric(P1); P2n <- as.numeric(P2)
#   P1c <- as.character(P1); P2c <- as.character(P2)
#   combination <- factor(apply(cbind(P1n*10 + P2n, P2n * 10 + P1n), 1, min))
#   combLev <- factor( ifelse(P1c > P2c, paste(P1c, P2c, sep = ":"), paste(P2c, P1c, sep = ":") ) )
#   mating <- factor(P1n*10 + P2n)
#   p <- length(levels(factor(c(levels(P1), levels(P2)) )))
#   n <- length(combination)
#
#   # Step 1. gets the parameters to be estimated, removing
#   # the unnecessary combinations
#   last <- seq(10, p*10, 10) + p
#   #levs <- as.numeric(levels(combination))
#   levs <- levels(combination)
#   idx <- c() # Identifica la posizione degli ultimi
#    for(i in 1:length(last)){
#        #i <- 2
#        y <- which(levs == last[i])
#        idx[i] <- y
#   }
#   levs <- as.character(levs[-idx])
#   SCA <- matrix(0, nrow = n, ncol = length(levs))
#   # colnames(SCA) <- levs
#   colnames(SCA) <- as.character(combLev)[-idx]
#
#   # Step 2. Insert 1s for all levels, but the last one
#   for(i in 1:length(levs)){
#     cond <- (combination == colnames(SCA)[i])*1
#     SCA[, i] <- cond
#    }
#   # Step 3. Insert the -1s for the last level. The last level of
#   # Par2, within each level of Par1. The last level of Par1
#   # requires another step
#   for(i in 1:(length(last) - 1)){
#       start <- ceiling(ifelse(i == 1, 1, last[i-1])/10) * 10 +1
#       arrival <- last[i]
#       sel <- seq(start, arrival, 1)
#       tmp <- as.character( sel[1:i] ) # Se necessario, inverte i reciproci
#       splits <- strsplit(tmp, "")
#       reversed <- lapply(splits, rev)
#       tmp <- as.character(lapply(reversed, paste, collapse = ""))
#       sel[1:i] <- as.numeric(tmp)
#       sel
#       idx <- c() # Identifica la posizione di quelli da scrivere
#       for(j in 1:length(sel)){
#            #i <- 7
#            y <- which(levs == sel[j])
#            if(length(y) > 0) idx[j] <- y
#       }
#       idx
#       SCA[,idx]
#       idx1 <- last[i]
#       SCA[combination == idx1, idx] <- -1
#   }
#   # Scrive il self dell'ultimo livello
#   #SCA[combination == last[p],] <- - apply(SCA, 2, sum)
#   SCA[combination == last[p],] <- 2
#   for(i in 1:p) {SCA[combination == last[p], colnames(SCA) == paste(i, i, sep = "")] <- 1}
#   colnames(SCA) <- paste("ts_", colnames(SCA), sep = "")
#   SCA
# }

SCA <- function(P1, P2, data = NULL){
    if(!is.null(data)){
    P1Name <- deparse(substitute(P1))
    P2Name <- deparse(substitute(P2))
    P1 <- data[[P1Name]]
    P2 <- data[[P2Name]]
  }
  # Matrix for SCA in heterosis model Hayman2
  # P1 <- df$Par1; P2 <- df$Par2
  P1 <- factor(as.character(P1))
  P2 <- factor(as.character(P2))
  P1n <- as.numeric(P1); P2n <- as.numeric(P2)
  P1c <- as.character(P1); P2c <- as.character(P2)
  combination <- factor(apply(cbind(P1n*10 + P2n, P2n * 10 + P1n), 1, min))
  combLev <- factor( ifelse(P1c < P2c, paste(P1c, P2c, sep = ":"), paste(P2c, P1c, sep = ":") ) )
  mating <- factor(P1n*10 + P2n)
  p <- length(levels(factor(c(levels(P1), levels(P2)) )))
  n <- length(combination)
  last <- seq(10, p*10, 10) + p
  levs <- as.numeric(levels(combination))
  idx1 <- c()
  for(i in 1:length(last)){ #Indica l'ultimo
         y <- which(levs == last[i])
         idx1[i] <- y
  }
  idx2 <- c()
  for(i in 1:length(last)){ #Indica il penultimo
          y <- which(levs == (last[i] - 1))
          if(length(y) > 0) idx2[i] <- y
  }
  selflist <- as.numeric(levels(factor(combination[P1n == P2n])))
  idx3 <- c()
  for(i in 1:length(selflist)){ #Indica il penultimo
          #i <- 2
          y <- which(levs == selflist[i])
          if(length(y) > 0) idx3[i] <- y
  }

  idx <- c(idx1, idx3)
  rimossi <- levs[idx]
  levs <- as.character(levs[-idx])
  rimossi <- c(rimossi, levs[length(levs)])
  levs <- levs[-length(levs)]
  SCA <- matrix(0, nrow = n, ncol = length(levs))
  colnames(SCA) <- paste(levs)
  #colnames(SCA)
  colNamsOrd <- levels(combLev)[-idx][-length(levels(combLev)[-idx])]

    # Step 2. Insert 1s for all the levels, which are
    # in the SCA matrix
    for(i in 1:length(levs)){
          cond <- (combination == colnames(SCA)[i])*1
          SCA[, i] <- cond
    }
    # Step 3. Insert the -1s for the last level. The last level of
    # Par2, within each level of Par1. The last level of Par1
    # requires another step
    for(i in 1:(length(last) - 1)){
        start <- ceiling(ifelse(i == 1, 1, last[i-1])/10) * 10 +1
        arrival <- last[i]
        sel <- seq(start, arrival, 1)
        tmp <- as.character( sel[1:i] ) # Se necessario, inverte i reciproci
        splits <- strsplit(tmp, "")
        reversed <- lapply(splits, rev)
        tmp <- as.character(lapply(reversed, paste, collapse = ""))
        sel[1:i] <- as.numeric(tmp)
        #sel
        idx <- c() # Identifica la posizione di quelli da scrivere
        for(j in 1:length(sel)){
             #i <- 7
             y <- which(levs == sel[j])
             if(length(y) > 0) idx[j] <- y
        }
        idx
        SCA[,idx]
        idx1 <- last[i]
        SCA[combination == idx1, idx] <- -1
    }
     # SCA
     # Mancano le combinazioni degli ultimi 3 ibridi
     # 67, 68, 76, 78, 86, 87
     tmp <- seq(p-2, p, 1)
     tmp <- as.data.frame(combn(tmp, 2))
     tmp <- apply(tmp, 2, function(x) paste(x[1], x[2], sep = ""))
     SCA[combination == tmp[1], ] <- -1
     SCA[combination == tmp[2], ] <- SCA[combination == tmp[2], ] + 1
     SCA[combination == tmp[3], ] <- SCA[combination == tmp[3], ] + 1
     #colnames(SCA) <- paste("s_", colnames(SCA), sep = "")
     colnames(SCA) <- paste("s_", colNamsOrd, sep = "")
     SCA
}

SCA.GE <- function(P1, P2, data = NULL){
    if(!is.null(data)){
    P1Name <- deparse(substitute(P1))
    P2Name <- deparse(substitute(P2))
    P1 <- data[[P1Name]]
    P2 <- data[[P2Name]]
  }
  # SCA for GE models
  P1 <- factor(as.character(P1))
  P2 <- factor(as.character(P2))
  P1n <- as.numeric(P1); P2n <- as.numeric(P2)
  combination <- factor(apply(cbind(P1n*10 + P2n, P2n * 10 + P1n), 1, min))
  mating <- factor(P1n*10 + P2n)
  p <- length(levels(factor(c(levels(P1), levels(P2)) )))
  n <- length(combination)
  last <- seq(10, p*10, 10) + p
  levs <- as.numeric(levels(combination))
  idx1 <- c()
   for(i in 1:length(last)){ #Indica l'ultimo
       #i <- 2
       y <- which(levs == last[i])
       idx1[i] <- y
   }
  idx2 <- c()
   for(i in 1:length(last)){ #Indica il penultimo
       #i <- 2
       y <- which(levs == (last[i] - 1))
       if(length(y) > 0) idx2[i] <- y
   }
  selflist <- as.numeric(levels(factor(combination[P1n == P2n])))
  idx3 <- c()
   for(i in 1:length(selflist)){ #Indica il penultimo
       #i <- 2
       y <- which(levs == selflist[i])
       if(length(y) > 0) idx3[i] <- y
   }

  idx <- c(idx1, idx3)
  rimossi <- levs[idx]
  levs <- as.character(levs[-idx])
  rimossi <- c(rimossi, levs[length(levs)])
  levs <- levs[-length(levs)]
  SCA <- matrix(0, nrow = n, ncol = length(levs))
  colnames(SCA) <- paste(levs)
  # Step 2. Insert 1s for all the levels, which are
  # in the SCA matrix
  for(i in 1:length(levs)){
       cond <- (combination == colnames(SCA)[i])*1
       SCA[, i] <- cond
  }
  # Step 3. Insert the -1s for the last level. The last level of
  # Par2, within each level of Par1. The last level of Par1
  # requires another step
  for(i in 1:(length(last) - 1)){
      start <- ceiling(ifelse(i == 1, 1, last[i-1])/10) * 10 +1
      arrival <- last[i]
      sel <- seq(start, arrival, 1)
      tmp <- as.character( sel[1:i] ) # Se necessario, inverte i reciproci
      splits <- strsplit(tmp, "")
      reversed <- lapply(splits, rev)
      tmp <- as.character(lapply(reversed, paste, collapse = ""))
      sel[1:i] <- as.numeric(tmp)
      #sel
      idx <- c() # Identifica la posizione di quelli da scrivere
      for(j in 1:length(sel)){
           #i <- 7
           y <- which(levs == sel[j])
           if(length(y) > 0) idx[j] <- y
      }
      #idx
      SCA[,idx]
      idx1 <- last[i]
      SCA[combination == idx1, idx] <- -1
  }
  tmp <- seq(p-2, p, 1)
  tmp <- as.data.frame(combn(tmp, 2))
  tmp <- apply(tmp, 2, function(x) paste(x[1], x[2], sep = ""))
  SCA[combination == tmp[1], ] <- -1
  SCA[combination == tmp[2], ] <- SCA[combination == tmp[2], ] + 1
  SCA[combination == tmp[3], ] <- SCA[combination == tmp[3], ] + 1
  SCA
}


# SCA <- function(combination, p){
#   # Step 1. gets the parameters to be estimated, removing
#   # the unnecessary combinations
#   n <- length(combination)
#   last <- seq(10, p*10, 10) + p
#   levs <- as.numeric(levels(combination))
#   idx <- c() # Identifica la posizione degli ultimi
#   for(i in 1:length(last)){
#       #i <- 2
#       y <- which(levs == last[i])
#       idx[i] <- y
#     }
#   levs <- as.character(levs[-idx])
#   SCA <- matrix(0, nrow = n, ncol = length(levs))
#   colnames(SCA) <- levs
#   # Step 2. Insert 1s for all levels, but the last one
#   for(i in 1:length(levs)){
#       cond <- (combination == colnames(SCA)[i])*1
#       SCA[, i] <- cond
#       }
#   # Step 3. Insert the -1s for the last level of
#   # Par2, within each level of Par1. The last level of Par1
#   # requires another step
#   for(i in 1:(length(last) - 1)){
#       start <- ceiling(ifelse(i == 1, 1, last[i-1])/10) * 10 +1
#       arrival <- last[i]
#       sel <- seq(start, arrival, 1)
#       tmp <- as.character( sel[1:i] ) # Se necessario, inverte i reciproci
#       splits <- strsplit(tmp, "")
#       reversed <- lapply(splits, rev)
#       tmp <- as.character(lapply(reversed, paste, collapse = ""))
#       sel[1:i] <- as.numeric(tmp)
#       sel
#       idx <- c() # Identifica la posizione di quelli da scrivere
#       for(j in 1:length(sel)){
#         #i <- 7
#         y <- which(levs == sel[j])
#         if(length(y) > 0) idx[j] <- y
#         }
#       idx
#       SCA[,idx]
#       idx1 <- last[i]
#       SCA[combination == idx1, idx] <- -1
#   }
#   #SCA[combination == last[p],] <- - apply(SCA, 2, sum)/length(SCA[combination == last[p],1])
#   SCA[combination == last[p],] <- 2
#   for(i in 1:p) {SCA[combination == last[p], colnames(SCA) == paste(i, i, sep = "")] <- 1}
#   #SCA[combination == last[p], ]
#   SCA
# }

RSCA <- function(P1, P2, data = NULL){
    if(!is.null(data)){
    P1Name <- deparse(substitute(P1))
    P2Name <- deparse(substitute(P2))
    P1 <- data[[P1Name]]
    P2 <- data[[P2Name]]
  }
  # Derive the dummies and other infos
  # P1 <- df$Par1; P2 <- df$Par2
  P1 <- factor(as.character(P1))
  P2 <- factor(as.character(P2))
  n <- length(P1)
  p <- length(levels(P1))
  P1n <- as.numeric(P1); P2n <- as.numeric(P2)
  mate <- factor(P1n*10 + P2n)
  combination <- factor(apply(cbind(P1n*10 + P2n, P2n * 10 + P1n), 1, min))
  P1c <- as.character(P1); P2c <- as.character(P2)
  combLev <- factor( ifelse(P1c < P2c, paste(P1c, P2c, sep = ":"), paste(P2c, P1c, sep = ":") ) )

  # Empty matrix
  rec <- matrix(0, nrow = n, ncol = (p - 1)*(p - 2)/2 )
  cont <- 0
  nams <- c(); nams2 <- c()
  for(i in 1:(p-2)) { for(j in (i+1):(p-1)){
     cont <- cont + 1
     nams[cont] <- paste(i, j, sep="")
     nams2[cont] <- paste(levels(P1)[i], levels(P2)[j], sep = ":")
     }}
  colnames(rec) <- nams

  for(i in 1:(p-2)) { for(j in (i+1):(p-1)){
    cond <- paste(i,j, sep="")
     rec[mate == cond, colnames(rec) == cond] <- 1
     rec[combination == cond & mate != cond, colnames(rec) == cond] <- -1
  }}

  leftr <- substr(mate, 1, 1)
  rightr <- substr(mate, 2, 2)
  leftc <- substr(nams, 1, 1)
  rightc <- substr(nams, 2, 2)

  for(i in 1:7){
    #i <- 1
     rec[rightr==p & leftr == i, leftc == i] <- -1
     rec[rightr==p & leftr == i, rightc == i] <- 1
     rec[leftr==p & rightr == i, leftc == i] <- 1
     rec[leftr==p & rightr == i, rightc == i] <- -1
    }
  colnames(rec) <- paste("rs_", nams2, sep = "")
  rec
  }


REC <- function(P1, P2, data = NULL){
    if(!is.null(data)){
    P1Name <- deparse(substitute(P1))
    P2Name <- deparse(substitute(P2))
    P1 <- data[[P1Name]]
    P2 <- data[[P2Name]]
  }
  # P1 <- df$Par1; P2 <- df$Par2
  P1 <- factor(as.character(P1))
  P2 <- factor(as.character(P2))
  n <- length(P1)
  p <- length(levels(P1))
  P1n <- as.numeric(P1); P2n <- as.numeric(P2)
  P1c <- as.character(P1); P2c <- as.character(P2)
  mate <- factor(P1n*10 + P2n)
  combination <- factor(apply(cbind(P1n*10 + P2n, P2n * 10 + P1n), 1, min))
  dr <- ifelse(P1c == P2c, 0, ifelse(P1c < P2c, -1, 1))
  combLev <- factor( paste(P1c, P2c, sep = ":") )

  last <- c(); cont = 1
  for(i in 1:p){ for(j in 1:i) { last[cont] <- paste(i, j, sep=""); cont = cont + 1 } }
  last <- as.numeric(last)
  levs <- as.numeric(levels(mate))
  idx <- c()
    for(i in 1:length(last)){
    y <- which(levs == last[i])
    if(length(y) > 0) idx[i] <- y
    }
  levs <- as.character(levs[-idx])
  # levs
  rec <- matrix(0, nrow = n, ncol = length(levs))
  colnames(rec) <- paste(levs)
  colNamsOrd <- levels(combLev)[-idx]

    for(i in 1:length(levs)){
        cond <- (combination == colnames(rec)[i] ) * 1
        rec[, i] <- cond
    }
  rec <- rec*dr
  # colnames(rec) <- paste("r_", colnames(rec), sep = "")
  colnames(rec) <- paste("r_", colNamsOrd, sep = "")
  rec
}

H.BAR <- function(P1, P2, data = NULL){
    if(!is.null(data)){
    P1Name <- deparse(substitute(P1))
    P2Name <- deparse(substitute(P2))
    P1 <- data[[P1Name]]
    P2 <- data[[P2Name]]
  }
  P1 <- factor(as.character(P1))
  P2 <- factor(as.character(P2))
  P1c <- as.character(P1)
  P2c <- as.character(P2)
  cr <- ifelse(P1c == P2c, 0, 1)
  cr <- factor(cr)
  n <- length(cr)
  contrasts(cr) <- c("contr.treatment")
  crM <- model.matrix(~cr)
  crM <- crM[,-1]
  if(is.vector(crM) == T) crM <- matrix(crM, n, 1)
  colnames(crM) <- "h.bar"
  crM
}

MDD <- function(P1, P2, data = NULL){
    if(!is.null(data)){
    P1Name <- deparse(substitute(P1))
    P2Name <- deparse(substitute(P2))
    P1 <- data[[P1Name]]
    P2 <- data[[P2Name]]
  }
  P1 <- factor(as.character(P1))
  P2 <- factor(as.character(P2))
  p <- length(levels(P1))
  P1c <- as.character(P1)
  P2c <- as.character(P2)
  cr <- ifelse(P1c == P2c, 0, 1)
  cr <- factor(cr)
  n <- length(cr)
  contrasts(cr) <- c("contr.sum")
  crM <- model.matrix(~cr)
  crM <- crM[,-1]
  crM <- ifelse(crM == 1, - (p - 1), 1)
  if(is.vector(crM) == T) crM <- matrix(crM, n, 1)
  colnames(crM) <- "m"
  crM
}

# matHi <- function(P1, P2){
#   # For GE models ??
#   P1c <- as.character(P1)
#   P2c <- as.character(P2)
#   selfs <- ifelse(P1c == P2c, 1, 0)
#   contrasts(P1) <- c("contr.sum")
#   contrasts(P2) <- c("contr.sum")
#   Z1 <- model.matrix(~P1)
#   Z2 <- model.matrix(~P2)
#   H <- (Z1 + Z2) * selfs
#   H <- H[,-1]
# }

Hi <- function(P1, P2, data = NULL){
    if(!is.null(data)){
    P1Name <- deparse(substitute(P1))
    P2Name <- deparse(substitute(P2))
    P1 <- data[[P1Name]]
    P2 <- data[[P2Name]]
  }
  # For GE2 and GE3 models
  P1 <- factor(as.character(P1))
  P2 <- factor(as.character(P2))
  P1c <- as.character(P1)
  P2c <- as.character(P2)
  #selfs <- ifelse(P1c == P2c, 1, 0)
  crosses <- ifelse(P1c == P2c, 0, 1)
  contrasts(P1) <- c("contr.sum")
  contrasts(P2) <- c("contr.sum")
  Z1 <- model.matrix(~P1)
  Z2 <- model.matrix(~P2)
  H <- (Z1 + Z2) * crosses
  H <- H[,-1]
  nams <- paste("h_", levels(P1)[1:(length(levels(P1))-1)], sep="")
  colnames(H) <- c(nams)
  H
}

GCAC <- function(P1, P2, data = NULL){
    if(!is.null(data)){
    P1Name <- deparse(substitute(P1))
    P2Name <- deparse(substitute(P2))
    P1 <- data[[P1Name]]
    P2 <- data[[P2Name]]
  }
  # For GE2 and GE3 models
  P1 <- factor(as.character(P1))
  P2 <- factor(as.character(P2))
  P1c <- as.character(P1)
  P2c <- as.character(P2)
  #selfs <- ifelse(P1c == P2c, 1, 0)
  crosses <- ifelse(P1c == P2c, 0, 1)
  contrasts(P1) <- c("contr.sum")
  contrasts(P2) <- c("contr.sum")
  Z1 <- model.matrix(~P1)
  Z2 <- model.matrix(~P2)
  H <- (Z1 + Z2) * crosses
  H <- H[,-1]
  nams <- paste("gc_", levels(P1)[1:(length(levels(P1))-1)], sep="")
  colnames(H) <- c(nams)
  H
}


DD <- function(P1, P2, data = NULL){
    if(!is.null(data)){
    P1Name <- deparse(substitute(P1))
    P2Name <- deparse(substitute(P2))
    P1 <- data[[P1Name]]
    P2 <- data[[P2Name]]
  }
  # For Hyman model 2
  P1 <- factor(as.character(P1))
  P2 <- factor(as.character(P2))
  p <- length(levels(P1))
  contrasts(P1) <- c("contr.sum")
  contrasts(P2) <- c("contr.sum")
  Z1 <- model.matrix(~P1)
  Z2 <- model.matrix(~P2)
  H <- (Z1 + Z2)
  H[H == 2] <- -(p - 2)
  H[H == -2] <- (p - 2)
  H <- H[,-1]
  nams <- paste("d_", levels(P1)[1:(length(levels(P1))-1)], sep="")
  colnames(H) <- c(nams)
  H
}

