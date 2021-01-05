# Wrapper for mmer to fit random diallel models
# Last edited: 4/1/2021
mmer.diallel <- function(formula, Block = NULL, Env = NULL,
                       fct = "GRIFFING2", data){
  # formula <- Yield ~ Par1 + Par2
  # data <- Griffing56
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
  # Ci sono repliche ?
  cond <- tapply(Y, list(Par1, Par2), length)
  cond[is.na(cond)] <- 1
  withRep <- any(cond > 1)

  # Dummy variables
  if(!is.null(Block)){
    df <- data.frame(Y, Block, Par1, Par2)
  } else {
    df <- data.frame(Y, Par1, Par2)
  }

  df$dr <- ifelse(as.character(Par1) < as.character(Par2), -1,
              ifelse(as.character(Par1) == as.character(Par2), 0, 1))
  df$combination = factor(ifelse(as.character(Par1) <= as.character(Par2),
                              paste(Par1, Par2, sep =""),
                              paste(Par2, Par1, sep ="")))
  df$selfs = ifelse(as.character(Par1) == as.character(Par2), 1, 0)
  df$crosses <- ifelse(as.character(Par1) == as.character(Par2), 0, 1)
  # print(df)

  # HAYMAN1 ##############
  if(fct == "HAYMAN1"){
    if(!is.null(Block)){
      form <- ~ Block + overlay(Par1, Par2) + overlay(Par1, Par2):dr +
        combination + combination:dr # GCA + RGCA + tSCA + RSCA
      rnam <- c("Block", "GCA", "RGCA", "tSCA", "RSCA", "Residuals")
    } else if(is.null(Block) & withRep == T){
      form <- ~ overlay(Par1, Par2) + overlay(Par1, Par2):dr +
        combination + combination:dr # GCA + RGCA + tSCA + RSCA
      rnam <- c("GCA", "RGCA", "tSCA", "RSCA", "Residuals")
    } else if(is.null(Block) & withRep == F){
      form <- ~ overlay(Par1, Par2) + overlay(Par1, Par2):dr +
        combination # GCA + RGCA + tSCA + RSCA
      rnam <- c("GCA", "RGCA", "tSCA", "RSCA")
    }

    modh <- sommer::mmer(Y ~ 1,
             random = form,
             verbose = F,
             data = df)
    res <- summary(modh)$varcomp[,-c(3:4)]
    row.names(res) <- rnam

  # HAYMAN2 #########
  }else if(fct == "HAYMAN2"){
    if(!is.null(Block)){
      form <- ~ Block + crosses + overlay(Par1, Par2) + overlay(Par1, Par2):dr +
        overlay(Par1, Par2):selfs + combination + combination:dr
      rnam <- c("Block", "MDD", "GCA", "RGCA", "DD", "SCA", "RSCA", "Residuals")
    } else if(is.null(Block) & withRep == T){
      form <- ~ crosses + overlay(Par1, Par2) + overlay(Par1, Par2):dr +
        overlay(Par1, Par2):selfs + combination + combination:dr
      rnam <- c("MDD", "GCA", "RGCA", "DD", "SCA", "RSCA", "Residuals")
    } else if(is.null(Block) & withRep == F){
      form <- ~ crosses + overlay(Par1, Par2) + overlay(Par1, Par2):dr +
        overlay(Par1, Par2):selfs + combination
      rnam <- c("MDD", "GCA", "RGCA", "DD", "SCA", "RSCA")
    }
    modh <- sommer::mmer(Y ~ 1,
             random = form,
             verbose = F,
             data = df)
    res <- summary(modh)$varcomp[,-c(3:4)]
    row.names(res) <- rnam

  # GRIFFING1 and 3  ####
  }else if(fct == "GRIFFING1" | fct == "GRIFFING3"){
    if(!is.null(Block)){
      form <- ~ Block + overlay(Par1, Par2) + combination + combination:dr
      rnam <- c("Block", "GCA", "tSCA", "Reciprocals", "Residuals")
    } else if(is.null(Block) & withRep == T){
      form <- ~ overlay(Par1, Par2) + combination + combination:dr
      rnam <- c("GCA", "tSCA", "Reciprocals", "Residuals")
    } else if(is.null(Block) & withRep == F){
      form <- ~ overlay(Par1, Par2) + combination
      rnam <- c("GCA", "tSCA", "Reciprocals")
    }
    modh <- sommer::mmer(Y ~ 1,
             random = form,
             verbose = F,
             data = df)
    res <- summary(modh)$varcomp[,-c(3:4)]
    row.names(res) <- rnam

  # GRIFFING2 and 4 ####
  } else if(fct == "GRIFFING2" | fct == "GRIFFING4"){
    if(!is.null(Block)){
      form <- ~ Block + overlay(Par1, Par2) + combination
      rnam <- c("Block", "GCA", "tSCA", "Residuals")
    } else if(is.null(Block) & withRep == T){
      form <- ~ overlay(Par1, Par2) + combination
      rnam <- c("GCA", "tSCA", "Residuals")
    } else if(is.null(Block) & withRep == F){
      form <- ~ overlay(Par1, Par2)
      rnam <- c("GCA", "tSCA")
    }
    modh <- sommer::mmer(Y ~ 1,
             random = form,
             verbose = F,
             data = df)
    res <- summary(modh)$varcomp[,-c(3:4)]
    row.names(res) <- rnam

  # GE2r ############
  } else if(fct == "GE2r"){
    if(!is.null(Block)){
      form <- ~ Block + overlay(Par1, Par2) + overlay(Par1, Par2):crosses +
        combination:crosses + combination:crosses:dr
      rnam <- c("Block", "Variety", "h.i", "SCA", "Reciprocals", "Residuals")
    } else if(is.null(Block) & withRep == T){
      form <- ~ overlay(Par1, Par2) + overlay(Par1, Par2):crosses +
        combination:crosses + combination:crosses:dr
      rnam <- c("Variety", "h.i", "SCA", "Reciprocals", "Residuals")
    } else if(is.null(Block) & withRep == F){
      form <- ~ overlay(Par1, Par2) + overlay(Par1, Par2):crosses +
        combination:crosses
      rnam <- c("Variety", "h.i", "SCA", "Reciprocals")
    }
    modh <- sommer::mmer(Y ~ crosses,
             random = form,
             verbose = F,
             data = df)
    res <- summary(modh)$varcomp[,-c(3:4)]
    row.names(res) <- rnam

  # GE2 ######
  } else if(fct == "GE2"){
    if(!is.null(Block)){
      form <- ~ Block + overlay(Par1, Par2) + overlay(Par1, Par2):crosses +
        combination:crosses
      rnam <- c("Block", "Variety", "h.i", "SCA", "Residuals")
    } else if(is.null(Block) & withRep == T){
      form <- ~ overlay(Par1, Par2) + overlay(Par1, Par2):crosses +
        combination:crosses
      rnam <- c("Variety", "h.i", "SCA", "Residuals")
    } else if(is.null(Block) & withRep == F){
      form <- ~ overlay(Par1, Par2) + overlay(Par1, Par2):crosses
      rnam <- c("Variety", "h.i", "SCA")
    }
    modh <- sommer::mmer(Y ~ crosses,
             random = form,
             verbose = F,
             data = df)
    res <- summary(modh)$varcomp[,-c(3:4)]
    row.names(res) <- rnam

  # GE3r ####
  } else if(fct == "GE3r"){
    if(!is.null(Block)){
      form <- ~ Block + overlay(Par1, Par2):crosses + Par1:selfs +
        combination:crosses + combination:crosses:dr
      rnam <- c("Block", "GCAC", "Selfed par.", "SCA", "Reciprocals", "Residuals")
    } else if(is.null(Block) & withRep == T){
      form <- ~ overlay(Par1, Par2):crosses + Par1:selfs +
        combination:crosses + combination:crosses:dr
      rnam <- c("GCAC", "Selfed par.", "SCA", "Reciprocals", "Residuals")
    } else if(is.null(Block) & withRep == F){
      form <- ~ overlay(Par1, Par2):crosses + Par1:selfs +
        combination:crosses
      rnam <- c("GCAC", "Selfed par.", "SCA", "Reciprocals")
    }
    modh <- sommer::mmer(Y ~ crosses,
             random = form,
             verbose = F,
             data = df)
    res <- summary(modh)$varcomp[,-c(3:4)]
    row.names(res) <- rnam

  # GE3 ####
  } else if(fct == "GE3"){
    if(!is.null(Block)){
      form <- ~ Block + overlay(Par1, Par2):crosses + Par1:selfs +
        combination:crosses
      rnam <- c("Block", "GCAC", "Selfed par.", "SCA", "Residuals")
    } else if(is.null(Block) & withRep == T){
      form <- ~ overlay(Par1, Par2):crosses + Par1:selfs +
        combination:crosses
      rnam <- c("GCAC", "Selfed par.", "SCA", "Residuals")
    } else if(is.null(Block) & withRep == F){
      form <- ~ overlay(Par1, Par2):crosses + Par1:selfs
      rnam <- c("GCAC", "Selfed par.", "SCA")
    }
    modh <- sommer::mmer(Y ~ crosses,
             random = form,
             verbose = F,
             data = df)
    res <- summary(modh)$varcomp[,-c(3:4)]
    row.names(res) <- rnam
  }

 return(res)
}

