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
  withRep <- any(tapply(Y, list(Par1, Par2), length) > 1)

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
  df$selfs = ifelse(Par1 == Par2, 1, 0)
  df$crosses <- ifelse(df$Par1 == df$Par2, 0, 1)
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

  # GRIFFING1
  }else if(fct == "GRIFFING1" & !is.null(Block)){
    modh <- sommer::mmer(Y ~ 1,
             random = ~ Block
               + overlay(Par1, Par2) # GCA
               + combination # tSCA
               + combination:dr, verbose = F, # REC
             data = df)
    res <- summary(modh)$varcomp[,-c(3:4)]
    row.names(res) <- c("Block", "GCA", "tSCA", "Reciprocals", "Residuals")
  } else if(fct == "GRIFFING2" & !is.null(Block)){
    modh <- sommer::mmer(Y ~ 1,
             random = ~ Block
               + overlay(Par1, Par2) # GCA
               + combination, verbose = F, # tSCA
             data = df)
    res <- summary(modh)$varcomp[,-c(3:4)]
    row.names(res) <- c("Block", "GCA", "tSCA", "Residuals")
  } else if(fct == "GRIFFING3" & !is.null(Block)){
    modh <- sommer::mmer(Y ~ 1,
             random = ~ Block
               + overlay(Par1, Par2) # GCA
               + combination # SCA
               + combination:dr, verbose = F, # REC
             data = df)
    res <- summary(modh)$varcomp[,-c(3:4)]
    row.names(res) <- c("Block", "GCA", "SCA", "Reciprocals", "Residuals")
  } else if(fct == "GRIFFING4" & !is.null(Block)){
    modh <- sommer::mmer(Y ~ 1,
             random = ~ Block
               + overlay(Par1, Par2) # GCA
               + combination, # SCA
               verbose = F,
             data = df)
    res <- summary(modh)$varcomp[,-c(3:4)]
    row.names(res) <- c("Block", "GCA", "SCA", "Residuals")
  }

 return(res)
}

