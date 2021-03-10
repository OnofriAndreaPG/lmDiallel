# Wrapper for mmer to fit random diallel models
# Last edited: 23/2/2021
# Dummy variables for selfs, crosses, combinations
# crosses <- ifelse(Par1 == Par2, 0, 1)
# selfs <- ifelse(Par1 == Par2, 1, 0)
# dr <- ifelse(as.character(Par1) < as.character(Par2), -1,
#                 ifelse(as.character(Par1) == as.character(Par2), 0, 1))
# combination <- factor( ifelse(as.character(Par1) <= as.character(Par2),
#                                  paste(Par1, Par2, sep =""),
#                                  paste(Par2, Par1, sep ="")) )

mmer.diallel <- function(formula, Block = NULL, Env = NULL,
                       fct = "GRIFFING2", data, type = "all"){
  # formula <- Yield ~ Par1 + Par2
  # data <- Griffing56

  if(type != "all" & type != "environment"){
    print("The argument type only accepts 'all' or 'environment' values")
    stop()
  }

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
  # Type of design. Env? Block? Reps?


  # Dummy variables
  if(!is.null(Env)){
    # Do we have replicates?
    cond <- tapply(Y, list(Env, Par1, Par2), length)
    cond[is.na(cond)] <- 1
    withRep <- any(cond > 1)
    if(!is.null(Block)){
    df <- data.frame(Y, Env, Block, Par1, Par2)
    } else {
    df <- data.frame(Y, Env, Par1, Par2)
    }
  } else {
    cond <- tapply(Y, list(Par1, Par2), length)
    cond[is.na(cond)] <- 1
    withRep <- any(cond > 1)
    if(!is.null(Block)){
    df <- data.frame(Y, Block, Par1, Par2)
    } else {
    df <- data.frame(Y, Par1, Par2)
    }}
  df$dr <- ifelse(as.character(df$Par1) < as.character(df$Par2), -1,
              ifelse(as.character(df$Par1) == as.character(df$Par2), 0, 1))
  df$combination <- factor(ifelse(as.character(df$Par1) <= as.character(df$Par2),
                              paste(df$Par1, df$Par2, sep =""),
                              paste(df$Par2, df$Par1, sep ="")))
  df$selfs <- ifelse(as.character(df$Par1) == as.character(df$Par2), 1, 0)
  df$crosses <- ifelse(as.character(df$Par1) == as.character(df$Par2), 0, 1)

  # Model specifications ########
  if(!is.null(Env) & type == "environment"){
  ## Environment random and genetical effects fixed #####
    # HAYMAN1 ####
  if(fct == "HAYMAN1"){
    formFix <- Y ~ GCA(Par1, Par2) + tSCA(Par1, Par2) + RGCA(Par1, Par2) +
      RSCA(Par1, Par2)
    if(!is.null(Block)){
      form <- ~ Env + Env:Block + GCA(Par1:Env, Par2:Env, type = "random") +
        tSCA(Par1:Env, Par2:Env, type = "random") +
        RGCA(Par1:Env, Par2:Env, type = "random") +
        RSCA(Par1:Env, Par2:Env, type = "random")
      rnam <- c("Env", "Env:Block", "GCA:Env", "tSCA:Env", "RGCA:Env", "RSCA:Env", "Residuals")
    } else if(is.null(Block) & withRep == T){
      form <- ~ Env + GCA(Par1:Env, Par2:Env, type = "random") +
        tSCA(Par1:Env, Par2:Env, type = "random") +
        RGCA(Par1:Env, Par2:Env, type = "random") +
        RSCA(Par1:Env, Par2:Env, type = "random")
      rnam <- c("Env", "GCA:Env", "tSCA:Env", "RGCA:Env", "RSCA:Env", "Residuals")
    } else if(is.null(Block) & withRep == F){
      form <- ~ Env + GCA(Par1:Env, Par2:Env, type = "random") +
        tSCA(Par1:Env, Par2:Env, type = "random") +
        RGCA(Par1:Env, Par2:Env, type = "random")
      rnam <- c("Env", "GCA:Env", "tSCA:Env", "RGCA:Env", "RSCA:Env")
    }
  # HAYMAN2 #########
  }else if(fct == "HAYMAN2"){
    print("Yet to be implemented. Sorry!")
    stop()

  # GRIFFING 1   ####
  }else if(fct == "GRIFFING1"){
   formFix <- Y ~ GCA(Par1, Par2) + tSCA(Par1, Par2) + REC(Par1, Par2)
    if(!is.null(Block)){
      form <- ~ Env + Env:Block + GCA(Par1:Env, Par2:Env, type = "random") +
        tSCA(Par1:Env, Par2:Env, type = "random") +
        REC(Par1:Env, Par2:Env, type = "random")
      rnam <- c("Env", "Env:Block", "GCA:Env", "tSCA:Env", "REC:Env", "Residuals")
    } else if(is.null(Block) & withRep == T){
      form <- ~ Env + GCA(Par1:Env, Par2:Env, type = "random") +
        tSCA(Par1:Env, Par2:Env, type = "random") +
        REC(Par1:Env, Par2:Env, type = "random")
      rnam <- c("Env", "GCA:Env", "tSCA:Env", "REC:Env", "Residuals")
    } else if(is.null(Block) & withRep == F){
      form <- ~ Env + GCA(Par1:Env, Par2:Env, type = "random") +
        tSCA(Par1:Env, Par2:Env, type = "random")
      rnam <- c("Env", "GCA:Env", "tSCA:Env", "REC:Env")
    }

   # GRIFFING 3  ####
  } else if(fct == "GRIFFING3"){
   formFix <- Y ~ GCA(Par1, Par2) + SCA.G3(Par1, Par2) + REC.G3(Par1, Par2)
    if(!is.null(Block)){
      form <- ~ Env + Env:Block + GCA(Par1:Env, Par2:Env, type = "random") +
        tSCA(Par1:Env, Par2:Env, type = "random") +
        REC(Par1:Env, Par2:Env, type = "random")
      rnam <- c("Env", "Env:Block", "GCA:Env", "SCA:Env", "REC:Env", "Residuals")
    } else if(is.null(Block) & withRep == T){
      form <- ~ Env + GCA(Par1:Env, Par2:Env, type = "random") +
        tSCA(Par1:Env, Par2:Env, type = "random") +
        REC(Par1:Env, Par2:Env, type = "random")
      rnam <- c("Env", "GCA:Env", "SCA:Env", "REC:Env", "Residuals")
    } else if(is.null(Block) & withRep == F){
      form <- ~ Env + GCA(Par1:Env, Par2:Env, type = "random") +
        tSCA(Par1:Env, Par2:Env, type = "random")
      rnam <- c("Env", "GCA:Env", "SCA:Env", "REC:Env")
    }

  # GRIFFING 2 ####
  } else if(fct == "GRIFFING2"){
      formFix <- Y ~ GCA(Par1, Par2) + tSCA(Par1, Par2)
    if(!is.null(Block)){
      form <- ~ Env + Env:Block + GCA(Par1:Env, Par2:Env, type = "random") +
        tSCA(Par1:Env, Par2:Env, type = "random")
      rnam <- c("Env", "Env:Block", "GCA:Env", "tSCA:Env","Residuals")
    } else if(is.null(Block) & withRep == T){
      form <- ~ Env + GCA(Par1:Env, Par2:Env, type = "random") +
        tSCA(Par1:Env, Par2:Env, type = "random")
      rnam <- c("Env", "GCA:Env", "tSCA:Env", "Residuals")
    } else if(is.null(Block) & withRep == F){
      form <- ~ Env + GCA(Par1:Env, Par2:Env, type = "random")
      rnam <- c("Env", "GCA:Env", "tSCA:Env")
    }

    # GRIFFING 4 ####
  } else if(fct == "GRIFFING4"){
    formFix <- Y ~ GCA(Par1, Par2) + SCA.G3(Par1, Par2)
    if(!is.null(Block)){
      form <- ~ Env + Env:Block + GCA(Par1:Env, Par2:Env, type = "random") +
        tSCA(Par1:Env, Par2:Env, type = "random")
      rnam <- c("Env", "Env:Block", "GCA:Env", "SCA:Env","Residuals")
    } else if(is.null(Block) & withRep == T){
      form <- ~ Env + GCA(Par1:Env, Par2:Env, type = "random") +
        tSCA(Par1:Env, Par2:Env, type = "random")
      rnam <- c("Env", "GCA:Env", "SCA:Env", "Residuals")
    } else if(is.null(Block) & withRep == F){
      form <- ~ Env + GCA(Par1:Env, Par2:Env, type = "random")
      rnam <- c("Env", "GCA:Env", "SCA:Env")
    }

  # GE2r ############
  }else if(fct == "GE2r"){
     formFix <- Y ~ H.BAR(Par1, Par2) + VEi(Par1, Par2) + Hi(Par1, Par2) +
       SCA(Par1, Par2) + REC(Par1, Par2)
    if(!is.null(Block)){
      form <- ~ Env + Env:Block +
        crosses:Env +
        VEi(Par1:Env, Par2:Env, type = "random") +
        Hi(Par1:Env, Par2:Env, type = "random") +
        SCA(Par1:Env, Par2:Env, type = "random") +
        REC(Par1:Env, Par2:Env, type = "random")
      rnam <- c("Env", "Env:Block", "H:BAR:Env", "VEi:Env", "Hi:Env",
                "SCA:Env", "REC:Env", "Residuals")
    } else if(is.null(Block) & withRep == T){
      form <- ~ Env +
        crosses:Env +
        VEi(Par1:Env, Par2:Env, type = "random") +
        Hi(Par1:Env, Par2:Env, type = "random") +
        SCA(Par1:Env, Par2:Env, type = "random") +
        REC(Par1:Env, Par2:Env, type = "random")
      rnam <- c("Env", "H:BAR:Env", "VEi:Env", "Hi:Env",
                "SCA:Env", "REC:Env", "Residuals")
    } else if(is.null(Block) & withRep == F){
      form <- ~ Env +
        crosses:Env +
        VEi(Par1:Env, Par2:Env, type = "random") +
        Hi(Par1:Env, Par2:Env, type = "random") +
        SCA(Par1:Env, Par2:Env, type = "random")
      rnam <- c("Env", "H:BAR:Env", "VEi:Env", "Hi:Env",
                "SCA:Env", "REC:Env")
    }

  # GE2 ######
  } else if(fct == "GE2"){
     formFix <- Y ~ H.BAR(Par1, Par2) + VEi(Par1, Par2) + Hi(Par1, Par2) +
       SCA(Par1, Par2)
    if(!is.null(Block)){
      form <- ~ Env + Env:Block + crosses:Env +
        VEi(Par1:Env, Par2:Env, type = "random") +
        Hi(Par1:Env, Par2:Env, type = "random") +
        SCA(Par1:Env, Par2:Env, type = "random")
      rnam <- c("Env", "Env:Block", "H:BAR:Env", "VEi:Env", "Hi:Env",
                "SCA:Env", "Residuals")
    } else if(is.null(Block) & withRep == T){
      form <- ~ Env + crosses:Env +
        VEi(Par1:Env, Par2:Env, type = "random") +
        Hi(Par1:Env, Par2:Env, type = "random") +
        SCA(Par1:Env, Par2:Env, type = "random")
      rnam <- c("Env", "H:BAR:Env", "VEi:Env", "Hi:Env",
                "SCA:Env", "Residuals")
    }  else if(is.null(Block) & withRep == F){
      form <- ~ Env + crosses:Env +
        VEi(Par1:Env, Par2:Env, type = "random") +
        Hi(Par1:Env, Par2:Env, type = "random")
      rnam <- c("Env", "H:BAR:Env", "VEi:Env", "Hi:Env",
                "SCA:Env")
    }


  # GE3r ####
  } else if(fct == "GE3r"){
     formFix <- Y ~ H.BAR(Par1, Par2) + SP(Par1, Par2) + GCAC(Par1, Par2) +
       SCA(Par1, Par2) + REC(Par1, Par2)
    if(!is.null(Block)){
      form <- ~ Env + Env:Block +
        crosses:Env +
        SP(Par1:Env, Par2:Env, type = "random") +
        GCAC(Par1:Env, Par2:Env, type = "random") +
        SCA(Par1:Env, Par2:Env, type = "random") +
        REC(Par1:Env, Par2:Env, type = "random")
      rnam <- c("Env", "Env:Block", "H:BAR:Env", "Selfs:Env", "GCAC:Env",
                "SCA:Env", "REC:Env", "Residuals")
    } else if(is.null(Block) & withRep == T){
      form <- ~ Env +
        crosses:Env +
        SP(Par1:Env, Par2:Env, type = "random") +
        GCAC(Par1:Env, Par2:Env, type = "random") +
        SCA(Par1:Env, Par2:Env, type = "random") +
        REC(Par1:Env, Par2:Env, type = "random")
      rnam <- c("Env", "H:BAR:Env", "Selfs:Env", "GCAC:Env",
                "SCA:Env", "REC:Env", "Residuals")
    } else if(is.null(Block) & withRep == F){
      form <- ~ Env +
        crosses:Env +
        VEi(Par1:Env, Par2:Env, type = "random") +
        Hi(Par1:Env, Par2:Env, type = "random") +
        SCA(Par1:Env, Par2:Env, type = "random")
      rnam <- c("Env", "H:BAR:Env", "Selfs:Env", "GCAC:Env",
                "SCA:Env", "REC:Env")
    }

  # GE3 ####
  } else if(fct == "GE3"){
    formFix <- Y ~ H.BAR(Par1, Par2) + SP(Par1, Par2) + GCAC(Par1, Par2) +
       SCA(Par1, Par2)
    if(!is.null(Block)){
      form <- ~ Env + Env:Block +
        crosses:Env +
        SP(Par1:Env, Par2:Env, type = "random") +
        GCAC(Par1:Env, Par2:Env, type = "random") +
        SCA(Par1:Env, Par2:Env, type = "random")
      rnam <- c("Env", "Env:Block", "H:BAR:Env", "Selfs:Env", "GCAC:Env",
                "SCA:Env", "Residuals")
    } else if(is.null(Block) & withRep == T){
      form <- ~ Env +
        crosses:Env +
        SP(Par1:Env, Par2:Env, type = "random") +
        GCAC(Par1:Env, Par2:Env, type = "random") +
        SCA(Par1:Env, Par2:Env, type = "random")
      rnam <- c("Env", "H:BAR:Env", "Selfs:Env", "GCAC:Env",
                "SCA:Env", "Residuals")
    } else if(is.null(Block) & withRep == F){
      form <- ~ Env +
        crosses:Env +
        VEi(Par1:Env, Par2:Env, type = "random") +
        Hi(Par1:Env, Par2:Env, type = "random")
      rnam <- c("Env", "H:BAR:Env", "Selfs:Env", "GCAC:Env",
                "SCA:Env")
    }}

  }else if(!is.null(Env) & type == "all"){
  # Environment random and genetical effects random ######
 # HAYMAN1 ####
  if(fct == "HAYMAN1"){
    formFix <- Y ~ 1
    if(!is.null(Block)){
      form <- ~ Env + Env:Block +
        GCA(Par1, Par2, type = "random") +
        tSCA(Par1, Par2, type = "random") +
        RGCA(Par1, Par2, type = "random") +
        RSCA(Par1, Par2, type = "random") +
        GCA(Par1:Env, Par2:Env, type = "random") +
        tSCA(Par1:Env, Par2:Env, type = "random") +
        RGCA(Par1:Env, Par2:Env, type = "random") +
        RSCA(Par1:Env, Par2:Env, type = "random")
      rnam <- c("Env", "Env:Block", "GCA", "tSCA", "RGCA", "RSCA", "GCA:Env", "tSCA:Env", "RGCA:Env", "RSCA:Env", "Residuals")
    } else if(is.null(Block) & withRep == T){
      form <- ~ Env +
        GCA(Par1, Par2, type = "random") +
        tSCA(Par1, Par2, type = "random") +
        RGCA(Par1, Par2, type = "random") +
        RSCA(Par1, Par2, type = "random") +
        GCA(Par1:Env, Par2:Env, type = "random") +
        tSCA(Par1:Env, Par2:Env, type = "random") +
        RGCA(Par1:Env, Par2:Env, type = "random") +
        RSCA(Par1:Env, Par2:Env, type = "random")
      rnam <- c("Env", "GCA", "tSCA", "RGCA", "RSCA", "GCA:Env", "tSCA:Env", "RGCA:Env", "RSCA:Env", "Residuals")
    } else if(is.null(Block) & withRep == F){
      form <- ~ Env +
        GCA(Par1, Par2, type = "random") +
        tSCA(Par1, Par2, type = "random") +
        RGCA(Par1, Par2, type = "random") +
        RSCA(Par1, Par2, type = "random") +
        GCA(Par1:Env, Par2:Env, type = "random") +
        tSCA(Par1:Env, Par2:Env, type = "random") +
        RGCA(Par1:Env, Par2:Env, type = "random")
      rnam <- c("Env", "GCA", "tSCA", "RGCA", "RSCA", "GCA:Env", "tSCA:Env", "RGCA:Env", "RSCA:Env")
    }

  # HAYMAN2 #########
  }else if(fct == "HAYMAN2"){
    print("Yet to be implemented. Sorry!")
    stop()

  # GRIFFING 1   ####
  }else if(fct == "GRIFFING1"){
   formFix <- Y ~ 1
    if(!is.null(Block)){
      form <- ~ Env + Env:Block +
        GCA(Par1, Par2, type = "random") +
        tSCA(Par1, Par2, type = "random") +
        REC(Par1, Par2, type = "random") +
        GCA(Par1:Env, Par2:Env, type = "random") +
        tSCA(Par1:Env, Par2:Env, type = "random") +
        REC(Par1:Env, Par2:Env, type = "random")
      rnam <- c("Env", "Env:Block", "GCA", "tSCA", "REC", "GCA:Env", "tSCA:Env", "REC:Env", "Residuals")
    } else if(is.null(Block) & withRep == T){
      form <- ~ Env +
        GCA(Par1, Par2, type = "random") +
        tSCA(Par1, Par2, type = "random") +
        REC(Par1, Par2, type = "random") +
        GCA(Par1:Env, Par2:Env, type = "random") +
        tSCA(Par1:Env, Par2:Env, type = "random") +
        REC(Par1:Env, Par2:Env, type = "random")
      rnam <- c("Env", "GCA", "tSCA", "REC", "GCA:Env", "tSCA:Env", "REC:Env", "Residuals")
    } else if(is.null(Block) & withRep == F){
      form <- ~ Env +
        GCA(Par1, Par2, type = "random") +
        tSCA(Par1, Par2, type = "random") +
        REC(Par1, Par2, type = "random") +
        GCA(Par1:Env, Par2:Env, type = "random") +
        tSCA(Par1:Env, Par2:Env, type = "random")
      rnam <- c("Env", "GCA", "tSCA", "REC",  "GCA:Env", "tSCA:Env", "REC:Env")
    }

   # GRIFFING 3  ####
  }else if(fct == "GRIFFING3"){
   formFix <- Y ~ 1
    if(!is.null(Block)){
      form <- ~ Env + Env:Block +
        GCA(Par1, Par2, type = "random") +
        SCA(Par1, Par2, type = "random") +
        REC(Par1, Par2, type = "random") +
        GCA(Par1:Env, Par2:Env, type = "random") +
        tSCA(Par1:Env, Par2:Env, type = "random") +
        REC(Par1:Env, Par2:Env, type = "random")
      rnam <- c("Env", "Env:Block", "GCA", "SCA", "REC", "GCA:Env", "SCA:Env", "REC:Env", "Residuals")
    } else if(is.null(Block) & withRep == T){
      form <- ~ Env +
        GCA(Par1, Par2, type = "random") +
        SCA(Par1, Par2, type = "random") +
        REC(Par1, Par2, type = "random") +
        GCA(Par1:Env, Par2:Env, type = "random") +
        tSCA(Par1:Env, Par2:Env, type = "random") +
        REC(Par1:Env, Par2:Env, type = "random")
      rnam <- c("Env", "GCA", "SCA", "REC", "GCA:Env", "SCA:Env", "REC:Env", "Residuals")
    } else if(is.null(Block) & withRep == F){
      form <- ~ Env +
        GCA(Par1, Par2, type = "random") +
        SCA(Par1, Par2, type = "random") +
        REC(Par1, Par2, type = "random") +
        GCA(Par1:Env, Par2:Env, type = "random") +
        tSCA(Par1:Env, Par2:Env, type = "random")
      rnam <- c("Env", "GCA", "SCA", "REC", "GCA:Env", "tSCA:Env", "REC:Env")
    }

  # GRIFFING 2 ####
  } else if(fct == "GRIFFING2"){

    formFix <- Y ~ 1
    if(!is.null(Block)){
      form <- ~ Env + Env:Block +
        GCA(Par1, Par2, type = "random") + tSCA(Par1, Par2, type = "random") +
        GCA(Par1:Env, Par2:Env, type = "random") +
        tSCA(Par1:Env, Par2:Env, type = "random")
      rnam <- c("Env", "Env:Block", "GCA", "tSCA", "GCA:Env", "tSCA:Env","Residuals")
    } else if(is.null(Block) & withRep == T){
      form <- ~ Env + GCA(Par1:Env, Par2:Env, type = "random") +
        GCA(Par1, Par2, type = "random") + tSCA(Par1, Par2, type = "random") +
        GCA(Par1:Env, Par2:Env, type = "random") +
        tSCA(Par1:Env, Par2:Env, type = "random")
      rnam <- c("Env", "GCA", "tSCA", "GCA:Env", "tSCA:Env", "Residuals")
    } else if(is.null(Block) & withRep == F){
      form <- ~ Env +
        GCA(Par1, Par2, type = "random") + tSCA(Par1, Par2, type = "random") +
        GCA(Par1:Env, Par2:Env, type = "random")
      rnam <- c("Env", "GCA", "tSCA", "GCA:Env", "tSCA:Env")
     }

    # GRIFFING 4 ####
  } else if(fct == "GRIFFING4"){
    formFix <- Y ~ 1
    if(!is.null(Block)){
      form <- ~ Env + Env:Block +
        GCA(Par1, Par2, type = "random") + SCA(Par1, Par2, type = "random") +
      GCA(Par1:Env, Par2:Env, type = "random") +
        tSCA(Par1:Env, Par2:Env, type = "random")
      rnam <- c("Env", "Env:Block", "GCA", "SCA", "GCA:Env", "tSCA:Env","Residuals")
    } else if(is.null(Block) & withRep == T){
      form <- ~ Env +
        GCA(Par1, Par2, type = "random") + SCA(Par1, Par2, type = "random") +
        GCA(Par1:Env, Par2:Env, type = "random") +
        tSCA(Par1:Env, Par2:Env, type = "random")
      rnam <- c("Env", "GCA", "SCA", "GCA:Env", "tSCA:Env", "Residuals")
    } else if(is.null(Block) & withRep == F){
      form <- ~ Env +
        GCA(Par1, Par2, type = "random") + SCA(Par1, Par2, type = "random") +
        GCA(Par1:Env, Par2:Env, type = "random")
      rnam <- c("Env", "GCA", "SCA", "GCA:Env", "tSCA:Env")
    }

  # GE2r ############
  } else if(fct == "GE2r"){
     formFix <- Y ~ H.BAR(Par1, Par2)
    if(!is.null(Block)){
      form <- ~ Env + Env:Block +
        crosses:Env +
        VEi(Par1, Par2, type = "random") +
        Hi(Par1, Par2, type = "random") +
        SCA(Par1, Par2, type = "random") +
        REC(Par1, Par2, type = "random") +
        VEi(Par1:Env, Par2:Env, type = "random") +
        Hi(Par1:Env, Par2:Env, type = "random") +
        SCA(Par1:Env, Par2:Env, type = "random") +
        REC(Par1:Env, Par2:Env, type = "random")
      rnam <- c("Env", "Env:Block", "VEi", "Hi", "SCA", "REC", "H:BAR:Env", "VEi:Env", "Hi:Env",
                "SCA:Env", "REC:Env", "Residuals")
    } else if(is.null(Block) & withRep == T){
      form <- ~ Env +
        crosses:Env +
        VEi(Par1, Par2, type = "random") +
        Hi(Par1, Par2, type = "random") +
        SCA(Par1, Par2, type = "random") +
        REC(Par1, Par2, type = "random") +
        VEi(Par1:Env, Par2:Env, type = "random") +
        Hi(Par1:Env, Par2:Env, type = "random") +
        SCA(Par1:Env, Par2:Env, type = "random") +
        REC(Par1:Env, Par2:Env, type = "random")
      rnam <- c("Env","VEi", "Hi", "SCA", "REC", "H:BAR:Env", "VEi:Env", "Hi:Env",
                "SCA:Env", "REC:Env", "Residuals")
    } else if(is.null(Block) & withRep == F){
      form <- ~ Env +
        crosses:Env +
        VEi(Par1, Par2, type = "random") +
        Hi(Par1, Par2, type = "random") +
        SCA(Par1, Par2, type = "random") +
        REC(Par1, Par2, type = "random") +
        VEi(Par1:Env, Par2:Env, type = "random") +
        Hi(Par1:Env, Par2:Env, type = "random") +
        SCA(Par1:Env, Par2:Env, type = "random")
      rnam <- c("Env", "VEi", "Hi", "SCA", "REC", "H:BAR:Env", "VEi:Env", "Hi:Env",
                "SCA:Env", "REC:Env")
    }

  # GE2 ######
  } else if(fct == "GE2"){
      formFix <- Y ~ H.BAR(Par1, Par2)
    if(!is.null(Block)){
      form <- ~ Env + Env:Block +
        crosses:Env +
        VEi(Par1, Par2, type = "random") +
        Hi(Par1, Par2, type = "random") +
        SCA(Par1, Par2, type = "random") +
        VEi(Par1:Env, Par2:Env, type = "random") +
        Hi(Par1:Env, Par2:Env, type = "random") +
        SCA(Par1:Env, Par2:Env, type = "random")
      rnam <- c("Env", "Env:Block", "VEi", "Hi", "SCA", "H:BAR:Env", "VEi:Env", "Hi:Env",
                "SCA:Env", "Residuals")
    } else if(is.null(Block) & withRep == T){
      form <- ~ Env +
        crosses:Env +
        VEi(Par1, Par2, type = "random") +
        Hi(Par1, Par2, type = "random") +
        SCA(Par1, Par2, type = "random") +
        VEi(Par1:Env, Par2:Env, type = "random") +
        Hi(Par1:Env, Par2:Env, type = "random") +
        SCA(Par1:Env, Par2:Env, type = "random")
      rnam <- c("Env","VEi", "Hi", "SCA", "H:BAR:Env", "VEi:Env", "Hi:Env",
                "SCA:Env", "Residuals")
    } else if(is.null(Block) & withRep == F){
      form <- ~ Env +
        crosses:Env +
        VEi(Par1, Par2, type = "random") +
        Hi(Par1, Par2, type = "random") +
        SCA(Par1, Par2, type = "random") +
        VEi(Par1:Env, Par2:Env, type = "random") +
        Hi(Par1:Env, Par2:Env, type = "random")
      rnam <- c("Env", "VEi", "Hi", "SCA", "H:BAR:Env", "VEi:Env", "Hi:Env",
                "SCA:Env")
    }

  # GE3r ####
  } else if(fct == "GE3r"){
     formFix <- Y ~ H.BAR(Par1, Par2)
    if(!is.null(Block)){
      form <- ~ Env + Env:Block +
        SP(Par1, Par2, type = "random") +
        GCAC(Par1, Par2, type = "random") +
        SCA(Par1, Par2, type = "random") +
        REC(Par1, Par2, type = "random") +
        crosses:Env +
        SP(Par1:Env, Par2:Env, type = "random") +
        GCAC(Par1:Env, Par2:Env, type = "random") +
        SCA(Par1:Env, Par2:Env, type = "random") +
        REC(Par1:Env, Par2:Env, type = "random")
      rnam <- c("Env", "Env:Block", "SP", "GCAC", "SCA", "REC", "H:BAR:Env", "Selfs:Env", "GCAC:Env",
                "SCA:Env", "REC:Env", "Residuals")
    } else if(is.null(Block) & withRep == T){
      form <- ~ Env +
        SP(Par1, Par2, type = "random") +
        GCAC(Par1, Par2, type = "random") +
        SCA(Par1, Par2, type = "random") +
        REC(Par1, Par2, type = "random") +
        crosses:Env +
        SP(Par1:Env, Par2:Env, type = "random") +
        GCAC(Par1:Env, Par2:Env, type = "random") +
        SCA(Par1:Env, Par2:Env, type = "random") +
        REC(Par1:Env, Par2:Env, type = "random")
      rnam <- c("Env",  "SP", "GCAC", "SCA", "REC", "H:BAR:Env", "Selfs:Env", "GCAC:Env",
                "SCA:Env", "REC:Env", "Residuals")
    } else if(is.null(Block) & withRep == F){
      form <- ~ Env +
        SP(Par1, Par2, type = "random") +
        GCAC(Par1, Par2, type = "random") +
        SCA(Par1, Par2, type = "random") +
        REC(Par1, Par2, type = "random") +
        crosses:Env +
        VEi(Par1:Env, Par2:Env, type = "random") +
        Hi(Par1:Env, Par2:Env, type = "random") +
        SCA(Par1:Env, Par2:Env, type = "random")
      rnam <- c("Env",  "SP", "GCAC", "SCA", "REC", "H:BAR:Env", "Selfs:Env", "GCAC:Env",
                "SCA:Env", "REC:Env")
    }


  # GE3 ####
  } else if(fct == "GE3"){
      formFix <- Y ~ H.BAR(Par1, Par2)
    if(!is.null(Block)){
      form <- ~ Env + Env:Block +
        SP(Par1, Par2, type = "random") +
        GCAC(Par1, Par2, type = "random") +
        SCA(Par1, Par2, type = "random") +
        crosses:Env +
        SP(Par1:Env, Par2:Env, type = "random") +
        GCAC(Par1:Env, Par2:Env, type = "random") +
        SCA(Par1:Env, Par2:Env, type = "random")
      rnam <- c("Env", "Env:Block", "SP", "GCAC", "SCA", "H:BAR:Env", "Selfs:Env", "GCAC:Env",
                "SCA:Env", "Residuals")
    } else if(is.null(Block) & withRep == T){
      form <- ~ Env +
        SP(Par1, Par2, type = "random") +
        GCAC(Par1, Par2, type = "random") +
        SCA(Par1, Par2, type = "random") +
        crosses:Env +
        SP(Par1:Env, Par2:Env, type = "random") +
        GCAC(Par1:Env, Par2:Env, type = "random") +
        SCA(Par1:Env, Par2:Env, type = "random")
      rnam <- c("Env",  "SP", "GCAC", "SCA", "H:BAR:Env", "Selfs:Env", "GCAC:Env",
                "SCA:Env", "Residuals")
    } else if(is.null(Block) & withRep == F){
      form <- ~ Env +
        SP(Par1, Par2, type = "random") +
        GCAC(Par1, Par2, type = "random") +
        SCA(Par1, Par2, type = "random") +
        crosses:Env +
        VEi(Par1:Env, Par2:Env, type = "random") +
        Hi(Par1:Env, Par2:Env, type = "random")
      rnam <- c("Env",  "SP", "GCAC", "SCA", "H:BAR:Env", "Selfs:Env", "GCAC:Env",
                "SCA:Env")
    }}
  } else {
  ## No environment and genetical effects random ####

  # HAYMAN1 ######
  if(fct == "HAYMAN1"){
    formFix <- Y ~ 1
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

  # HAYMAN2 #########
  }else if(fct == "HAYMAN2"){
    formFix <- Y ~ 1
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

  # GRIFFING1 and 3  ####
  }else if(fct == "GRIFFING1"){
    formFix <- Y ~ 1
    if(!is.null(Block)){
      form <- ~ Block + overlay(Par1, Par2) + combination + combination:dr
      rnam <- c("Block", "GCA", "tSCA", "REC", "Residuals")
    } else if(is.null(Block) & withRep == T){
      form <- ~ overlay(Par1, Par2) + combination + combination:dr
      rnam <- c("GCA", "tSCA", "REC", "Residuals")
    } else if(is.null(Block) & withRep == F){
      form <- ~ overlay(Par1, Par2) + combination
      rnam <- c("GCA", "tSCA", "REC")
    }

  }else if(fct == "GRIFFING3"){
    formFix <- Y ~ 1
    if(!is.null(Block)){
      form <- ~ Block + overlay(Par1, Par2) + combination + combination:dr
      rnam <- c("Block", "GCA", "SCA", "REC", "Residuals")
    } else if(is.null(Block) & withRep == T){
      form <- ~ overlay(Par1, Par2) + combination + combination:dr
      rnam <- c("GCA", "SCA", "REC", "Residuals")
    } else if(is.null(Block) & withRep == F){
      form <- ~ overlay(Par1, Par2) + combination
      rnam <- c("GCA", "SCA", "REC")
    }

  # GRIFFING2 and 4 ####
  } else if(fct == "GRIFFING2"){
    formFix <- Y ~ 1
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

    # GRIFFING2 and 4 ####
  } else if(fct == "GRIFFING4"){
    formFix <- Y ~ 1
    if(!is.null(Block)){
      form <- ~ Block + overlay(Par1, Par2) + combination
      rnam <- c("Block", "GCA", "SCA", "Residuals")
    } else if(is.null(Block) & withRep == T){
      form <- ~ overlay(Par1, Par2) + combination
      rnam <- c("GCA", "SCA", "Residuals")
    } else if(is.null(Block) & withRep == F){
      form <- ~ overlay(Par1, Par2)
      rnam <- c("GCA", "SCA")
    }

  # GE2r ############
  } else if(fct == "GE2r"){
    formFix <- Y ~ crosses
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

  # GE2 ######
  } else if(fct == "GE2"){
    formFix <- Y ~ crosses
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

  # GE3r ####
  } else if(fct == "GE3r"){
    formFix <- Y ~ crosses
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

  # GE3 ####
  } else if(fct == "GE3"){
    formFix <- Y ~ crosses
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
  }
  }
 modh <- sommer::mmer(formFix,
                       random = form,
             verbose = F,
             data = df)

 vc <- summary(modh)$varcomp[,-c(3:4)]
 beta <- summary(modh)$beta
 row.names(vc) <- rnam
 returnList <- list(mod = modh, beta = beta, varcomp = vc)
 return(returnList)
}

