# This functions create model matrices for diallel models
# Date of last edit: 19/07/2021
# Removing limitation for less than 10 parentals
model.matrix.diallel <- function(object, ...){
  return(object$modMatrix)
}

model.matrixDiallel <- function(formula, Block = NULL, Env = NULL,
                                 fct = NULL, data = NULL,
                                type = "nested"){

  if(is.null(fct)){
    # This is a formula based output
    # This is used by one who knows what he is doing, and makes
    # direct use of the functions GCA, SCA and so on
    X <- model.matrix(formula, data)
    # X <- X1[,-1]
    # attr(X, "assign") <- attr(X1, "assign")[-1] - 1
    # attr(X, "assign") <- attr(X1, "assign")[-1] - 1
    } else {
    # fct based output.
    # It creates the incidence matrices for specific models
    mf <- match.call(expand.dots = FALSE) # Riprende la chiamata, con i nomi
    m <- match(c("formula", "Block", "Env", "data"), names(mf), 0L) # Trova nella chiamata la formula. m rappresenta la posizione della formula nella chiamata
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
      # Experiments in one environment
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

        # Building incidence matrix (0:5)
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
          namEffs <- c(namEffs, "GCA", "tSCA", "Reciprocals",
                                     "Residuals")
          attr(X, "namEff") <- namEffs

        } else if(fct == "GRIFFING2"){
          # GRIFFING 2 - No reciprocals #######################
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
          namEffs <- c(namEffs, "GCA", "tSCA",
                                     "Residuals")
          attr(X, "namEff") <- namEffs
        } else if(fct == "GRIFFING3"){
          # GRIFFING 3 - Reciprocal effects, no selfs ###################
          Z <- GCA(P1, P2)
          # SCA <- SCA.G3(P1, P2)
          SCA <- SCA(P1, P2) # Corrected on 1/7/21
          # rec <- REC.G3(P1, P2)
          rec <- REC(P1, P2) # Corrected on 2/3/21
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

        } else if(fct == "GRIFFING4"){
          # GRIFFING 4 - No reciprocals, no selfs #########
          # Z <- GCA(P1, P2) # Original
          Z <- GCAmis(P1, P2) # Edited on 20/03/23
          # SCA <- SCA.G3(P1, P2) # Original
          # SCA <- SCA(P1, P2) # Edited on 1/7/21
          SCA <- SCAmis(P1, P2) # Edited on 20/03/23
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
          crM <- H.BAR(P1, P2)
          H <- SP(P1, P2)
          Z <- GCAC(P1, P2)
          SCA <- SCA(P1, P2)

          # Building incidence matrix (0:5)
          X <- cbind(X, crM, H, Z, SCA)
          groups <- c(groups, seq(nGroups+1, nGroups+4, 1))
          reps <- c(reps,  1 ,length(H[1,])
                             ,length(Z[1, ])
                             ,length(SCA[1,]))
          levs <- rep(groups, reps)
          attr(X, "assign") <- levs
          namEffs <- c(namEffs, "h.bar",
                                 "Selfed parents",
                       "gcac", "SCA", "Residuals")
          attr(X, "namEff") <- namEffs

        } else if(fct == "GE3r"){
          # GE3r - Con reciproci ###############################################
          # 23/03/2020

          # Matrix for crosses
          #crM <- matrix(crosses, n, 1)
          crM <- H.BAR(P1, P2)
          # colnames(crM) <- "h.bar"

          H <- SP(P1, P2)
          # nams <- paste("sp_", levels(P1)[1:(length(levels(P1))-1)], sep="")
          # colnames(H) <- c(nams)

          Z <- GCAC(P1, P2)
          # nams <- paste("gcac_", levels(P1)[1:(length(levels(P1))-1)], sep="")
          # colnames(Z) <- c(nams)

          # Matrix for sca
          SCA <- SCA(P1, P2)
          # colnames(SCA) <- paste("sca_", colnames(SCA), sep = "")

          rec <- REC(P1, P2)
          #  colnames(rec) <- paste("rec_", colnames(rec), sep = "")

          # Building incidence matrix (0:6)
          X <- cbind(X, crM, H, Z, SCA, rec)
          groups <- c(groups, seq(nGroups+1, nGroups+5, 1))
          reps <- c(reps,  1 ,length(H[1,])
                             ,length(Z[1, ])
                             ,length(SCA[1,]),
                    length(rec[1,]))
          levs <- rep(groups, reps)
          attr(X, "assign") <- levs
          namEffs <- c(namEffs, "h.bar",
                                 "Selfed parents", "gcac",
                                 "SCA", "Reciprocals",
                                     "Residuals")
          attr(X, "namEff") <- namEffs

        }else{
             stop("Model not yet implemented")
        }
      } else {
        # GE Data ####################################
        # Creazione matrice incidenza
        # effetti genetici can be crossed or nested
        # Nested is the rule for fitting, but crossed are
        # necessary for ANOVA
        X <- model.matrixDiallel.MET(Par1, Par2, Block, Env ,
                                    fct, type = type)

      }
    }
  return(X)
}

model.matrixDiallel.MET <- function(Par1, Par2, Block, Env ,
                                    fct, type = "crossed"){

  if(type != "crossed" & type != "nested")
    stop("'Incidence matrices 'Type' can only be 'nested' ore 'crossed'")

  if(type == "crossed"){
      # Effetti genetici crossed
      n <- length(Par1)
      P1 <- factor(as.character(Par1))
      P2 <- factor(as.character(Par2))
      nGroups <- 0
      groups <- c(nGroups)
      Env <- factor(Env)
      contrasts(Env) <- c("contr.sum")
      EnvMat <- model.matrix(~ Env)

      if(!is.null(Block)) {
        Block <- factor(Block)
        contrasts(Block) <- c("contr.sum")
        X <- model.matrix(~ Env/Block)
        asgn1 <- attr(X, "assign")
        nGroups <- 2
        groups <- 0:2
        namEffs <- c("Intercept", "Env", "Env/Block")
      } else {
        X <- EnvMat
        asgn1 <- attr(X, "assign")
        nGroups <- 1
        groups <- 0:1
        namEffs <- c("Intercept", "Env")
      }
      EnvMat <- as.matrix(EnvMat[,-1])

      if(fct == "HAYMAN1"){
        # HAYMAN1 - con tSCA #############################
        # 6/04/2020
        Z <- GCA(P1, P2)
        SCA <- tSCA(P1, P2)
        RGCA <- RGCA(P1, P2)
        rec <- RSCA(P1, P2)
        Zenv <- int.matrix(Z, EnvMat)
        SCAenv <- int.matrix(SCA, EnvMat)
        RGCAenv <- int.matrix(RGCA, EnvMat)
        RSCAenv <- int.matrix(rec, EnvMat)

        X <- cbind(X, Z, SCA, RGCA, rec, Zenv, SCAenv, RGCAenv, RSCAenv)
        groups <- c(seq(nGroups+1, nGroups+8, 1))
        reps <- c(length(Z[1,]), length(SCA[1,])
                                ,length(RGCA[1,])
                                ,length(rec[1,])
                                ,length(Zenv[1,])
                                ,length(SCAenv[1,])
                  ,length(RGCAenv[1,])
                  ,length(RSCAenv[1,]))
        levs <- rep(groups, reps)
        attr(X, "assign") <- c(asgn1, levs)
        namEffs <- c(namEffs, "GCA","tSCA", "RGCA", "RSCA",
                              "GCA:Env", "tSCA:Env", "RGCA:Env", "RSCA:Env",
                                   "Residuals")
        attr(X, "namEff") <- namEffs

      } else if(fct == "HAYMAN2"){
        # HAYMAN2 - SCA decomposta ###########
        # 23/03/2020
       crM <- MDD(P1, P2)
       Z <- GCA(P1, P2)
       H <- DD(P1, P2)
       SCA <- SCA(P1, P2)
       RGCA <- RGCA(P1, P2)
       rec <- RSCA(P1, P2)
       crMenv <- int.matrix(crM, EnvMat)
       Zenv <- int.matrix(Z, EnvMat)
       Henv <- int.matrix(H, EnvMat)
       SCAenv <- int.matrix(SCA, EnvMat)
       RGCAenv <- int.matrix(RGCA, EnvMat)
       RSCAenv <- int.matrix(rec, EnvMat)

       # Building incidence matrix (0:7)
       X <- cbind(X, crM, Z, H, SCA, RGCA, rec, crMenv, Zenv, Henv,
                  SCAenv, RGCAenv, RSCAenv)
       groups <- c(seq(nGroups+1, nGroups+12, 1))
       reps <- c(length(crM[1,]), length(Z[1,]),
                 length(H[1,]), length(SCA[1,]),
                 length(RGCA[1,]), length(rec[1,]),
                 length(crMenv[1,]), length(Zenv[1,]),
                 length(Henv[1,]), length(SCAenv[1,]),
                 length(RGCAenv[1,]),length(RSCAenv[1,]))
       levs <- rep(groups, reps)
       attr(X, "assign") <- c(asgn1, levs)
       namEffs <- c(namEffs, "MDD", "GCA", "DD", "SCA",
                     "RGCA", "RSCA", "MDD:Env", "GCA:Env", "DD:Env",
                     "SCA:Env", "RGCA:Env", "RSCA:Env", "Residuals")
       attr(X, "namEff") <- namEffs

      } else if(fct == "GRIFFING1"){
        # GRIFFING 1 - Reciprocal effects #######################
        Z <- GCA(P1, P2)
        SCA <- tSCA(P1, P2)
        rec <- REC(P1, P2)
        Zenv <- int.matrix(Z, EnvMat)
        SCAenv <- int.matrix(SCA, EnvMat)
        recEnv <- int.matrix(rec, EnvMat)

        X <- cbind(X, Z, SCA, rec, Zenv, SCAenv, recEnv)
        groups <- c(seq(nGroups+1, nGroups+6, 1))
        reps <- c(length(Z[1,]), length(SCA[1,]),
                  length(rec[1,]), length(Zenv[1,]),
                  length(SCAenv[1,]), length(recEnv[1,]))
        levs <- rep(groups, reps)
        attr(X, "assign") <- c(asgn1, levs)
        namEffs <- c(namEffs, "GCA", "tSCA", "Reciprocals",
                     "GCA:Env", "tSCA:Env", "Reciprocals:Env",
                                   "Residuals")
        attr(X, "namEff") <- namEffs

      } else if(fct == "GRIFFING2"){
        # GRIFFING 2 - No reciprocals #######################
        # 23/03/2020
        Z <- GCA(P1, P2)
        SCA <- tSCA(P1, P2)
        Zenv <- int.matrix(Z, EnvMat)
        SCAenv <- int.matrix(SCA, EnvMat)
        X <- cbind(X, Z, SCA, Zenv, SCAenv)
        groups <- c(seq(nGroups+1, nGroups+4, 1))
        reps <- c(length(Z[1,]), length(SCA[1,]),
                  length(Zenv[1,]), length(SCAenv[1,]))
        levs <- rep(groups, reps)
        attr(X, "assign") <- c(asgn1, levs)
        namEffs <- c(namEffs, "GCA", "tSCA", "GCA:Env", "tSCA:Env",
                                   "Residuals")
        attr(X, "namEff") <- namEffs
      } else if(fct == "GRIFFING3"){
        # GRIFFING 3 - Reciprocal effects, no selfs ###################
        Z <- GCA(P1, P2)
        SCA <- SCA(P1, P2) # Corrected on 1/7/21
        rec <- REC(P1, P2) # Corrected on 2/3/21
        Zenv <- int.matrix(Z, EnvMat)
        SCAenv <- int.matrix(SCA, EnvMat)
        recEnv <- int.matrix(rec, EnvMat)
        X <- cbind(X, Z, SCA, rec, Zenv, SCAenv, recEnv)
        groups <- c(seq(nGroups+1, nGroups+6, 1))
        reps <- c(length(Z[1,]), length(SCA[1,]), length(rec[1,]),
                  length(Zenv[1,]), length(SCAenv[1,]), length(recEnv[1,]))
        levs <- rep(groups, reps)
        attr(X, "assign") <- c(asgn1, levs)
        namEffs <- c(namEffs, "GCA", "SCA", "Reciprocals",
                     "GCA:Env", "SCA:Env", "Reciprocals:Env",
                     "Residuals")
        attr(X, "namEff") <- namEffs

      } else if(fct == "GRIFFING4"){
        # GRIFFING 4 - No reciprocals, no selfs #########
        Z <- GCAmis(P1, P2) # Edited on 20/03/23
        SCA <- SCAmis(P1, P2) # Edited on 20/03/23
        Zenv <- int.matrix(Z, EnvMat)
        SCAenv <- int.matrix(SCA, EnvMat)
        X <- cbind(X, Z, SCA, Zenv, SCAenv)
        groups <- c(seq(nGroups+1, nGroups+4, 1))
        reps <- c(length(Z[1,]), length(SCA[1,]),
                  length(Zenv[1,]), length(SCAenv[1,]))
        levs <- rep(groups, reps)
        attr(X, "assign") <- c(asgn1, levs)
        namEffs <- c(namEffs, "GCA", "SCA", "GCA:Env", "SCA:Env",
                     "Residuals")
        attr(X, "namEff") <- namEffs

      } else if(fct == "GE2"){
        # GE2 - Senza reciproci ####################
        # 23/03/2020
        crM <- H.BAR(P1, P2)
        Z <- VEi(P1, P2)
        H <- Hi(P1, P2)
        SCA <- SCA(P1, P2)
        crMenv <- int.matrix(crM, EnvMat)
        Zenv <- int.matrix(Z, EnvMat)
        Henv <- int.matrix(H, EnvMat)
        SCAenv <- int.matrix(SCA, EnvMat)
        X <- cbind(X, crM, Z, H, SCA, crMenv, Zenv, Henv, SCAenv)
        groups <- c(seq(nGroups+1, nGroups+8, 1))
        reps <- c(1, length(Z[1,]) ,length(H[1, ]),length(SCA[1,]),
                  length(crMenv[1,]) ,length(Zenv[1,]),
                  length(Henv[1, ]),length(SCAenv[1,]))
        levs <- rep(groups, reps)
        attr(X, "assign") <- c(asgn1, levs)
        namEffs <- c(namEffs, "h.bar", "Variety", "h.i", "SCA",
                     "h.bar:Env", "Variety:Env", "h.i:Env", "SCA:Env",
                                   "Residuals")
        attr(X, "namEff") <- namEffs

      } else if(fct == "GE2r"){
        # GE2r - With reciprocals ####################
        # 23/03/2020
        crM <- H.BAR(P1, P2)
        Z <- VEi(P1, P2)
        H <- Hi(P1, P2)
        SCA <- SCA(P1, P2)
        rec <- REC(P1, P2)
        crMenv <- int.matrix(crM, EnvMat)
        Zenv <- int.matrix(Z, EnvMat)
        Henv <- int.matrix(H, EnvMat)
        SCAenv <- int.matrix(SCA, EnvMat)
        recEnv <- int.matrix(rec, EnvMat)
        X <- cbind(X, crM, Z, H, SCA, rec, crMenv, Zenv, Henv, SCAenv,
                   recEnv)
        groups <- c(seq(nGroups+1, nGroups+10, 1))
        reps <- c(1, length(Z[1,]) ,length(H[1, ]),length(SCA[1,]),
                  length(rec[1,]),
                  length(crMenv[1,]) ,length(Zenv[1,]),
                  length(Henv[1, ]),length(SCAenv[1,]),
                  length(recEnv[1,]))
        levs <- rep(groups, reps)
        attr(X, "assign") <- c(asgn1, levs)
        namEffs <- c(namEffs, "h.bar", "Variety", "h.i", "SCA", "Reciprocal",
                     "h.bar:Env", "Variety:Env", "h.i:Env", "SCA:Env",
                     "Reciprocal:Env", "Residuals")
        attr(X, "namEff") <- namEffs

      } else if(fct == "GE3"){
        # GE3 - Senza reciproci ##############################
        # 23/03/2020
        crM <- H.BAR(P1, P2)
        H <- SP(P1, P2)
        Z <- GCAC(P1, P2)
        SCA <- SCA(P1, P2)
        crMenv <- int.matrix(crM, EnvMat)
        Henv <- int.matrix(H, EnvMat)
        Zenv <- int.matrix(Z, EnvMat)
        SCAenv <- int.matrix(SCA, EnvMat)
        X <- cbind(X, crM, H, Z, SCA, crMenv, Henv, Zenv, SCAenv)
        groups <- c(seq(nGroups+1, nGroups+8, 1))
        reps <- c(1 ,length(H[1,]), length(Z[1, ]), length(SCA[1,]),
                  length(crMenv[1,]), length(Henv[1,]),
                  length(Zenv[1, ]), length(SCAenv[1,]))
        levs <- rep(groups, reps)
        attr(X, "assign") <- c(asgn1, levs)
        namEffs <- c(namEffs, "h.bar", "Selfed parents", "gcac", "SCA",
                     "h.bar:Env", "Selfed parents:Env", "gcac:Env",
                     "SCA:Env", "Residuals")
        attr(X, "namEff") <- namEffs

      } else if(fct == "GE3r"){
      # GE3r - Con reciproci ###############################################
      # 23/03/2020
      crM <- H.BAR(P1, P2)
      H <- SP(P1, P2)
      Z <- GCAC(P1, P2)
      SCA <- SCA(P1, P2)
      rec <- REC(P1, P2)
      crMenv <- int.matrix(crM, EnvMat)
      Henv <- int.matrix(H, EnvMat)
      Zenv <- int.matrix(Z, EnvMat)
      SCAenv <- int.matrix(SCA, EnvMat)
      recEnv <- int.matrix(rec, EnvMat)
      X <- cbind(X, crM, H, Z, SCA, rec,
                 crMenv, Henv, Zenv, SCAenv, recEnv)
      groups <- c(seq(nGroups+1, nGroups+10, 1))
      reps <- c(1, length(H[1,]), length(Z[1, ]), length(SCA[1,]),
                length(rec[1,]),
                length(crMenv[1,]), length(Henv[1,]),
                length(Zenv[1, ]), length(SCAenv[1,]),
                length(recEnv[1,]))
      levs <- rep(groups, reps)
      attr(X, "assign") <- c(asgn1, levs)
      namEffs <- c(namEffs, "h.bar", "Selfed parents", "gcac", "SCA",
                   "Reciprocal",
                   "h.bar:Env", "Selfed parents:Env", "gcac:Env",
                   "SCA:Env", "Reciprocal:Env", "Residuals")
      attr(X, "namEff") <- namEffs

      }else{
        stop("Model not yet implemented")
      }
      } else {
      # Effetti genetici nested (normal)
      # Par1 <- P1
      # Par2 <- P2
      Env <- factor(Env)
      if(!is.null(Block)){
        Block <- factor(Block)
        datasetS <- data.frame(Id = 1:length(Par1), Env, Block, Par1, Par2)
        datasetS <- datasetS[order(datasetS$Env, datasetS$Par1,
                     datasetS$Par2, datasetS$Block), ]
        matsOr <- plyr::dlply(datasetS, c("Env"), function(df){
                model.matrixDiallel(~ df$Par1 + df$Par2, df$Block,
                    fct = fct)})
      } else {
        datasetS <- data.frame(Id = 1:length(Par1), Env, Par1, Par2)
        datasetS <- datasetS[order(datasetS$Env, datasetS$Par1,
                     datasetS$Par2), ]
        matsOr <- plyr::dlply(datasetS, c("Env"), function(df){
                model.matrixDiallel(~ df$Par1 + df$Par2, fct = fct)})
      }

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
      # asgn2 <- c(unlist(asgnList[1]),unlist(lapply(asgnList[-1], function(x) unlist(x))) )
      # asgn2 <- as.numeric(c(rep(0, length(levels(datasetS$Env))), asgn2))
      # attr(X, "assign2") <- asgn2
      }
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


GCA <- function(P1, P2, type = "fix", data = NULL){
  # This is modified to work with mating design 4
  if(!is.null(data)){
    P1Name <- deparse(substitute(P1))
    P2Name <- deparse(substitute(P2))
    P1 <- data[[P1Name]]
    P2 <- data[[P2Name]]
  }
  if(type == "random"){
      Z <- sommer::overlay(P1, P2, sparse = F)
  } else {
  P1 <- factor(as.character(P1))
  P2 <- factor(as.character(P2))
  levs <- c(levels(P1), levels(P2))
  levs <- levels(factor(levs, levels = unique(levs)))
  Z1n <- factor(P1, levels = levs, ordered = T)
  Z2n <- factor(P2, levels = levs)
  contrasts(Z1n) <- c("contr.sum")
  contrasts(Z2n) <- c("contr.sum")
  Z1 <- model.matrix(~ Z1n)
  Z2 <- model.matrix(~ Z2n)
  Z <- (Z1 + Z2)
  Z <- Z[,-1]
  nams <- paste("g_", levs[1:(length(levs)-1)], sep="")
  colnames(Z) <- c(nams)
  }
  Z
}

VEi <- function(P1, P2, type = "fix", data = NULL){
    if(!is.null(data)){
    P1Name <- deparse(substitute(P1))
    P2Name <- deparse(substitute(P2))
    P1 <- data[[P1Name]]
    P2 <- data[[P2Name]]
    }
  if(type == "random"){
    Z <- GCA(P1, P2, type = "random")
    # Z <- Z/2
    return(Z)
  } else {

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
  Z }
}

SP <- function(P1, P2, type = "fix", data = NULL){
    if(!is.null(data)){
    P1Name <- deparse(substitute(P1))
    P2Name <- deparse(substitute(P2))
    P1 <- data[[P1Name]]
    P2 <- data[[P2Name]]
    }
  if(type == "random"){
    crosses <- ifelse(P1 == P2, 0, 1)
    Z <- GCA(P1, P2, type = "random") * crosses
    #colnames(Z) <- sub("combination", "", colnames(Z))
    Z <- Z[, apply(Z, 2, function(x) !all(x==0))]
    return(Z)
  } else {
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
  Z}
}

RGCA <- function(P1, P2, type = "fix", data = NULL){
  if(!is.null(data)){
    P1Name <- deparse(substitute(P1))
    P2Name <- deparse(substitute(P2))
    P1 <- data[[P1Name]]
    P2 <- data[[P2Name]]
  }
    if(type == "random"){
    dr <- ifelse(as.character(P1) < as.character(P2), -1,
              ifelse(as.character(P1) == as.character(P2), 0, 1))
    Z <- GCA(P1, P2, type = "random") * dr
    # colnames(Z) <- sub("combination", "", colnames(Z))
    # Z <- Z[, apply(Z, 2, function(x) !all(x==0))]
    return(Z)
  } else {
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
  RGCA }
}

tSCA <- function(P1, P2, type = "fix", data = NULL){
    if(!is.null(data)){
    P1Name <- deparse(substitute(P1))
    P2Name <- deparse(substitute(P2))
    P1 <- data[[P1Name]]
    P2 <- data[[P2Name]]
    }
  if(type == "random"){
    combination <- factor( ifelse(as.character(P1) <= as.character(P2),
                                 paste(P1, P2, sep =""),
                                 paste(P2, P1, sep ="")) )
    Z <- model.matrix(~ combination - 1)
    colnames(Z) <- sub("combination", "", colnames(Z))
    return(Z)
  } else {
  # Matrix tSCA: final version: 30/6/2020
  P1 <- factor(as.character(P1))
  P2 <- factor(as.character(P2))
  P1c <- as.character(P1); P2c <- as.character(P2)

  # Combination
  tmp <- ifelse(P1c < P2c, paste(P1c, P2c, sep =":"),
         paste(P2c, P1c, sep = ":"))
  combination <- factor(tmp) #, levels = unique(tmp))
  combLev <- NA
  mating <- P1:P2
  p <- length(levels(factor(c(levels(P1), levels(P2)) )))
  n <- length(combination)

  # Step 1. gets the parameters to be estimated, removing
  # the unnecessary combinations
  tmp <- sapply(by(P2, P1, function(x) levels(x)), function(x) max(as.character(x)))
  tmp <- ifelse(names(tmp) < tmp, paste(names(tmp), tmp, sep = ":"), paste(tmp, names(tmp), sep = ":"))
  last <- levels(factor(tmp, levels = unique(tmp)))

  levs <- levels(combination)
  idx <- c() # Identifica la posizione degli ultimi livelli
   for(i in 1:length(last)){
       #i <- 2
       y <- which(levs == last[i])
       idx[i] <- y
  }
  levs <- as.character(levs[-idx]) # Esclude gli ultimi livelli per ogni Par1
  SCA <- matrix(0, nrow = n, ncol = length(levs))
  colnames(SCA) <- levs

  # Step 2. Insert 1s for all levels, but the last one
  for(i in 1:length(levs)){
    # i <- 1
    cond <- (combination == colnames(SCA)[i])*1
    SCA[, i] <- cond
   }

  # Step 3. Insert the -1s for the last level. The last level of
  # Par2, within each level of Par1. The last level of Par1
  # requires another step
  for(i in 1:(length(last) - 1)){
     arrival <- last[i]
    tmp <- strsplit(arrival, ":")[[1]]
    revArrival <- paste(rev(tmp), collapse = ":")
    lastEl <- c(arrival, revArrival)
    sel <- sapply(strsplit(colnames(SCA), ":"), function(i) any(i == tmp[1]))
    idx <- sapply(1:length(combination), function(i) any(lastEl == mating[i]))
    SCA[idx, sel] <- -1
  }
  SCA
  # Scrive il self dell'ultimo livello
  SCA[combination == last[p], ] <- 2
  for(i in 1:p) {
    SCA[combination == last[p], colnames(SCA) == paste(levels(P1)[i], levels(P2)[i], sep = ":")] <- 1
    }
  colnames(SCA) <- paste("ts_", colnames(SCA), sep = "")
  row.names(SCA) <- c(1:length(SCA[,1]))
  return(SCA) }
}

SCA <- function(P1, P2, type = "fix", data = NULL){
  # Edited on 10/3/2023
  if(!is.null(data)){
    P1Name <- deparse(substitute(P1))
    P2Name <- deparse(substitute(P2))
    P1 <- data[[P1Name]]
    P2 <- data[[P2Name]]
  }
  if(type == "random"){
    crosses <- ifelse(as.character(P1) == as.character(P2), 0, 1)
    combination <- factor( ifelse(as.character(P1) <= as.character(P2),
                                 paste(P1, P2, sep =""),
                                 paste(P2, P1, sep ="")) )
    Z <- model.matrix(~ combination - 1) * crosses
    colnames(Z) <- sub("combination", "", colnames(Z))
    Z <- Z[, apply(Z, 2, function(x) !all(x==0))]
    return(Z)
  } else {

  # Matrix for SCA in heterosis model Hayman2
  # It is also used where the selfs are not includede
  P1 <- factor(as.character(P1)) #, levels = unique(P1))
  P2 <- factor(as.character(P2)) #, levels = unique(P1)) # Livelli uguali?
  P1c <- as.character(P1); P2c <- as.character(P2)

  # create the combination levels
  tmp <- ifelse(P1c < P2c, paste(P1c, P2c, sep =":"),
         paste(P2c, P1c, sep = ":"))
  combination <- factor(tmp) #, levels = unique(tmp))
  combLev <- NA

  # Create the matings (considers the reciprocals)
  mating <- P1:P2
  p <- length(levels(factor(c(levels(P1), levels(P2)) )))
  n <- length(combination)

  # See whether selfs are included and find the last level for each
  # combination
  selflist <- levels(factor(combination[P1c == P2c]))
  levs <- sort(unique(c(levels(P1), levels(P2))))
  tmp <- paste(levs, max(levs), sep = ":")
  parLevs <- levs
  last <- levels(factor(tmp, levels = unique(tmp)))
  levs <- levels(combination)

  # Cerca la posizione dell'ultimo in ogni gruppo
  idx1 <- c()
  for(i in 1:length(last)){
         y <- which(levs == last[i])
         # print(i); print(length(y))
         if(length(y) != 0) idx1[i] <- y else next
  }

  # cerca la posizione dei selfs, se esistenti
  idx3 <- c()
  for(i in 1:length(selflist)){
          #i <- 2
          y <- which(levs == selflist[i])
          if(length(y) > 0) idx3[i] <- y
  }

  # Definisce i parametr da stimare. Deve rimuovere gli ultimi livelli
  # e i selfs, se esistenti. Deeve anche rimuover l'elemento con parentali
  # (p - 2) e (p - 1)
  idx <- c(idx1, idx3)
  rimossi <- levs[idx]
  levs <- as.character(levs[-idx])
  rimossi <- c(rimossi, levs[length(levs)])
  levs <- levs[-length(levs)] # rimuove l'ultimo

  # Step 1. Create an empty matrix
  SCA <- matrix(0, nrow = n, ncol = length(levs))
  colnames(SCA) <- paste(levs)

  # Step 2. Insert 1s for all the levels, which correspond
  # to estimands
  for(i in 1:length(levs)){
          cond <- (combination == colnames(SCA)[i]) * 1
          SCA[, i] <- cond
  }

  # Step 3. Insert the -1s for the last level of
  # Par2, within each level of Par1. This is only done for Par 1
  # going from 1 to (p - 3), per gli altri ci vuole un altro step
  for(i in 1:(length(last) - 3)){
    arrival <- last[i]
    tmp <- strsplit(arrival, ":")[[1]]
    revArrival <- paste(rev(tmp), collapse = ":")
    lastEl <- c(arrival, revArrival)
    sel <- sapply(strsplit(colnames(SCA), ":"), function(i) any(i == tmp[1]))
    idx <- sapply(1:length(combination), function(i) any(lastEl == mating[i]))
    SCA[idx, sel] <- -1
    # SCA
    }

  # Mancano le combinazioni degli ultimi 3 ibridi
    revParents <- function(x){
        tmp <- strsplit(x, ":")[[1]]
        paste(rev(tmp), collapse = ":")}

    tmp <- seq(p - 2, p, 1)
    tmp <- as.data.frame(combn(tmp, 2))
    # Edited on 10/3/2023, to avoid an error for mating schemes without selfs
    # tmp <- apply(tmp, 2, function(x) paste(levels(P1)[x[1]], levels(P2)[x[2]], sep = ":"))
    tmp <- apply(tmp, 2, function(x) paste(parLevs[x[1]], parLevs[x[2]], sep = ":"))
    tmp2 <- mapply(revParents, tmp)

     # Si occupa del livello P1 = p -2 e P2 = P-1 che deve essere pari
     # all'opposto della somma di tutti gli altri parametri
     SCA[combination == tmp[1], ] <- -1
     SCA[combination == tmp2[1], ] <- -1

     # Si occupa del parametro per (p-2, p), che si ottiene come somma
     # di tutti i parametri i cui parentali non sono uguali a (p - 2)
     tmp3 <- strsplit(tmp[2], ":")[[1]]
     sel <- sapply(strsplit(colnames(SCA), ":"), function(i) !any(i == tmp3[1]))
     SCA[combination == tmp[2], sel] <- 1
     SCA[combination == tmp2[2], sel] <- 1

     # Si occupa del parametro per (p-1, p), che si ottiene come somma
     # di tutti i parametri i cui parentali non sono uguali a (p - 1)
     tmp3 <- strsplit(tmp[3], ":")[[1]]
     sel <- sapply(strsplit(colnames(SCA), ":"), function(i) !any(i == tmp3[1]))
     SCA[combination == tmp[3], sel] <- 1
     SCA[combination == tmp2[3], sel] <- 1
     colnames(SCA) <- paste("s_", colnames(SCA), sep = "")
     SCA
     }
}


# SCA <- function(P1, P2, type = "fix", data = NULL){
#     if(!is.null(data)){
#     P1Name <- deparse(substitute(P1))
#     P2Name <- deparse(substitute(P2))
#     P1 <- data[[P1Name]]
#     P2 <- data[[P2Name]]
#     }
#     if(type == "random"){
#     crosses <- ifelse(P1 == P2, 0, 1)
#     combination <- factor( ifelse(as.character(P1) <= as.character(P2),
#                                  paste(P1, P2, sep =""),
#                                  paste(P2, P1, sep ="")) )
#     Z <- model.matrix(~ combination - 1) * crosses
#     colnames(Z) <- sub("combination", "", colnames(Z))
#     Z <- Z[, apply(Z, 2, function(x) !all(x==0))]
#     return(Z)
#   } else {
#
#   # Matrix for SCA in heterosis model Hayman2
#   P1 <- factor(as.character(P1)) #, levels = unique(P1))
#   P2 <- factor(as.character(P2)) #, levels = unique(P1)) # Livelli uguali?
#   P1c <- as.character(P1); P2c <- as.character(P2)
#
#   # combination
#   tmp <- ifelse(P1c < P2c, paste(P1c, P2c, sep =":"),
#          paste(P2c, P1c, sep = ":"))
#   combination <- factor(tmp) #, levels = unique(tmp))
#
#   combLev <- NA
#
#   mating <- P1:P2
#
#   p <- length(levels(factor(c(levels(P1), levels(P2)) )))
#   n <- length(combination)
#
#   selflist <- levels(factor(combination[P1c == P2c]))
#
#   # See whether selfs are included and find the last level for each P1
#   levs <- sort(unique(c(levels(P1), levels(P2))))
#   # tmp <- sapply(by(P1, P2, function(x) levels(x)), function(x) max(as.character(x)))
#   # tmp <- ifelse(names(tmp) < tmp, paste(names(tmp), tmp, sep = ":"), paste(tmp, names(tmp), sep = ":"))
#   tmp <- paste(levs, max(levs), sep = ":")
#   last <- levels(factor(tmp, levels = unique(tmp)))
#
#   # tmp <- sapply(by(P2, P1, function(x) levels(x)), function(x) sort(as.character(x))[(length(x) - 1)])
#   # tmp <- ifelse(names(tmp) < tmp, paste(names(tmp), tmp, sep = ":"), paste(tmp, names(tmp), sep = ":"))
#   # lastButOne <- levels(factor(tmp, levels = unique(tmp)))
#   #
#   # tmp <- sapply(by(P2, P1, function(x) levels(x)), function(x) sort(as.character(x))[(length(x) - 2)])
#   # tmp <- ifelse(names(tmp) < tmp, paste(names(tmp), tmp, sep = ":"), paste(tmp, names(tmp), sep = ":"))
#   # lastButTwo <- levels(factor(tmp, levels = unique(tmp)))
#
#   levs <- levels(combination)
#   idx1 <- c()
#   for(i in 1:length(last)){ #Indica la posizione dell'ultimo in ogni gruppo
#          y <- which(levs == last[i])
#          # print(i); print(length(y))
#          if(length(y) != 0) idx1[i] <- y else next
#   }
#   # idx2 <- c()
#   # for(i in 1:length(last)){ #Indica il penultimo
#   #         y <- which(levs == lastButOne[i])
#   #         if(length(y) > 0) idx2[i] <- y
#   # }
#   #
#   # idx2b <- c()
#   # for(i in 1:length(last)){ #Indica il terzultimo
#   #         y <- which(levs == lastButTwo[i])
#   #         if(length(y) > 0) idx2b[i] <- y
#   # }
#
#   idx3 <- c()
#   for(i in 1:length(selflist)){ #Indica i selfs
#           #i <- 2
#           y <- which(levs == selflist[i])
#           if(length(y) > 0) idx3[i] <- y
#   }
#   idx <- c(idx1, idx3)
#   rimossi <- levs[idx]
#   levs <- as.character(levs[-idx])
#   rimossi <- c(rimossi, levs[length(levs)])
#   levs <- levs[-length(levs)]
#   SCA <- matrix(0, nrow = n, ncol = length(levs))
#   colnames(SCA) <- paste(levs)
#
#   # Step 2. Insert 1s for all the levels, which are
#   # in the SCA matrix
#   for(i in 1:length(levs)){
#           cond <- (combination == colnames(SCA)[i]) * 1
#           SCA[, i] <- cond
#     }
#
#   # Step 3. Insert the -1s for the last level. The last level of
#   # Par2, within each level of Par1. The last level of Par1
#   # requires another step
#   for(i in 1:(length(last) - 3)){
#     # i <- 1
#     arrival <- last[i]
#     tmp <- strsplit(arrival, ":")[[1]]
#     revArrival <- paste(rev(tmp), collapse = ":")
#     lastEl <- c(arrival, revArrival)
#     sel <- sapply(strsplit(colnames(SCA), ":"), function(i) any(i == tmp[1]))
#     idx <- sapply(1:length(combination), function(i) any(lastEl == mating[i]))
#     SCA[idx, sel] <- -1
#     SCA
#     }
#
#   # Mancano le combinazioni degli ultimi 3 ibridi
#     revParents <- function(x){
#         tmp <- strsplit(x, ":")[[1]]
#         paste(rev(tmp), collapse = ":")}
#
#     tmp <- seq(p - 2, p, 1)
#     tmp <- as.data.frame(combn(tmp, 2))
#     tmp <- apply(tmp, 2, function(x) paste(levels(P1)[x[1]], levels(P2)[x[2]], sep = ":"))
#     tmp2 <- mapply(revParents, tmp)
#
#
#      SCA[combination == tmp[1], ] <- -1
#      SCA[combination == tmp2[1], ] <- -1
#
#      tmp3 <- strsplit(tmp[2], ":")[[1]]
#      sel <- sapply(strsplit(colnames(SCA), ":"), function(i) !any(i == tmp3[1]))
#      SCA[combination == tmp[2], sel] <- 1
#      SCA[combination == tmp2[2], sel] <- 1
#
#      tmp3 <- strsplit(tmp[3], ":")[[1]]
#      sel <- sapply(strsplit(colnames(SCA), ":"), function(i) !any(i == tmp3[1]))
#      SCA[combination == tmp[3], sel] <- 1
#      SCA[combination == tmp2[3], sel] <- 1
#      #colnames(SCA) <- paste("s_", colnames(SCA), sep = "")
#      colnames(SCA) <- paste("s_", colnames(SCA), sep = "")
#      SCA
#      }
# }

# SCA.old <- function(P1, P2, type = "fix", data = NULL){
#     if(!is.null(data)){
#     P1Name <- deparse(substitute(P1))
#     P2Name <- deparse(substitute(P2))
#     P1 <- data[[P1Name]]
#     P2 <- data[[P2Name]]
#     }
#     if(type == "random"){
#     crosses <- ifelse(P1 == P2, 0, 1)
#     combination <- factor( ifelse(as.character(P1) <= as.character(P2),
#                                  paste(P1, P2, sep =""),
#                                  paste(P2, P1, sep ="")) )
#     Z <- model.matrix(~ combination - 1) * crosses
#     colnames(Z) <- sub("combination", "", colnames(Z))
#     Z <- Z[, apply(Z, 2, function(x) !all(x==0))]
#     return(Z)
#   } else {
#
#   # Matrix for SCA in heterosis model Hayman2
#   # P1 <- df$Par1; P2 <- df$Par2
#   P1 <- factor(as.character(P1))
#   P2 <- factor(as.character(P2))
#   P1n <- as.numeric(P1); P2n <- as.numeric(P2)
#   P1c <- as.character(P1); P2c <- as.character(P2)
#   combination <- factor(apply(cbind(P1n*10 + P2n, P2n * 10 + P1n), 1, min))
#   combLev <- factor( ifelse(P1c < P2c, paste(P1c, P2c, sep = ":"), paste(P2c, P1c, sep = ":") ) )
#   mating <- factor(P1n*10 + P2n)
#   p <- length(levels(factor(c(levels(P1), levels(P2)) )))
#   n <- length(combination)
#   last <- seq(10, p*10, 10) + p
#   levs <- as.numeric(levels(combination))
#   idx1 <- c()
#   for(i in 1:length(last)){ #Indica l'ultimo
#          y <- which(levs == last[i])
#          idx1[i] <- y
#   }
#
#   idx2 <- c()
#   for(i in 1:length(last)){ #Indica il penultimo
#           y <- which(levs == (last[i] - 1))
#           if(length(y) > 0) idx2[i] <- y
#   }
#   selflist <- as.numeric(levels(factor(combination[P1n == P2n])))
#   idx3 <- c()
#   for(i in 1:length(selflist)){ #Indica il penultimo
#           #i <- 2
#           y <- which(levs == selflist[i])
#           if(length(y) > 0) idx3[i] <- y
#   }
#   idx <- c(idx1, idx3)
#   rimossi <- levs[idx]
#   levs <- as.character(levs[-idx])
#   rimossi <- c(rimossi, levs[length(levs)])
#   levs <- levs[-length(levs)]
#   SCA <- matrix(0, nrow = n, ncol = length(levs))
#   colnames(SCA) <- paste(levs)
#   #colnames(SCA)
#   colNamsOrd <- levels(combLev)[-idx][-length(levels(combLev)[-idx])]
#
#     # Step 2. Insert 1s for all the levels, which are
#     # in the SCA matrix
#     for(i in 1:length(levs)){
#           cond <- (combination == colnames(SCA)[i])*1
#           SCA[, i] <- cond
#     }
#     # Step 3. Insert the -1s for the last level. The last level of
#     # Par2, within each level of Par1. The last level of Par1
#     # requires another step
#     for(i in 1:(length(last) - 1)){
#         start <- ceiling(ifelse(i == 1, 1, last[i-1])/10) * 10 +1
#         arrival <- last[i]
#         sel <- seq(start, arrival, 1)
#         tmp <- as.character( sel[1:i] ) # Se necessario, inverte i reciproci
#         splits <- strsplit(tmp, "")
#         reversed <- lapply(splits, rev)
#         tmp <- as.character(lapply(reversed, paste, collapse = ""))
#         sel[1:i] <- as.numeric(tmp)
#         #sel
#         idx <- c() # Identifica la posizione di quelli da scrivere
#         for(j in 1:length(sel)){
#              #i <- 7
#              y <- which(levs == sel[j])
#              if(length(y) > 0) idx[j] <- y
#         }
#         idx
#         SCA[,idx]
#         idx1 <- last[i]
#         SCA[combination == idx1, idx] <- -1
#     }
#      # SCA
#      # Mancano le combinazioni degli ultimi 3 ibridi
#      # 67, 68, 76, 78, 86, 87
#      tmp <- seq(p-2, p, 1)
#      tmp <- as.data.frame(combn(tmp, 2))
#      tmp <- apply(tmp, 2, function(x) paste(x[1], x[2], sep = ""))
#      SCA[combination == tmp[1], ] <- -1
#      SCA[combination == tmp[2], ] <- SCA[combination == tmp[2], ] + 1
#      SCA[combination == tmp[3], ] <- SCA[combination == tmp[3], ] + 1
#      #colnames(SCA) <- paste("s_", colnames(SCA), sep = "")
#      colnames(SCA) <- paste("s_", colNamsOrd, sep = "")
#      SCA }
# }

SCA.G3 <- function(P1, P2, type = "fix", data = NULL){
  # tSCA effect in absence of selfed parents (only crosses)
  # and reciprocals
  # It is superseeded!!!!!!!!!
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
  last <- seq(10, (p-1)*10, 10) + p
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

SCA.GE <- function(P1, P2, type = "fix", data = NULL){
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
  colnames(SCA) <- paste(levs) #paste("s_", levs, sep = "")
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
  colnames(SCA) <- paste("s_", colnames(SCA), sep = "")
  SCA

}

RSCA <- function(P1, P2, type = "fix", data = NULL){
    if(!is.null(data)){
    P1Name <- deparse(substitute(P1))
    P2Name <- deparse(substitute(P2))
    P1 <- data[[P1Name]]
    P2 <- data[[P2Name]]
  }
    if(type == "random"){
    dr <- ifelse(as.character(P1) < as.character(P2), -1,
              ifelse(as.character(P1) == as.character(P2), 0, 1))

    Z <- tSCA(P1, P2, type = "random") * dr
    # colnames(Z) <- sub("combination", "", colnames(Z))
    Z <- Z[, apply(Z, 2, function(x) !all(x==0))]
    return(Z)
  } else {

  # Derive the dummies and other infos
  # P1 <- df$Par1; P2 <- df$Par2
  P1 <- factor(as.character(P1))
  P2 <- factor(as.character(P2))
  P1c <- as.character(P1); P2c <- as.character(P2)
  n <- length(P1)
  p <- length(levels(P1))

  # P1n <- as.numeric(P1); P2n <- as.numeric(P2) # Problematico
  # mate <- factor(P1n*10 + P2n)
  mate <- P1:P2
  # combination <- factor(apply(cbind(P1n*10 + P2n, P2n * 10 + P1n), 1, min))
  # combLev <- factor( ifelse(P1c < P2c, paste(P1c, P2c, sep = ":"), paste(P2c, P1c, sep = ":") ) )
  tmp <- ifelse(P1c < P2c, paste(P1c, P2c, sep =":"),
         paste(P2c, P1c, sep = ":"))
  combination <- factor(tmp) #, levels = unique(tmp))
  combLev <- NA

  # Empty matrix
  rec <- matrix(0, nrow = n, ncol = (p - 1)*(p - 2)/2 )
  cont <- 0
  nams <- c(); nams2 <- c()

  # Select the names of columns
  for(i in 1:(p-2)) { for(j in (i + 1):(p-1)){
     cont <- cont + 1
     nams[cont] <- paste(i, j, sep="")
     nams2[cont] <- paste(levels(P1)[i], levels(P2)[j], sep = ":")
     }}
  colnames(rec) <- nams2

  # Step 1. Add 1 for the crosses corresponding to column name
  for(i in 1:(p - 2)) { for(j in (i + 1):(p - 1)){
    #i <- 1; j <- 2
    cond <- paste(levels(P1)[i], levels(P2)[j], sep=":")
    rec[mate == cond, colnames(rec) == cond] <- 1
    rec[combination == cond & mate != cond, colnames(rec) == cond] <- -1
  }}

  # Step 2. Work on the last level
  leftr <- P1
  rightr <- P2
  leftc <- as.character(do.call(rbind, mapply(strsplit, nams2, split = ":"))[,1])
  rightc <- as.character(do.call(rbind, mapply(strsplit, nams2, split = ":"))[,2])

  for(i in 1:length(rec[1,])){
     rec[rightr == levels(P2)[p] & leftr == levels(P1)[i], leftc == levels(P1)[i]] <- -1
     rec[rightr == levels(P2)[p] & leftr == levels(P1)[i], rightc == levels(P1)[i]] <- 1
     rec[leftr == levels(P2)[p] & rightr == levels(P2)[i], leftc == levels(P1)[i]] <- 1
     rec[leftr == levels(P2)[p] & rightr == levels(P2)[i], rightc == levels(P1)[i]] <- -1
    }
  colnames(rec) <- paste("rs_", nams2, sep = "")
  rec
  }
  }

REC <- function(P1, P2, type = "fix", data = NULL){
    if(!is.null(data)){
    P1Name <- deparse(substitute(P1))
    P2Name <- deparse(substitute(P2))
    P1 <- data[[P1Name]]
    P2 <- data[[P2Name]]
    }
   if(type == "random"){
    combination <- factor( ifelse(as.character(P1) <= as.character(P2),
                                 paste(P1, P2, sep =""),
                                 paste(P2, P1, sep ="")) )
    dr <- ifelse(as.character(P1) < as.character(P2), -1,
              ifelse(as.character(P1) == as.character(P2), 0, 1))

    Z <- model.matrix(~ combination - 1) * dr
    colnames(Z) <- sub("combination", "", colnames(Z))
    Z <- Z[, apply(Z, 2, function(x) !all(x==0))]
    return(Z)
  } else {
  # P1 <- factor(as.character(P1))
  # P2 <- factor(as.character(P2))
  P1 <- factor(as.character(P1))
  P2 <- factor(as.character(P2))
  n <- length(P1)
  p <- length(levels(P1))
  # P1n <- as.numeric(P1); P2n <- as.numeric(P2)
  P1c <- as.character(P1); P2c <- as.character(P2)

  # mate <- factor(P1n*10 + P2n)
  # tmp <- as.numeric(paste(P1n, P2n, sep = ""))
  # mate <- factor(tmp, levels = unique(tmp))
  # as.character(P1:P2)

  # combination <- factor(apply(cbind(P1n*10 + P2n, P2n * 10 + P1n), 1, min))
  tmp <- ifelse(P1c < P2c, paste(P1c, P2c, sep =":"),
         paste(P2c, P1c, sep = ":"))
  combination <- factor(tmp) #, levels = unique(tmp))

  # dr <- ifelse(P1c == P2c, 0, ifelse(P1c > P2c, -1, 1))
  dr <- ifelse(P1c == P2c, 0, ifelse(P1c > P2c, -1, 1))
  # combLev <- factor( paste(P1c, P2c, sep = ":") )
  combLev <- P1:P2
  # tmp <- ifelse(P1n < P2n, paste(P1c, P2c, sep =":"),
  #        paste(P2c, P1c, sep =":"))
  # combLev <- factor(tmp, levels = unique(tmp))

  last <- c(); cont = 1
  for(i in 1:p){ for(j in 1:i) {
    last[cont] <- paste(i, j, sep=":")
    last[cont] <- paste(levels(P1)[i], levels(P2)[j], sep=":")
    cont = cont + 1
    } }
  # last <- as.numeric(last) # self + reciprocals ?
  levs <- levels(combLev) # All levels
  idx <- c()
    for(i in 1:length(last)){
    y <- which(levs == last[i])
    if(length(y) > 0) idx[i] <- y
    }
  idx <- idx[is.na(idx) == F]  # Added on 2/3/21
  levs <- as.character(levs[-idx]) # only crosses, without reciprocals

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
  rec }
}

REC.G3 <- function(P1, P2, type = "fix", data = NULL){
  # Reciprocal effects for designs with
  # no selfed parents
  # P1 <- df$Par1;P2 <- df$Par2
  # IT IS SUPERSEEDED
  if(!is.null(data)){
    P1Name <- deparse(substitute(P1))
    P2Name <- deparse(substitute(P2))
    P1 <- data[[P1Name]]
    P2 <- data[[P2Name]]
  }
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

  last <- c(); cont <- 1
  for(i in 1:p){ for(j in 1:i) {
    if(i != j) {
      last[cont] <- paste(i, j, sep="")
      cont = cont + 1 } } }
  # last
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


H.BAR <- function(P1, P2, type = "fix", data = NULL){
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

MDD <- function(P1, P2, type = "fix", data = NULL){
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

Hi <- function(P1, P2, type = "fix", data = NULL){
    if(!is.null(data)){
    P1Name <- deparse(substitute(P1))
    P2Name <- deparse(substitute(P2))
    P1 <- data[[P1Name]]
    P2 <- data[[P2Name]]
    }
  if(type == "random"){
    crosses <- ifelse(P1 == P2, 0, 1)
    Z <- GCA(P1, P2, type = "random") * crosses
    # colnames(Z) <- sub("combination", "", colnames(Z))
    # Z <- Z[, apply(Z, 2, function(x) !all(x==0))]
    return(Z)
  } else {

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
  H }
}

GCAC <- function(P1, P2, type = "fix", data = NULL){
  if(!is.null(data)){
    P1Name <- deparse(substitute(P1))
    P2Name <- deparse(substitute(P2))
    P1 <- data[[P1Name]]
    P2 <- data[[P2Name]]
  }
  if(type == "random"){
    selfs <- ifelse(P1 == P2, 1, 0)
    Z <- model.matrix(~ P1 - 1) * selfs
    # colnames(Z) <- sub("combination", "", colnames(Z))
    # Z <- Z[, apply(Z, 2, function(x) !all(x==0))]
    return(Z)
  } else {

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
  H}
}

DD <- function(P1, P2, type = "fix", data = NULL){
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

GCAmis <- function(P1, P2, type = "fix", data = NULL){
  # This is modified to work with mating design 4
  # in case of missing crosses, but it is supposed
  # to work always (to be tested)
  # Edited on 18/3/2023
  if(!is.null(data)){
    P1Name <- deparse(substitute(P1))
    P2Name <- deparse(substitute(P2))
    P1 <- data[[P1Name]]
    P2 <- data[[P2Name]]
  }
  if(type == "random"){
    Z <- sommer::overlay(P1, P2, sparse = F)
  } else {
    # tmp <- data.frame(P2, P1)
    # arrange(tmp, P1)
    P1 <- factor(as.character(P1))
    P2 <- factor(as.character(P2))
    levs <- c(levels(P1), levels(P2))
    levs <- levels(factor(levs, levels = unique(levs)))
    Z1n <- factor(P1, levels = levs, ordered = T)
    Z2n <- factor(P2, levels = levs)
    contrasts(Z1n) <- c("contr.sum")
    contrasts(Z2n) <- c("contr.sum")

    # Da qui cambia in caso di missing crosses
    Z1 <- model.matrix(~ Z1n - 1)
    Z2 <- model.matrix(~ Z2n - 1)
    Z <- (Z1 + Z2)
    nams <- paste("g_", levs[1:(length(levs))], sep="")
    colnames(Z) <- c(nams)
    tab <- checkScheme(P1, P2)

    # Individua quali parents non hanno missing crosses
    # Edited on 18/3/2023
    misCros <- nams[unique(as.vector(tab$missingCrosses))]
    if(!is.null(misCros)){
      # misCros <- paste("g_", misCros, sep ="")
      sel <- !(nams %in% misCros)
      wsel <- max(which(sel == TRUE))
      sel <- rep(TRUE, length(nams))
      sel[wsel] <- FALSE
    } else {
      sel <- rep(TRUE, length(nams))
      sel[length(sel)] <- FALSE
    }
    # tmp1 <- names(which(apply(tab$tab, 1, sum, na.rm = T) * apply(tab$tab, 2, sum, na.rm = T) == length(levs) - 1))
    # tmp1 <- tmp1[length(tmp1)]
    # tmp1 <- paste("g_", tmp1, sep = "")
    # sel <- colnames(Z) != tmp1
    # if(all(sel == FALSE)) sel[length(sel)] <- TRUE
    # Z[,!sel] == 1
    sel
    Z <- Z[,sel] - Z[,!sel]
  }
  Z
}

SCAmis <- function(P1, P2, type = "fix", data = NULL){
  # This is modified to work with mating design 4
  # in case of missing crosses, but it is supposed
  # to work always (to be tested)
  # Edited on 18/3/2023
  if(!is.null(data)){
    P1Name <- deparse(substitute(P1))
    P2Name <- deparse(substitute(P2))
    P1 <- data[[P1Name]]
    P2 <- data[[P2Name]]
  }
  if(type == "random"){
    crosses <- ifelse(P1 == P2, 0, 1)
    combination <- factor( ifelse(as.character(P1) <= as.character(P2),
                                  paste(P1, P2, sep =""),
                                  paste(P2, P1, sep ="")) )
    Z <- model.matrix(~ combination - 1) * crosses
    colnames(Z) <- sub("combination", "", colnames(Z))
    Z <- Z[, apply(Z, 2, function(x) !all(x==0))]
    return(Z)
  } else {

    # It is also used where the selfs are not includede
    P1 <- factor(as.character(P1)) #, levels = unique(P1))
    P2 <- factor(as.character(P2)) #, levels = unique(P1)) # Livelli uguali?
    P1c <- as.character(P1); P2c <- as.character(P2)

    # create the combination levels
    tmp <- ifelse(P1c < P2c, paste(P1c, P2c, sep =":"),
                  paste(P2c, P1c, sep = ":"))
    combination <- factor(tmp) #, levels = unique(tmp))
    combLev <- NA

    # Create the matings (considers the reciprocals)
    mating <- P1:P2
    p <- length(levels(factor(c(levels(P1), levels(P2)) )))
    n <- length(combination)

    # See whether selfs are included and find the last level for each
    # combination
    selflist <- levels(factor(combination[P1c == P2c]))
    levs <- sort(unique(c(levels(P1), levels(P2))))
    tmp <- paste(levs, max(levs), sep = ":")
    parLevs <- levs

    # Step 1. Create an empty matrix
    # Da qui cambia per la matrice sbilanciata
    # Populate the matrix from GCAmat
    gcamat <- GCAmis(P1, P2)
    np <- ncol(gcamat)
    n <- nrow(gcamat)
    tab <- checkScheme(P1, P2)$missingCrosses
    numMissing <- nrow(tab)

    # Create empty scamat
    scamat <- matrix(0, n, np * (np - 1)/2)
    nams <- rep(0, np * (np - 1)/2)
    cont <- 1
    # Cross multiplication of matrices
    for(i in 1:(np-1)){
      for(j in (i+1):np){
        scamat[,cont] <- gcamat[,i] * gcamat[,j]
        n1 <- substr(colnames(gcamat)[i], 3, nchar(colnames(gcamat)[j]))
        n2 <- substr(colnames(gcamat)[j], 3, nchar(colnames(gcamat)[j]))
        nams[cont] <- paste(n1,n2, sep = ":")
        cont <- cont + 1
      }
    }
    colnames(scamat) <- nams

    # Rimuove le colonne dei missing crosses
    if(!is.null(numMissing)){
      toRem <- c()
      for(i in 1:numMissing){
        # i <- 3
        sel <- paste(parLevs[tab[i,1]], parLevs[tab[i,2]], sep = ":")
        sel <- which(colnames(scamat)==sel)
        # print(sel)
        toRem <- c(toRem, sel)
      }
      scamat <- scamat[,-toRem]
    }

    # Rimuove l'ultima colonna e la sottrae dalle altre
    last <- length(scamat[1,])
    scamat <- scamat - scamat[,last]
    scamat <- scamat[,-last]
    # row.names()
    scamat
  }
}

