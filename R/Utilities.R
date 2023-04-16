emm.diallel <- function(obj){
  # This function operates on a summary.glht object, obtained from a
  # MET lmDiallel object. It calculates the marginal means across
  # environments.
  nam <- names(obj$test$coefficient)
  df <- data.frame("name" = nam)
  nam <- tidyr::separate(df, col="name", into = c("eff", "env"), sep = ":")
  df <- data.frame(nam, coef = obj$test$coefficient,
                   SEs = obj$test$sigm,
                   row.names = 1:length(nam[,1]))
  dw <- tidyr::pivot_wider(df, names_from = env, values_from = coef,
                  id_cols = eff, values_fn = mean)
  dw2 <- tidyr::pivot_wider(df, names_from = env, values_from = SEs,
                  id_cols = eff, values_fn = mean)
  dw <- as.data.frame(dw)
  effs <- apply(dw[,-1], 1, mean)
  dw2 <- as.data.frame(dw2)
  SEs <- apply(dw2[,-1], 1, function(row) sqrt(sum(row^2)/length(row)^2))
  data.frame("Estimate" = effs, "Std. Error" = SEs,
             "t value" = effs/SEs,
             "Pr(>|t|)" = 2 * pt(effs/SEs, obj$df),
             row.names = dw[,1],
             check.names = F)
}

checkScheme <- function(P1, P2){
  # This function checks the mating scheme and
  # checks whether there are missing crosses
  # Used internally and not exposed
  # Edited 15/4/2023
  P1 <- factor(P1)
  P2 <- factor(P2)
  P12 <- factor(paste(P1, P2, sep = ":"))
  tab <- table(P1, P2)
  parents <- unique(c(levels(P1), levels(P2)))
  crosses <- levels(P12) # selfs are included, if available

  # Buolds all possible (hypothetical) combinations
  selfsHyp <- paste(parents, parents, sep = ":")
  crossesHyp <- paste(combn(parents, 2)[1,], combn(parents, 2)[2,], sep = ":")
  recHyp <- paste(combn(parents, 2)[2,], combn(parents, 2)[1,], sep = ":")


  # finds which selfs/crosses/reciprocals are included
  # and checks for possible exchanges of parental roles
  ss <- selfsHyp %in% crosses
  cc <- crossesHyp %in% crosses
  rr <- recHyp %in% crosses
  rc <- c(matrix(cc,1) + matrix(rr,1))

  # Based on the above the scheme and balance are
  # detected
  if(any(rc > 1) & any(ss)){
    matingScheme <- 1
    if(all(rc == 2) & all(ss)){
      # balanced
      mis <- NULL
    } else {
      # missing crosses
      fulList <- c(selfsHyp, crossesHyp, recHyp)
      mis <- fulList[!(fulList %in% crosses)]
    }

  } else if(!any(rc > 1) & any(ss)){
    matingScheme <- 2
    if(all(rc == 1) & all(ss)){
      # balanced
      mis <- NULL
    } else {
      # missing crosses
      mis <- crossesHyp[!rc]
      }
  } else if(any(rc > 1) & !any(ss)){
    matingScheme <- 3
    if(all(rc == 2)){
      # balanced
      mis <- NULL
    } else {
      # missing crosses
      fulList <- c(crossesHyp, recHyp)
      mis <- fulList[!(fulList %in% crosses)]
    }
  } else if(!any(rc > 1) & !any(ss)){
    matingScheme <- 4
    if(all(rc == 1)){
      # balanced
      mis <- NULL
    } else {
      # missing crosses
      mis <- crossesHyp[!rc]
      }
  }

  if(!is.null(mis)){
    mis <- matrix(mis, length(mis), 1)
    mis <- data.frame(do.call(rbind, strsplit(as.character(mis[,1]), ":", fixed = T)))
    colnames(mis) <- c("Par1", "Par2")
    parMis <- sort(unique(unlist(mis)))
    parNoMis <- parents[!(parents %in% parMis)]
  } else {
    parNoMis <- parents
  }

  list(matingScheme = matingScheme, missingCrosses = mis,
       # missingSelfs = t(misSelfs),
       parNoMis = parNoMis, tab = tab)
}

int.matrix <- function(Xa, Xb){
  nca <- length(Xa[1,])
  ncb <- length(Xb[1,])
  nra <- length(Xa[,1])
  nrb <- length(Xb[,1])
  if(nra != nrb) stop ("The two matrices need to have the same number of rows")

  # Create an empty matrix
  intmat <- matrix(0, nra, nca*ncb)
  nams <- rep(0, nca*ncb)
  cont <- 1

  # Populate the new matrices
  for(j in 1:ncb){
    for(i in 1:nca){
      intmat[,cont] <- Xa[,i] * Xb[,j]
      n1 <- colnames(Xa)[i]
      n2 <- colnames(Xb)[j]
      nams[cont] <- paste(n1,n2, sep = ":")
      cont <- cont + 1
    }
  }
  colnames(intmat) <- nams
  intmat
}

