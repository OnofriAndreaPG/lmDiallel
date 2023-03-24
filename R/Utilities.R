checkScheme <- function(P1, P2){
  # This checks the mating scheme and checks whether there are missing crosses
  # Edited 10/3/2023
  # P1 <- dfb$Par1; P2 <- dfb$Par2
  P1 <- factor(P1)
  P2 <- factor(P2)
  levs <- unique(c(levels(P1), levels(P2)))
  P1 <- factor(P1, levels = levs)
  P2 <- factor(P2, levels = levs)
  tab <- table(P1, P2)
  selfs <- !all(diag(tab) == 0)
  fullUpper <- any(tab[upper.tri(tab, diag = FALSE)] != 0)
  fullLower <- any(tab[lower.tri(tab, diag = FALSE)] != 0)
  missingUpper <- any(tab[upper.tri(tab, diag = FALSE)] == 0)
  missingLower <- any(tab[lower.tri(tab, diag = FALSE)] == 0)
  missingSelfs <- any(diag(tab) == 0)

  if(selfs & fullUpper & fullLower) {
    # Full Diallel
    matingScheme <- 1
    if(missingUpper | missingLower | missingSelfs) {
      mis <- which(tab == 0, arr.ind = T)
      # misSelfs <- which(diag(tab) == 0, arr.ind = T)
      # tab[tab == 0] <- NA
    } else {
        mis <- NULL
        }
  } else if (!selfs & fullUpper & fullLower)  {
    # no selfs
    matingScheme <- 3
    diag(tab) <- NA
    if(missingUpper | missingLower) {
      mis <- which(tab == 0, arr.ind = T)
      mis <- mis[mis[,1] != mis[,2],]
      tab[tab == 0] <- NA
      } else { mis <- NULL }
  } else if (selfs & (!fullUpper | !fullLower)){
      # no reciprocals
      matingScheme <- 2
      if(fullUpper & missingUpper) {
        tab[lower.tri(tab, diag = F)] <- NA
        mis <- which(tab == 0, arr.ind = T)
        tab[tab == 0] <- NA
      } else if(fullLower & missingLower){
        tab[upper.tri(tab == 0, diag = F)] <- NA
        mis <- which(tab == 0, arr.ind = T)
        tab[tab == 0] <- NA
      } else { mis <- NULL
      if(fullLower) tab[upper.tri(tab, diag = F)] <- NA else tab[lower.tri(tab, diag = F)] <- NA
      }
  } else if (!selfs & (!fullUpper | !fullLower)){
      # no selfs, no reciprocals
      matingScheme <- 4
      diag(tab) <- NA
      if(fullUpper & missingUpper) {
        tab[lower.tri(tab, diag = F)] <- NA
        mis <- which(tab == 0, arr.ind = T)
        tab[tab == 0] <- NA
      } else if(fullLower & missingLower){
        tab[upper.tri(tab == 0, diag = F)] <- NA
        mis <- which(tab == 0, arr.ind = T)
        tab[tab == 0] <- NA
      } else { mis <- NULL
      if(fullLower) tab[upper.tri(tab, diag = F)] <- NA else tab[lower.tri(tab, diag = F)] <- NA
      }

  }
  # Find parents with no missing crosses (in any)
  # levs <- unique(c(levels(P1), levels(P2)))
  parMis <- sort(unique(c(mis)))
  parNoMis <- which(!(levs %in% levs[parMis]))
  list(matingScheme = matingScheme, missingCrosses = mis,
       # missingSelfs = t(misSelfs),
       parNoMis = levs[parNoMis], tab = tab)
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

