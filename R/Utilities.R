checkScheme <- function(P1, P2){
  # This checks the mating scheme and checks whether there are missing crosses
  # Edited 10/3/2023
  levs <- unique(c(levels(P1), levels(P2)))
  P1 <- factor(P1, levels = levs)
  P2 <- factor(P2, levels = levs)
  tab <- table(P1, P2)
  selfs <- !all(diag(tab) == 0)
  fullUpper <- any(tab[upper.tri(tab, diag = FALSE)] != 0)
  fullLower <- any(tab[lower.tri(tab, diag = FALSE)] != 0)
  missingUpper <- any(tab[upper.tri(tab, diag = FALSE)] == 0)
  missingLower <- any(tab[lower.tri(tab, diag = FALSE)] == 0)

  if(selfs & fullUpper & fullLower) {
    matingScheme <- 1
    if(missingUpper | missingLower) {
      mis <- which(tab == 0, arr.ind = T)
      tab[tab == 0] <- NA
      } else {mis <- NULL}
  } else if (!selfs & fullUpper & fullLower)  {
    matingScheme <- 3
    diag(tab) <- NA
    if(missingUpper | missingLower) {
      mis <- which(tab == 0, arr.ind = T)
      mis <- mis[mis[,1] != mis[,2],]
      tab[tab == 0] <- NA
      } else { mis <- NULL }
  } else if (selfs & (!fullUpper | !fullLower)){
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
 list(matingScheme = matingScheme, missingCrosses = mis, tab = tab)
}
