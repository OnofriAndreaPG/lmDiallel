gen.eff <- function(obj){
  # Get the data
    assign <- attr(model.matrix(obj), "assign")
    P1 <- obj$model[,2]
    P2 <- obj$model[,3]
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
# Yield <- c(12, 13, 14, 11, 15, 21, 17, 16, 19)
# Par1 <- c("A", "A", "A", "B", "B", "B", "C", "C", "C") #father
# Par2 <- c("A", "B", "C", "A", "B", "C","A", "B", "C") #mother
# df <- data.frame(Par1, Par2, Yield)
# rm(Yield, Par1, Par2)
# df
# library(lmDiallel)
# dMod2 <- lm.diallel(Yield ~ Par1 + Par2, data = df, fct = "HAYMAN1")
# obj <- dMod2
# summary(dMod2)
#
# GEeff(dMod2)
#
# str(dMod2)
# GCAeff(df$Par1, df$Par2)
# P1 <- c("A", "A", "B") #father
# P2 <- c("B", "C", "C") #mother

diallel.eff <- function(obj, MSE = NULL, dfr = NULL) {

    k <- gen.eff(obj)
    linfct.list <- list(linfct = k, MSE = MSE, dfr = dfr)
    class(linfct.list) <- "diallelMod"
    # attr(linfct, "MSE") <- MSE
    # attr(linfct, "dfr") <- dfr
    # class(linfct) <- "diallelMod"
    # print(linfct)
    return(linfct.list)
    # linfct <- list(...)
    #
    # linfct <- lapply(linfct, function(x) {
    #     if (is.numeric(x) && !is.matrix(x)) {
    #         return(matrix(x, nrow = 1))
    #     } else {
    #         return(x)
    #     }})
    #
    # if (is.null(names(linfct)))
    #     stop(sQuote("linfct"), " doesn't have a ", sQuote("names"),
    #          " attribute")
    #
    # classes <- sapply(linfct, function(x) inherits(x, "matrix") ||
    #                                       inherits(x, "character"))
    #
    # if (length(linfct) == 1 && linfct[[1]] == "Means") {
    #     class(linfct) <- "means"
    #     return(linfct)
    # }
    #
    # attr(linfct, "interaction_average") <- interaction_average
    # attr(linfct, "covariate_average") <- covariate_average
    #
    # if (all(classes)) {
    #     class(linfct) <- "mcp"
    #     return(linfct)
    # }
    #
    # stop("Arguments don't consist of either matrices or characters")
}


### multiple comparison procedures
glht.diallelMod <- function(obj, linfct, ...) {

    ### extract factors and contrast matrices from `model'
    MSE <- linfct$MSE
    dfr <- ifelse(is.null(linfct$dfr), obj$df.residual, linfct$dfr)
    k <- linfct$linfct
    coefMod <- coef(obj)
    vcovMod <- vcov(obj, MSE = MSE)
    args <- list(coef = coefMod, vcov = vcovMod, df = 26)
    class(args) <- "parm"
    # k <- matrix(linfct, nrow(linfct), ncol(linfct))
    # row.names(k) <- row.names(linfct)
    # print(args)
    # args <- list(model = args, linfct = k)

    ret <- multcomp::glht(args, k)
    # summary(ret)
    return(ret)
    # # print(args)
    # ret <- glht(args, linfct)

    # if (ia || ca) {
    #     ### experimental version
    #     tmp <- mcp2matrix2(model, linfct = linfct, interaction_average = ia,
    #                        covariate_average = ca)
    # } else {
    #     ### use old version
    #     tmp <- mcp2matrix(model, linfct = linfct)
    # }
    # args <- list(model = model, linfct = tmp$K)
    # if (!is.null(tmp$alternative))
    #     args$alternative <- tmp$alternative
    # if (any(tmp$m != 0))
    #     args$rhs <- tmp$m
    # args <- c(args, list(...))

    # ret <- do.call("glht", args)
    # ret$type <- tmp$type
    # ret$focus <- names(linfct)
    # return(ret)
}
