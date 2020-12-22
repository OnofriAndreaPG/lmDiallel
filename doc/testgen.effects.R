library(devtools)
install_github("OnofriAndreaPG/lmDiallel")
library(lmDiallel)

data(hayman54)
fit <- lm.diallel(Ftime ~ Par1 + Par2, Block = Block,
                  data = hayman54,
                  fct = "HAYMAN1")
summary(fit)
anova(fit)
source("R/glht.diallel_v1.1.R")

diallel.eff(fit)
GCA.eff(hayman54$Par1, hayman54$Par2)

data(hayman54)
head(hayman54)
fit <- lm.diallel(Ftime ~ Par1 + Par2, Block = Block,
                  data = hayman54,
                  fct = "HAYMAN2")
summary(fit)
anova(fit)
obj <- fit
hayman2.eff(obj)
library(multcomp)
gh <- glht(linfct = diallel.eff(fit))
summary(gh, test = adjusted(type = "none"))

