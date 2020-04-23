library(devtools)
install_github("OnofriAndreaPG/lmDiallel")
library(lmDiallel)
data(hayman54)
head(hayman54)

fit <- lm.diallel(Ftime ~ Par1 + Par2, Block = Block,
                  data = hayman54,
                  fct = "HAYMAN2")
summary(fit)
anova(fit)

