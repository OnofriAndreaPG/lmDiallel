# WORKING WITH THE MEANS ####################################
# Example 1. From Hayman 1954
# MSE: 416.8273; DF: 63; 2 Blocks

rm(list = ls())
library(sommer)
library(lmDiallel)

data("hayman54")
df <- hayman54
df$Block <- factor(df$Block)
df$Par1 <- factor(df$Par1)
df$Par2 <- factor(df$Par2)

# mod. ANOVA tradizionale
mod.aov <- lm(Ftime ~ Par1 + Par2, data = df)
anova(mod.aov)
gh <- glht(mod.aov, linfct = mcp(Par1 = "Tukey"))
summary(gh)

# Hayman - Model 1 ################
fit <- lm.diallel(Ftime ~ Par1 + Par2, data = df,
                  fct = "HAYMAN1")
summary(fit)
anova(fit)
gh <- glht(linfct = diallel.eff(fit), adjust = "none")
summary(gh)

fit2 <- lm(Ftime ~ GCA(Par1, Par2) +
             tSCA(Par1, Par2) + RGCA(Par1, Par2) +
             RSCA(Par1, Par2), data = df)
summary(fit2)
anova(fit2)
summary.diallel(fit2, MSE = 208.4136, dfr = 63)
anova.diallel(fit2, MSE = 208.4136, dfr = 63)

# Hayman - Model 2 ################
source("lm.diallel.R")
source("model.matrix.diallel.r")
fit <- lm.diallel(Ftime ~ Par1 + Par2, data = df,
                  fct = "HAYMAN2")
summary(fit, MSE = 208.4136, dfr = 63)
anova(fit, MSE = 208.4136, dfr = 63)
anova(fit)[,2] * 2

fit2 <- lm(Ftime ~ GCA(Par1, Par2) +
             MDD(Par1, Par2) + DD(Par1, Par2) +
             SCA(Par1, Par2) +
             RGCA(Par1, Par2) +
             RSCA(Par1, Par2), data = df)
summary.diallel(fit2, MSE = 208.4136, dfr = 63)
aovTab <- anova.diallel(fit2, MSE = 208.4136, dfr = 63)
aovTab
aovTab[,2] * 2

# Griffing - Method 1 ###########################
# Including reciprocal effects
fit <- lm.diallel(Ftime ~ Par1 + Par2, data = df,
                  fct = "GRIFFING1")
round(summary(fit, MSE = 208.4136, dfr = 63), 3)
anova(fit, MSE = 208.4136, dfr = 63)

fit2 <- lm(Ftime ~ GCA(Par1, Par2) +
             tSCA(Par1, Par2) +
             REC(Par1, Par2), data = df)
summary(fit2)
summary.diallel(fit2, MSE = 208.4136, dfr = 63)
anova.diallel(fit2, MSE = 208.4136, dfr = 63)

# Griffing - Method 2 #################
# Solo GCA e tSCA - no reciproci solo per
# disegni 2 e 4
source("lm.diallel.R")
source("model.matrix.diallel.r")
fit <- lm.diallel(Ftime ~ Par1 + Par2, data = df,
                  fct = "GRIFFING2")
summary(fit)
anova(fit, MSE = 208.4136, dfr = 63)
anova(fit)[,2] * 2

fit2 <- lm(Ftime ~ GCA(Par1, Par2) +
             tSCA(Par1, Par2), data = df)
summary(fit2)
summary.diallel(fit2, MSE = 208.4136, dfr = 63)
anova.diallel(fit2, MSE = 208.4136, dfr = 63)

# GE2 ########################################
source("lm.diallel.R")
source("model.matrix.diallel.r")
fit <- lm.diallel(Ftime ~ Par1 + Par2, data = df,
                  fct = "GE2")
summary(fit)
tabaov <- anova(fit, MSE = 208.4136, dfr = 63)
tabaov
tabaov[,2]*2

fit2 <- lm(Ftime ~ H.BAR(Par1, Par2) +
             VEi(Par1, Par2) + Hi(Par1, Par2) +
             SCA(Par1, Par2), data = df)
anova(fit2)

fit <- lm.diallel(Ftime ~ Par1 + Par2, data = df,
                  fct = "GE2r")
round(summary(fit, MSE = 208.4136, dfr = 63), 4)
tabaov <- anova(fit, MSE = 208.4136, dfr = 63)
tabaov
tabaov[,2]*2

fit2 <- lm(Ftime ~ H.BAR(Par1, Par2) +
             VEi(Par1, Par2) + Hi(Par1, Par2) +
             SCA(Par1, Par2) +
             REC(Par1, Par2), data = df)
anova(fit2)

# GE3 ############################################
# GCA + SCA scomposta (senza hi) + selfed parents
source("lm.diallel.R")
source("model.matrix.diallel.r")
fit <- lm.diallel(Ftime ~ Par1 + Par2, data = df,
                  fct = "GE3")
summary(fit, MSE = 208.4136, dfr = 63)
anova(fit, MSE = 208.4136, dfr = 63)
anova(fit)[,2] * 2

fit2 <- lm(Ftime ~ H.BAR(Par1, Par2) +
             SP(Par1, Par2) + GCAC(Par1, Par2) +
             SCA(Par1, Par2), data = df)
anova.diallel(fit2, MSE = 208.4136, dfr = 63)


fit <- lm.diallel(Ftime ~ Par1 + Par2, data = df,
                  fct = "GE3r")
round(summary(fit, MSE = 208.4136, dfr = 63), 4)
anova(fit, MSE = 208.4136, dfr = 63)
anova(fit)[,2] * 2

fit2 <- lm(Ftime ~ H.BAR(Par1, Par2) +
             SP(Par1, Par2) + GCAC(Par1, Par2) +
             SCA(Par1, Par2) + REC(Par1, Par2), data = df)
anova.diallel(fit2, MSE = 208.4136, dfr = 63)


# Example 2 Lonnquist ###########################
# Method = 2
rm(list=ls())
library(agridat)
data(lonnquist.maize)
dat <- lonnquist.maize
dat <- transform(dat,
                 p1=factor(p1,
                   levels=c("C","L","M","H","G","P","B","RM","N","K","R2","K2")),
                 p2=factor(p2,
                   levels=c("C","L","M","H","G","P","B","RM","N","K","R2","K2")))
d2 <- subset(dat, is.element(p1, c("M","H","G","B","K","K2")) &
                      is.element(p2, c("M","H","G","B","K","K2")))
d2 <- droplevels(d2)
df <- d2
names(df) <- c("Par1", "Par2", "Yield")
df
# lonnquist61 <- df
# save(lonnquist61, file = "lonnquist61.rda")

# Traditional
mod.aov <- lm(Yield ~ Par1 + Par2, data = df)
anova(mod.aov)
residuals(mod.aov)

# Hayman analysis (in Gardner and Eberhart)
source("lm.diallel.R")
source("model.matrix.diallel.R")

# Model 1 = Griffing - Effects are not orthogonal to each other !!!!!
fit2 <- lm(Yield ~ GCA(Par1, Par2), data = df)
summary(fit2)
cbind(df, residuals(fit2))

96.157 + 1.4875 - 2.1625 + 5.81785714 #(media_BG)
fit3 <- lm(Yield ~ GCA(Par1, Par2) +
             tSCA(Par1, Par2), data = df)

96.7750 + 1.2917 - 2.1250 + 5.3583 #(media_BG)


summary(fit3)
summary.diallel(fit3, MSE = 7.1, dfr = 60)
anova.diallel(fit3, MSE = 7.1, dfr = 60)

# Model 2
fit2 <- lm(Yield ~ GCA(Par1, Par2)
           + MDD(Par1, Par2) +
             DD(Par1, Par2) +
             SCA(Par1, Par2), data = df)
summary.diallel(fit2, MSE = 7.1, dfr = 60)
anova.diallel(fit2, MSE = 7.1, dfr = 60)

coefs[1] + 2 * coefs[2] - 5 * coefs[7] - 4 * coefs[8]
96.157 + coefs[2] + 1.483 + coefs[8]
coefs[1] - 5 * coefs[7] - 4 * coefs[8]
coefs[1] + coefs[7]

# Griffing - Method 2 ###############################
# Solo GCA e SCAt - no reciproci
source("lm.diallel.R")
source("model.matrix.diallel.R")
fit <- lm.diallel(Yield ~ Par1 + Par2, data = df,
                  fct = "GRIFFING2")
summary(fit)
summary(fit, MSE = 7.1, dfr = 60)
anova(fit)
anova(fit, MSE = 7.1, dfr = 60)

fit2 <- lm(Yield ~ GCA(Par1, Par2) + tSCA(Par1, Par2), data = df)
summary(fit2)
summary.diallel(fit2, MSE = 7.1, dfr = 60)
anova(fit2)
anova.diallel(fit2, MSE = 7.1, dfr = 60)

# GE 2 ###########################################
# La GCA viene decomposta e le medie sono diverse
# per selfed e crossed
#GE2
source("lm.diallel.R")
source("model.matrix.diallel.R")
fit <- lm.diallel(Yield ~ Par1 + Par2, data = df,
                  fct = "GE2")

summary(fit, MSE = 7.1, dfr = 60)
summary(fit)
anova(fit, MSE = 7.1, dfr = 60)
residuals(fit)

fit2 <- lm(Yield ~ H.BAR(Par1, Par2) +
             VEi(Par1, Par2) + Hi(Par1, Par2) +
             SCA(Par1, Par2), data = df)
summary.diallel(fit2, MSE = 7.1, dfr = 60)
summary(fit2)
anova(fit2)
anova.diallel(fit2, MSE = 7.1, dfr = 60)
residuals(fit2)


source("heterosis_fun.R")
res <- heterosis(df$Par1, df$Par2, Value = df$Yield)
res

92.450 - 1.450 #MM: 91
92.450 - 0.75 #HH: 91.7
92.450 - 0.5*1.450 - 0.5*0.75 #Media = 91.35
98.8 - 91.35 # Eterosi abs: 7.45
7.45/91.35*100 #Eterosi perc.
98.8 - 91.7 #BPH absol.: 7.1
5.190 - 1.475 + 0.35 + 3.385 #MPH: From model parameters
5.190 - 1.475 + 3.385 # BPH: from model parameters (abs)
# From the best


# GE3 ###############################################
# GCA + SCA scomposta (senza hi) + selfed parents
source("lm.diallel.R")
source("model.matrix.diallel.R")
fit <- lm.diallel(Yield ~ Par1 + Par2, data = df,
                  fct = "GE3")
summary(fit, MSE = 7.1, dfr = 60)
anova(fit, MSE = 7.1, dfr = 60)
residuals(fit)

fit2 <- lm(Yield ~ H.BAR(Par1, Par2) +
             SP(Par1, Par2) + GCAC(Par1, Par2) +
             SCA(Par1, Par2), data = df)
summary.diallel(fit2, MSE = 7.1, dfr = 60)
anova.diallel(fit2, MSE = 7.1, dfr = 60)

# Example of mating method 4 ###########################
# Balanced design ######################################
rm(list =ls())
fileName <- "https://www.casaonofri.it/_datasets/Griffing.csv"
df <- read_csv(fileName)

df$Par1 <- factor(df$Par1)
df$Par2 <- factor(df$Par2)

library(lmDiallel)
# source("lm.diallel.R")
# source("model.matrix.diallel.R")

# fit <- lm.diallel(Yield ~ Par1 + Par2, data = df,
#                    fct = "GRIFFING2") # Non funge!
# summary(fit)
# anova(fit, MSE = 21.05, dfr = 2558)
# residuals(fit)

fit <- lm(Yield ~ GCA(Par1, Par2) + SCA(Par1, Par2), data = df)
summary.diallel(fit, MSE = 21.05, dfr = 2558)
anova.diallel(fit, MSE = 21.05, dfr = 2558)

