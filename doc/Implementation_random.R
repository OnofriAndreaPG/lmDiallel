# Reproducing the results in Piepho e Mohring, 2011 ##########
# Example 1, Hyman 1954. MODEL 1: full diallel ###############
rm(list = ls())
library(sommer)
library(lmDiallel)

data("hayman54")
df <- hayman54 # Working with row data
df$Block <- factor(df$Block)
df$Par1 <- factor(df$Par1)
df$Par2 <- factor(df$Par2)

# Dummy variables for selfs, crosses, combinations
df$crosses <- ifelse(df$Par1 == df$Par2, 0, 1)
df$selfs <- ifelse(df$Par1 == df$Par2, 1, 0)
df$dr <- ifelse(as.character(df$Par1) < as.character(df$Par2), -1,
                ifelse(as.character(df$Par1) == as.character(df$Par2), 0, 1))
df$combination <- factor( ifelse(as.character(df$Par1) <= as.character(df$Par2),
                                 paste(df$Par1, df$Par2, sep =""),
                                 paste(df$Par2, df$Par1, sep ="")) )

# Hayman 1 ##############################################
# Equation
# mod1r <- asreml(Ftime ~ 1, data=df,
#                 random = ~ Block + Par1 + and(Par2) # GCA
#                 + Par1:dr + and(Par2:dr) # RGCA: exclude selfs
#                 + combination # SCA
#                 + combination:dr) #RSCAA: exclude selfs
# summary(mod1r)$varcomp

mod1h <- mmer(Ftime ~ 1, data=df,
              random = ~ Block + overlay(Par1, Par2)
              + overlay(Par1, Par2):dr
              + combination
              + combination:dr)
summary(mod1h)$varcomp
mod1m <- mmer.diallel(Ftime ~ Par1 + Par2, Block = Block, data = df,
                      fct = "HAYMAN1")
mod1m

mod2m <- mmer.diallel(Ftime ~ Par1 + Par2, data = df,
                      fct = "HAYMAN1")
mod2m


# Hayman 2 ##############################################
# Parting RSCA e RGCA and adding heterosis
# mod1r <- asreml(Ftime ~ Block + crosses, data=df,
#                 random = ~ Block + crosses
#                 + Par1 + and(Par2) # GCA
#                 + Par1:dr + and(Par2:dr) #RGCA
#                 + Par1:selfs + and(Par2:selfs) #h.i
#                 + combination #SCA
#                 + combination:dr) # RSCA
# summary(mod1r)$varcomp

mod1h <- mmer(Ftime ~ 1, data=df,
              random = ~ Block + overlay(Par1, Par2) # GCA
              + crosses #MDD
              + overlay(Par1, Par2):dr #RGCA
              + overlay(Par1, Par2):selfs #h.i
              + combination #SCA
              + combination:dr) # RSCA
summary(mod1h)$varcomp

mod1m <- mmer.diallel(Ftime ~ Par1 + Par2, Block = Block,
                      data = df, fct = "HAYMAN2")
mod1m

mod1m <- mmer.diallel(Ftime ~ Par1 + Par2,
                      data = df, fct = "HAYMAN2")
mod1m


# GRIFFING MODEL 1 ##############################
# GCA + tSCA + Reciprocals
# mod1r <- asreml(Ftime ~ 1, data=df,
#                 random = ~ Block
#                 + Par1 + and(Par2) #GCA
#                 + combination #SCA
#                 + combination:dr) #exclude selfs, add reciprocals
# summary(mod1r)$varcomp
mod1h <- mmer(Ftime ~ 1, data=df,
              random = ~ Block
              + overlay(Par1, Par2) # GCA
              + combination #tSCA
              + combination:dr) # REC
summary(mod1h)$varcomp

mod1m <- mmer.diallel(Ftime ~ Par1 + Par2, Block = Block,
                      data = df, fct = "GRIFFING1")
mod1m

mod1m <- mmer.diallel(Ftime ~ Par1 + Par2,
                      data = df, fct = "GRIFFING1")
mod1m

# GRIFFING MODEL 2 ##############################
# GCA + tSCA, no Reciprocals
# mod1r <- asreml(Ftime ~ 1, data=df,
#                 random = ~ Block
#                 + Par1 + and(Par2) #GCA
#                 + combination) #SCA
# summary(mod1r)$varcomp
dfNoRec <- df %>%
  filter(as.character(Par1) <= as.character(Par2)) #Removing reciprocals

mod1h <- mmer(Ftime ~ 1, data=dfNoRec,
              random = ~ Block
              + overlay(Par1, Par2) # GCA
              + combination) #tSCA
summary(mod1h)$varcomp

mod1m <- mmer.diallel(Ftime ~ Par1 + Par2, Block = Block,
                      data = dfNoRec, fct = "GRIFFING2")
mod1m

mod1m <- mmer.diallel(Ftime ~ Par1 + Par2,
                      data = dfNoRec, fct = "GRIFFING2")
mod1m

# GRIFFING MODEL 3 ##############################
# GCA + tSCA + Rec, selfs
dfNoSel <- df %>%
  filter(as.character(Par1) != as.character(Par2)) # removing selfs

mod1h <- mmer(Ftime ~ 1, data=dfNoSel,
              random = ~ Block
              + overlay(Par1, Par2) # GCA
              + combination #tSCA
              + combination:dr) # REC
summary(mod1h)$varcomp

mod1m <- mmer.diallel(Ftime ~ Par1 + Par2, Block = Block,
                      data = dfNoSel, fct = "GRIFFING3")
mod1m

mod1m <- mmer.diallel(Ftime ~ Par1 + Par2,
                      data = dfNoSel, fct = "GRIFFING3")
mod1m

# GRIFFING MODEL 4 ##############################
# GCA + SCA no Rec no selfs
dfNone <- df %>%
  filter(as.character(Par1) != as.character(Par2)) %>% # removing selfs
  filter(as.character(Par1) <= as.character(Par2))

mod1h <- mmer(Ftime ~ 1, data = dfNone,
              random = ~ Block
              + overlay(Par1, Par2) # GCA
              + combination) #SCA
summary(mod1h)$varcomp

mod1m <- mmer.diallel(Ftime ~ Par1 + Par2, Block = Block,
                      data = dfNone, fct = "GRIFFING4")
mod1m

mod1m <- mmer.diallel(Ftime ~ Par1 + Par2,
                      data = dfNone, fct = "GRIFFING4")
mod1m


# GE2r ###########################################
# with reciprocals
# mod1r <- asreml(Ftime ~ Block + crosses, data=df,
#                 random = ~ Par1 + and(Par2) # GCA
#                 + Par1:crosses + and(Par2:crosses)  #h.i
#                 + combination:crosses) # SCA
# summary(mod1r)$varcomp

mod1h <- mmer(Ftime ~ crosses, data=df,
              random = ~ Block + overlay(Par1, Par2) # GCA
              + overlay(Par1, Par2):crosses  #h.i
              + combination:crosses # SCA
              + combination:crosses:dr)
summary(mod1h)$varcomp

mod1m <- mmer.diallel(Ftime ~ Par1 + Par2, Block = Block,
                      data = df, fct = "GE2r")
mod1m

mod1m <- mmer.diallel(Ftime ~ Par1 + Par2, Block = Block,
                      data = df, fct = "GE2r")
mod1m

mod1m <- mmer.diallel(Ftime ~ Par1 + Par2,
                      data = df, fct = "GE2r")
mod1m

# GE2 ###########################################
# without reciprocals
# mod1r <- asreml(Ftime ~ Block + crosses, data=df,
#                 random = ~ Par1 + and(Par2) # GCA
#                 + Par1:crosses + and(Par2:crosses)  #h.i
#                 + combination:crosses) # SCA
# summary(mod1r)$varcomp

mod1h <- mmer(Ftime ~ crosses, data=df,
              random = ~ Block + overlay(Par1, Par2) # GCA
              + overlay(Par1, Par2):crosses  #h.i
              + combination:crosses)
summary(mod1h)$varcomp

mod1m <- mmer.diallel(Ftime ~ Par1 + Par2, Block = Block,
                      data = df, fct = "GE2")
mod1m

mod1m <- mmer.diallel(Ftime ~ Par1 + Par2,
                      data = df, fct = "GE2")
mod1m

# GE3r ################################
# with reciprocals. In Mohring: mixed, model 3 reduced
# mod1r <- asreml(Ftime ~ Block + crosses, data=df,
#                 random = ~ Par1:crosses + and(Par2:crosses) #GCA solo crosses
#                 + Par1:selfs # Selfed parents
#                 + combination:crosses # SCA solo crosses
#                 + combination:dr) # reciprocals
# summary(mod1r)$varcomp
mod1h <- mmer(Ftime ~ crosses, data=df,
                random = ~ Block + overlay(Par1, Par2):crosses #GCA, only crosses
                + Par1:selfs # Selfed parents
                + combination:crosses # SCA solo crosses
                + combination:dr) # reciprocals
summary(mod1h)$varcomp

mod1m <- mmer.diallel(Ftime ~ Par1 + Par2, Block = Block,
                      data = df, fct = "GE3r")
mod1m

mod1m <- mmer.diallel(Ftime ~ Par1 + Par2,
                      data = df, fct = "GE3r")
mod1m


# GE3 ############################################
# no reciprocals
# mod1r <- asreml(Ftime ~ Block + crosses, data=df,
#                 random = ~ Par1:crosses + and(Par2:crosses) #GCA solo crosses
#                 + Par1:selfs # Selfed parents
#                 + combination:crosses) # SCA solo crosses
# summary(mod1r)$varcomp

mod1h <- mmer(Ftime ~ crosses, data=df,
              random = ~ Block + overlay(Par1, Par2):crosses #GCA solo crosses
              + Par1:selfs # Selfed parents
              + combination:crosses) # SCA solo crosses
summary(mod1h)$varcomp

mod1m <- mmer.diallel(Ftime ~ Par1 + Par2, Block = Block,
                      data = df, fct = "GE3")
mod1m

mod1m <- mmer.diallel(Ftime ~ Par1 + Par2,
                      data = df, fct = "GE3")
mod1m

# GE3 with RGCA + RSCA ##########################
#In Mohring: mixed, model 3 completo
# mod1r <- asreml(Ftime ~ Block + crosses, data=df,
#                 random = ~ Par1:crosses + and(Par2:crosses) #GCA solo crosses
#                 + Par1:selfs # Selfed parents
#                 + combination:crosses # SCA solo crosses
#                 + Par1:dr + and(Par2:dr) # RGCA: exclude selfs
#                 + combination:dr) #RSCA: exclude selfs
# summary(mod1r)$varcomp

# Not included in lmDiallel
mod1h <- mmer(Ftime ~ Block + crosses, data=df,
              random = ~ overlay(Par1, Par2):crosses #GCA solo crosses
              + Par1:selfs # Selfed parents
              + combination:crosses # SCA solo crosses
              + overlay(Par1,Par2):dr # RGCA: exclude selfs
              + combination:dr) #RSCA: exclude selfs

summary(mod1h)$varcomp


# WORKING WITH THE MEANS ####################################
# Example 1. From Hayman 1954
rm(list = ls())
library(sommer)
library(lmDiallel)
library(tidyverse)

data("hayman54")

df <- hayman54 %>%  # Working with means
  group_by(Par1, Par2) %>%
  summarise(Ftime = mean(Ftime))
df$Par1 <- factor(df$Par1)
df$Par2 <- factor(df$Par2)


# Dummy variables for selfs, crosses, combinations
df$crosses <- ifelse(df$Par1 == df$Par2, 0, 1)
df$selfs <- ifelse(df$Par1 == df$Par2, 1, 0)
df$dr <- ifelse(as.character(df$Par1) < as.character(df$Par2), -1,
                ifelse(as.character(df$Par1) == as.character(df$Par2), 0, 1))
df$combination <- factor( ifelse(as.character(df$Par1) <= as.character(df$Par2),
                                 paste(df$Par1, df$Par2, sep =""),
                                 paste(df$Par2, df$Par1, sep ="")) )


# Hayman - Model 1 ################
# GCA + tSCA + RGCA + RSCA
# # Mate design 1 or 3
# mod1 <- asreml(Ftime ~ 1, data = df,
#                random = ~ Par1 + and(Par2) # GCA
#                + Par1:dr + and(Par1:dr) # RGCA
#                + combination) #SCA. RSCA è residuo + pure error
# summary(mod1)$varcomp
modh <- mmer(Ftime ~ 1,
             random = ~ overlay(Par1, Par2)
             + overlay(Par1, Par2):dr
             + combination, verbose = F,
             data=df)
summary(modh)$varcomp

mod1m <- mmer.diallel(Ftime ~ Par1 + Par2, data = df,
                      fct = "HAYMAN1")
mod1m

# Hayman - Model 2 ################
# with heterosis
# mod1 <- asreml(Ftime ~ crosses, data = df,
#                random = ~ Par1 + and(Par2) # GCA
#                + Par1:dr + and(Par1:dr) # RGCA
#                + Par1:selfs + and(Par2:selfs) #h.i
#                + combination) #SCA. RSCA è residuo + pure error
# summary(mod1)$varcomp

modh <- mmer(Ftime ~ 1,
             random = ~ crosses + overlay(Par1, Par2)
             + overlay(Par1, Par2):dr
             + overlay(Par1, Par2):selfs
             + combination,
             data=df)
summary(modh)$varcomp #Sembra più stabile

mod1m <- mmer.diallel(Ftime ~ Par1 + Par2, data = df,
                      fct = "HAYMAN2")
mod1m

# GRIFFING MODEL1 #######################
# GCA + SCA + Reciprocals (come residuo)
# mod1 <- asreml(Ftime ~ crosses, data = df,
#                random = ~ Par1 + and(Par2) # GCA
#                + combination) #SCA. Reciprocals è residuo + pure error
# summary(mod1)$varcomp

modh <- mmer(Ftime ~ 1,
             random = ~ overlay(Par1, Par2)
             + combination,
             data=df)
summary(modh)$varcomp #Sembra più stabile

mod1m <- mmer.diallel(Ftime ~ Par1 + Par2, data = df,
                      fct = "GRIFFING1")
mod1m

# GRIFFING MODEL 2 ##########################
# No reciprocals
dfNoRec <- df %>%
  filter(as.character(Par1) <= as.character(Par2)) #Removing reciprocals

modh <- mmer(Ftime ~ 1,
             random = ~ overlay(Par1, Par2),
             data=dfNoRec)
summary(modh)$varcomp #Sembra più stabile

mod1m <- mmer.diallel(Ftime ~ Par1 + Par2, data = dfNoRec,
                      fct = "GRIFFING2")
mod1m

# GRIFFING MODEL 3 ##########################
# No selfs
dfNoSel <- df %>%
  filter(as.character(Par1) != as.character(Par2)) # removing selfs

modh <- mmer(Ftime ~ 1,
             random = ~ overlay(Par1, Par2) +
               combination,
             data=dfNoSel)
summary(modh)$varcomp #Sembra più stabile

mod1m <- mmer.diallel(Ftime ~ Par1 + Par2, data = dfNoSel,
                      fct = "GRIFFING3")
mod1m

# GRIFFING MODEL 4 ##########################
# No selfs
dfNone <- df %>%
  filter(as.character(Par1) != as.character(Par2)) %>% # removing selfs
  filter(as.character(Par1) <= as.character(Par2))

modh <- mmer(Ftime ~ 1,
             random = ~ overlay(Par1, Par2),
             data=dfNone)
summary(modh)$varcomp #Sembra più stabile

mod1m <- mmer.diallel(Ftime ~ Par1 + Par2, data = dfNone,
                      fct = "GRIFFING4")
mod1m

# GE2r ###########################################
# with reciprocals
# mod1r <- asreml(Ftime ~ Block + crosses, data=df,
#                 random = ~ Par1 + and(Par2) # GCA
#                 + Par1:crosses + and(Par2:crosses)  #h.i
#                 + combination:crosses) # SCA
# summary(mod1r)$varcomp

mod1h <- mmer(Ftime ~ crosses, data=df,
              random = ~ overlay(Par1, Par2) # GCA
              + overlay(Par1, Par2):crosses  #h.i
              + combination:crosses)
summary(mod1h)$varcomp

mod1m <- mmer.diallel(Ftime ~ Par1 + Par2,
                      data = df, fct = "GE2r")
mod1m

# GE2 ###########################################
# without reciprocals
# mod1r <- asreml(Ftime ~ Block + crosses, data=df,
#                 random = ~ Par1 + and(Par2) # GCA
#                 + Par1:crosses + and(Par2:crosses)  #h.i
#                 + combination:crosses) # SCA
# summary(mod1r)$varcomp

mod1h <- mmer(Ftime ~ crosses, data=df,
              random = ~ overlay(Par1, Par2) # GCA
              + overlay(Par1, Par2):crosses)
summary(mod1h)$varcomp

mod1m <- mmer.diallel(Ftime ~ Par1 + Par2,
                      data = df, fct = "GE2")
mod1m


# GE3r ################################
# with reciprocals. In Mohring: mixed, model 3 reduced
# mod1r <- asreml(Ftime ~ Block + crosses, data=df,
#                 random = ~ Par1:crosses + and(Par2:crosses) #GCA solo crosses
#                 + Par1:selfs # Selfed parents
#                 + combination:crosses # SCA solo crosses
#                 + combination:dr) # reciprocals
# summary(mod1r)$varcomp
mod1h <- mmer(Ftime ~ crosses, data=df,
                random = ~ overlay(Par1, Par2):crosses #GCA, only crosses
                + Par1:selfs # Selfed parents
                + combination:crosses) # reciprocals
summary(mod1h)$varcomp

mod1m <- mmer.diallel(Ftime ~ Par1 + Par2,
                      data = df, fct = "GE3r")
mod1m

# GE3 ############################################
# no reciprocals
# mod1r <- asreml(Ftime ~ Block + crosses, data=df,
#                 random = ~ Par1:crosses + and(Par2:crosses) #GCA solo crosses
#                 + Par1:selfs # Selfed parents
#                 + combination:crosses) # SCA solo crosses
# summary(mod1r)$varcomp

mod1h <- mmer(Ftime ~ crosses, data=df,
              random = ~ overlay(Par1, Par2):crosses #GCA solo crosses
              + Par1:selfs) # SCA solo crosses
summary(mod1h)$varcomp

mod1m <- mmer.diallel(Ftime ~ Par1 + Par2,
                      data = df, fct = "GE3")
mod1m

# Example 2 Lonnquist #####################
# Reproducing the results in Mohring, 2011
# Method = 2
rm(list=ls())
library(lmDiallel)
data(lonnquist61)
df <- lonnquist61

# Dummy variables for selfs, crosses, combinations
df$crosses <- ifelse(df$Par1 == df$Par2, 0, 1)
df$selfs <- ifelse(df$Par1 == df$Par2, 1, 0)
df$dr <- ifelse(as.character(df$Par1) < as.character(df$Par2), -1,
                ifelse(as.character(df$Par1) == as.character(df$Par2), 0, 1))
df$combination <- factor( ifelse(as.character(df$Par1) <= as.character(df$Par2),
                                 paste(df$Par1, df$Par2, sep =""),
                                 paste(df$Par2, df$Par1, sep ="")) )


# Griffing MODEL 2 ###############################
# Solo GCA e SCAt - no reciproci
# mod1r <- asreml(Yield ~ 1, data=df,
#                random = ~ Par1 + and(Par2))
# summary(mod1r)$varcomp

mod1h <- mmer(Yield ~ 1, data=df,
                random = ~ overlay(Par1, Par2))
summary(mod1h)$varcomp

mod1m <- mmer.diallel(Yield ~ Par1 + Par2,
                      data = df, fct = "GRIFFING2")
mod1m

# GE 2 ###########################################
# modr <- asreml(Yield ~ crosses, data = df,
#                random = ~ Par1 + and(Par2) #GCA
#                + Par1:crosses + and(Par2:crosses) ) #h.i
# summary(modr)$varcomp

modh <- mmer(Yield ~ crosses, data = df,
               random = ~ overlay(Par1, Par2) #GCA
               + overlay(Par1, Par2):crosses) #h.i
summary(modh)$varcomp

mod1m <- mmer.diallel(Yield ~ Par1 + Par2,
                      data = df, fct = "GE2")
mod1m


# GE3 ###############################################
# GCA + SCA scomposta (senza hi) + selfed parents
# mod3r <- asreml(Yield ~ crosses,
#                 random = ~ Par1:crosses + and(Par2:crosses)
#                 + Par1:selfs,
#                data=df)
# summary(mod3r)$varcomp

mod1m <- mmer.diallel(Yield ~ Par1 + Par2,
                      data = df, fct = "GE3")
mod1m


# Example of mating method 4 ###########################
# Balanced design ######################################
rm(list =ls())
data(Griffing56)
df <- Griffing56
head(df)

# GRIFFING 2 ###############
# asreml-r not OK!!!!!!! the overlay is not good
# mod1r <- asreml(Yield ~ 1,
#                 random = ~Par1 + and(Par2), data=df)
# summary(mod1r)$varcomp

mod1h <- mmer(Yield ~ 1,
                random = ~overlay(Par1, Par2), data=df)
summary(mod1h)$varcomp
mod1m <- mmer.diallel(Yield ~ Par1 + Par2,
                      data = df, fct = "GRIFFING4")
mod1m

# Repeated experiment ###########################################
data("diallelMET")
contrasts(diallelMET$Block) <- c("contr.sum")
contrasts(diallelMET$Env) <- c("contr.sum")

# GE2 ###########################################
# no reciprocals
mod1r <- asreml(Ftime ~ Block + crosses, data=df,
                random = ~ Par1 + and(Par2) # GCA
                + Par1:crosses + and(Par2:crosses)  #h.i
                + combination:crosses) # SCA
summary(mod1r)$varcomp

mod1h <- mmer(Ftime ~ Block + crosses, data=df,
              random = ~ overlay(Par1, Par2) # GCA
              + overlay(Par1, Par2):crosses  #h.i
              + combination:crosses) # SCA
summary(mod1h)$varcomp
