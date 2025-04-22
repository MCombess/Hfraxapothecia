

rm(list = ls())

library(dplyr)
library(lmerTest)
library(lme4)
library(MuMIn)
library(qpcR)
library(emmeans)

#Reading in the data with the mean Temperature and Humidity over different periods of time prior to sampling, the Gap Fraction (VisSky), vegetation cover,
#the Diameter, radius (r), and area of fruit bodies on a rachis sampled from a given site quadrat on a given date in um^2
#the sum of these areas per quadrat per location at each sampling point is then presented as either mm (sumareaqmm) or um (sumareaqum)

df <- read.csv("fruitbodydataformodelling.csv")
#construct maximum models, reduce models to only significant effects, rank models

#want sqrt of area so the dependent variable (sumareaq=the total apothecia area measured at a sampling location on a given date)
#is linear to meet model assumption. Then log this to meet model assumption that residuals have constant variance and they are normally distributed
df <- df %>% mutate(sumareaqsqrt=sqrt(sumareaqum))
df <- df %>% mutate(lsumareaqsqrt=log10(sumareaqsqrt))
colnames(df)
#drop redundant column, and classify factor variables as factors
df <- df[,-1]
df$Site <- as.factor(df$Site)
df$Quadrat <- as.factor(df$Quadrat)
df$date <- as.Date(df$date, format = "%d/%m/%Y")
df$Vegetation <- as.factor(df$Vegetation)

#remove duplicated data so have one observation per sampling point of log10(sqrt(total apothecia area))
colnames(df)
df2 <- df[,-c(18,19,23,24,25,26,27)]
df2 <- df2[!duplicated(df2),]
df <- df2
rm(df)
#This gives 399 observations to work with#

####################################################################################################

#Construct maximum models for those with ground cover

RHlmers1 <- lmer(lsumareaqsqrt~scale(meanRH7)*scale(meanT7)*Vegetation+(1|Site/Quadrat) + (1|Year), data=df2)

RHlmers2 <- lmer(lsumareaqsqrt~scale(meanRH7)*scale(meanT14)*Vegetation+(1|Site/Quadrat) + (1|Year), data=df2)

RHlmers3 <- lmer(lsumareaqsqrt~scale(meanRH7)*scale(meanT21)*Vegetation+(1|Site/Quadrat) + (1|Year), data=df2)

RHlmers4 <- lmer(lsumareaqsqrt~scale(meanRH7)*scale(meanT28)*Vegetation+(1|Site/Quadrat) + (1|Year), data=df2)

RHlmers5 <- lmer(lsumareaqsqrt~scale(meanRH7)*scale(meanT35)*Vegetation+(1|Site/Quadrat) + (1|Year), data=df2)

RHlmers6 <- lmer(lsumareaqsqrt~scale(meanRH7)*scale(meanT42)*Vegetation+(1|Site/Quadrat) + (1|Year), data=df2)



RHlmers7 <- lmer(lsumareaqsqrt~scale(meanRH14)*scale(meanT7)*Vegetation+(1|Site/Quadrat) + (1|Year), data=df2)

RHlmers8 <- lmer(lsumareaqsqrt~scale(meanRH14)*scale(meanT14)*Vegetation+(1|Site/Quadrat) + (1|Year), data=df2)

RHlmers9 <- lmer(lsumareaqsqrt~scale(meanRH14)*scale(meanT21)*Vegetation+(1|Site/Quadrat) + (1|Year), data=df2)

RHlmers10 <- lmer(lsumareaqsqrt~scale(meanRH14)*scale(meanT28)*Vegetation+(1|Site/Quadrat) + (1|Year), data=df2)

RHlmers11 <- lmer(lsumareaqsqrt~scale(meanRH14)*scale(meanT35)*Vegetation+(1|Site/Quadrat) + (1|Year), data=df2)

RHlmers12 <- lmer(lsumareaqsqrt~scale(meanRH14)*scale(meanT42)*Vegetation+(1|Site/Quadrat) + (1|Year), data=df2)


RHlmers13 <- lmer(lsumareaqsqrt~scale(meanRH21)*scale(meanT7)*Vegetation+(1|Site/Quadrat) + (1|Year), data=df2)

RHlmers14 <- lmer(lsumareaqsqrt~scale(meanRH21)*scale(meanT14)*Vegetation+(1|Site/Quadrat) + (1|Year), data=df2)

RHlmers15 <- lmer(lsumareaqsqrt~scale(meanRH21)*scale(meanT21)*Vegetation+(1|Site/Quadrat) + (1|Year), data=df2)

RHlmers16 <- lmer(lsumareaqsqrt~scale(meanRH21)*scale(meanT28)*Vegetation+(1|Site/Quadrat) + (1|Year), data=df2)

RHlmers17 <- lmer(lsumareaqsqrt~scale(meanRH21)*scale(meanT35)*Vegetation+(1|Site/Quadrat) + (1|Year), data=df2)

RHlmers18 <- lmer(lsumareaqsqrt~scale(meanRH21)*scale(meanT42)*Vegetation+(1|Site/Quadrat) + (1|Year), data=df2)


RHlmers19 <- lmer(lsumareaqsqrt~scale(meanRH28)*scale(meanT7)*Vegetation+(1|Site/Quadrat) + (1|Year), data=df2)

RHlmers20 <- lmer(lsumareaqsqrt~scale(meanRH28)*scale(meanT14)*Vegetation+(1|Site/Quadrat) + (1|Year), data=df2)

RHlmers21 <- lmer(lsumareaqsqrt~scale(meanRH28)*scale(meanT21)*Vegetation+(1|Site/Quadrat) + (1|Year), data=df2)

RHlmers22 <- lmer(lsumareaqsqrt~scale(meanRH28)*scale(meanT28)*Vegetation+(1|Site/Quadrat) + (1|Year), data=df2)

RHlmers23 <- lmer(lsumareaqsqrt~scale(meanRH28)*scale(meanT35)*Vegetation+(1|Site/Quadrat) + (1|Year), data=df2)

RHlmers24 <- lmer(lsumareaqsqrt~scale(meanRH28)*scale(meanT42)*Vegetation+(1|Site/Quadrat) + (1|Year), data=df2)


RHlmers25 <- lmer(lsumareaqsqrt~scale(meanRH35)*scale(meanT7)*Vegetation+(1|Site/Quadrat) + (1|Year), data=df2)

RHlmers26 <- lmer(lsumareaqsqrt~scale(meanRH35)*scale(meanT14)*Vegetation+(1|Site/Quadrat) + (1|Year), data=df2)

RHlmers27 <- lmer(lsumareaqsqrt~scale(meanRH35)*scale(meanT21)*Vegetation+(1|Site/Quadrat) + (1|Year), data=df2)

RHlmers28 <- lmer(lsumareaqsqrt~scale(meanRH35)*scale(meanT28)*Vegetation+(1|Site/Quadrat) + (1|Year), data=df2)

RHlmers29 <- lmer(lsumareaqsqrt~scale(meanRH35)*scale(meanT35)*Vegetation+(1|Site/Quadrat) + (1|Year), data=df2)

RHlmers30 <- lmer(lsumareaqsqrt~scale(meanRH35)*scale(meanT42)*Vegetation+(1|Site/Quadrat) + (1|Year), data=df2)


RHlmers31 <- lmer(lsumareaqsqrt~scale(meanRH42)*scale(meanT7)*Vegetation+(1|Site/Quadrat) + (1|Year), data=df2)

RHlmers32 <- lmer(lsumareaqsqrt~scale(meanRH42)*scale(meanT14)*Vegetation+(1|Site/Quadrat) + (1|Year), data=df2)

RHlmers33 <- lmer(lsumareaqsqrt~scale(meanRH42)*scale(meanT21)*Vegetation+(1|Site/Quadrat) + (1|Year), data=df2)

RHlmers34 <- lmer(lsumareaqsqrt~scale(meanRH42)*scale(meanT28)*Vegetation+(1|Site/Quadrat) + (1|Year), data=df2)

RHlmers35 <- lmer(lsumareaqsqrt~scale(meanRH42)*scale(meanT35)*Vegetation+(1|Site/Quadrat) + (1|Year), data=df2)

RHlmers36 <- lmer(lsumareaqsqrt~scale(meanRH42)*scale(meanT42)*Vegetation+(1|Site/Quadrat) + (1|Year), data=df2)


#stepwise reduction of non-significant effects from maximum models

RHs1 <-  step(RHlmers1,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs2 <-  step(RHlmers2,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs3 <-  step(RHlmers3,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs4 <-  step(RHlmers4,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs5 <-  step(RHlmers5,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs6 <-  step(RHlmers6,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs7 <-  step(RHlmers7,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs8 <-  step(RHlmers8,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs9 <-  step(RHlmers9,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs10 <-  step(RHlmers10,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs11 <-  step(RHlmers11,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs12 <-  step(RHlmers12,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs13 <-  step(RHlmers13,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs14 <-  step(RHlmers14,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs15 <-  step(RHlmers15,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs16 <-  step(RHlmers16,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs17 <-  step(RHlmers17,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs18 <-  step(RHlmers18,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs19 <-  step(RHlmers19,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs20 <-  step(RHlmers20,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs21 <-  step(RHlmers21,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs22 <-  step(RHlmers22,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs23 <-  step(RHlmers23,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs24 <-  step(RHlmers24,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs25 <-  step(RHlmers25,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs26 <-  step(RHlmers26,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs27 <-  step(RHlmers27,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs28 <-  step(RHlmers28,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs29 <-  step(RHlmers29,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs30 <-  step(RHlmers30,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs31 <-  step(RHlmers31,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs32 <-  step(RHlmers32,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs33 <-  step(RHlmers33,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs34 <-  step(RHlmers34,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs35 <-  step(RHlmers35,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs36 <-  step(RHlmers36,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)

#AICc values of the models to then perform model selection
RHs1
RHflmers1 <- AICc(lmer(lsumareaqsqrt ~ scale(meanRH7) + scale(meanT7) + Vegetation + (1 | Site/Quadrat) + (1 | Year) + scale(meanRH7):scale(meanT7) + scale(meanT7):Vegetation, data=df2, REML=FALSE))
RHs2
RHflmers2 <- AICc(lmer(lsumareaqsqrt ~ scale(meanRH7) + scale(meanT14) + Vegetation + (1 | Site/Quadrat) + (1 | Year) + scale(meanRH7):scale(meanT14) + scale(meanT14):Vegetation, data=df2, REML=FALSE))
RHs3
RHflmers3 <- AICc(lmer(lsumareaqsqrt ~ scale(meanRH7) + scale(meanT21) + Vegetation + (1 | Site/Quadrat) + (1 | Year) + scale(meanRH7):scale(meanT21) + scale(meanT21):Vegetation, data=df2, REML=FALSE))
RHs4
RHflmers4 <- AICc(lmer(lsumareaqsqrt ~ scale(meanRH7) + scale(meanT28) + Vegetation + (1 | Site/Quadrat) + (1 | Year) + scale(meanRH7):scale(meanT28) + scale(meanT28):Vegetation, data=df2, REML=FALSE))
RHs5
RHflmers5 <- AICc(lmer(lsumareaqsqrt ~ scale(meanT35) + (1 | Site/Quadrat) + (1 | Year), data=df2, REML=FALSE))
RHs6
RHflmers6 <- AICc(lmer(lsumareaqsqrt ~ scale(meanT42) + (1 | Site/Quadrat) + (1 | Year), data=df2, REML=FALSE))
RHs7
RHflmers7 <- AICc(lmer(lsumareaqsqrt ~ scale(meanRH14) + scale(meanT7) + Vegetation + (1 | Site/Quadrat) + (1 | Year) + scale(meanRH14):scale(meanT7) + scale(meanT7):Vegetation, data=df2, REML=FALSE))
RHs8
RHflmers8 <- AICc(lmer(lsumareaqsqrt ~ scale(meanRH14) + scale(meanT14) + Vegetation + (1 | Site/Quadrat) + (1 | Year) + scale(meanRH14):scale(meanT14) + scale(meanT14):Vegetation, data=df2, REML=FALSE))
RHs9
RHflmers9 <- AICc(lmer(lsumareaqsqrt ~ scale(meanT21) + Vegetation + (1 | Site/Quadrat) + (1 | Year) + scale(meanT21):Vegetation, data=df2, REML=FALSE))

#wont't converge without REML=FALSE - use bobyqa optimiser
RHs10
#RHflmers10 <- AICc(lmer(lsumareaqsqrt ~ scale(meanT28) + (1 | Site/Quadrat) + (1 | Year), data=df2, REML=FALSE))
RHflmers10 <- AICc(lmer(lsumareaqsqrt ~ scale(meanT28) + (1 | Site/Quadrat) + (1 | Year), lmerControl(optimizer = "bobyqa"), data=df2, REML=FALSE))


RHs11
RHflmers11 <- AICc(lmer(lsumareaqsqrt ~ scale(meanT35) + (1 | Site/Quadrat) + (1 | Year), data=df2, REML=FALSE))
RHs12
RHflmers12 <- AICc(lmer(lsumareaqsqrt ~ scale(meanT42) + (1 | Site/Quadrat) + (1 | Year), data=df2, REML=FALSE))
RHs13
RHflmers13 <- AICc(lmer(lsumareaqsqrt ~ scale(meanRH21) + scale(meanT7) + Vegetation + (1 | Site/Quadrat) + (1 | Year) + scale(meanRH21):scale(meanT7) + scale(meanT7):Vegetation, data=df2, REML=FALSE))
RHs14
RHflmers14 <- AICc(lmer(lsumareaqsqrt ~ scale(meanRH21) + scale(meanT14) + Vegetation + (1 | Site/Quadrat) + (1 | Year) + scale(meanRH21):scale(meanT14) + scale(meanT14):Vegetation, data=df2, REML=FALSE))
RHs15
RHflmers15 <- AICc(lmer(lsumareaqsqrt ~ scale(meanRH21) + scale(meanT21) + Vegetation + (1 | Site/Quadrat) + (1 | Year) + scale(meanRH21):Vegetation, data=df2, REML=FALSE))


#wont't converge without REML=FALSE - use bobyqa optimiser
RHs16
#RHflmers16 <- AICc(lmer(lsumareaqsqrt ~ scale(meanT28) + (1 | Site/Quadrat) + (1 | Year), data=df2, REML=FALSE))
RHflmers16 <- AICc(lmer(lsumareaqsqrt ~ scale(meanT28) + (1 | Site/Quadrat) + (1 | Year), lmerControl(optimizer = "bobyqa"), data=df2, REML=FALSE))


RHs17
RHflmers17 <- AICc(lmer(lsumareaqsqrt ~ scale(meanT35) + (1 | Site/Quadrat) + (1 | Year), data=df2, REML=FALSE))
RHs18
RHflmers18 <- AICc(lmer(lsumareaqsqrt ~ scale(meanT42) + (1 | Site/Quadrat) + (1 | Year), data=df2, REML=FALSE))
RHs19
RHflmers19 <- AICc(lmer(lsumareaqsqrt ~ scale(meanRH28) + scale(meanT7) + Vegetation + (1 | Site/Quadrat) + (1 | Year) + scale(meanRH28):scale(meanT7) + scale(meanT7):Vegetation, data=df2, REML=FALSE))
RHs20
RHflmers20 <- AICc(lmer(lsumareaqsqrt ~ scale(meanRH28) + scale(meanT14) + Vegetation + (1 | Site/Quadrat) + (1 | Year) + scale(meanRH28):scale(meanT14) + scale(meanT14):Vegetation, data=df2, REML=FALSE))
RHs21
RHflmers21 <- AICc(lmer(lsumareaqsqrt ~ scale(meanRH28) + scale(meanT21) + Vegetation + (1 | Site/Quadrat) + (1 | Year) + scale(meanRH28):Vegetation, data=df2, REML=FALSE))
RHs22
RHflmers22 <- AICc(lmer(lsumareaqsqrt ~ scale(meanRH28) + scale(meanT28) + Vegetation + (1 | Site/Quadrat) + (1 | Year) + scale(meanRH28):Vegetation, data=df2, REML=FALSE))
RHs23
RHflmers23 <- AICc(lmer(lsumareaqsqrt ~ scale(meanT35) + (1 | Site/Quadrat) + (1 | Year), data=df2, REML=FALSE))
RHs24
RHflmers24 <- AICc(lmer(lsumareaqsqrt ~ scale(meanT42) + (1 | Site/Quadrat) + (1 | Year), data=df2, REML=FALSE))
RHs25
RHflmers25 <- AICc(lmer(lsumareaqsqrt ~ scale(meanRH35) + scale(meanT7) + Vegetation + (1 | Site/Quadrat) + (1 | Year) + scale(meanRH35):scale(meanT7) + scale(meanT7):Vegetation, data=df2, REML=FALSE))
RHs26
RHflmers26 <- AICc(lmer(lsumareaqsqrt ~ scale(meanRH35) * scale(meanT14) * Vegetation + (1 | Site/Quadrat) + (1 | Year), data=df2, REML=FALSE))
RHs27
RHflmers27 <- AICc(lmer(lsumareaqsqrt ~ scale(meanRH35) * scale(meanT21) * Vegetation + (1 | Site/Quadrat) + (1 | Year), data=df2, REML=FALSE))

#wont't converge without REML=FALSE - use bobyqa optimiser
RHs28
#RHflmers28 <- AICc(lmer(lsumareaqsqrt ~ scale(meanT28) + (1 | Site/Quadrat) + (1 | Year), data=df2, REML=FALSE))
RHflmers28 <- AICc(lmer(lsumareaqsqrt ~ scale(meanT28) + (1 | Site/Quadrat) + (1 | Year), lmerControl(optimizer = "bobyqa"), data=df2, REML=FALSE))


RHs29
RHflmers29 <- AICc(lmer(lsumareaqsqrt ~ scale(meanT35) + (1 | Site/Quadrat) + (1 | Year), data=df2, REML=FALSE))
RHs30
RHflmers30 <- AICc(lmer(lsumareaqsqrt ~ scale(meanRH35) + scale(meanT42) + (1 | Site/Quadrat) + (1 | Year), data=df2, REML=FALSE))
RHs31
RHflmers31 <- AICc(lmer(lsumareaqsqrt ~ scale(meanRH42) + scale(meanT7) + Vegetation + (1 | Site/Quadrat) + (1 | Year) + scale(meanRH42):scale(meanT7) + scale(meanT7):Vegetation, data=df2, REML=FALSE))
RHs32
RHflmers32 <- AICc(lmer(lsumareaqsqrt ~ scale(meanRH42) * scale(meanT14) * Vegetation + (1 | Site/Quadrat) + (1 | Year), data=df2, REML=FALSE))
RHs33
RHflmers33 <- AICc(lmer(lsumareaqsqrt ~ scale(meanRH42) * scale(meanT21) * Vegetation + (1 | Site/Quadrat) + (1 | Year), data=df2, REML=FALSE))
RHs34
RHflmers34 <- AICc(lmer(lsumareaqsqrt ~ scale(meanRH42) + scale(meanT28) + (1 | Site/Quadrat) + (1 | Year) + scale(meanRH42):scale(meanT28), data=df2, REML=FALSE))
RHs35
RHflmers35 <- AICc(lmer(lsumareaqsqrt ~ scale(meanRH42) + scale(meanT35) + (1 | Site/Quadrat) + (1 | Year), data=df2, REML=FALSE))
RHs36
RHflmers36 <- AICc(lmer(lsumareaqsqrt ~ scale(meanRH42) + scale(meanT42) + (1 | Site/Quadrat) + (1 | Year), data=df2, REML=FALSE))


############################################################################################################################################

#Construct maxmimal Gap fraction (VisSky) models - could not fit Site/Quadrat as a random effect so do Site:quadrat,
#this defines that Site does not have a random effect, but each Quadrat Site combination does#

RHlmers1Vis <- lmer(lsumareaqsqrt~Site+scale(meanRH7)*scale(meanT7)*scale(VisSky)+(1|Site:Quadrat)+(1|Year), data=df2)

RHlmers2Vis <- lmer(lsumareaqsqrt~Site+scale(meanRH7)*scale(meanT14)*scale(VisSky)+(1|Site:Quadrat) + (1|Year), data=df2)

RHlmers3Vis <- lmer(lsumareaqsqrt~Site+scale(meanRH7)*scale(meanT21)*scale(VisSky)+(1|Site:Quadrat) + (1|Year), data=df2)

RHlmers4Vis <- lmer(lsumareaqsqrt~Site+scale(meanRH7)*scale(meanT28)*scale(VisSky)+(1|Site:Quadrat) + (1|Year), data=df2)

RHlmers5Vis <- lmer(lsumareaqsqrt~Site+scale(meanRH7)*scale(meanT35)*scale(VisSky)+(1|Site:Quadrat) + (1|Year), data=df2)

RHlmers6Vis <- lmer(lsumareaqsqrt~Site+scale(meanRH7)*scale(meanT42)*scale(VisSky)+(1|Site:Quadrat) + (1|Year), data=df2)



RHlmers7Vis <- lmer(lsumareaqsqrt~Site+scale(meanRH14)*scale(meanT7)*scale(VisSky)+(1|Site:Quadrat) + (1|Year), data=df2)

RHlmers8Vis <- lmer(lsumareaqsqrt~Site+scale(meanRH14)*scale(meanT14)*scale(VisSky)+(1|Site:Quadrat) + (1|Year), control = lmerControl("bobyqa"), data=df2)

RHlmers9Vis <- lmer(lsumareaqsqrt~Site+scale(meanRH14)*scale(meanT21)*scale(VisSky)+(1|Site:Quadrat) + (1|Year), data=df2)

RHlmers10Vis <- lmer(lsumareaqsqrt~Site+scale(meanRH14)*scale(meanT28)*scale(VisSky)+(1|Site:Quadrat) + (1|Year), control=lmerControl("bobyqa"), data=df2)

RHlmers11Vis <- lmer(lsumareaqsqrt~Site+scale(meanRH14)*scale(meanT35)*scale(VisSky)+(1|Site:Quadrat) + (1|Year), control=lmerControl("bobyqa"), data=df2)

RHlmers12Vis <- lmer(lsumareaqsqrt~Site+scale(meanRH14)*scale(meanT42)*scale(VisSky)+(1|Site:Quadrat) + (1|Year), data=df2)


RHlmers13Vis <- lmer(lsumareaqsqrt~Site+scale(meanRH21)*scale(meanT7)*scale(VisSky)+(1|Site:Quadrat) + (1|Year), data=df2)

RHlmers14Vis <- lmer(lsumareaqsqrt~Site+scale(meanRH21)*scale(meanT14)*scale(VisSky)+(1|Site:Quadrat) + (1|Year), data=df2)

RHlmers15Vis <- lmer(lsumareaqsqrt~Site+scale(meanRH21)*scale(meanT21)*scale(VisSky)+(1|Site:Quadrat) + (1|Year), data=df2)

RHlmers16Vis <- lmer(lsumareaqsqrt~Site+scale(meanRH21)*scale(meanT28)*scale(VisSky)+(1|Site:Quadrat) + (1|Year), data=df2)

RHlmers17Vis <- lmer(lsumareaqsqrt~Site+scale(meanRH21)*scale(meanT35)*scale(VisSky)+(1|Site:Quadrat) + (1|Year), data=df2)

RHlmers18Vis <- lmer(lsumareaqsqrt~Site+scale(meanRH21)*scale(meanT42)*scale(VisSky)+(1|Site:Quadrat) + (1|Year), data=df2)


RHlmers19Vis <- lmer(lsumareaqsqrt~Site+scale(meanRH28)*scale(meanT7)*scale(VisSky)+(1|Site:Quadrat) + (1|Year), data=df2)

RHlmers20Vis <- lmer(lsumareaqsqrt~Site+scale(meanRH28)*scale(meanT14)*scale(VisSky)+(1|Site:Quadrat) + (1|Year), data=df2)

RHlmers21Vis <- lmer(lsumareaqsqrt~Site+scale(meanRH28)*scale(meanT21)*scale(VisSky)+(1|Site:Quadrat) + (1|Year), data=df2)

RHlmers22Vis <- lmer(lsumareaqsqrt~Site+scale(meanRH28)*scale(meanT28)*scale(VisSky)+(1|Site:Quadrat) + (1|Year), data=df2)

RHlmers23Vis <- lmer(lsumareaqsqrt~Site+scale(meanRH28)*scale(meanT35)*scale(VisSky)+(1|Site:Quadrat) + (1|Year), data=df2)

RHlmers24Vis <- lmer(lsumareaqsqrt~Site+scale(meanRH28)*scale(meanT42)*scale(VisSky)+(1|Site:Quadrat) + (1|Year), data=df2)


RHlmers25Vis <- lmer(lsumareaqsqrt~Site+scale(meanRH35)*scale(meanT7)*scale(VisSky)+(1|Site:Quadrat) + (1|Year), data=df2)

RHlmers26Vis <- lmer(lsumareaqsqrt~Site+scale(meanRH35)*scale(meanT14)*scale(VisSky)+(1|Site:Quadrat) + (1|Year), data=df2)

RHlmers27Vis <- lmer(lsumareaqsqrt~Site+scale(meanRH35)*scale(meanT21)*scale(VisSky)+(1|Site:Quadrat) + (1|Year), control=lmerControl("bobyqa"), data=df2)

RHlmers28Vis <- lmer(lsumareaqsqrt~Site+scale(meanRH35)*scale(meanT28)*scale(VisSky)+(1|Site:Quadrat) + (1|Year), control=lmerControl("bobyqa"), data=df2)

RHlmers29Vis <- lmer(lsumareaqsqrt~Site+scale(meanRH35)*scale(meanT35)*scale(VisSky)+(1|Site:Quadrat) + (1|Year), control=lmerControl("bobyqa"), data=df2)

RHlmers30Vis <- lmer(lsumareaqsqrt~Site+scale(meanRH35)*scale(meanT42)*scale(VisSky)+(1|Site:Quadrat) + (1|Year), data=df2)


RHlmers31Vis <- lmer(lsumareaqsqrt~Site+scale(meanRH42)*scale(meanT7)*scale(VisSky)+(1|Site:Quadrat) + (1|Year), data=df2)

RHlmers32Vis <- lmer(lsumareaqsqrt~Site+scale(meanRH42)*scale(meanT14)*scale(VisSky)+(1|Site:Quadrat) + (1|Year), data=df2)

RHlmers33Vis <- lmer(lsumareaqsqrt~Site+scale(meanRH42)*scale(meanT21)*scale(VisSky)+(1|Site:Quadrat) + (1|Year), data=df2)

RHlmers34Vis <- lmer(lsumareaqsqrt~Site+scale(meanRH42)*scale(meanT28)*scale(VisSky)+(1|Site:Quadrat) + (1|Year), data=df2)

RHlmers35Vis <- lmer(lsumareaqsqrt~Site+scale(meanRH42)*scale(meanT35)*scale(VisSky)+(1|Site:Quadrat) + (1|Year), data=df2)

RHlmers36Vis <- lmer(lsumareaqsqrt~Site+scale(meanRH42)*scale(meanT42)*scale(VisSky)+(1|Site:Quadrat) + (1|Year), data=df2)

#stepwise reduction of non-significant effects
RHs1Vis <-  step(RHlmers1Vis,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs2Vis <-  step(RHlmers2Vis,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs3Vis <-  step(RHlmers3Vis,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs4Vis <-  step(RHlmers4Vis,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs5Vis <-  step(RHlmers5Vis,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs6Vis <-  step(RHlmers6Vis,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs7Vis <-  step(RHlmers7Vis,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs8Vis <-  step(RHlmers8Vis,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs9Vis <-  step(RHlmers9Vis,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs10Vis <-  step(RHlmers10Vis,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs11Vis <-  step(RHlmers11Vis,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs12Vis <-  step(RHlmers12Vis,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs13Vis <-  step(RHlmers13Vis,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs14Vis <-  step(RHlmers14Vis,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs15Vis <-  step(RHlmers15Vis,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs16Vis <-  step(RHlmers16Vis,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs17Vis <-  step(RHlmers17Vis,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs18Vis <-  step(RHlmers18Vis,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs19Vis <-  step(RHlmers19Vis,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs20Vis <-  step(RHlmers20Vis,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs21Vis <-  step(RHlmers21Vis,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs22Vis <-  step(RHlmers22Vis,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs23Vis <-  step(RHlmers23Vis,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs24Vis <-  step(RHlmers24Vis,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs25Vis <-  step(RHlmers25Vis,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs26Vis <-  step(RHlmers26Vis,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs27Vis <-  step(RHlmers27Vis,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs28Vis <-  step(RHlmers28Vis,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs29Vis <-  step(RHlmers29Vis,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs30Vis <-  step(RHlmers30Vis,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs31Vis <-  step(RHlmers31Vis,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs32Vis <-  step(RHlmers32Vis,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs33Vis <-  step(RHlmers33Vis,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs34Vis <-  step(RHlmers34Vis,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs35Vis <-  step(RHlmers35Vis,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)
RHs36Vis <-  step(RHlmers36Vis,  ddf = "Kenward-Roger",alpha.fixed=0.01,reduce.random=FALSE)

#AICc values of the models to then perform model selection
RHs1Vis
RHflmers1Vis <- AICc(lmer(lsumareaqsqrt ~ scale(meanRH7) + scale(meanT7) + scale(VisSky) + (1 | Site:Quadrat) + (1 | Year) + scale(meanRH7):scale(meanT7) + scale(meanT7):scale(VisSky), data=df2, REML=FALSE))
RHs2Vis
RHflmers2Vis <- AICc(lmer(lsumareaqsqrt ~ scale(meanRH7) + scale(meanT14) + scale(VisSky) + (1 | Site:Quadrat) + (1 | Year) + scale(meanRH7):scale(meanT14) + scale(meanT14):scale(VisSky), data=df2, REML=FALSE))
RHs3Vis
RHflmers3Vis <- AICc(lmer(lsumareaqsqrt ~ scale(meanRH7) + scale(meanT21) + scale(VisSky) + (1 | Site:Quadrat) + (1 | Year) + scale(meanRH7):scale(meanT21) + scale(meanT21):scale(VisSky), data=df2, REML=FALSE))
RHs4Vis
RHflmers4Vis <- AICc(lmer(lsumareaqsqrt ~ scale(meanRH7) + scale(meanT28) + scale(VisSky) + (1 | Site:Quadrat) + (1 | Year) + scale(meanRH7):scale(meanT28) + scale(meanT28):scale(VisSky), data=df2, REML=FALSE))
RHs5Vis
RHflmers5Vis <- AICc(lmer(lsumareaqsqrt ~ scale(meanT35) + scale(VisSky) + (1 | Site:Quadrat) + (1 | Year), data=df2, REML=FALSE))
RHs6Vis
RHflmers6Vis <- AICc(lmer(lsumareaqsqrt ~ scale(meanT42) + scale(VisSky) + (1 | Site:Quadrat) + (1 | Year), data=df2, REML=FALSE))
RHs7Vis
RHflmers7Vis <- AICc(lmer(lsumareaqsqrt ~ scale(meanRH14) + scale(meanT7) + scale(VisSky) + (1 | Site:Quadrat) + (1 | Year) + scale(meanRH14):scale(meanT7) + scale(meanT7):scale(VisSky), data=df2, REML=FALSE))
RHs8Vis
RHflmers8Vis <- AICc(lmer(lsumareaqsqrt ~ scale(meanRH14) + scale(meanT14) + scale(VisSky) + (1 | Site:Quadrat) + (1 | Year) + scale(meanRH14):scale(meanT14) + scale(meanT14):scale(VisSky), data=df2, REML=FALSE))
RHs9Vis
RHflmers9Vis <- AICc(lmer(lsumareaqsqrt ~ scale(meanT21) + scale(VisSky) + (1 | Site:Quadrat) + (1 | Year) + scale(meanT21):scale(VisSky), data=df2, REML=FALSE))
RHs10Vis
RHflmers10Vis <- AICc(lmer(lsumareaqsqrt ~ scale(meanT28) + scale(VisSky) + (1 | Site:Quadrat) + (1 | Year) + scale(meanT28):scale(VisSky), data=df2, REML=FALSE))
RHs11Vis
RHflmers11Vis <- AICc(lmer(lsumareaqsqrt ~ scale(meanT35) + scale(VisSky) + (1 | Site:Quadrat) + (1 | Year), data=df2, REML=FALSE))
RHs12Vis
RHflmers12Vis <- AICc(lmer(lsumareaqsqrt ~ scale(meanT42) + scale(VisSky) + (1 | Site:Quadrat) + (1 | Year), data=df2, REML=FALSE))
RHs13Vis
RHflmers13Vis <- AICc(lmer(lsumareaqsqrt ~ scale(meanRH21) + scale(meanT7) + scale(VisSky) + (1 | Site:Quadrat) + (1 | Year) + scale(meanRH21):scale(meanT7) + scale(meanT7):scale(VisSky), data=df2, REML=FALSE))
RHs14Vis
RHflmers14Vis <- AICc(lmer(lsumareaqsqrt ~ scale(meanRH21) + scale(meanT14) + scale(VisSky) + (1 | Site:Quadrat) + (1 | Year) + scale(meanRH21):scale(meanT14) + scale(meanT14):scale(VisSky), data=df2, REML=FALSE))
RHs15Vis
RHflmers15Vis <- AICc(lmer(lsumareaqsqrt ~ scale(meanT21) + scale(VisSky) + (1 | Site:Quadrat) + (1 | Year) + scale(meanT21):scale(VisSky), data=df2, REML=FALSE))
RHs16Vis
RHflmers16Vis <- AICc(lmer(lsumareaqsqrt ~ scale(meanT28) + scale(VisSky) + (1 | Site:Quadrat) + (1 | Year) + scale(meanT28):scale(VisSky), data=df2, REML=FALSE))
RHs17Vis
RHflmers17Vis <- AICc(lmer(lsumareaqsqrt ~ scale(meanT35) + scale(VisSky) + (1 | Site:Quadrat) + (1 | Year), data=df2, REML=FALSE))
RHs18Vis
RHflmers18Vis <- AICc(lmer(lsumareaqsqrt ~ scale(meanT42) + scale(VisSky) + (1 | Site:Quadrat) + (1 | Year), data=df2, REML=FALSE))
RHs19Vis
RHflmers19Vis <- AICc(lmer(lsumareaqsqrt ~ scale(meanRH28) + scale(meanT7) + scale(VisSky) + (1 | Site:Quadrat) + (1 | Year) + scale(meanRH28):scale(meanT7) + scale(meanT7):scale(VisSky), data=df2, REML=FALSE))
RHs20Vis
RHflmers20Vis <- AICc(lmer(lsumareaqsqrt ~ scale(meanRH28) + scale(meanT14) + scale(VisSky) + (1 | Site:Quadrat) + (1 | Year) + scale(meanRH28):scale(meanT14) + scale(meanT14):scale(VisSky), data=df2, REML=FALSE))
RHs21Vis
RHflmers21Vis <- AICc(lmer(lsumareaqsqrt ~ scale(meanT21) + scale(VisSky) + (1 | Site:Quadrat) + (1 | Year) + scale(meanT21):scale(VisSky), data=df2, REML=FALSE))
RHs22Vis
RHflmers22Vis <- AICc(lmer(lsumareaqsqrt ~ scale(meanT28) + scale(VisSky) + (1 | Site:Quadrat) + (1 | Year) + scale(meanT28):scale(VisSky), data=df2, REML=FALSE))
RHs23Vis
RHflmers23Vis <- AICc(lmer(lsumareaqsqrt ~ scale(meanT35) + scale(VisSky) + (1 | Site:Quadrat) + (1 | Year), data=df2, REML=FALSE))
RHs24Vis
RHflmers24Vis <- AICc(lmer(lsumareaqsqrt ~ scale(meanT42) + scale(VisSky) + (1 | Site:Quadrat) + (1 | Year), data=df2, REML=FALSE))
RHs25Vis
RHflmers25Vis <- AICc(lmer(lsumareaqsqrt ~ scale(meanRH35) + scale(meanT7) + scale(VisSky) + (1 | Site:Quadrat) + (1 | Year) + scale(meanRH35):scale(meanT7) + scale(meanT7):scale(VisSky), data=df2, REML=FALSE))
RHs26Vis
RHflmers26Vis <- AICc(lmer(lsumareaqsqrt ~ scale(meanRH35) + scale(meanT14) + scale(VisSky) + (1 | Site:Quadrat) + (1 | Year) + scale(meanRH35):scale(meanT14) + scale(meanT14):scale(VisSky), data=df2, REML=FALSE))
RHs27Vis
RHflmers27Vis <- AICc(lmer(lsumareaqsqrt ~ scale(meanRH35) + scale(meanT21) + scale(VisSky) + (1 | Site:Quadrat) + (1 | Year) + scale(meanRH35):scale(meanT21) + scale(meanT21):scale(VisSky), data=df2, REML=FALSE))
RHs28Vis
RHflmers28Vis <- AICc(lmer(lsumareaqsqrt ~ scale(meanT28) + scale(VisSky) + (1 | Site:Quadrat) + (1 | Year) + scale(meanT28):scale(VisSky), data=df2, REML=FALSE))
RHs29Vis
RHflmers29Vis <- AICc(lmer(lsumareaqsqrt ~ scale(meanT35) + scale(VisSky) + (1 | Site:Quadrat) + (1 | Year), data=df2, REML=FALSE))
RHs30Vis
RHflmers30Vis <- AICc(lmer(lsumareaqsqrt ~ scale(meanRH35) + scale(meanT42) + scale(VisSky) + (1 | Site:Quadrat) + (1 | Year), data=df2, REML=FALSE))
RHs31Vis
RHflmers31Vis <- AICc(lmer(lsumareaqsqrt ~ scale(meanRH42) + scale(meanT7) + scale(VisSky) + (1 | Site:Quadrat) + (1 | Year) + scale(meanRH42):scale(meanT7) + scale(meanT7):scale(VisSky), data=df2, REML=FALSE))
RHs32Vis
RHflmers32Vis <- AICc(lmer(lsumareaqsqrt ~ scale(meanRH42) + scale(meanT14) + scale(VisSky) + (1 | Site:Quadrat) + (1 | Year) + scale(meanRH42):scale(meanT14) + scale(meanT14):scale(VisSky), data=df2, REML=FALSE))
RHs33Vis
RHflmers33Vis <- AICc(lmer(lsumareaqsqrt ~ scale(meanRH42) + scale(meanT21) + scale(VisSky) + (1 | Site:Quadrat) + (1 | Year) + scale(meanRH42):scale(meanT21) + scale(meanT21):scale(VisSky), data=df2, REML=FALSE))
RHs34Vis
RHflmers34Vis <- AICc(lmer(lsumareaqsqrt ~ scale(meanRH42) + scale(meanT28) + scale(VisSky) + (1 | Site:Quadrat) + (1 | Year) + scale(meanRH42):scale(meanT28) + scale(meanT28):scale(VisSky), data=df2, REML=FALSE))
RHs35Vis
RHflmers35Vis <- AICc(lmer(lsumareaqsqrt ~ scale(meanRH42) + scale(meanT35) + scale(VisSky) + (1 | Site:Quadrat) + (1 | Year), data=df2, REML=FALSE))
RHs36Vis
RHflmers36Vis <- AICc(lmer(lsumareaqsqrt ~ scale(meanRH42) + scale(meanT42) + scale(VisSky) + (1 | Site:Quadrat) + (1 | Year), data=df2, REML=FALSE))


###########################################################################################################################################################################################################################################################################################################################################

#now put AICc values into a dataframe and save as a file to then calculate delta AICc and rank models

a <- data.frame(RHflmers1, RHflmers2)
a <- data.frame(a, RHflmers3)
a <- data.frame(a, RHflmers4)
a <- data.frame(a, RHflmers5)
a <- data.frame(a, RHflmers6)
a <- data.frame(a, RHflmers7)
a <- data.frame(a, RHflmers8)
a <- data.frame(a, RHflmers9)
a <- data.frame(a, RHflmers10)
a <- data.frame(a, RHflmers11)
a <- data.frame(a, RHflmers12)
a <- data.frame(a, RHflmers13)
a <- data.frame(a, RHflmers14)
a <- data.frame(a, RHflmers15)
a <- data.frame(a, RHflmers16)
a <- data.frame(a, RHflmers17)
a <- data.frame(a, RHflmers18)
a <- data.frame(a, RHflmers19)
a <- data.frame(a, RHflmers20)
a <- data.frame(a, RHflmers21)
a <- data.frame(a, RHflmers22)
a <- data.frame(a, RHflmers23)
a <- data.frame(a, RHflmers24)
a <- data.frame(a, RHflmers25)
a <- data.frame(a, RHflmers26)
a <- data.frame(a, RHflmers27)
a <- data.frame(a, RHflmers28)
a <- data.frame(a, RHflmers29)
a <- data.frame(a, RHflmers30)
a <- data.frame(a, RHflmers31)
a <- data.frame(a, RHflmers32)
a <- data.frame(a, RHflmers33)
a <- data.frame(a, RHflmers34)
a <- data.frame(a, RHflmers35)
a <- data.frame(a, RHflmers36)

a <- t(a)
a
#write.csv() #Write a file of the vegetation model AICc outputs
rm(a)
a <- data.frame(RHflmers1Vis, RHflmers2Vis)
a <- data.frame(a, RHflmers3Vis)
a <- data.frame(a, RHflmers4Vis)
a <- data.frame(a, RHflmers5Vis)
a <- data.frame(a, RHflmers6Vis)
a <- data.frame(a, RHflmers7Vis)
a <- data.frame(a, RHflmers8Vis)
a <- data.frame(a, RHflmers9Vis)
a <- data.frame(a, RHflmers10Vis)
a <- data.frame(a, RHflmers11Vis)
a <- data.frame(a, RHflmers12Vis)
a <- data.frame(a, RHflmers13Vis)
a <- data.frame(a, RHflmers14Vis)
a <- data.frame(a, RHflmers15Vis)
a <- data.frame(a, RHflmers16Vis)
a <- data.frame(a, RHflmers17Vis)
a <- data.frame(a, RHflmers18Vis)
a <- data.frame(a, RHflmers19Vis)
a <- data.frame(a, RHflmers20Vis)
a <- data.frame(a, RHflmers21Vis)
a <- data.frame(a, RHflmers22Vis)
a <- data.frame(a, RHflmers23Vis)
a <- data.frame(a, RHflmers24Vis)
a <- data.frame(a, RHflmers25Vis)
a <- data.frame(a, RHflmers26Vis)
a <- data.frame(a, RHflmers27Vis)
a <- data.frame(a, RHflmers28Vis)
a <- data.frame(a, RHflmers29Vis)
a <- data.frame(a, RHflmers30Vis)
a <- data.frame(a, RHflmers31Vis)
a <- data.frame(a, RHflmers32Vis)
a <- data.frame(a, RHflmers33Vis)
a <- data.frame(a, RHflmers34Vis)
a <- data.frame(a, RHflmers35Vis)
a <- data.frame(a, RHflmers36Vis)


a <- t(a)
#write.csv() #Write a file of the gap fraction model AICc outputs
rm(a)
a

#######################################################################################################################

#Load dataframes, calculate delta AICc and rank models for models with ground cover
#hm1 <- read.csv() #Read the file of the vegetation model AICc outputs

colnames(hm1)[1] <- "Model"
colnames(hm1)[2] <- "AICc"
AICcW <- akaike.weights(hm1$AICc)
final <- data.frame(AICcW, hm1)
#Models RHflmers2 and RHflmers32 with delta AICc less than 7

#Calculate delta AICc and rank models for models with gap fraction
#hm2 <- read.csv() #Read the file of the gap model AICc outputs

colnames(hm2)[1] <- "Model"
colnames(hm2)[2] <- "AICc"

AICcVW <- akaike.weights(hm2$AICc)
finalV <- data.frame(AICcVW, hm2)
#Model RHflmers2Vis with delta AICc less than 7

#Models within AICc of less than seven for using ground cover
#test whether including (Site/Quadrat), or (Site:Quadrat)
#significantly affects models - if it doesn't then I can compare
#these models against the gap fraction model, because gap fraction models
#would not converge with (Site/Quadrat)
veg1 <- lmer(lsumareaqsqrt ~ scale(meanRH7) + scale(meanT14) + Vegetation + (1 | Site/Quadrat) + (1 | Year) + scale(meanRH7):scale(meanT14) + scale(meanT14):Vegetation, data=df2, REML=FALSE)
veg2 <- lmer(lsumareaqsqrt ~ scale(meanRH42) * scale(meanT14) * Vegetation + (1 | Site/Quadrat) + (1 | Year), data=df2, REML=FALSE)

veg1ns <- lmer(lsumareaqsqrt ~ scale(meanRH7) + scale(meanT14) + Vegetation + (1 | Site:Quadrat) + (1 | Year) + scale(meanRH7):scale(meanT14) + scale(meanT14):Vegetation, data=df2, REML=FALSE)
veg2ns <- lmer(lsumareaqsqrt ~ scale(meanRH42) * scale(meanT14) * Vegetation + (1 | Site:Quadrat) + (1 | Year), data=df2, REML=FALSE)

anova(veg1,veg1ns)
anova(veg2,veg2ns)
#Does not significantly affect models if constructed with (Site/Quadrat) or (Site:Quadrat)
#Therefore construct models with AICc less than 7 (i.e. plausible models) with (Site:Quadrat)
#and compare model with gap fraction, to those including ground cover
#get AICc values of these models - merge into a dataframe to then
#rank based on delta AICc values

veg1ns <- AICc(lmer(lsumareaqsqrt ~ scale(meanRH7) + scale(meanT14) + Vegetation + (1 | Site:Quadrat) + (1 | Year) + scale(meanRH7):scale(meanT14) + scale(meanT14):Vegetation, data=df2, REML=FALSE))
veg2ns <- AICc(lmer(lsumareaqsqrt ~ scale(meanRH42) * scale(meanT14) * Vegetation + (1 | Site:Quadrat) + (1 | Year), data=df2, REML=FALSE))

vis1 <- AICc(lmer(lsumareaqsqrt ~ scale(meanRH7) + scale(meanT14) + scale(VisSky) + (1 | Site:Quadrat) + (1 | Year) + scale(meanRH7):scale(meanT14) + scale(meanT14):scale(VisSky), data=df2, REML=FALSE))

a <- data.frame(veg1ns, veg2ns)
a <- data.frame(a, vis1)

a <- t(a)

#write.csv() #Write a file for the AICc values to compare the above most supported models
rm(a)
#######################################################################################################################

#Load dataframe and calculate the delta AICc of the plausible models with gap fraction and
#ground cover

#hm1 <- read.csv() #Read the file for the AICc values to compare the above most supported models

colnames(hm1)[1] <- "Model"
colnames(hm1)[2] <- "AICc"

AICcW <- akaike.weights(hm1$AICc)
final <- data.frame(AICcW, hm1)


#Plausible models (delta AICc less than 7)
#1 - lmer(lsumareaqsqrt ~ scale(meanRH7) + scale(meanT14) + Vegetation + (1 | Site:Quadrat) + (1 | Year) + scale(meanRH7):scale(meanT14) + scale(meanT14):Vegetation, data=df2, REML=FALSE)
#2 - lmer(lsumareaqsqrt ~ scale(meanRH7) + scale(meanT14) + scale(VisSky) + (1 | Site:Quadrat) + (1 | Year) + scale(meanRH7):scale(meanT14) + scale(meanT14):scale(VisSky), data=df2, REML=FALSE)
#3 - lmer(lsumareaqsqrt ~ scale(meanRH42) * scale(meanT14) * Vegetation + (1 | Site:Quadrat) + (1 | Year), data=df2, REML=FALSE)
###################################################################################################################################################################################################################################

#construct models fitted with REML and (Site/Quadrat) as random effect where possible
#, calculate the variance explained by model (marginal and conditional R squared)
a <- lmer(lsumareaqsqrt ~ scale(meanRH7) + scale(meanT14) + Vegetation + (1 | Site/Quadrat) + (1 | Year) + scale(meanRH7):scale(meanT14) + scale(meanT14):Vegetation, data=df2, REML=TRUE)
r.squaredGLMM(a)
#> r.squaredGLMM(a)
#R2m      R2c
#[1,] 0.289441 0.607931

summary(a)
library(lmerTest)
anova(a, ddf = "Kenward-Roger")
#Type III Analysis of Variance Table with Kenward-Roger's method
#                               Sum Sq Mean Sq NumDF  DenDF  F value    Pr(>F)    
#scale(meanRH7)                 0.0278  0.0278     1 370.49   0.3562    0.5510    
#scale(meanT14)                10.0623 10.0623     1 367.37 128.8054 < 2.2e-16 ***
#Vegetation                     0.4025  0.1342     3  23.11   1.6908    0.1967    
#scale(meanRH7):scale(meanT14)  2.4289  2.4289     1 368.19  31.0916  4.77e-08 ***
#scale(meanT14):Vegetation      2.0488  0.6829     3 370.54   8.7422  1.29e-05 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


b <- lmer(lsumareaqsqrt ~ scale(meanRH7) + scale(meanT14) + scale(VisSky) + (1 | Site:Quadrat) + (1 | Year) + scale(meanRH7):scale(meanT14) + scale(meanT14):scale(VisSky), data=df2, REML=TRUE)
r.squaredGLMM(b)
#R2m       R2c
#[1,] 0.2977217 0.5674551

summary(b)
anova(b, ddf = "Kenward-Roger")
#Type III Analysis of Variance Table with Kenward-Roger's method
#                               Sum Sq Mean Sq NumDF  DenDF  F value    Pr(>F)    
#scale(meanRH7)                 0.0074  0.0074     1 383.85   0.0919  0.761968    
#scale(meanT14)                12.8813 12.8813     1 376.95 159.7119 < 2.2e-16 ***
#scale(VisSky)                  0.9400  0.9400     1  25.22  11.6547  0.002171 ** 
#scale(meanRH7):scale(meanT14)  2.2461  2.2461     1 371.49  27.8491 2.235e-07 ***
#scale(meanT14):scale(VisSky)   1.2667  1.2667     1 377.10  15.7053 8.853e-05 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



c <- lmer(lsumareaqsqrt ~ scale(meanRH42) * scale(meanT14) * Vegetation + (1 | Site/Quadrat) + (1 | Year), data=df2, REML=TRUE)
r.squaredGLMM(c)
#R2m       R2c
#[1,] 0.2989616 0.6460298
summary(c)
anova(c, ddf = "Kenward-Roger")
#Type III Analysis of Variance Table with Kenward-Roger's method
#                                          Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
#scale(meanRH42)                           0.4294  0.4294     1 366.51  5.5669  0.018827 *  
#scale(meanT14)                            6.2351  6.2351     1 359.67 80.8415 < 2.2e-16 ***
#Vegetation                                0.4751  0.1584     3  26.95  2.0248  0.134152    
#scale(meanRH42):scale(meanT14)            0.5389  0.5389     1 362.61  6.9873  0.008565 ** 
#scale(meanRH42):Vegetation                0.1107  0.0369     3 373.66  0.4784  0.697528    
#scale(meanT14):Vegetation                 1.5993  0.5331     3 360.95  6.9120  0.000155 ***
#scale(meanRH42):scale(meanT14):Vegetation 1.2214  0.4071     3 360.36  5.2786  0.001423 ** 
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Check model assumptions
#test for multicolinearity
library(car)
vif(a)
vif(b)
vif(c)

#residual plots

#model a
plot(a)
qqnorm(resid(a))
qqline(resid(a))
plot(fitted(a),resid(a))
plot(df2$Vegetation,resid(a))
plot(df2$meanT14,resid(a))
plot(df2$meanRH7,resid(a))
plot(df2$Site,resid(a))
df2$SiteQuadrat <- paste0(df2$Site, df2$Quadrat)
df2$SiteQuadrat <- as.factor(df2$SiteQuadrat)
plot(df2$SiteQuadrat,resid(a))
plot(df2$Year,resid(a))

#model b
plot(b)
qqnorm(resid(b))
qqline(resid(b))
plot(fitted(b),resid(b))
plot(df2$Vegetation,resid(b))
plot(df2$meanT14,resid(b))
plot(df2$meanRH7,resid(b))
plot(df2$Site,resid(b))
df2$SiteQuadrat <- paste0(df2$Site, df2$Quadrat)
df2$SiteQuadrat <- as.factor(df2$SiteQuadrat)
plot(df2$SiteQuadrat,resid(b))
plot(df2$Year,resid(b))

#model c
plot(c)
qqnorm(resid(c))
qqline(resid(c))
plot(fitted(c),resid(c))
plot(df2$Vegetation,resid(c))
plot(df2$meanT14,resid(c))
plot(df2$meanRH7,resid(c))
plot(df2$Site,resid(c))
df2$SiteQuadrat <- paste0(df2$Site, df2$Quadrat)
df2$SiteQuadrat <- as.factor(df2$SiteQuadrat)
plot(df2$SiteQuadrat,resid(c))
plot(df2$Year,resid(c))


#Now check model residuals, overdispersion and outliers using DHARMa package
library(DHARMa)
#have to have variables and then fit model when using this package
dfscaled <- df2 %>% mutate(smeanRH7=scale(meanRH7),smeanT14=scale(meanT14),smeanRH42=scale(meanRH42),sVisSky=scale(VisSky))

a2 <- lmer(lsumareaqsqrt ~ smeanRH7 + smeanT14 + Vegetation + (1 | Site/Quadrat) + (1 | Year) + smeanRH7:smeanT14 + smeanT14:Vegetation, data=dfscaled, REML=TRUE)
b2 <- lmer(lsumareaqsqrt ~ smeanRH7 + smeanT14 + sVisSky + (1 | Site:Quadrat) + (1 | Year) + smeanRH7:smeanT14 + smeanT14:sVisSky, data=dfscaled, REML=TRUE)
c2 <- lmer(lsumareaqsqrt ~ smeanRH42 * smeanT14 * Vegetation + (1 | Site/Quadrat) + (1 | Year), data=dfscaled, REML=TRUE)

#model a 
#with random effects set on the fitted values
simulationOutput <- simulateResiduals(fittedModel=a2, re.form=NULL,plot=T,n=1000)
plotResiduals(simulationOutput, dfscaled$smeanT14)
plotResiduals(simulationOutput, dfscaled$smeanRH7)
plotResiduals(simulationOutput, dfscaled$Vegetation)
plotResiduals(simulationOutput, dfscaled$Site)
plotResiduals(simulationOutput, dfscaled$Year)
plotResiduals(simulationOutput, dfscaled$SiteQuadrat)
plotResiduals(simulationOutput, dfscaled$SiteQuadrat)
#all look fine


#model b
#with random effects set on the fitted values
simulationOutput <- simulateResiduals(fittedModel=b2, re.form=NULL,plot=T,n=1000)
plotResiduals(simulationOutput, dfscaled$smeanT14)
plotResiduals(simulationOutput, dfscaled$smeanRH7)
plotResiduals(simulationOutput, dfscaled$sVisSky)
plotResiduals(simulationOutput, dfscaled$Site)
plotResiduals(simulationOutput, dfscaled$Year)
plotResiduals(simulationOutput, dfscaled$SiteQuadrat)

#all look fine

#model c
#with random effects set on the fitted values
simulationOutput <- simulateResiduals(fittedModel=c2, re.form=NULL,plot=T,n=1000)
plotResiduals(simulationOutput, dfscaled$smeanT14)
plotResiduals(simulationOutput, dfscaled$smeanRH42)
plotResiduals(simulationOutput, dfscaled$Vegetation)
plotResiduals(simulationOutput, dfscaled$Site)
plotResiduals(simulationOutput, dfscaled$Year)
plotResiduals(simulationOutput, dfscaled$SiteQuadrat)

#all look fine
