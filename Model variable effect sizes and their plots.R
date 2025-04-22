
#plot models and calculate estimated marginal means for effect sizes

rm(list = ls())
library(dplyr)
library(lmerTest)
library(lme4)
library(MuMIn)
library(qpcR)
library(emmeans)
library(merTools)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)

#setwd("~/Final Write Up Folder/Papers post thesiss/Fungal Ecology/Apothecia paper/R to submit with paper")

df <- read.csv("fruitbodydataformodelling.csv")
colnames(df)
#construct maximum models, reduce models to only significant effects, rank models

#want sqrt of area so the dependent variable (sumareaq=the total apothecia area measured at a sampling location on a given date)
#is linear to meet model assumption.
#Log to meet model assumption that residuals have constant variance and they are normally distributed
df <- df %>% mutate(sumareaqsqrt=sqrt(sumareaqum))
df <- df %>% mutate(lsumareaqsqrt=log10(sumareaqsqrt))

#drop redundant column, and classify factor variables as factors
df <- df[,-1]
df$Site <- as.factor(df$Site)
df$Quadrat <- as.factor(df$Quadrat)
df$date <- as.Date(df$date, format = "%d/%m/%Y")
df$Vegetation <- as.factor(df$Vegetation)

#remove duplicated data so have one observation per sampling point of log10(sqrt(total apothecia area))
colnames(df)
df2 <- df[,-c(18,19,23,24,25,26)]
df2 <- df2[!duplicated(df2),]
df <- df2
rm(df)


#construct plausible models to then plot

a <- lmer(lsumareaqsqrt ~ scale(meanRH7) + scale(meanT14) + Vegetation + (1 | Site/Quadrat) + (1 | Year) + scale(meanRH7):scale(meanT14) + scale(meanT14):Vegetation, data=df2, REML=TRUE)

b <- lmer(lsumareaqsqrt ~ scale(meanRH7) + scale(meanT14) + scale(VisSky) + (1 | Site:Quadrat) + (1 | Year) + scale(meanRH7):scale(meanT14) + scale(meanT14):scale(VisSky), data=df2, REML=TRUE)

c <- lmer(lsumareaqsqrt ~ scale(meanRH42) * scale(meanT14) * Vegetation + (1 | Site/Quadrat) + (1 | Year), data=df2, REML=TRUE)

#So I need to make a column called 'Site2' that combines Location and Quadrat e.g. 'Carmarthenshire A',
#another one called 'Site3' that adds on the gap fraction e.g. 'Carmarthenshire A (Gf = 0.019)'
#Also need a new column vegetation where Moss = Moss, Litter = Litter, Mud = Soil, Herb = Vegetation
df2$Site2 <- paste(df2$Site, df2$Quadrat, sep = " ")
df2$vegetation <- df2$Vegetation

levels(df2$vegetation)
levels(df2$vegetation) <- c("Vegetation", "Litter", "Moss", "Soil")

df2$Site3 <- paste(df2$Site2, df2$VisSky, sep = "; ")

df3 <- df2

#################################################################################################################################################################################################################


#emmeans plot
#estimated margin mean on each ground cover type at high and low relative humidity (95 and 80) and temperature (14 and 16.5)
#Original
a3a <- emmeans(a, ~scale(meanRH7) + scale(meanT14) + Vegetation +scale(meanRH7):scale(meanT14) + scale(meanT14):Vegetation, at = list(meanRH7 = c(95,80), meanT14 =c(14,16.5)), lmer.df = "Kenward-Roger")

a3a
a3aplot <- data.frame(a3a)
a3aplot$Vegetation <- as.factor(a3aplot$Vegetation)
library(plyr)
a3aplote <- mapvalues(a3aplot$Vegetation, from = c("Herb", "Mud"), to = c("Vegetation", "Soil"))
a3aplot <- cbind(a3aplot, a3aplote)
colnames(a3aplot)[9] <- "GroundCover"
rm(a3aplote)
a3aplot$GroundCover <- ordered(a3aplot$GroundCover, levels = c("Vegetation", "Moss", "Soil", "Litter"))
a3aplot$meanT14 <- as.factor(a3aplot$meanT14)
a3aplot$meanRH7 <- as.factor(a3aplot$meanRH7)

a3a <- ggplot(a3aplot,aes(meanRH7,emmean))+
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE, color=factor(meanT14)), position=position_dodge(width=0.8), width=.1, linewidth=0.5) +
  geom_point(aes(color=meanT14, shape=meanT14), position=position_dodge(width=0.8), size=3)+
  geom_vline(xintercept = 1.5, colour = c('black'), linetype = "dashed")+
  scale_color_manual(values = c("14" = "#808080", "16.5" = "black"), name = "14 d mean\ntemperature (°C)")+
  scale_shape_manual(values=c(9,15), name = "14 d mean\ntemperature (°C)")+
  ylim(2.5,3.5)+
  xlab("7 d mean relative humidity (%)")+
  ylab(expression(Log[10](sqrt('apothecia area'(µm^2)))))+
  facet_wrap(~GroundCover, nrow=1)+
  theme_bw()+
  theme(axis.ticks.x = element_blank(), text=element_text(size=12),legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),strip.background =element_rect(fill="white"))
a3a



#########################################################################################################################################################################################

#Effect sizes of significant variables and interactions for model b and plots
#model b - significant interaction between gap fraction and mean temperature 14
#days prior to sample collection, and between mean temperature 14 days prior to sample
#collection and mean relative humidity 7 days prior to sample collection

#Gap fraction is a fraction - i.e. values less than 1. Therefore for effects size per unit
#variable I need to alter the gap fraction variable. I want increase per 0.1 gap fraction
#thus creat new dataframe and model (b4) for emtrends with VisSky multiplied by 10

df4 <- df2 %>% mutate(VisSky=VisSky*10)
b4 <- lmer(lsumareaqsqrt ~ scale(meanRH7) + scale(meanT14) + scale(VisSky) + (1 | Site:Quadrat) + (1 | Year) + scale(meanRH7):scale(meanT14) + scale(meanT14):scale(VisSky), data=df4, REML=TRUE)

#now the estimated marginal mean at low and high humidity (80 and 95), low and high gap fraction (0.05 and 0.15)
#and low and high temperature (14 and 16.5)

b4b <- emmeans(b, ~ scale(meanRH7) + scale(meanT14) + scale(VisSky) +scale(meanRH7):scale(meanT14) + scale(meanT14):scale(VisSky), at = list(meanRH7=c(80,95), VisSky=c(0.05,0.15), meanT14=c(14,16.5)))
b4b

b4bplot <- data.frame(b4b)
b4bplot$VisSky <- as.factor(b4bplot$VisSky)
b4bplot$meanRH7 <- as.factor(b4bplot$meanRH7)
b4bplot$meanT14 <- as.factor(b4bplot$meanT14)
b4bplot
vis_names <- c(
  `0.05` = "0.05 Gap fraction",
  `0.15` = "0.15 Gap fraction"
)


b4b <- ggplot(b4bplot,aes(meanRH7,emmean))+
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE, color=factor(meanT14)), position=position_dodge(width=0.8), width=.1, linewidth=0.2) +
  geom_point(aes(color=meanT14,shape=meanT14), position=position_dodge(width=0.8), size=2)+
  geom_vline(xintercept = 1.5, colour = c('black'), linetype = "dashed")+
  scale_color_manual(values = c("14" = "#808080", "16.5" = "black"), name = "14 d mean\ntemperature (°C)")+
  scale_shape_manual(values=c(9,15), name = "14 d mean\ntemperature (°C)")+
  ylim(2.5,3.58)+
  xlab("7 d mean relative humidity (%)")+
  ylab(expression(Log[10](sqrt('apothecia area'(µm^2)))))+
  scale_y_continuous(limits = c(2.5, 3.6), breaks = seq(2.5, 3.5, 0.5))+
  facet_wrap(~VisSky, scales='free', nrow=1,labeller = as_labeller(vis_names))+
  theme_bw()+
  theme(axis.ticks.x = element_blank(), text=element_text(size=11),legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),strip.background =element_rect(fill="white"))

b4b


#######################################################################################################################################


#now the estimated marginal mean at lower and higher humidity (85 and 92), lower and higher temperature (14 and 16.5)
#on all ground cover types
cm1 <- emmeans(c, ~scale(meanRH42)*scale(meanT14)*Vegetation, at = list(meanRH42 = c(85,92), meanT14 =c(14,16.5)), lmer.df = "Kenward-Roger")
cm1
cmplot <- data.frame(cm1)
cmplot$Vegetation <- as.factor(cmplot$Vegetation)
library(plyr)
cmplote <- mapvalues(cmplot$Vegetation, from = c("Herb", "Mud"), to = c("Vegetation", "Soil"))
cmplot <- cbind(cmplot, cmplote)
colnames(cmplot)[9] <- "GroundCover"
rm(cmplote)
cmplot$GroundCover <- ordered(cmplot$GroundCover, levels = c("Vegetation", "Moss", "Soil", "Litter"))
cmplot$meanT14 <- as.factor(cmplot$meanT14)
cmplot$meanRH42 <- as.factor(cmplot$meanRH42)


cm2 <- ggplot(cmplot,aes(meanRH42,emmean))+
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE, color=factor(meanT14)), position=position_dodge(width=0.8), width=.1, size=0.5) +
  geom_point(aes(color=meanT14,shape=meanT14), position=position_dodge(width=0.8), size=3)+
  geom_vline(xintercept = 1.5, colour = c('black'), linetype = "dashed")+
  scale_color_manual(values = c("14" = "#808080", "16.5" = "black"), name = "14 d mean\ntemperature (°C)")+
  scale_shape_manual(values=c(9,15), name = "14 d mean\ntemperature (°C)")+
  ylim(2.5,3.6)+
  xlab("42 d mean relative humidity (%)")+
  ylab(expression(Log[10](sqrt('apothecia area'(µm^2)))))+
  facet_wrap(~GroundCover, nrow=1)+
  scale_y_continuous(limits = c(2.5, 3.6), breaks = seq(2.5, 3.5, 0.5))+
  theme_bw()+
  theme(axis.ticks.x = element_blank(), text=element_text(size=12),legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),strip.background =element_rect(fill="white"))
cm2





