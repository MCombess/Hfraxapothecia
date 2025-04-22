
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

###################################################################################################################################

#model a plots

df3$vegetation <- as.factor(df3$vegetation)
df3$Vegetation <- as.factor(df3$Vegetation)
df3$Quadrat <- as.factor(df3$Quadrat)
df3$Site <- as.factor(df3$Site)
df3$Site2 <- as.factor(df3$Site2)
df3$Site3 <- as.factor(df3$Site3)

#get 95% confidence interval of model to plot with the observed datapoints

preds <- predictInterval(merMod = a, newdata = df3,
                         level = 0.95, n.sims = 1000,
                         stat = "median", type="linear.prediction",
                         include.resid.var = FALSE)

#divide the predictions into 2018 and 2019 data so I can plot separately
dfpredsa <- cbind(df3,preds)
dfpredsa$date <- as.Date(dfpredsa$date, "%d/%m/%Y")

dfpredsa$Vegetation <- as.factor(dfpredsa$Vegetation)
dfpredsa$Quadrat <- as.factor(dfpredsa$Quadrat)
dfpredsa$Site <- as.factor(dfpredsa$Site)
dfpredsa$Site2 <- as.factor(dfpredsa$Site2)

dfpreds18a <- dfpredsa %>% filter(Year==18)
dfpreds19a <- dfpredsa %>% filter(Year==19)

ma18 <- ggplot(dfpreds18a,aes(date,lsumareaqsqrt))+
  geom_point(aes(color=vegetation, shape=vegetation), size=2,alpha=10)+
  geom_errorbar(aes(group=Quadrat,ymin=lwr,ymax=upr,color=vegetation),alpha=2)+
  scale_colour_manual(values=c(Vegetation="darkolivegreen",Soil="gold3", Litter="brown", Moss="blue4"))+
  scale_fill_manual(values=c(Vegetation="darkolivegreen3",Soil="gold3", Litter="brown2", Moss="blue1"))+
  scale_shape_manual(values=c(15,8,19,17))+
  ylim(0,4)+
  scale_y_continuous(limits = c(1.5, NA), breaks = pretty_breaks(n = 3))+
  xlim(as.Date(c('2018/05/01', '2018/10/10')))+
  #ylab(expression(Log[10](sqrt('apothecia area'(µm^2)))))+
  ggtitle("(A)")+
  facet_wrap(~Site2, nrow=6)+
  theme_bw()+
  theme(text=element_text(size=14),legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),strip.background =element_rect(fill="white"),axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 60, vjust = 0.5),axis.title.y = element_blank())
ma18



ma19 <- ggplot(dfpreds19a,aes(date,lsumareaqsqrt))+
  geom_point(aes(color=vegetation, shape=vegetation), size=2, alpha = 10)+
  geom_errorbar(aes(group=Quadrat,ymin=lwr,ymax=upr,color=vegetation),alpha=2)+
  scale_colour_manual(values=c(Vegetation="darkolivegreen",Soil="gold3", Litter="brown", Moss="blue4"))+
  scale_fill_manual(values=c(Vegetation="darkolivegreen3",Soil="gold3", Litter="brown2", Moss="blue1"))+
  scale_shape_manual(values=c(15,8,19,17))+
  ylim(0,4)+
  scale_y_continuous(limits = c(1.5, NA), breaks = pretty_breaks(n = 3))+
  xlim(as.Date(c('2019/05/01', '2019/10/10')))+
  #ylab(expression(Log[10](sqrt('apothecia area'(µm^2)))))+
  ggtitle("(B)")+
  facet_wrap(~Site2, nrow=7)+
  theme_bw()+
  theme(text=element_text(size=14),legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),strip.background =element_rect(fill="white"),axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 60, vjust = 0.5),axis.title.y = element_blank())
ma19

library(cowplot)
library(ggpubr)

combined <- plot_grid(ma18, ma19, ncol = 1, align = "v", rel_heights = c(1, 1))
y_label <- text_grob(expression(Log[10](sqrt('apothecia area'(µm^2)))), rot = 90, vjust = 0.5, hjust = 0.5)
final_plot <- plot_grid(y_label, combined, ncol = 2, rel_widths = c(0.05, 1))


final_plot
#############################################################################################################################################

#model b plots
preds <- predictInterval(merMod = b, newdata = df3,
                         level = 0.95, n.sims = 1000,
                         stat = "median", type="linear.prediction",
                         include.resid.var = FALSE)

#divide the predictions into 2018 and 2019 data so I can plot separately
dfpredsb <- cbind(df3,preds)
dfpredsb$date <- as.Date(dfpredsb$date, "%d/%m/%Y")

dfpredsb$Vegetation <- as.factor(dfpredsb$Vegetation)
dfpredsb$Quadrat <- as.factor(dfpredsb$Quadrat)
dfpredsb$Site <- as.factor(dfpredsb$Site)
dfpredsb$Site3 <- as.factor(dfpredsb$Site3)

dfpreds18b <- dfpredsb %>% filter(Year==18)
dfpreds19b <- dfpredsb %>% filter(Year==19)


mb18 <- ggplot(dfpreds18b,aes(date,lsumareaqsqrt))+
  geom_point(aes(color=VisSky), size=2)+
  geom_errorbar(aes(group=Quadrat,ymin=lwr,ymax=upr,color=VisSky))+
  scale_fill_gradient(low = "grey90", high = "black", na.value = NA)+
  scale_color_gradient(low = "grey90", high = "black", na.value = NA)+
  labs(fill="Gap fraction") +
  labs(color="Gap fraction") +
  ylim(0,4)+
  scale_y_continuous(limits = c(1.5, NA), breaks = pretty_breaks(n = 3))+
  xlim(as.Date(c('2018/05/01', '2018/10/10')))+
  #ylab(expression(Log[10](sqrt('apothecia area'(µm^2)))))+
  ggtitle("(A)")+
  facet_wrap(~Site3, nrow=6)+
  theme_bw()+
  theme(text=element_text(size=14),legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),strip.background =element_rect(fill="white"),axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 60, vjust = 0.5),strip.text = element_text(size = 10),axis.title.y = element_blank())
mb18



mb19 <- ggplot(dfpreds19b,aes(date,lsumareaqsqrt))+
  geom_point(aes(color=VisSky), size=2)+
  geom_errorbar(aes(group=Quadrat,ymin=lwr,ymax=upr,color=VisSky))+
  scale_fill_gradient(low = "grey90", high = "black", na.value = NA)+
  scale_color_gradient(low = "grey90", high = "black", na.value = NA)+
  labs(fill="Gap fraction") +
  labs(color="Gap fraction") +
  ylim(0,4)+
  scale_y_continuous(limits = c(1.5, NA), breaks = pretty_breaks(n = 3))+
  xlim(as.Date(c('2019/05/01', '2019/10/10')))+
  #ylab(expression(Log[10](sqrt('apothecia area'(µm^2)))))+
  xlab("Date (2019)")+
  ggtitle("(B)")+
  facet_wrap(~Site3, nrow=7)+
  theme_bw()+
  theme(text=element_text(size=14),legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),strip.background =element_rect(fill="white"),axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 60, vjust = 0.5),strip.text = element_text(size = 10),axis.title.y = element_blank())
mb19

combined <- plot_grid(mb18, mb19, ncol = 1, align = "v", rel_heights = c(1, 1))
y_label <- text_grob(expression(Log[10](sqrt('apothecia area'(µm^2)))), rot = 90, vjust = 0.5, hjust = 0.5)
final_plot <- plot_grid(y_label, combined, ncol = 2, rel_widths = c(0.05, 1))


final_plot


############################################################################################################################################

#model c plots
#get 95% confidence interval of model to plot with the observed datapoints
preds <- predictInterval(merMod = c, newdata = df3,
                         level = 0.95, n.sims = 1000,
                         stat = "median", type="linear.prediction",
                         include.resid.var = FALSE)

#divide the predictions into 2018 and 2019 data so I can plot separately
dfpredsc <- cbind(df3,preds)
dfpredsc$date <- as.Date(dfpredsc$date, "%d/%m/%Y")

dfpredsc$Vegetation <- as.factor(dfpredsc$Vegetation)
dfpredsc$Quadrat <- as.factor(dfpredsc$Quadrat)
dfpredsc$Site <- as.factor(dfpredsc$Site)
dfpredsc$Site2 <- as.factor(dfpredsc$Site2)

dfpreds18c <- dfpredsc %>% filter(Year==18)
dfpreds19c <- dfpredsc %>% filter(Year==19)

mc18 <- ggplot(dfpreds18c,aes(date,lsumareaqsqrt))+
  geom_point(aes(color=vegetation, shape=vegetation),size=2, alpha = 10)+
  geom_errorbar(aes(group=Quadrat,ymin=lwr,ymax=upr,color=vegetation),alpha=2)+
  scale_colour_manual(values=c(Vegetation="darkolivegreen",Soil="gold3", Litter="brown", Moss="blue4"))+
  scale_fill_manual(values=c(Vegetation="darkolivegreen3",Soil="gold3", Litter="brown2", Moss="blue1"))+
  scale_shape_manual(values=c(15,8,19,17))+
  ylim(0,4)+
  scale_y_continuous(limits = c(1.5, NA), breaks = pretty_breaks(n = 3))+
  xlim(as.Date(c('2018-05-01', '2018-10-10')))+
  ylab(expression(Log[10](sqrt('apothecia area'(µm^2)))))+
  ggtitle("(A)")+
  facet_wrap(~Site2, nrow=6)+
  theme_bw()+
  theme(text=element_text(size=14),legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),strip.background =element_rect(fill="white"),axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 60, vjust = 0.5),axis.title.y = element_blank())
mc18



mc19 <- ggplot(dfpreds19c,aes(date,lsumareaqsqrt))+
  geom_point(aes(color=vegetation, shape=vegetation),size=2, alpha = 10)+
  geom_errorbar(aes(group=Quadrat,ymin=lwr,ymax=upr,color=vegetation))+#,alpha=0.5)+
  scale_colour_manual(values=c(Vegetation="darkolivegreen",Soil="gold3", Litter="brown", Moss="blue4"))+
  scale_fill_manual(values=c(Vegetation="darkolivegreen3",Soil="gold3", Litter="brown2", Moss="blue1"))+
  scale_shape_manual(values=c(15,8,19,17))+
  ylim(0,4)+
  scale_y_continuous(limits = c(1.5, NA), breaks = pretty_breaks(n = 3))+
  xlim(as.Date(c('2019-05-01', '2019-10-10')))+
  ylab(expression(Log[10](sqrt('apothecia area'(µm^2)))))+
  ggtitle("(B)")+
  facet_wrap(~Site2, nrow=7)+
  theme_bw()+
  theme(text=element_text(size=14),legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),strip.background =element_rect(fill="white"),axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 60, vjust = 0.5),axis.title.y = element_blank())
mc19

combined <- plot_grid(mc18, mc19, ncol = 1, align = "v", rel_heights = c(1, 1))
y_label <- text_grob(expression(Log[10](sqrt('apothecia area'(µm^2)))), rot = 90, vjust = 0.5, hjust = 0.5)
final_plot <- plot_grid(y_label, combined, ncol = 2, rel_widths = c(0.05, 1))


final_plot



