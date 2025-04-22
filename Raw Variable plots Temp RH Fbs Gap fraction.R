

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
library(ggpubr)
library(scales)

#plot raw variables
#first plot the raw fruit body data

site_names <-c("A"="A","B"="B","C"="C","D"="D","Carmarthenshire" = "C'thenshire","Devon" = "Devon", "Durham" = "C' Durham", "Hampshire" = "Hampshire", "Shropshire"="Shropshire", "Wiltshire"="Wiltshire", "Hexham"="N'berland")

library(ggh4x)


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
df2 <- df[,-c(18,19,23,24)]
df2 <- df2[!duplicated(df2),]
df <- df2
rm(df)


#So I need to make a column called 'Site2' that combines Location and Quadrat e.g. 'Carmarthenshire A',
#another one called 'Site3' that adds on the gap fraction e.g. 'Carmarthenshire A (Gf = 0.019)'
#Also need a new column vegetation where Moss = Moss, Litter = Litter, Mud = Soil, Herb = Vegetation
df2$Site2 <- paste(df2$Site, df2$Quadrat, sep = " ")
df2$vegetation <- df2$Vegetation

levels(df2$vegetation)
levels(df2$vegetation) <- c("Vegetation", "Litter", "Moss", "Soil")

df2$Site3 <- paste(df2$Site2, df2$VisSky, sep = "; ")

df3 <- df2


df318 <- df3 %>% filter(date<"2019-01-01")
df319 <- df3 %>% filter(date>"2019-01-01")

rawfbplot19 <- ggplot(df319,aes(date,sumareaqmm))+
  geom_point(aes(color=factor(vegetation, levels = c("Vegetation", "Moss", "Soil", "Litter")), shape=factor(vegetation, levels = c("Vegetation", "Moss", "Soil", "Litter"))), size=1.5,alpha=10)+
  #geom_line(aes(group=Quadrat,colour=vegetation))+
  scale_color_manual(values=c(Vegetation="darkolivegreen", Moss="blue4", Soil="gold3", Litter="brown"), name="Ground\ncover")+
  scale_shape_manual(values=c(15,19,17,8), name="Ground\ncover")+
  scale_y_continuous(limits = c(0, NA), breaks = pretty_breaks(n = 2.8),expand = expansion(mult = c(0.1, 0.15)))+
  xlim(as.Date(c('2019-06-01', '2019-10-10')))+
  ylab(expression(Log[10](sqrt('apothecia area'(µm^2)))))+
  ggtitle("(B)")+
  facet_wrap(~Site2, scales="free_y",nrow=7)+
  theme_bw()+
  theme(text=element_text(size=14),legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),strip.background =element_rect(fill="white"),axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 60, vjust = 0.5),axis.title.y = element_blank(),axis.text.y = element_text(size=10))
rawfbplot19


rawfbplot18 <- ggplot(df318,aes(date,sumareaqmm, colour=vegetation))+
  geom_point(aes(color=factor(vegetation, levels = c("Vegetation", "Moss", "Soil", "Litter")), shape=factor(vegetation, levels = c("Vegetation", "Moss", "Soil", "Litter"))), size=1.5,alpha=10)+
  scale_color_manual(values=c(Vegetation="darkolivegreen", Moss="blue4", Soil="gold3", Litter="brown"), name="Ground\ncover")+
  scale_shape_manual(values=c(15,19,17,8), name="Ground\ncover")+
  xlim(as.Date(c('2018-05-15', '2018-10-10')))+
  ylab(expression(Log[10](sqrt('apothecia area'(µm^2)))))+
  ggtitle("(A)")+
  facet_wrap(~Site2, scales="free_y",nrow=6)+
  scale_y_continuous(limits = c(0, NA), breaks = pretty_breaks(n = 3),expand = expansion(mult = c(0.1, 0.15)))+
  theme_bw()+
  theme(text=element_text(size=14),legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),strip.background =element_rect(fill="white"),axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 60, vjust = 0.5),axis.title.y = element_blank(),axis.text.y = element_text(size=10))
rawfbplot18

combined <- plot_grid(rawfbplot18, rawfbplot19, ncol = 1, align = "v", rel_heights = c(1, 1))
y_label <- text_grob(expression('Apothecia area ' (mm^2)), rot = 90, vjust = 0.5, hjust = 0.5)
final_plot <- plot_grid(y_label, combined, ncol = 2, rel_widths = c(0.05, 1))


final_plot



############################################################################################################################################################################################################################

#ViskSky

library(plyr)
df3e <- mapvalues(df3$Vegetation, from = c("Herb", "Mud"), to = c("Vegetation", "Soil"))
df3a <- cbind(df3e, df3)
df3a <- df3a[,-2]
colnames(df3a)[1] <- "Site"
df3plot <- df3
df3plot$vegetation <- ordered(df3plot$vegetation, levels = c("Vegetation", "Moss", "Soil", "Litter"))

df3plote <- mapvalues(df3plot$Site, from = c("Durham", "Hexham"), to = c("County Durham", "Northumberland"))
df3plot <- cbind(df3plote, df3plot)
df3plot <- df3plot[,-c(2,4)]
colnames(df3plot)[1] <- "Site"


df3plot$Site <- ordered(df3plot$Site, levels = c("Carmarthenshire", "Devon", "County Durham", "Hampshire", "Shropshire", "Wiltshire", "Northumberland"))

library(ggpattern)


vsdf3plot <- df3plot[c(1, 7, 14, 21, 27, 34, 41, 51, 56, 68, 70, 78, 85, 93, 98, 102, 108, 116, 125, 130, 135, 147, 151, 153, 310, 313, 320, 327),]
library(ggrepel)
vsplot2 <- ggplot(vsdf3plot,aes(Site,VisSky, color=vegetation))+
  geom_point(aes(shape=vegetation),size=1)+
  geom_text_repel(data = vsdf3plot,
                  aes(label = Quadrat, x=Site, y=VisSky), show.legend = FALSE)+
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5), colour = c('black'), linetype = "dashed")+
  scale_color_manual(values=c(Vegetation="darkolivegreen", Moss="blue4", Soil="gold3", Litter="brown"), name="Ground cover")+
  scale_shape_manual(values=c(Vegetation=15, Moss=19, Soil=17, Litter=8), name="Ground cover")+
  labs(x = NULL) + 
  ylab("Gap fraction")+
  scale_x_discrete(labels = c("Carmarthenshire" = "C'marthenshire","Devon" = "Devon", "County Durham" = "C' Durham", "Hampshire" = "Hampshire", "Shropshire"="Shropshire", "Wiltshire"="Wiltshire", "Northumberland"="N'mberland"))+
  labs(fill="Ground cover")+
  ggtitle("(A)")+
  theme_bw()+
  theme(text=element_text(size=12),legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),strip.background =element_rect(fill="white"),axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 60, vjust = 0.5),axis.title.y = element_blank())
vsplot2

vplot <- ggplot(vsdf3plot,aes(vegetation,VisSky))+
  geom_boxplot(notch=FALSE)+
  geom_point(position=position_jitter(w=0.05),size=3,shape=1, colour="black")+
  labs(x = NULL) + 
  ylab("Gap fraction")+
  labs(fill="Ground cover")+
  ggtitle("(B)")+
  theme_classic()+
  theme(text=element_text(size=10),legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.ticks.x = element_blank(),axis.title.y = element_blank())
vplot

combined <- plot_grid(vsplot2, vplot, ncol = 2, align = "h", rel_heights = c(1, 1))
y_label <- text_grob("Gap fraction", rot = 90, vjust = 0.5, hjust = 0.5)
final_plot <- plot_grid(y_label, combined, ncol = 2, rel_widths = c(0.05, 1))
final_plot


#####################################################################################################################################################################################################################################################################
rm(list = ls())
#Mean T and RH plots
alldata <- read.csv("TandRHvalues.csv")

alldata$date <- as.Date(alldata$date, "%Y-%m-%d")

library(plyr)
alldatae <- mapvalues(alldata$Site, from = c("Hexham", "Durham"), to = c("Northumberland", "County Durham"))
alldata <- cbind(alldatae, alldata)
alldata <- alldata[,-c(2,3)]
colnames(alldata)[1] <- "Site"
detach("package:plyr", unload = TRUE)
alldata18 <- alldata %>% filter(date<"2019-01-01")
alldata19 <- alldata %>% filter(date>"2019-01-01")
alldata18 <- alldata18 %>% filter(Site!="Northumblerand")
#split T and RH data by month
alldata18June <- alldata18 %>% filter(date>"2018-05-31", date<"2018-07-01")
alldata18July <- alldata18 %>% filter(date>"2018-06-30", date<"2018-08-01")
alldata18August <- alldata18 %>% filter(date>"2018-07-31", date<"2018-09-01")
alldata18June["Month"] <- "June"
alldata18July["Month"] <- "July"
alldata18August["Month"] <- "August"
alldata18month <- rbind(alldata18June, alldata18July)
alldata18month <- rbind(alldata18month, alldata18August)
rm(alldata18June, alldata18July, alldata18August)

alldata19June <- alldata19 %>% filter(date>"2019-05-31", date<"2019-07-01")
alldata19July <- alldata19 %>% filter(date>"2019-06-30", date<"2019-08-01")
alldata19August <- alldata19 %>% filter(date>"2019-07-31", date<"2019-09-01")
alldata19June["Month"] <- "June"
alldata19July["Month"] <- "July"
alldata19August["Month"] <- "August"
alldata19month <- rbind(alldata19June, alldata19July)
alldata19month <- rbind(alldata19month, alldata19August)
rm(alldata19June, alldata19July, alldata19August)


library(ggplot2)
RHplot18 <- ggplot(alldata18,aes(date,meanRH))+
  geom_line()+
  xlim(as.Date(c("2018-05-15","2018-09-30")))+
  scale_y_continuous(limits = c(54, 100), breaks = seq(60, 100, by = 10))+
  ylab("Daily mean relative humidity (%)")+
  ggtitle("(A)")+
  facet_wrap(~Site, nrow=3,scales="free_x")+
  theme_bw()+
  theme(text=element_text(size=12),legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),strip.background =element_rect(fill="white"),
        axis.title.x = element_blank(),axis.title.y = element_blank())
RHplot18


RHplot19 <- ggplot(alldata19,aes(date,meanRH))+
  geom_line()+
  xlim(as.Date(c("2019-05-15","2019-09-30")))+
  scale_y_continuous(limits = c(54, 100), breaks = seq(60, 100, by = 10))+
  ylab("Daily mean relative humidity (%)")+
  ggtitle("(B)")+
  facet_wrap(~Site, nrow=4,scales="free_x")+
  theme_bw()+
  theme(text=element_text(size=12),legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),strip.background =element_rect(fill="white"),
        axis.title.x = element_blank(),axis.title.y = element_blank())
RHplot19


Tplot18 <- ggplot(alldata18,aes(date,meanT))+
  geom_line()+
  xlim(as.Date(c("2018-05-15","2018-09-30")))+
  scale_y_continuous(breaks = seq(8, 21, by = 2))+
  ylab("Daily mean temperature (°C)")+
  scale_y_continuous(limits = c(8, 23), breaks = seq(10, 22, by = 2))+
  ggtitle("(A)")+
  facet_wrap(~Site, nrow=3,scales="free_x")+
  theme_bw()+
  theme(text=element_text(size=12),legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),strip.background =element_rect(fill="white"),
        axis.title.x = element_blank(),axis.title.y = element_blank())
Tplot18


Tplot19 <- ggplot(alldata19,aes(date,meanT))+
  geom_line()+
  xlim(as.Date(c("2019-05-15","2019-09-30")))+
  scale_y_continuous(breaks = seq(8, 22, by = 2))+
  ylab("Daily mean temperature (°C)")+
  scale_y_continuous(limits = c(8, 23), breaks = seq(10, 22, by = 2))+
  facet_wrap(~Site, nrow=4,scales="free_x")+
  ggtitle("(B)")+
  theme_bw()+
  theme(text=element_text(size=12),legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),strip.background =element_rect(fill="white"),
        axis.title.x = element_blank(),axis.title.y = element_blank())
Tplot19

####################################################################################################

#Boxplot the monthly T and RH by site, facet grids month
alldata18month$Y <- 2018
alldata19month$Y <- 2019
alldatamonth <- rbind(alldata18month,alldata19month)
alldatamonth$Site <- as.factor(alldatamonth$Site)
alldatamonth$Month <- as.factor(alldatamonth$Month)
alldatamonth$Y <- as.factor(alldatamonth$Y)
alldatamonth$MonthY <- paste(alldatamonth$Month, alldatamonth$Y, sep = " ")
alldatamonth$MonthY <- as.factor(alldatamonth$MonthY)
alldatamonth$MonthY <- ordered(alldatamonth$MonthY, levels = c("June 2018", "June 2019", "July 2018", "July 2019", "August 2018", "August 2019"))
alldatamonth$Site <- ordered(alldatamonth$Site, levels = c("Carmarthenshire", "Devon", "County Durham", "Hampshire", "Shropshire", "Wiltshire", "Northumberland"))
alldatamonth$Month <- ordered(alldatamonth$Month, levels = c("June", "July", "August"))
library(ggh4x)
library(ggpattern)
Tbplot <- ggplot(alldatamonth,aes(x=Site, y=meanT))+
  geom_boxplot(notch=FALSE)+
  ylab("Daily mean temperature (°C)")+
  scale_y_continuous(limits = c(8, 23), breaks = seq(10, 22, by = 2))+
  scale_x_discrete(labels = c("Carmarthenshire" = "C'shire","Devon" = "Devon", "County Durham" = "C' Dur'", "Hampshire" = "H'shire", "Shropshire"="S'shire", "Wiltshire"="W'shire", "Northumberland"="N'land"))+
  facet_wrap(Month~Y,scales="free_x",nrow=3)+
  ggtitle("(C)")+
  theme_bw()+
  theme(text=element_text(size=12),legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),strip.background =element_rect(fill="white"),
        axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_text(angle = 60, vjust = 0.5))
Tbplot


RHbplot <- ggplot(alldatamonth,aes(x=Site, y=meanRH))+
  geom_boxplot(notch=FALSE)+
  labs(x = NULL) + 
  ylab("Daily mean relative humidity (%)")+
  scale_y_continuous(limits = c(54, 100), breaks = seq(60, 100, by = 10))+
  scale_x_discrete(labels = c("Carmarthenshire" = "C'shire","Devon" = "Devon", "County Durham" = "C' Dur'", "Hampshire" = "H'shire", "Shropshire"="S'shire", "Wiltshire"="W'shire", "Northumberland"="N'land"))+
  facet_wrap(Month~Y,scales="free_x",nrow=3)+
  ggtitle("(C)")+
  theme_bw()+
  theme(text=element_text(size=12),legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),strip.background =element_rect(fill="white"),
        axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_text(angle = 60, vjust = 0.5))
RHbplot


combined <- plot_grid(Tplot18,Tplot19, ncol = 1, rel_heights = c(1, 1.333333))
combinedb <- plot_grid(Tbplot, nrow = 1, rel_heights = c(1, 1))
combinedf <- plot_grid(combined, combinedb,
                       ncol=2, rel_heights = c(1, 1))

combinedf

y_label <- text_grob("Daily mean temperature (°C)", rot = 90, vjust = 0.5, hjust = 0.5)
final_plotT <- plot_grid(y_label, combinedf, ncol = 2, rel_widths = c(0.05, 1))
final_plotT




combined <- plot_grid(RHplot18,RHplot19, ncol = 1, rel_heights = c(1, 1.333333))
combinedb <- plot_grid(RHbplot, nrow = 1, rel_heights = c(1, 1))
combinedf <- plot_grid(combined, combinedb,
                       ncol=2, rel_heights = c(1, 1))

combinedf

y_label <- text_grob("Daily mean relative humidity (%)", rot = 90, vjust = 0.5, hjust = 0.5)
final_plotRH <- plot_grid(y_label, combinedf, ncol = 2, rel_widths = c(0.05, 1))
final_plotRH

