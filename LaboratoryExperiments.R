
rm(list = ls())

library(car)
library(dplyr)
library(lubridate)
library(ggplot2)

#Experiment 2.2.1 analyses - incubation of rachises in dark incubators at 4, 10, 15, and 20 degrees
e1 <- read.csv("exp221data.csv")


#Chronologically number the sampling dates
e1$Date <- as.Date(e1$Date, "%d/%m/%Y",  tz="GMT") #Date represents the assessment date
e1$Idate <- as.Date(e1$Idate, "%d/%m/%Y",  tz="GMT") #Idate represents the date boxes were incubated
e1 <- e1 %>% mutate(aday=yday(Date))
e1 <- e1 %>% mutate(Iday=yday(Idate))
e1 <- e1 %>% mutate(dy=aday-Iday) #Calculate the number of days after incubation
e1 <- e1[,-c(9)] #Drop redundant column
dy <- e1$dy
e1 <- e1[,-c(11)]

#Chronologically ranking the number of days after incubation, and recombining with the dataset 
dy <- data.frame(dy)
dy <- dy %>% group_by(dy) %>% mutate(n = cur_group_id())
e1 <- cbind(dy,e1)
rm(dy)

#Incubation temperature as a factor
e1$T <- as.factor(e1$T)

#Only keeping data for boxes with the maximum proportion of rachises with immature apothecia (Stipes) for each box
maxe1 <- e1 %>% 
  group_by(T, Box) %>%
  filter(P == max(P))


#Find the earliest date this maximum was reached for each box
df <- maxe1 %>% group_by(T, Box) %>% filter(dy==min(dy))
rm(maxe1)

#Calculating the summary statistics
df %>% group_by(T) %>% summarise(res=quantile(P,probs=0.5))
#median 15 degrees = 62.5, 20 degrees = 100

df %>% group_by(T) %>% summarise(res=quantile(P,probs=0.75))
#upper quantile 15 degrees = 83.3, 20 degrees = 100

df %>% group_by(T) %>% summarise(res=quantile(P,probs=0.25))
#lower quantile 15 degrees = 56.7, 20 degrees = 92.6

#IQR
df %>% group_by(T) %>% summarise(res=IQR(P))

########Construct binomial glm to examine the effect of temperature on the probability of rachises developing immature apothecia (stipes)#############

df <- df %>% mutate(NS=Total-Stipe) #Calculating the number of rachises that didn't develop immature apothecia (stipes)#

mone <- glm(cbind(Stipe,NS)~T,family=binomial(link="logit"), data=df)

Anova(mone)

#Analysis of Deviance Table (Type II tests)

#Response: cbind(Stipe, NS)
#LR Chisq Df Pr(>Chisq)    
#T   35.608  1  2.413e-09 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Check model assumptions
library(DHARMa)
simulationOutput <- simulateResiduals(fittedModel = mone, plot = T)
#Difficulty calculating with the number of datapoints, but eyeballing the residuals they look fine#


#Now perform a Wilcoxon rank-sum test to examine whether the timing of maximum number of rachises that developed iummature apothecia
#differed between the two temperatures

#Make dataframe with the chronologically ordered sampling date number of maximum proportion
#of rachises with immature apothecia (stipes) in each box for 15 degree incubated (15D) and 20 degree incubated (20D)
#then do the wilcoxon rank sum test with this data, because the wilcoxon rank sum test takes ordinal data
df1 <- df
df1w15 <- df1 %>% filter(T==15)
df1w20 <- df1 %>% filter(T==20)

n <- df1w15$n
dfw <- data.frame(n)

colnames(dfw)[1] <- "15D"

n <- df1w20$n
dfw <- cbind(n,dfw)
colnames(dfw)[1] <- "20D"
rm(n)

#Converted the time to an ordinal scale for analyses using Wilcoxon rank sum test
e1w <- wilcox.test(n ~ T, data = df1, exact=FALSE)
e1w
??wilcox.test
#Wilcoxon rank sum test with continuity correction

#data:  n by T
#W = 47.5, p-value = 0.003103
#alternative hypothesis: true location shift is not equal to 0


#Checking Mann Whitney U test (another term for Wilcoxon rank sum test) assumption that the distributions are not normal, but similar.
boxplot(n ~ T, data = df1)
hist(dfw$'15D')
hist(dfw$'20D')
#Seems fine given the sample size#

#Summary statistics
median(df1w15$dy)
#74

quantile(df1w15$dy,probs=0.25)
#67

quantile(df1w15$dy,probs=0.75)
#74

median(df1w20$dy)
#53

quantile(df1w20$dy,probs=0.25)
#53

quantile(df1w20$dy,probs=0.75)
#60

####################   Plot the outputs####################

#The proportion of infected rachises over time, until the last day of the maximum for all boxes.
e1plot <- e1 %>% mutate(P=P/100)

#Creating labels for plot
Ts <- c(
  `15` = "15\u00B0C",
  `20` = "20\u00B0C")


library(plotrix)
#overlay the mean and se for each sampling point,
#Only plot up until the day after incubation where each box at a given temperature has reached the maximum proportion of rachises with 
#immature apothecia development (stipes)

e1plot <- e1plot %>% group_by(T,dy) %>% mutate(avP=median(P),) #Calculating the mean proportion of rachises with
#immature apothecia development (stipes) at each timepoint

e1plot15 <- e1plot %>% filter(T=="15") %>% filter(dy<75) #it's day 74 for 20 degrees
e1plot20 <- e1plot %>% filter(T=="20") %>% filter(dy<61) #it's day 60 for 20 degrees
e1plota <- merge(e1plot15, e1plot20,all="TRUE")

e1rawplot <- ggplot(e1plota,aes(x= dy,y = P))+
  geom_boxplot(aes(x= dy,y = P, group=dy),outliers=FALSE) +
  geom_point(position=position_jitter(w=0.1),shape = 1, alpha = 0.3, size = 3)+
  ylab("Proportion of rachises with immature apothecia")+
  xlab("Days after the start of incubation")+
  scale_x_continuous(limits = c(-5, 80), breaks = seq(0, 80, by = 10))+
  facet_grid(~T, labeller = as_labeller(Ts))+
  theme_bw()+
  theme(text=element_text(size=10),legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),strip.background =element_rect(fill="white"))
e1rawplot



#Visualising the fit of the glm model
df1 <- df1 %>% mutate(P=P/100)
predsone <- data.frame(predict(mone, list(T = df1$T),se.fit = TRUE,type="response"))
pone <- cbind(df1,predsone)
pone <- pone %>% mutate(CI=se.fit*1.96)
pone <- pone %>% mutate(UCL=fit+CI)
pone <- pone %>% mutate(LCL=fit-CI)

pmone <- ggplot(pone,aes(T,P))+
  geom_errorbar(aes(ymin=LCL, ymax=UCL), position=position_dodge(width=0.8), width=.15, size=1) +
  geom_point(aes(x=T,y=fit), size=3,data=pone)+
  geom_point(position=position_jitter(w=0.05),size=3,shape=1, colour="black")+
  ylim(0,1)+
  ylab("Probability of immature apothecia on rachis")+
  xlab("Incubation temperature (\u00B0C)")+
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, by = 0.2))+
  theme_bw()+
  theme(text=element_text(size=10),legend.position = "none")
pmone


###########   Section 2.2.2 Analyses    #########################################

#Experiment 2 - effect of previous incubation temperature (from experiment 1)
#on apothecia development when rachises are moved to ambient coniditions in the
#laboratory
#binomial glm examining previous incubation temperature effects
#on the probability that rachises will develop apothecia - use estimated margins means and Tukey
#method for calculation of p values for post hoc analysis examination of significant differences
#between groups
#Examine difference between treatments in the time taken to reach the maximum proportion of rachises
#with apothecia per box using Kruskal-Wallis rank sum test, followed by Wilcoxon rank sum test with continuity
#correction for post hoc pairwise comparisons between treatments
rm(list = ls())
#read dataframe
e3 <- read.csv("exp222data.csv")

#chronologically number the sampling dates as previously
e3$Date <- as.Date(e3$Date, "%d/%m/%Y",  tz="GMT")
e3$Idate <- as.Date(e3$Idate, "%d/%m/%Y",  tz="GMT")
e3 <- e3 %>% mutate(aday=yday(Date)) #date of assessment
e3 <- e3 %>% mutate(Iday=yday(Idate)) #date of incubation
e3 <- e3 %>% mutate(dy=aday-Iday)
dy <- e3$dy
e3 <- e3[,-c(10)]

dy <- data.frame(dy)
dy <- dy %>% group_by(dy) %>% mutate(n = cur_group_id())
e3 <- cbind(dy,e3)
rm(dy)

e3$T <- as.factor(e3$T)

#Now filter the maximum proportion of rachises with apothecia for each box as above

maxe3 <- e3 %>% 
  group_by(T, Box) %>%
  filter(P == max(P))

#Filter the earliest maximum for each box at each temperature
df3 <- maxe3 %>% group_by(T, Box) %>% filter(dy==min(dy))
rm(maxe3)

#Calculating the summary statistics
df3 %>% group_by(T) %>% summarise(res=quantile(P,probs=0.5))
#median 4 degrees = 96.2; median 10 degrees = 83.9; median 15 degrees = 75; 20 degrees = 63.6

df3 %>% group_by(T) %>% summarise(res=quantile(P,probs=0.75))
#upper quantile 4 degrees = 100; median 10 degrees = 98.2; median 15 degrees = 75; 20 degrees = 73.6

df3 %>% group_by(T) %>% summarise(res=quantile(P,probs=0.25))
#lower quantile 4 degrees = 92.3; median 10 degrees = 73.3; median 15 degrees = 75; 20 degrees = 60.9

df3 %>% group_by(T) %>% summarise(res=IQR(P))

#construct binomial glm to examine the effect of previous incubation temperature on the probability of rachises
#developing apothecia

df3 <- df3 %>% mutate(NS=Total-Apo) #Rachises with no apothecia

#glm
mtwo <- glm(cbind(Apo,NS)~T,family=binomial(link="logit"), data=df3)
Anova(mtwo)
#Analysis of Deviance Table (Type II tests)

#Response: cbind(Apo, NS)
#LR Chisq Df Pr(>Chisq)    
#T   18.429  3  0.0003588 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Check model assumptions with the DHARMa package
library(DHARMa)
simulationOutput <- simulateResiduals(fittedModel = mtwo, plot = T)
#Fine

#estimated margins means and Tukey method for calculation of p values for post hoc analysis examination of significant differences
#between groups
e2n <- emmeans(mtwo, pairwise~T, at = list(Type =c("4","10","15","20"), type="response", lmer.df = "Kenward-Roger"))
summary(e2n, type="response")

#$emmeans
#T   prob     SE  df asymp.LCL asymp.UCL
#4  0.931 0.0253 Inf     0.862     0.967
#10 0.851 0.0414 Inf     0.751     0.916
#15 0.732 0.0592 Inf     0.602     0.832
#20 0.708 0.0564 Inf     0.587     0.805

#Confidence level used: 0.95 
#Intervals are back-transformed from the logit scale 

#$contrasts
#contrast  odds.ratio    SE  df null z.ratio p.value
#T4 / T10        2.34 1.200 Inf    1   1.670  0.3395
#T4 / T15        4.91 2.430 Inf    1   3.219  0.0070
#T4 / T20        5.55 2.650 Inf    1   3.589  0.0019
#T10 / T15       2.10 0.932 Inf    1   1.663  0.3434
#T10 / T20       2.37 1.010 Inf    1   2.023  0.1794
#T15 / T20       1.13 0.459 Inf    1   0.298  0.9908

#P value adjustment: tukey method for comparing a family of 4 estimates 
#Tests are performed on the log odds ratio scale 


#Kruskal-Wallis rank sum Test to examine differences in timing of maximum apothecia development between treatments
#Assumptions are the same as the Wilcoxon. E.g. samples come from the same
#distribution shape


n4 <- df3 %>% filter(T==4)
n10 <- df3 %>% filter(T==10)
n15 <- df3 %>% filter(T==15)
n20 <- df3 %>% filter(T==20)


boxplot(n ~ T, data = df3)
hist(n4$n)
hist(n10$n)
hist(n15$n)
hist(n20$n)
#looks fine for this amount of data

K3 <- kruskal.test(n ~ T, data = df3)
K3

#Kruskal-Wallis rank sum test

#data:  n by T
#Kruskal-Wallis chi-squared = 16.563, df = 3, p-value = 0.0008692

pairwise.wilcox.test(df3$n, df3$T, p.adjust.method = "bonf", exact=FALSE)

#Pairwise comparisons using Wilcoxon rank sum test with continuity correction 

#data:  df3$n and df3$T 

        #4     10    15   
#10   0.335     -   -    
#15   0.020   0.742  -    
#20   0.012   0.126 1.000

#P value adjustment method: bonferroni 
#bonferroni is the conservative method#

#summary statistics for treatments
median(n4$dy)
#66

IQR(n4$dy)
#18


median(n10$dy)
#53.5

IQR(n10$dy)
#7

median(n15$dy)
#50

IQR(n15$dy)
#6

median(n20$dy)
#44

IQR(n20$dy)
#11.25


####################   Plot


#Plot 1: The proportion of infected rachises over time, until the last day of the maximum for all boxes.
e2plot <- e3 %>% mutate(P=P/100)
Ts2 <- c(
  `4` = "4\u00B0C",`10` = "10\u00B0C",`15` = "15\u00B0C",
  `20` = "20\u00B0C")

e2plot4 <- e2plot %>% filter(T=="4") %>% filter(dy<85)
e2plot10 <- e2plot %>% filter(T=="10") %>% filter(dy<91)
e2plot15 <- e2plot %>% filter(T=="15") %>% filter(dy<51)
e2plot20 <- e2plot %>% filter(T=="20") %>% filter(dy<51)
e2plota <- merge(e2plot15, e2plot20,all="TRUE")
e2plota <- merge(e2plota, e2plot10,all="TRUE")
e2plota <- merge(e2plota, e2plot4,all="TRUE")


library(plotrix)
e2rawplot <- ggplot(e2plota,aes(dy,P))+
  geom_boxplot(aes(x= dy,y = P, group=dy),outliers=FALSE) +
  geom_point(position=position_jitter(w=0.1),shape = 1, alpha = 0.2, size = 3)+
  ylab("Proportion of rachises with mature apothecia")+
  xlab("Days after removal from incubators")+
  facet_grid(~T, labeller = as_labeller(Ts2))+
  theme_bw()+
  theme(text=element_text(size=10),legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),strip.background =element_rect(fill="white"))
e2rawplot



#Now examining the glm model fit
df3 <- df3 %>% mutate(P=P/100)
predstwo <- data.frame(predict(mtwo, list(T = df3$T),se.fit = TRUE,type="response"))
ptwo <- cbind(df3,predstwo)
ptwo <- ptwo %>% mutate(CI=se.fit*1.96)
ptwo <- ptwo %>% mutate(UCL=fit+CI)
ptwo <- ptwo %>% mutate(LCL=fit-CI)

pmtwo <- ggplot(ptwo,aes(T,P))+
  geom_errorbar(aes(ymin=LCL, ymax=UCL), position=position_dodge(width=0.8), width=.15, size=1) +
  geom_point(aes(x=T,y=fit), size=3,data=ptwo)+
  geom_point(position=position_jitter(w=0.05),size=3,shape=1, colour="black")+
  ylab("Probability of mature apothecia on rachis")+
  xlab("Pre-laboratory incubation temperature (\u00B0C)")+
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, by = 0.2))+
  theme_bw()+
  theme(text=element_text(size=10),legend.position = "none")
pmtwo

#Done#

