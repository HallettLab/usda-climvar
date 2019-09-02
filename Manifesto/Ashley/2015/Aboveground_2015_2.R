##The purpose of this script is to test treatment effects on ANPP and cover for 2015
#To compare with Lina's belowground results

library(tidyverse)
library(nlme)
library(ggplot2)
library(multcomp)
library(lsmeans)
library(MuMIn)
library(gridExtra)

setwd("~/Dropbox/ClimVar/DATA/")
May_ANPP_FG<-read.csv("Plant_composition_data/ANPP/ANPP_CleanedData/ClimVar_ANPP-func-groups.csv")
#remove compost subplots for later
May_ANPP_FG_noC <- May_ANPP_FG %>% filter(subplot!="C")
#create subset with 2015 data only, remove compost subplot
May_ANPPfg_2015<-filter(May_ANPP_FG, year=='2015', subplot !="C")

###################
## look at composition plots only
###################
#3. compensation (H2)
#H2 will be confirmed if mixed plots have greater ANPP compared to grass only or forb only
May_ANPPfg_2015_noXC<-filter(May_ANPPfg_2015, subplot !="XC", harvest=="Second")
mfg<-lme(log(weight_g_m+1) ~treatment*subplot*func, random=~1|shelterBlock, May_ANPPfg_2015_noXC, na.action=na.exclude)
summary(mfg)
anova(mfg) #subplot significant
r.squaredGLMM(mfg)#32% of variation explained by fixed effects, 32% by whole model
qqnorm(residuals(mfg))
qqline(residuals(mfg))
shapiro.test(residuals(mfg)) #normal
LSfg<-lsmeans(mfg, ~subplot*func)
contrast(LSfg, "pairwise") #mixed plots have greater biomass than forb, but not grass

ggplot(d=May_ANPPfg_2015_noXC, aes(x=subplot, y=weight_g_m, fill=func)) +
  theme_linedraw()+
  facet_wrap(~treatment)+
  labs(x="composition treatment", y="ANPP g/m2")+
  geom_boxplot(aes(y=weight_g_m), shape=16)

ggplot(d=May_ANPPfg_2015_noXC, aes(x=subplot, y=weight_g_m, fill = treatment)) +
  theme_linedraw()+
  facet_wrap(~func)+
  labs(x="composition treatment", y="ANPP g/m2")+
  geom_boxplot(aes(y=weight_g_m), shape=16)

ggplot(d=traits_2015_comp, aes(x=CWM.DiamC, y=CWM.SRLC, color=subplot, shape=treatment))+
  geom_point()
  #geom_smooth(method="lm", se=F)

ggplot(d=traits_2015_comp, aes(x=CWM.DiamC, y=CWM.SRLC, color=treatment))+
  geom_point()

ggplot(d=traits_2015_comp, aes(x=CWM.DiamC, y=CWM.SRLF, color=subplot, shape=treatment))+
  geom_point()

ggplot(d=traits_2015_comp, aes(x=CWM.DiamC, y=CWM.PropF, color=subplot, shape=treatment))+
  geom_point()

ggplot(d=traits_2015_comp, aes(x=CWM.DiamC, y=CWM.Dens, color=subplot))+
  geom_point()+
  geom_smooth(method="lm")

