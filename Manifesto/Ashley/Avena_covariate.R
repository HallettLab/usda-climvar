library(tidyverse)
library(nlme)
library(ggplot2)
library(multcomp)
library(lsmeans)
library(MuMIn)
setwd("~/Dropbox/ClimVar/DATA/")
May_ANPP<-read.csv("Plant_composition_data/ANPP/ANPP_CleanedData/ClimVar_ANPP-peak.csv")
cover<-read.csv("Plant_composition_data/Cover/Cover_CleanedData/ClimVar_species-cover.csv")
data<- cover %>% dplyr::select(-X, -species, -genus, -status, -func, -func2) %>% spread(species_name, cover)
names(data)<-str_replace_all(names(data), c(" " = "." , "," = "" ))

df <- data %>% mutate(AvDom = if_else(Avena.barbata >= 40, 1, 0))
May_ANPP<- May_ANPP %>% dplyr::select(-date)
May_ANPP<-merge(May_ANPP, df)      

May_ANPP[,'AvDom'] <- as.factor(as.character(May_ANPP[,'AvDom']))
May_ANPP[,'year'] <- as.factor(as.character(May_ANPP[,'year']))
May_ANPP[,'shelter'] <- as.factor(as.character(May_ANPP[,'shelter']))

##1. Does variability of rainfall affect forage production (H1)?
#Using a mixed model to examine the effects of rainfall timing (fixed effect) with shelterblock nested within year as random effect

#create subset with no species manipulations (control community) only
May_XC<-filter(May_ANPP, subplot=='XC')
May_noAvDom <- filter(May_ANPP, Avena.barbata  <= 30, subplot =="XC")
May_XC_2017 <-filter(May_XC, year=="2017")

m1<-lme(weight_g_m ~shelter*AvDom*year, random=~1|shelterBlock, May_XC, na.action=na.exclude)
summary(m1)
anova(m1)
r.squaredGLMM(m1) #30% of variation explained by fixed effects, 39% by whole model (interannual variation?)
qqnorm(residuals(m1))
qqline(residuals(m1))
shapiro.test(residuals(m1))
#normally distributed, continue
LS1<-lsmeans(m1, ~year*AvDom)
contrast(LS1, "pairwise")


ggplot(d=May_XC, aes(x=treatment, y=weight_g_m)) +
  theme_linedraw()+
  facet_grid(~AvDom)+
  labs(x="rainfall treatment", y="ANPP g/m2")+
  #annotate("text", x= c("consistentDry", "controlRain","fallDry","springDry"), y = c(900, 975, 900,900), label = c("a", "b", "ab", "ab"), color = "#999999") +
  geom_boxplot(aes(y=weight_g_m), shape=16)

ggplot(d=May_XC, aes(x=year, y=weight_g_m)) +
  theme_linedraw()+
  facet_grid(~treatment*AvDom)+
  labs(x="", y="ANPP g/m2")+
  #annotate("text", x= c("2015", "2016","2017"), y = c(900, 975, 975), label = c("a", "b", "b"), color = "#999999") +
  geom_boxplot(aes(y=weight_g_m), shape=16)

May_ANPP_noC<-filter(May_ANPP, subplot!='C', subplot!='XC')

m3<-lme(weight_g_m ~treatment*subplot*AvDom, random=~1|shelterBlock, May_ANPP_noC, na.action=na.exclude)
summary(m3)
anova(m3)
r.squaredGLMM(m3)
qqnorm(residuals(m3))
qqline(residuals(m3))
shapiro.test(residuals(m3)) #normal
LS4<-lsmeans(m3, ~AvDom)
contrast(LS4, "pairwise")

ggplot(d=May_ANPP_noC, aes(x=AvDom, y=weight_g_m, fill=subplot)) +
  #facet_wrap(~year)+
  theme_linedraw()+
  labs(x="composition treatment", y="ANPP g/m2")+
  geom_boxplot(aes(y=weight_g_m), shape=16)

May_ANPP_2015<- May_ANPP %>% filter(year=="2015", subplot!="XC", subplot != "C")
May_ANPP_2015_noF <-May_ANPP_2015 %>% filter(subplot != "F")

m4<-lme(weight_g_m ~subplot*treatment*AvDom, random=~1|shelterBlock, May_ANPP_2015_noF, na.action=na.exclude)
summary(m4)
anova(m4)
r.squaredGLMM(m4)
qqnorm(residuals(m4))
qqline(residuals(m4))
shapiro.test(residuals(m4)) #normal
LS4<-lsmeans(m4, ~subplot*treatment*AvDom)
contrast(LS4, "pairwise")

ggplot(d=May_ANPP_2015_noF, aes(x=AvDom, y=weight_g_m, fill=subplot)) +
  #facet_wrap(~treatment)+
  theme_linedraw()+
  labs(x="composition treatment", y="ANPP g/m2")+
  geom_boxplot(aes(y=weight_g_m), shape=16)


