library(tidyverse)
library(nlme)
library(ggplot2)
library(multcomp)
library(lsmeans)
library(MuMIn)
setwd("~/Dropbox/ClimVar/DATA/Plant_composition_data")
May_ANPP<-read.csv("ANPP/ANPP_CleanedData/ClimVar_ANPP-peak.csv")
levels(May_ANPP$plot)
levels(May_ANPP$year)
levels(May_ANPP$date)
levels(May_ANPP$subplot)
levels(May_ANPP$treatment)
levels(May_ANPP$shelter)
levels(May_ANPP$shelterBlock)

#change plots, years, shelter to factors
May_ANPP[,'plot'] <- as.factor(as.character(May_ANPP[,'plot']))
May_ANPP[,'year'] <- as.factor(as.character(May_ANPP[,'year']))
May_ANPP[,'shelter'] <- as.factor(as.character(May_ANPP[,'shelter']))

##1. Does variability of rainfall affect forage production (H1)?
#Using a mixed model to examine the effects of rainfall timing (fixed effect) with shelterbloc nested within year as random effect

#create subset with no species manipulations (control community) only
May_XC<-filter(May_ANPP, subplot=='XC')

m1<-lme(weight_g_m ~treatment, random=~1|year/shelterBlock, May_XC, na.action=na.exclude)
summary(m1)
anova(m1)
r.squaredGLMM(m1) #12% of variation explained by fixed effects, 57% by whole model (interannual variation?)
qqnorm(residuals(m1))
qqline(residuals(m1))
shapiro.test(residuals(m1))
#normally distributed, continue
LS1<-lsmeans(m1, ~treatment)
contrast(LS1, "pairwise")
#control rain is significantly greater than all the drought treatments, no surprise
#control rain is most similar to spring dry

##2. Does seasonality of rainfall affect forage production (H1)?
##Expect peak ANPP to be highest when rainfall occurs during peak season or consistently
##Expect peak ANPP to be lowest when rainfall occurs late in season
#try again with control rain removed to compare only treatments with same total rainfall
May_XC_drought<-filter(May_XC, treatment!='controlRain')

m2<-lme(weight_g_m ~treatment, random=~1|year/shelterBlock, May_XC_drought, na.action=na.exclude)
summary(m2)
anova(m2)
r.squaredGLMM(m2) #only 2% of variation explained by fixed effects, 45% explained by whole model (interannual variation?)
qqnorm(residuals(m2))
qqline(residuals(m2))
shapiro.test(residuals(m2))
#normally distributed, continue
LS2<-lsmeans(m2, ~treatment)
contrast(LS2, "pairwise")
#no differences in total ANPP among drought treatments

#does variability in soil moisture affect variability in ANPP (H1: variability)?
#NEED TO DO: regression model
#SOIL MOISTURE DATA? ?DECAGON?? note that the below coefficient of variation calculation won't run yet
CV <- function(x){(sd(x)/mean(x))*100}
moistCV<-aggregate(moisture ~ plot*year, data= X, FUN = CV)



#3. compensation (H2)
#H2 will be confirmed if mixed plots have greater ANPP compared to grass only or forb only
#first remove compost treatment 
May_ANPP_noC<-filter(May_ANPP, subplot!='C')

m3<-lme(weight_g_m ~treatment*subplot, random=~1|year/shelterBlock, May_ANPP_noC, na.action=na.exclude)
summary(m3)
anova(m3)
r.squaredGLMM(m3)
qqnorm(residuals(m3))
qqline(residuals(m3))
shapiro.test(residuals(m3))
#not normally distributed, try log transform

m4<-lme(log(weight_g_m+1) ~treatment*subplot, random=~1|year/shelterBlock, May_ANPP_noC, na.action=na.exclude)
summary(m4)
anova(m4)
r.squaredGLMM(m4)#7% of variation explained by fixed effects, 47% explained by entire model (lots of interannual variation?)
qqnorm(residuals(m4))
qqline(residuals(m4))
shapiro.test(residuals(m4))
#barely normal
LS4<-lsmeans(m4, ~subplot)
contrast(LS4, "pairwise")
#mixed plots have greater ANPP compared to forb-only, but not grass-only
#no sig difference in ANPP between mixed plots and XC (no manipulation)




