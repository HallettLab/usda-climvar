library(tidyverse)
library(nlme)
library(ggplot2)
library(multcomp)
library(lsmeans)
library(MuMIn)
library(RColorBrewer)
setwd("~/Dropbox/ClimVar/DATA/")
FG_ANPP<-read.csv("Plant_composition_data/ANPP/ANPP_CleanedData/ClimVar_ANPP-func-groups.csv")
levels(FG_ANPP$plot)
levels(FG_ANPP$year)
levels(FG_ANPP$func) #note FG names are not consistent
levels(FG_ANPP$subplot)
levels(FG_ANPP$treatment)
levels(FG_ANPP$shelter)
levels(FG_ANPP$shelterBlock)
head(FG_ANPP)

#change plots, years, shelter to factors
FG_ANPP[,'plot'] <- as.factor(as.character(FG_ANPP[,'plot']))
FG_ANPP[,'year'] <- as.factor(as.character(FG_ANPP[,'year']))
FG_ANPP[,'shelter'] <- as.factor(as.character(FG_ANPP[,'shelter']))

#rename func groups
FG_ANPP$func[FG_ANPP$func=="N"] <- "N-fixer"
FG_ANPP$func[FG_ANPP$func=="G"] <- "Grass"
FG_ANPP$func[FG_ANPP$func=="F"] <- "Forb"

##1. Does variability of rainfall affect forage production (H1)?
#Using a mixed model to examine the effects of rainfall timing (fixed effect) and subplot (fixed) with shelterblock nested within year as random effect
#create subset with only species manipulations (community treatments) only
FG_2015<-filter(FG_ANPP, year=='2015', harvest=="Second", subplot != "XC") %>% dplyr::select(-harvest)

m1<-lme(log(weight_g_m+1) ~treatment*func*subplot, random=~1|year/shelterBlock, FG_2015, na.action=na.exclude)
summary(m1)
anova(m1) #significant effects of subplot, func, treatment, and subplot*func
r.squaredGLMM(m1) #73% of variation explained by fixed effects, 77% by whole model
qqnorm(residuals(m1))
qqline(residuals(m1))
shapiro.test(residuals(m1))
#normally distributed, continue
LS1<-lsmeans(m1, ~treatment*subplot*func)
contrast(LS1, "pairwise")
#control rain ANPP is significantly greater than all the drought except spring dry, no surprise
#control rain is most similar to spring dry
#let's see it
ggplot(d=FG_2015, aes(x=treatment, y=weight_g_m, fill=func)) +
  theme_linedraw()+
  facet_wrap(~subplot)+
  labs(x="rainfall treatment", y="ANPP g/m2")+
  #annotate("text", x= c("consistentDry", "controlRain","fallDry","springDry"), y = c(900, 975, 900,900), label = c("a", "b", "ab", "ab"), color = "#999999") +
  geom_boxplot(aes(y=weight_g_m), shape=16)

FG_2015$treatment2 <- as.character(FG_2015$treatment)
#Then turn it back into a factor with the levels in the correct order
FG_2015$treatment2 <- factor(FG_2015$treatment2, levels = c("controlRain","fallDry", "springDry","consistentDry"))
ggplot(d=FG_2015, aes(x=treatment2, y=weight_g_m, fill=func)) +
  theme_linedraw()+
  labs(x="rainfall treatment", y="ANPP g/m2")+
  #annotate("text", x= c("controlRain","consistentDry", "fallDry","springDry"), y = c( 975, 900, 900,900), label = c("a", "b", "ab", "ab"), color = "#999999") +
  geom_boxplot(aes(y=weight_g_m), shape=16)

#let's do this again but test functional groups individually
##GRASS##
FG_2015_grass<-filter(FG_2015, func=="Grass", subplot!="F") 

m2<-lme(weight_g_m ~treatment, random=~1|shelterBlock/subplot, FG_2015_grass, na.action=na.exclude, contrasts=list(treatment=contr.treatment))
summary(m2)
anova(m2) #no subplot and no treatment effect
r.squaredGLMM(m2) #22% of variation explained by fixed effects, 22% by whole model (interannual variation?)
qqnorm(residuals(m2))
qqline(residuals(m2))
shapiro.test(residuals(m2))
#normally distributed, continue
LS2<-lsmeans(m2, ~treatment*subplot)
contrast(LS2, "pairwise")

FG_2015_grass$treatment2 <- as.character(FG_2015_grass$treatment)
#Then turn it back into a factor with the levels in the correct order
FG_2015_grass$treatment2 <- factor(FG_2015_grass$treatment2, levels = c("controlRain","springDry", "fallDry","consistentDry"))

ggplot(FG_2015_grass, aes(x = treatment2, y = weight_g_m)) + 
  facet_wrap(~subplot)+
  geom_boxplot(aes(fill = treatment2)) + 
  #scale_fill_manual(values = c("royalblue2","lightsteelblue1", "peachpuff", "sienna"), guide = guide_legend(title = "Treatment")) +
  geom_jitter(position=position_jitter(0.2), aes(shape=shelterBlock)) +
  #geom_jitter(position=position_jitter(0.2), aes(color=May_XC$shelterBlock, shape=as.factor(year))) +
  #scale_color_manual(values = c("gray80","gray50", "black"), guide = guide_legend(title = "Year")) +
  theme_bw(base_size = 14) +
  #facet_wrap(~year)+
  #annotate("text", x= c("controlRain","consistentDry", "fallDry","springDry"), y = c( 990, 900, 900,900), label = c("a", "b", "b", "ab"), color = "black") +
  xlab("Rainfall Treatment") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ylab("Grass ANPP g/m2")

##FORB##
#create subset with no species manipulations (control community) only
FG_2015_forb<-filter(FG_2015, func=='Forb', subplot!="G")

m3<-lme(weight_g_m ~treatment, random=~1|shelterBlock/subplot, FG_2015_forb, na.action=na.exclude)
summary(m3)
anova(m3)#no treatment effect
r.squaredGLMM(m3) #<1% of variation explained by fixed effects, 35% by whole model (interannual variation?)
qqnorm(residuals(m3))
qqline(residuals(m3))
shapiro.test(residuals(m3))
#normally distributed, continue
LS3<-lsmeans(m3, ~treatment*subplot)
contrast(LS3, "pairwise") #no differences

FG_2015_forb$treatment2 <- as.character(FG_2015_forb$treatment)
#Then turn it back into a factor with the levels in the correct order
FG_2015_forb$treatment2 <- factor(FG_2015_forb$treatment2, levels = c("controlRain","springDry", "fallDry","consistentDry"))

ggplot(FG_2015_forb, aes(x = treatment2, y = weight_g_m)) + 
  facet_wrap(~subplot)+
  geom_boxplot(aes(fill = treatment2)) + 
  #scale_fill_manual(values = c("royalblue2","lightsteelblue1", "peachpuff", "sienna"), guide = guide_legend(title = "Treatment")) +
  geom_jitter(position=position_jitter(0.2), aes(shape=shelterBlock)) +
  #geom_jitter(position=position_jitter(0.2), aes(color=May_XC$shelterBlock, shape=as.factor(year))) +
  #scale_color_manual(values = c("gray80","gray50", "black"), guide = guide_legend(title = "Year")) +
  theme_bw(base_size = 14) +
  #facet_wrap(~year)+
  #annotate("text", x= c("controlRain","consistentDry", "fallDry","springDry"), y = c( 175, 175, 175, 175), label = c("a", "a", "a", "a"), color = "black") +
  xlab("Rainfall Treatment") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ylab("Forb ANPP g/m2")


##N-FIXER##
FG_2015_nfix<-filter(FG_2015, func=="N-fixer", subplot!="G")

m4<-lme(log(weight_g_m+1) ~treatment, random=~1|subplot/shelterBlock, FG_2015_nfix, na.action=na.exclude)
summary(m4)
anova(m4) #treatment is significant
r.squaredGLMM(m4) #30% of variation explained by fixed effects, 39% by whole model (interannual variation?)
qqnorm(residuals(m4))
qqline(residuals(m4))
shapiro.test(residuals(m4))
#normally distributed, continue
LS4<-lsmeans(m4, ~treatment*subplot)
contrast(LS4, "pairwise") #consistent drought is different from fall and spring, but not control rain

FG_2015_nfix$treatment2 <- as.character(FG_2015_nfix$treatment)
#Then turn it back into a factor with the levels in the correct order
FG_2015_nfix$treatment2 <- factor(FG_2015_nfix$treatment2, levels = c("controlRain","springDry", "fallDry","consistentDry"))

ggplot(FG_2015_nfix, aes(x = treatment2, y = weight_g_m)) + 
  facet_wrap(~subplot)+
  geom_boxplot(aes(fill = treatment2)) + 
  #scale_fill_manual(values = c("royalblue2","lightsteelblue1", "peachpuff", "sienna"), guide = guide_legend(title = "Treatment")) +
  geom_jitter(position=position_jitter(0.2)) +
  #geom_jitter(position=position_jitter(0.2), aes(color=May_XC$shelterBlock, shape=as.factor(year))) +
 # scale_color_manual(values = c("gray80","gray50", "black"), guide = guide_legend(title = "Year")) +
  theme_bw(base_size = 14) +
  #facet_wrap(~year)+
 # annotate("text", x= c("controlRain","consistentDry", "fallDry","springDry"), y = c( 175, 175, 175, 175), label = c("ab", "b", "a", "a"), color = "black") +
  xlab("Rainfall Treatment") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ylab("N-Fixer ANPP g/m2")

