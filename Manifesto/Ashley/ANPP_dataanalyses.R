library(tidyverse)
library(nlme)
library(ggplot2)
library(multcomp)
library(lsmeans)
library(MuMIn)
setwd("~/Dropbox/ClimVar/DATA/")
May_ANPP<-read.csv("Plant_composition_data/ANPP/ANPP_CleanedData/ClimVar_ANPP-peak.csv")
levels(May_ANPP$plot)
levels(May_ANPP$year)
levels(May_ANPP$date)
levels(May_ANPP$subplot)
levels(May_ANPP$treatment)
levels(May_ANPP$shelter)
levels(May_ANPP$shelterBlock)
head(May_ANPP)

#change plots, years, shelter to factors
May_ANPP[,'plot'] <- as.factor(as.character(May_ANPP[,'plot']))
May_ANPP[,'year'] <- as.factor(as.character(May_ANPP[,'year']))
May_ANPP[,'shelter'] <- as.factor(as.character(May_ANPP[,'shelter']))


##1. Does variability of rainfall affect forage production (H1)?
#Using a mixed model to examine the effects of rainfall timing (fixed effect) with shelterblock nested within year as random effect

#create subset with no species manipulations (control community) only
May_XC<-filter(May_ANPP, subplot=='XC')

m1<-lme(weight_g_m ~treatment*year, random=~1|shelterBlock, May_XC, na.action=na.exclude)
summary(m1)
anova(m1)
r.squaredGLMM(m1) #12% of variation explained by fixed effects, 57% by whole model (interannual variation?)
qqnorm(residuals(m1))
qqline(residuals(m1))
shapiro.test(residuals(m1))
#normally distributed, continue
LS1<-lsmeans(m1, ~year)
contrast(LS1, "pairwise")
#control rain ANPP is significantly greater than all the drought except spring dry, no surprise
#control rain is most similar to spring dry
#let's see it
ggplot(d=May_XC, aes(x=treatment, y=weight_g_m)) +
  theme_linedraw()+
  labs(x="rainfall treatment", y="ANPP g/m2")+
  annotate("text", x= c("consistentDry", "controlRain","fallDry","springDry"), y = c(900, 975, 900,900), label = c("a", "b", "ab", "ab"), color = "#999999") +
  geom_boxplot(aes(y=weight_g_m), shape=16)

ggplot(d=May_XC, aes(x=year, y=weight_g_m)) +
  theme_linedraw()+
  labs(x="", y="ANPP g/m2")+
  annotate("text", x= c("2015", "2016","2017"), y = c(900, 975, 975), label = c("a", "b", "b"), color = "#999999") +
  geom_boxplot(aes(y=weight_g_m), shape=16)

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
#same plot as above but removes control rain
ggplot(d=May_XC_drought, aes(x=treatment, y=weight_g_m)) +
  theme_linedraw()+
  labs(x="rainfall treatment", y="ANPP g/m2")+
  geom_boxplot(aes(y=weight_g_m), shape=16)

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
LS4<-lsmeans(m4, ~treatment)
contrast(LS4, "pairwise")
#ANOVA: overall sign effect of subplot on ANPP
#lsmeans: B (mixed plots) have greater ANPP compared to F (forb-only), but not G (grass-only)
#no sig difference in ANPP between B (mixed plots) and XC (no manipulation)
#plot the effect of species composition on ANPP:
ggplot(d=May_ANPP_noC, aes(x=treatment, y=weight_g_m, fill=subplot)) +
  facet_wrap(~year)+
  theme_linedraw()+
  labs(x="composition treatment", y="ANPP g/m2")+
  geom_boxplot(aes(y=weight_g_m), shape=16)


#4. Does variability in soil moisture affect variability in ANPP (H1: variability)?
#first check interannual variability of ANPP for controls
ggplot(d=May_XC, aes(x=year, y=weight_g_m)) +
  theme_linedraw()+
  labs(x="year", y="ANPP g/m2")+
  geom_boxplot(aes(y=weight_g_m), shape=16)
#2015 has less ANPP than 2016 and 2017

#check block variability for controls
ggplot(d=May_XC, aes(x=shelterBlock, y=weight_g_m)) +
  theme_linedraw()+
  labs(x="Block", y="ANPP g/m2")+
  geom_boxplot(aes(y=weight_g_m), shape=16)
#lots of variability in ANPP by block

#bring in soil moisture data using Lauren's Decagon code
# create a header to pull in Andrews naming scheme
head <- read.csv("Decagon data/header_index.csv", header = F)[-3, -c(1:2)] %>%
  t() %>%
  tbl_df() 

names(head)=c("plot", "subplot")

head2 <- head %>%
  mutate(plot_subplot = paste("P", plot, subplot, sep = "_"))

head3 <- c("year", "doy", head2$plot_subplot)

# read in data
dat <- read.csv("Decagon data/Decagon_master_RESCALED.csv", header = F) 
names(dat) = head3

# format for normal use
dat2 <- dat %>%
  gather(plot, sm, P_1_1:P_16_5) %>%
  separate(plot, c("P", "plot", "spptrtNum")) %>%
  dplyr::select(-P) %>%
  mutate(sm = ifelse(sm == 99999, NA, sm)) %>%
  tbl_df() %>%
  mutate(doy2 = trunc(doy), plot = as.numeric(plot)) %>%
  group_by(year, plot, spptrtNum, doy2) %>%
  summarize(sm = mean(sm, na.rm=T))

#AS: updated code so spptreatment matches ANPP data, call it subplot
shelterkey <- read.csv("Decagon data/Shelter_key.csv")
#create a dummy dataframe that links treatments with a numeric code
dkey1<-expand.grid(seq(1,16,1), c("B", "C",  "F", "G", "XC"))
names(dkey1)=c("plot", "subplot")

dkey2 <- merge(dkey1, shelterkey[c("plot", "treatment", "shelterBlock")])%>%
  mutate(uniqueID=paste(plot, subplot, treatment, sep="_"))%>%
  mutate(plotNum=parse_number(plot), spptrtNum=as.numeric(subplot))%>%
  mutate(raintrtNum=1, raintrtNum=ifelse(treatment=="fallDry", 2, raintrtNum),
         raintrtNum=ifelse(treatment=="controlRain", 3, raintrtNum),
         raintrtNum=ifelse(treatment=="springDry", 4, raintrtNum)) %>%
  arrange(uniqueID) %>%
  tbl_df()

# merge with key
smdat <- merge(dkey2, dat2) %>%
  tbl_df() 


# pull out focal plots, an
#AS: updated subplot names
smdat2 <- smdat %>%
  # filter(subplot == "B" | subplot == "G" | subplot == "F" ) %>%
  # filter(subplot == "XC" | subplot == "C") %>%
  #  filter(subplot == "XC" | spptreatment == "B") %>%
  filter( subplot == "G" | subplot == "F" ) %>%
  tbl_df() %>%
  mutate(doy3= doy2/365,
         doy4 = year + doy3) 


smdat2 %>%
  group_by(treatment, subplot) %>%
  summarize(meansm = mean(sm, na.rm=T))

smdat2 %>%
  group_by(treatment, year) %>%
  summarize(meansm = mean(sm, na.rm=T))

ggplot(subset(smdat2), aes(x=doy4, y=sm,  color = subplot, group = interaction(subplot, plot))) +
  geom_line(size = .5) + theme_bw() + facet_wrap(~treatment) # facet_grid(shelterBlock~treatment) 

ggplot(subset(smdat2), aes(x=doy4, y=sm,  color = subplot)) + geom_line() + facet_wrap(~treatment)
ggplot(subset(smdat2, subplot != "both"), aes(x=doy4, y=sm,  color = treatment)) + geom_line() + facet_wrap(~subplot)

# aggregate across treatments
smdat3 <- smdat2 %>%
  group_by(subplot, treatment, year, doy4, doy3) %>%
  summarize(sm = mean(sm, na.rm=T))


ggplot(subset(smdat3), aes(x=doy4, y=sm,  color = subplot, group = (subplot))) +
  geom_line(size = .5) + theme_bw() + facet_wrap(~treatment) # facet_grid(shelterBlock~treatment) 

#create a plot showing sm data by treatment
ggplot(subset(smdat3), aes(x=doy4, y=sm,  color = treatment, group = (treatment))) +
  geom_line(size = .5) + theme_bw() # facet_wrap(~subplot) # facet_grid(shelterBlock~treatment) 

#were the treatments effective? 
#summarise by treatment
smdat %>%
  group_by(treatment) %>%
  summarize(meansm = mean(sm, na.rm=T))
#looks like consistent dry has lowest sm
ggplot(d=smdat, aes(x=treatment, y=sm)) +
  theme_linedraw()+
  labs(x="treatment", y="soil moisture")+
  geom_boxplot(aes(y=sm), shape=16)



#calculate coefficient of variation for soil moisture 
CV <- function(x){(sd(x)/mean(x))*100}
moistCV<-aggregate(sm ~ treatment*shelterBlock*subplot*year, data= smdat, FUN = CV)

colnames(moistCV)[colnames(moistCV)=="sm"] <- "sm_cv"
May_ANPP2<-May_ANPP %>% group_by(treatment,shelterBlock,subplot,year)%>%summarise(meanANPP=mean(weight_g_m))

May_ANPP2<- merge(May_ANPP2, moistCV)

ggplot(d=May_ANPP2, aes(x=treatment, y=sm_cv)) +
  theme_linedraw()+
  labs(x="treatment", y="soil moisture CV")+
  geom_boxplot(aes(y=sm_cv), shape=16)

sm_mean<-aggregate(sm~treatment*shelterBlock*subplot*year, data=smdat, FUN=mean)
May_ANPP2<- merge(May_ANPP2, sm_mean)

May_ANPP2_noC<-filter(May_ANPP2, subplot!='C')

#how does sm vary by treatment?
ggplot(d=May_ANPP2_noC, aes(x=treatment, y=sm)) +
  theme_linedraw()+
  labs(x="treatment", y="soil moisture")+
  geom_boxplot(aes(y=sm), shape=16)

#regression
#how does ANPP relate to soil moisture?
m5<-lme(meanANPP ~ sm, random=~1|shelterBlock/subplot, data=May_ANPP2_noC, na.action=na.exclude)
summary(m5)
anova(m5)
r.squaredGLMM(m5) #26% of variation explained by fixed effects, 37% by whole model 
qqnorm(residuals(m5))
qqline(residuals(m5))
shapiro.test(residuals(m5))
#normally distributed, continue


ggplot(May_ANPP2_noC, aes(x=sm, y=meanANPP))+
  geom_smooth(method='lm')+
  labs(x="Soil Moisture", y="ANPP g/m2")
#ANPP increases with soil moisture, no surprise there

#Does soil moisture vary by block?
m6<-lme(sm ~ shelterBlock, random=~1|treatment, data=May_ANPP2_noC, na.action=na.exclude)
summary(m6)
anova(m6) #block not significant
r.squaredGLMM(m6) #3% of variation explained by fixed effects, 12% by whole model 
qqnorm(residuals(m6))
qqline(residuals(m6))
shapiro.test(residuals(m6))#not normal
LS6<-lsmeans(m6, ~shelterBlock)
contrast(LS6, "pairwise")

ggplot(May_ANPP2_noC, aes(x=shelterBlock, y=sm))+
  geom_boxplot(aes(y=sm))+
  theme_linedraw()+
  labs(x="Block", y="Soil Moisture")
#block B has higher sm

#Does treatment affect sm?
m7<-lme(sm ~ treatment, random=~1|shelterBlock, data=May_ANPP2_noC, na.action=na.exclude)
summary(m7)
anova(m7) #treatment is significant
r.squaredGLMM(m7) #8% of variation explained by fixed effects, 9% by whole model 
qqnorm(residuals(m7))
qqline(residuals(m7))
shapiro.test(residuals(m7))
#normal

LS7<-lsmeans(m7, ~treatment)
contrast(LS7, "pairwise")
#sm sign different between all consistent dry and control rain

ggplot(May_ANPP2_noC, aes(x=treatment, y=sm))+
  geom_boxplot(aes(y=sm))+
  labs(x="Treatment", y="Soil Moisture")
#control rain is highest sm, consistent dry is lowest

#see how treatments affect sm by block:
ggplot(May_ANPP2_noC, aes(x=treatment, y=sm))+
  geom_boxplot(aes(y=sm))+
  facet_wrap(~shelterBlock)+
  labs(x="Block", y="Soil Moisture")
#treatments look consistent across blocks, maybetreatments not as effective in Block C?

#how does sm affect ANPP by treatment?
ggplot(May_ANPP2_noC, aes(x=sm, y=meanANPP, color = treatment, group= (treatment)))+
  geom_smooth(method='lm')+
  geom_point()+
  labs(x="Soil Moisture", y="ANPP g/m2")

ggplot(May_ANPP2_noC, aes(x=sm, y=meanANPP, color = subplot, group= (subplot)))+
  geom_smooth(method='lm')+
  geom_point()+
  labs(x="Soil Moisture", y="ANPP g/m2")

ggplot(May_ANPP2_noC, aes(x=sm, y=meanANPP, color = subplot, group= (subplot)))+
  geom_smooth(method='lm')+
  facet_wrap(~treatment)+
  labs(x="Soil Moisture", y="ANPP g/m2") 

#how does ANPP relate to variation in soil moisture
ggplot(May_ANPP2_noC, aes(x=sm_cv, y=meanANPP))+
  geom_smooth(method='lm')+
  labs(x="Soil Moisture CV", y="ANPP g/m2")

ggplot(May_ANPP2_noC, aes(x=sm_cv, y=meanANPP, color = treatment, group= (treatment)))+
  geom_smooth(method='lm')+
  labs(x="Soil Moisture CV", y="ANPP g/m2")

ggplot(May_ANPP2_noC, aes(x=sm_cv, y=meanANPP, color = subplot, group= (subplot)))+
  geom_smooth(method='lm')+
  labs(x="Soil Moisture CV", y="ANPP g/m2")
  
ggplot(May_ANPP2_noC, aes(x=sm_cv, y=meanANPP, color = subplot, group= (subplot)))+
  geom_smooth(method='lm')+
  facet_wrap(~treatment)+
  labs(x="Soil Moisture CV", y="ANPP g/m2") 

#how does variation in ANPP relate to variation in soil moisture?
#calculate coefficient of variation for ANPP
anppCV<-aggregate(weight_g_m ~ treatment*shelterBlock*subplot*year, data= May_ANPP, FUN = CV)
colnames(anppCV)[colnames(anppCV)=="weight_g_m"] <- "ANPP_cv"
anppCV_noC<-filter(anppCV, subplot!='C')
May_ANPP3_noC<- merge(May_ANPP2_noC, anppCV_noC)

ggplot(May_ANPP3_noC, aes(x=sm_cv, y=ANPP_cv))+
  geom_smooth(method='lm')+
  labs(x="Soil Moisture CV", y="ANPP CV")

ggplot(May_ANPP3_noC, aes(x=sm_cv, y=ANPP_cv, color=treatment, group=treatment))+
  geom_smooth(method='lm')+
  labs(x="Soil Moisture CV", y="ANPP CV")

ggplot(May_ANPP3_noC, aes(x=subplot, y=ANPP_cv, color=subplot, group=subplot))+
  geom_boxplot()+
  facet_grid(~treatment)+
  labs(x="treatment", y="ANPP CV")

ggplot(May_ANPP3_noC, aes(x=sm_cv, y=ANPP_cv, color = subplot, group= (subplot)))+
  geom_smooth(method='lm') +
  labs(x="Soil Moisture CV", y="ANPP CV")

ggplot(May_ANPP3_noC, aes(x=sm_cv, y=ANPP_cv, color = subplot, group= (subplot)))+
  geom_smooth(method='lm')+
  facet_wrap(~treatment)+
  labs(x="Soil Moisture CV", y="ANPP CV") 

###############################################
#############################################
######FUNCTIONAL GROUP ANPP################

Func_ANPP<-read.csv("Plant_composition_data/ANPP/ANPP_CleanedData/ClimVar_ANPP-func-groups.csv")
levels(Func_ANPP$plot)
levels(Func_ANPP$year)
levels(Func_ANPP$date)
levels(Func_ANPP$subplot)
levels(Func_ANPP$treatment)
levels(Func_ANPP$shelter)
levels(Func_ANPP$shelterBlock)
levels(Func_ANPP$func)
#need to fix levels of func


library(forcats)
Func_ANPP1<-Func_ANPP %>%
  mutate(func = fct_collapse(func, G=c("Grass", "G"), F=c("Forb", "F"), N =c("N-fixer","N")))
levels(Func_ANPP1$func)

#change plots, years, shelter to factors
Func_ANPP1[,'plot'] <- as.factor(as.character(Func_ANPP1[,'plot']))
Func_ANPP1[,'year'] <- as.factor(as.character(Func_ANPP1[,'year']))
Func_ANPP1[,'shelter'] <- as.factor(as.character(Func_ANPP1[,'shelter']))

#create a subset, second harvest only
Func_ANPP2<-filter(Func_ANPP1, harvest=='Second')

#create another subset, XC only 
Func_ANPP2_XC <-filter(Func_ANPP2, subplot=='XC')

#Does the biomass of functional groups change with drought treatment?
m8<-lme(weight_g_m ~func*treatment, random=~1|year/shelterBlock, Func_ANPP2_XC, na.action=na.exclude)
summary(m8)
anova(m8)
r.squaredGLMM(m8) #86% of variation explained by fixed effects, 88% by whole model 
qqnorm(residuals(m8))
qqline(residuals(m8))
shapiro.test(residuals(m8))
#def not normal, try log transformation
m9<-lme(log(weight_g_m+1) ~func*treatment, random=~1|year/shelterBlock, Func_ANPP2_XC, na.action=na.exclude)
summary(m9)
anova(m9)
r.squaredGLMM(m9) #73% of variation explained by fixed effects, 77% by whole model 
qqnorm(residuals(m9))
qqline(residuals(m9))
shapiro.test(residuals(m9))
#not great, but better
LS9<-lsmeans(m9, ~treatment*func)
contrast(LS9, "pairwise")

#total biomass by functional group 
ggplot(Func_ANPP2_XC, aes(x=func, y=weight_g_m, color = func, group= (func)))+
geom_boxplot()+
  facet_grid(~treatment)+
  labs(x="Functional Group", y="ANPP g/m2") 

ggplot(Func_ANPP2_XC, aes(x=func, y=weight_g_m, color = func, group= (func)))+
  geom_boxplot()+
  facet_grid(~treatment*year)+
  labs(x="Functional Group", y="ANPP g/m2")

#relative biomass by functional group
m10<-lme(propweight ~func*treatment, random=~1|year/shelterBlock, Func_ANPP2_XC, na.action=na.exclude)
summary(m10)
anova(m10)
r.squaredGLMM(m10) #94% of variation explained by fixed effects, 94% by whole model 
qqnorm(residuals(m10))
qqline(residuals(m10))
shapiro.test(residuals(m10))
#not normal, log transform
m11<-lme(log(propweight+1) ~func*treatment, random=~1|year/shelterBlock, Func_ANPP2_XC, na.action=na.exclude)
summary(m11)
anova(m11)
r.squaredGLMM(m11) #95% of variation explained by fixed effects, 95% by whole model 
qqnorm(residuals(m11))
qqline(residuals(m11))
shapiro.test(residuals(m11))
#no better, try sqrt
m12<-lme(sqrt(propweight) ~func*treatment, random=~1|year/shelterBlock, Func_ANPP2_XC, na.action=na.exclude)
summary(m12)
anova(m12)
r.squaredGLMM(m12) #95% of variation explained by fixed effects, 95% by whole model 
qqnorm(residuals(m12))
qqline(residuals(m12))
shapiro.test(residuals(m12))
LS12<-lsmeans(m12, ~treatment*func)
contrast(LS12, "pairwise")

ggplot(Func_ANPP2_XC, aes(x=func, y=propweight, color = func, group= (func)))+
  geom_boxplot()+
  facet_grid(~treatment*year)+
  labs(x="Functional Group", y="Relative ANPP")

#total biomass by functional group per block
ggplot(Func_ANPP2_XC, aes(x=func, y=weight_g_m, color = func, group= (func)))+
  geom_boxplot()+
  facet_grid(~shelterBlock)+
  labs(x="Functional Group", y="ANPP g/m2") 

#relative biomass by functional group per block
ggplot(Func_ANPP2_XC, aes(x=func, y=propweight, color = func, group= (func)))+
  geom_boxplot()+
  facet_grid(~shelterBlock)+
  labs(x="Functional Group", y="Relative ANPP")
#site B is grassier? A is 'forbier'?

#does grass ANPP differ by treatment?
Func_ANPP2_G <- Func_ANPP2 %>% filter(func=="G") %>% filter(subplot!="F")
m14<-lme(weight_g_m ~treatment, random=~1|year/shelterBlock, Func_ANPP2_G, na.action=na.exclude)
summary(m14)
anova(m14) #no treatment effect on grass ANPP
r.squaredGLMM(m14) #2% of variation explained by fixed effects, 42% by whole model 
qqnorm(residuals(m14))
qqline(residuals(m14))
shapiro.test(residuals(m14))#normal
LS14<-lsmeans(m14, ~treatment)
contrast(LS14, "pairwise")

ggplot(Func_ANPP2_G, aes(x=year, y=weight_g_m, color=treatment))+
  geom_boxplot()+
  labs(x="year",y="Grass ANPP")

Func_ANPP2_F <- Func_ANPP2 %>% filter(func=="F") %>% filter(subplot!="G")
ggplot(Func_ANPP2_F, aes(x=year, y=weight_g_m, color=treatment))+
  geom_boxplot()+
  labs(x="year",y="Forb ANPP")


# compensation (H2)
#H2 will be confirmed if mixed plots have greater ANPP compared to grass only or forb only
#Does forb biomass increase with fall drought?
Func_ANPP3 <-filter(Func_ANPP2, treatment!='controlRain')
Func_ANPP3_F <-filter(Func_ANPP3, func=='F')
m13<-lme(weight_g_m ~treatment*subplot, random=~1|year/shelterBlock, data=Func_ANPP3_F, na.action=na.exclude)
summary(m13)
anova(m13)
r.squaredGLMM(m13) #25% of variation explained by fixed effects
qqnorm(residuals(m13))
qqline(residuals(m13))
shapiro.test(residuals(m13))
#not normally distributed, try log transform
m14<-lme(log(weight_g_m+1) ~treatment*subplot, random=~1|year/shelterBlock, Func_ANPP3_F, na.action=na.exclude)
summary(m14)
anova(m14)#significant effects of subplot
r.squaredGLMM(m14)#15% of variation explained by fixed effects, 55% explained by entire model 
qqnorm(residuals(m14))
qqline(residuals(m14))
shapiro.test(residuals(m14))
#better
LS14<-lsmeans(m14, ~subplot)
contrast(LS14, "pairwise")
#ANOVA: overall sign effect of subplot on Forb ANPP
#no shock there... there  is sign. greater forb biomass in the forb treatment

#Does forb biomass increase with fall drought?
#try again with XC only 
Func_ANPP4_F <-filter(Func_ANPP3_F, subplot=='XC')
m15<-lme(weight_g_m ~treatment, random=~1|year/shelterBlock, data=Func_ANPP4_F, na.action=na.exclude)
summary(m15)
anova(m15)
r.squaredGLMM(m15)
qqnorm(residuals(m15))
qqline(residuals(m15))
shapiro.test(residuals(m15))
#not normally distributed, try log transform
m16<-lme(log(weight_g_m+1) ~treatment, random=~1|year/shelterBlock, Func_ANPP4_F, na.action=na.exclude)
summary(m16)
anova(m16)
r.squaredGLMM(m16)#<1% of variation explained by fixed effects, 23% explained by entire model 
qqnorm(residuals(m16))
qqline(residuals(m16))
shapiro.test(residuals(m16))
#normal
LS16<-lsmeans(m16, ~treatment)
contrast(LS16, "pairwise")

ggplot(Func_ANPP4_F, aes(x=treatment, y=weight_g_m, color = treatment, group= (treatment)))+
  geom_boxplot()+
  facet_grid(~year)+
  labs(x="Treatment", y="Forb ANPP")

#try with proportion
m17<-lme(propweight ~treatment, random=~1|year/shelterBlock, data=Func_ANPP4_F, na.action=na.exclude)
summary(m17)
anova(m17)
r.squaredGLMM(m17)
qqnorm(residuals(m17))
qqline(residuals(m17))
shapiro.test(residuals(m17))
#not normally distributed, try log transform
m18<-lme(log(propweight+1) ~treatment, random=~1|year/shelterBlock, Func_ANPP4_F, na.action=na.exclude)
summary(m18)
anova(m18)
r.squaredGLMM(m18)#2% of variation explained by fixed effects, 42% explained by entire model
qqnorm(residuals(m18))
qqline(residuals(m18))
shapiro.test(residuals(m18))
#not any better

ggplot(Func_ANPP4_F, aes(x=treatment, y=propweight, color = treatment, group= (treatment)))+
  geom_boxplot()+
  facet_grid(~year)+
  labs(x="Treatment", y="Prop Forb ANPP")

#does forb ANPP increase in the dry year (2015) for XC only?
mf1<-lme(weight_g_m ~year, random=~1|treatment/shelterBlock, data=Func_ANPP4_F, na.action=na.exclude)
summary(mf1)
anova(mf1)
r.squaredGLMM(mf1)
qqnorm(residuals(mf1))
qqline(residuals(mf1))
shapiro.test(residuals(mf1))
#not normally distributed, try log transform
mf2<-lme(log(weight_g_m+1) ~year, random=~1|treatment/shelterBlock, Func_ANPP4_F, na.action=na.exclude)
summary(mf2)
anova(mf2)#year highly significant
r.squaredGLMM(mf2)#21% of variation explained by fixed effects, 64% explained by entire model
qqnorm(residuals(mf2))
qqline(residuals(mf2))
shapiro.test(residuals(mf2))
#better
LSf2<-lsmeans(mf2, ~year)
contrast(LSf2, "pairwise")

ggplot(Func_ANPP4_F, aes(x=year,y=weight_g_m))+
  geom_boxplot()+
  labs(x="Year", y="Forb ANPP")

#does interaction of year*treatment affect forb ANPP?
mf3<-lme(weight_g_m ~year*treatment, random=~1|shelterBlock, data=Func_ANPP4_F, na.action=na.exclude)
summary(mf3)
anova(mf3)
r.squaredGLMM(mf3) #13% variation explained by fixed effects
qqnorm(residuals(mf3))
qqline(residuals(mf3))
shapiro.test(residuals(mf3))
#not normally distributed, try log transform
mf4<-lme(log(weight_g_m+1) ~year*treatment, random=~1|shelterBlock, Func_ANPP4_F, na.action=na.exclude)
summary(mf4)
anova(mf4)#year significant
r.squaredGLMM(mf4)#24% of variation explained by fixed effects, 38% explained by entire model
qqnorm(residuals(mf4))
qqline(residuals(mf4))
shapiro.test(residuals(mf4))
#better
LSf4<-lsmeans(mf4, ~year*treatment)
contrast(LSf4, "pairwise")

###################################################################
#now let's look at soil moisture effects on functional group ANPP
Func_ANPP2<-Func_ANPP1 %>% group_by(func,treatment,shelterBlock,subplot)%>%summarise(meanANPP=mean(weight_g_m)) 
Func_ANPP2<- merge(Func_ANPP2, moistCV)
Func_ANPP2<- merge(Func_ANPP2, sm_mean)

FanppCV<-aggregate(weight_g_m ~ treatment*shelterBlock*subplot*func, data= Func_ANPP1, FUN = CV)
colnames(FanppCV)[colnames(FanppCV)=="weight_g_m"] <- "ANPP_cv"
Func_ANPP4<- merge(Func_ANPP2, FanppCV)


#how does ANPP of functional groups relate to soil moisture?
MeanFunc<-Func_ANPP1 %>% group_by(func,treatment)%>%summarise(meanANPP=mean(weight_g_m))
SE <- function(x){(sd(x)/sqrt(length(x)))}
ANPPse<-aggregate(weight_g_m ~ func*treatment, data= Func_ANPP1, FUN = SE)
colnames(ANPPse)[colnames(ANPPse)=="weight_g_m"] <- "ANPP.SE"
MeanFunc<- merge(MeanFunc, ANPPse)
FanppCV2<-aggregate(weight_g_m ~ treatment*func, data= Func_ANPP1, FUN = CV)
colnames(FanppCV2)[colnames(FanppCV2)=="weight_g_m"] <- "ANPP_cv"
MeanFunc<- merge(MeanFunc, FanppCV2)
sm_mean2<-aggregate(sm~treatment, data=smdat, FUN=mean)
smCV2<-aggregate(sm ~ treatment, data=smdat, FUN = CV)
colnames(smCV2)[colnames(smCV2)=="sm"] <- "sm_cv"
sm2<- merge(sm_mean2, smCV2)
MeanFunc<- merge(MeanFunc, sm2)

pd <- position_dodge(0.1) # move points to the left and right

ggplot(MeanFunc, aes(x=sm_cv, y=meanANPP, color=func, shape=treatment))+
  geom_errorbar(aes(ymin=meanANPP-ANPP.SE, ymax=meanANPP+ANPP.SE), width=1, position=pd)+
  geom_point(position=pd)+
  xlim(55,65)+
  labs(x="Soil Moisture CV", y="ANPP g/m2")

ggplot(MeanFunc, aes(x=sm_cv, y=meanANPP, color=func))+
  geom_errorbar(aes(ymin=meanANPP-ANPP.SE, ymax=meanANPP+ANPP.SE), width=1, position=pd)+
  geom_point(position=pd)+
  geom_smooth(data=MeanFunc, method='lm', se=FALSE)+
  xlim(55,65)+
  labs(x="Soil Moisture CV", y="ANPP g/m2")

ggplot(MeanFunc, aes(x=sm_cv, y=ANPP_cv, color=func, shape=treatment))+
  #geom_errorbar(aes(ymin=meanANPP-ANPP.SE, ymax=meanANPP+ANPP.SE), width=1, position=pd)+
  geom_point(position=pd)+
  geom_smooth(method='lm')+
  labs(x="Soil Moisture CV", y="ANPP CV")

ggplot(Func_ANPP4, aes(x=sm_cv, y=ANPP_cv, color=func, group=func))+
  geom_point()+
  geom_smooth(method='lm')+
  labs(x="Soil Moisture CV", y="ANPP CV")

ggplot(Func_ANPP4, aes(x=sm_cv, y=meanANPP, color=func, group=func))+
  geom_point()+
  geom_smooth(method='lm')+
  labs(x="Soil Moisture CV", y="mean ANPP")

ggplot(Func_ANPP4, aes(x=sm, y=meanANPP, color=func, group=func))+
  geom_point()+
  geom_smooth(method='lm')+
  labs(x="mean soil moisture", y="mean ANPP")

ggplot(Func_ANPP4, aes(x=sm_cv, y=ANPP_cv, color=func, group=func))+
  geom_point()+
  geom_smooth(method='lm')+
  facet_wrap(~treatment)+
  labs(x="Soil Moisture CV", y="ANPP CV")

