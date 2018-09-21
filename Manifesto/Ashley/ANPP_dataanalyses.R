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
#control rain ANPP is significantly greater than all the drought except spring dry, no surprise
#control rain is most similar to spring dry
#let's see it
ggplot(d=May_XC, aes(x=treatment, y=weight_g_m)) +
  theme_linedraw()+
  labs(x="rainfall treatment", y="ANPP g/m2")+
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
LS4<-lsmeans(m4, ~subplot)
contrast(LS4, "pairwise")
#ANOVA: overall sign effect of subplot on ANPP
#lsmeans: B (mixed plots) have greater ANPP compared to F (forb-only), but not G (grass-only)
#no sig difference in ANPP between B (mixed plots) and XC (no manipulation)
#plot the effect of species composition on ANPP:
ggplot(d=May_ANPP_noC, aes(x=subplot, y=weight_g_m)) +
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



#were the treatments effective? 
#summarise by treatment
smdat %>%
  group_by(treatment) %>%
  summarize(meansm = mean(sm, na.rm=T))
#looks like consistent dry has lowest sm

#calculate coefficient of variation for soil moisture 
CV <- function(x){(sd(x)/mean(x))*100}
moistCV<-aggregate(sm ~ treatment*shelterBlock*subplot, data= smdat, FUN = CV)

May_ANPP2<- merge(May_ANPP, moistCV)
colnames(May_ANPP2)[colnames(May_ANPP2)=="sm"] <- "sm_cv"

sm_mean<-aggregate(sm~treatment*shelterBlock*subplot, data=smdat, FUN=mean)
May_ANPP2<- merge(May_ANPP2, sm_mean)

May_ANPP2_noC<-filter(May_ANPP2, subplot!='C')

#regression
#how does ANPP relate to soil moisture?
m5<-lme(weight_g_m ~ sm, random=~1|shelterBlock/subplot, data=May_ANPP2_noC, na.action=na.exclude)
summary(m5)
anova(m5)
r.squaredGLMM(m5) #2% of variation explained by fixed effects, 7% by whole model 
qqnorm(residuals(m5))
qqline(residuals(m5))
shapiro.test(residuals(m5))
#normally distributed, continue


ggplot(May_ANPP2_noC, aes(x=sm, y=weight_g_m))+
  geom_smooth(method='lm')+
  labs(x="Soil Moisture", y="ANPP g/m2")
#ANPP increases with soil moisture, no surprise there

#does soil moisture vary by block?
m6<-lme(sm ~ shelterBlock, random=~1|treatment, data=May_ANPP2_noC, na.action=na.exclude)
summary(m6)
anova(m6) #treatment is significant
r.squaredGLMM(m6) #43% of variation explained by fixed effects, 60% by whole model 
qqnorm(residuals(m6))
qqline(residuals(m6))
shapiro.test(residuals(m6))
#not normally distributed
LS6<-lsmeans(m6, ~shelterBlock)
contrast(LS6, "pairwise")
#B is different from the rest

ggplot(May_ANPP2_noC, aes(x=shelterBlock, y=sm))+
  geom_boxplot(aes(y=sm))+
  labs(x="Block", y="Soil Moisture")
#block B has higher sm

#does treatment affect sm?
m7<-lme(sm ~ treatment, random=~1|shelterBlock, data=May_ANPP2_noC, na.action=na.exclude)
summary(m7)
anova(m7) #treatment is significant
r.squaredGLMM(m7) #43% of variation explained by fixed effects, 60% by whole model 
qqnorm(residuals(m7))
qqline(residuals(m7))
shapiro.test(residuals(m7))
#not normally distributed
LS7<-lsmeans(m7, ~treatment)
contrast(LS7, "pairwise")
#sm sign different between all combinations except falldry-springdry

ggplot(May_ANPP2_noC, aes(x=treatment, y=sm))+
  geom_boxplot(aes(y=sm))+
  labs(x="Treatment", y="Soil Moisture")
#control rain is highest sm, consistent dry is lowest

#see how treatments affect sm by block:
ggplot(May_ANPP2_noC, aes(x=treatment, y=sm))+
  geom_boxplot(aes(y=sm))+
  facet_wrap(~shelterBlock)+
  labs(x="Block", y="Soil Moisture")
#treatments not as effective in Block C?

#how does sm affect ANPP by treatment?
ggplot(May_ANPP2_noC, aes(x=sm, y=weight_g_m, color = treatment, group= (treatment)))+
  geom_smooth(method='lm')+
  geom_point()+
  labs(x="Soil Moisture", y="ANPP g/m2")

ggplot(May_ANPP2_noC, aes(x=sm, y=weight_g_m, color = subplot, group= (subplot)))+
  geom_smooth(method='lm')+
  geom_point()+
  labs(x="Soil Moisture", y="ANPP g/m2")

ggplot(May_ANPP2_noC, aes(x=sm, y=weight_g_m, color = subplot, group= (subplot)))+
  geom_smooth(method='lm')+
  facet_wrap(~treatment)+
  labs(x="Soil Moisture", y="ANPP g/m2") 

#how does ANPP relate to variation in soil moisture
ggplot(May_ANPP2_noC, aes(x=sm_cv, y=weight_g_m))+
  geom_smooth(method='lm')+
  labs(x="Soil Moisture CV", y="ANPP g/m2")

ggplot(May_ANPP2_noC, aes(x=sm_cv, y=weight_g_m, color = treatment, group= (treatment)))+
  geom_smooth(method='lm')+
  labs(x="Soil Moisture CV", y="ANPP g/m2")

ggplot(May_ANPP2_noC, aes(x=sm_cv, y=weight_g_m, color = subplot, group= (subplot)))+
  geom_smooth(method='lm')+
  labs(x="Soil Moisture CV", y="ANPP g/m2")
  
ggplot(May_ANPP2_noC, aes(x=sm_cv, y=weight_g_m, color = subplot, group= (subplot)))+
  geom_smooth(method='lm')+
  facet_wrap(~treatment)+
  labs(x="Soil Moisture CV", y="ANPP g/m2") 

#how does variation in ANPP relate to variation in soil moisture?
#calculate coefficient of variation for ANPP
anppCV<-aggregate(weight_g_m ~ treatment*shelterBlock*subplot, data= May_ANPP2_noC, FUN = CV)
colnames(anppCV)[colnames(anppCV)=="weight_g_m"] <- "ANPP_cv"
May_ANPP3_noC<- merge(May_ANPP2_noC, anppCV)

ggplot(May_ANPP3_noC, aes(x=sm_cv, y=ANPP_cv))+
  geom_smooth(method='lm')+
  labs(x="Soil Moisture CV", y="ANPP CV")

ggplot(May_ANPP3_noC, aes(x=sm_cv, y=ANPP_cv, color=treatment, group=treatment))+
  geom_smooth(method='lm')+
  labs(x="Soil Moisture CV", y="ANPP CV")

ggplot(May_ANPP3_noC, aes(x=sm_cv, y=ANPP_cv, color = subplot, group= (subplot)))+
  geom_smooth(method='lm') +
  labs(x="Soil Moisture CV", y="ANPP CV")

ggplot(May_ANPP3_noC, aes(x=sm_cv, y=ANPP_cv, color = subplot, group= (subplot)))+
  geom_smooth(method='lm')+
  facet_wrap(~treatment)+
  labs(x="Soil Moisture CV", y="ANPP CV") 
