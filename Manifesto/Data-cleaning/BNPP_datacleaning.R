##R script for May 2015 BNPP data cleaning ##
## NOTE: This script sums the first 4 sampling intervals instead of projecting total ##
### LMH ##

library(tidyverse)
library(nlme)


#read in May 2015 BNPP .csv and shelter key#
rawBNPP <- read.csv("~/Dropbox/ClimVar/DATA/Plant_composition_data/BNPP/BNPP_EnteredData/ClimVar_BNPP_20150829.csv") %>%
  tbl_df()
names(rawBNPP)[8] = "plot"


sampleBNPP <- read.csv("~/Dropbox/ClimVar/DATA/Plant_composition_data/BNPP/BNPP_EnteredData/ClimVar_BNPP_20150520.csv") %>%
  tbl_df() %>%
  mutate(plot = Plot) %>%
  mutate(uniquePlot=paste(Plot, Subplot, sep="_")) %>%
  mutate(uniqueSample=paste(uniquePlot, Depth, sep="_"))
names(sampleBNPP)[14]="rmass"



key <- read.csv("~/Dropbox/ClimVar/DATA/Plant_composition_data/dd_Experiment/dd_EnteredData/Shelter_key.csv")%>%
  tbl_df()

## 
sampleExploreBNPP <- merge(sampleBNPP, key) %>%
  tbl_df() %>%
  select(plot, Subplot, Depth, Interval, rmass, Notes, treatment, shelterBlock, shelter)%>%
  tbl_df() %>%
  mutate(uniquePlot=paste(plot, Subplot, sep="_")) %>%
  mutate(uniqueSample=paste(uniquePlot, Depth, sep="_"))  %>%
  filter(rmass > 0) #to filter out missing values, incomplete dataset as of 8/28/15#
names(sampleExploreBNPP)=c("plot", "subplot", "depth", "interval", "rmass", "notes", "treatment", "shelterBlock", "shelter", "uniquePlot", "uniqueSample")


#Create exploratory dataset for graphing with rainfall treatment data added in#
exploreBNPP <- merge(rawBNPP, key)%>%
  select(plot, Subplot, Depth..cm., Interval, Root.weight..g., Notes, treatment, shelterBlock, shelter)%>%
  tbl_df()%>%
  mutate(uniquePlot=paste(plot, Subplot, sep="_")) %>%
  mutate(uniqueSample=paste(uniquePlot, Depth..cm., sep="_"))  %>%
  filter(Root.weight..g. > 0) #to filter out missing values, incomplete dataset as of 8/28/15#
names(exploreBNPP)=c("plot", "subplot", "depth", "interval", "rmass", "notes", "treatment", "shelterBlock", "shelter", "uniquePlot", "uniqueSample")

# Summarize the first 4 intervals
sumBNPP <- rbind(exploreBNPP, sampleExploreBNPP) %>%
  filter(interval < 4) %>%
  tbl_df() %>%
  group_by(plot, subplot, depth, treatment, shelterBlock, shelter) %>%
  summarize(totmass=sum(rmass)) %>%
  tbl_df() %>%
  # A coarse conversion for the times when depth could not be reached
  mutate(totmass = ifelse(depth == "10-19", totmass*1.11, totmass),
         totmass = ifelse(depth == "20-28", totmass*1.25, totmass),
         totmass = ifelse(depth == "0-6", totmass*(10/6), totmass),
         totmass = ifelse(depth == "10-14.5", totmass*(10/4.5), totmass),
         totmass = ifelse(depth == "20-23", totmass*(10/3), totmass)) %>%
  mutate(depth = as.character(depth),
         depth = ifelse(depth == "10-19", "10-20", depth),
         depth = ifelse(depth == "20-28", "20-30", depth),
         depth = ifelse(depth == "0-6", "0-10", depth),
         depth = ifelse(depth == "10-14.5", "10-20", depth),
         depth = ifelse(depth == "20-23", "20-30", depth)) %>%
  
  filter(depth=="0-10" | depth=="10-20"|depth=="20-30") %>%
  tbl_df() %>%
  mutate(fall=as.factor(ifelse(treatment=="springDry" | treatment== "controlRain", 1, 0))) %>%
  mutate(spring=as.factor(ifelse(treatment=="fallDry" | treatment == "controlRain", 1, 0))) %>%
  mutate(subplot=as.factor(subplot)) %>%
  #   mutate(forb = as.factor(ifelse(subplot == "G", 0, 1)),
  #          grass = as.factor(ifelse(subplot == "F", 0, 1))) %>%
  mutate(treatment=ordered(treatment, levels = c(consistentDry="consistentDry", springDry="springDry", fallDry="fallDry", controlRain="controlRain"))) %>%
  mutate(trtsub = as.factor(paste(treatment, subplot, sep = "_"))) %>%
  mutate(bmass_g_m2 = totmass * ((100*100)/(pi*(2.5^2)))) %>%
  select(-totmass)

#write.csv(sumBNPP, "BNPP_MayHarvest_2015.csv")




# Analysis musings: to clean up -------------------------------------------


ggplot(sumBNPP, aes(x=subplot, y=bmass_g_m2, color=spring)) + geom_boxplot() + facet_grid(depth~treatment)


ggplot(sumBNPP, aes(x=treatment, y=bmass_g_m2, color = subplot)) + geom_boxplot() + facet_wrap(~depth)

ggplot(sumBNPP, aes(x=treatment, y=bmass_g_m2)) + geom_boxplot() + facet_wrap(~depth)
m1<-lme(bmass_g_m2~treatment , random=~1|shelterBlock, 
        data = subset(sumBNPP, depth == "0-10" & (treatment == "controlRain" |
                                                    treatment == "consistentDry")),
        contrasts=list(treatment=contr.treatment)) #contrasts must be treatment for this
anova(m1)



ggplot(sumBNPP, aes(x=subplot, y=bmass_g_m2)) + geom_boxplot() + facet_wrap(~depth)
ggplot(sumBNPP, aes(x=treatment, y=bmass_g_m2)) + geom_boxplot() + facet_grid(depth~subplot)

ggplot(sumBNPP, aes(x=subplot, y=bmass_g_m2, color = fall)) + geom_boxplot() + facet_wrap(~depth)


library(nlme)
library(multcomp)
m1<-lme(bmass_g_m2~spring*fall  +subplot + depth, random=~1|shelterBlock, data = sumBNPP,
        contrasts=list(fall=contr.treatment, spring=contr.treatment, subplot = contr.treatment)) #contrasts must be treatment for this
anova(m1)

m1<-lme(bmass_g_m2~treatment*subplot + depth, random=~1|shelterBlock, data = sumBNPP,
        contrasts=list(treatment=contr.treatment, subplot = contr.treatment)) #contrasts must be treatment for this
anova(m1)


m1<-lme(bmass_g_m2~treatment*subplot, random=~1|shelterBlock, data = subset(sumBNPP, depth = "0-10"),
        contrasts=list(treatment=contr.treatment, subplot = contr.treatment)) #contrasts must be treatment for this
anova(m1)
summary(m1)
summary(glht(m1,linfct=mcp(subplot="Tukey", treatment="Tukey")), alternative="Bonferonni")


ggplot(subset(sumBNPP, depth = "0-10"), aes(x=subplot, y=bmass_g_m2)) + geom_boxplot()

sumBNPP %>%
  filter(depth == "0-10") %>%
  group_by(subplot) %>%
  summarize(bmass = mean(bmass_g_m2))

m1<-lme(bmass_g_m2~trtsub, random=~1|shelterBlock, data = subset(sumBNPP, depth = "0-10"),
        contrasts=list(trtsub=contr.treatment)) #contrasts must be treatment for this
anova(m1)

summary(glht(m1,linfct=mcp(trtsub="Tukey")), alternative="Bonferonni")


summary(glht(m1,linfct=mcp(subplot="Tukey")), alternative="Bonferonni")



m1<-lme(bmass_g_m2~trtsub + depth, random=~1|shelterBlock,data=sumBNPP,
        contrasts=list(trtsub = contr.treatment)) #contrasts must be treatment for this
anova(m1)
summary(glht(m1,linfct=mcp(trtsub="Tukey")), alternative="Bonferonni")


ggplot(sumBNPP, aes(x=subplot, y=totmass, color=treatment)) + geom_boxplot() + facet_grid(depth~treatment)
ggplot(sumBNPP, aes(x=treatment, y=totmass, color=treatment)) + geom_boxplot() + facet_grid(~subplot)
ggplot(sumBNPP, aes(x=depth, y=totmass, color=subplot)) + geom_boxplot() + facet_grid(~treatment)

ggplot(sumBNPP, aes(x=treatment, y=totmass)) + geom_boxplot() + facet_wrap(~depth)
ggplot(sumBNPP, aes(x=treatment, y=totmass, color=subplot)) + geom_boxplot() + facet_wrap(~depth)



#subset data for accumulation
dat<-exploreBNPP%>%
  select(uniqueSample, interval, rmass) %>%
  arrange(uniqueSample, interval)

#create dat2 with accumulated mass
dat2<-as.data.frame(cbind(uniqueSample=as.character(), interval=as.numeric(), rmass=as.numeric(),
                          allrmass=as.numeric()))

uniqueSamples<-unique(dat$uniqueSample)
for (i in 1: length(uniqueSamples)) {
  subber<-subset(dat, uniqueSample==uniqueSamples[i])
  subber$rmass2<-cumsum(subber$rmass)
  dat2<-rbind(dat2, subber)
}

