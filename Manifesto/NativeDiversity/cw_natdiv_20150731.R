##R script for spring 2015 Native Diversity data analysis##
library(tidyr)
library(dplyr)
library(ggplot2)
library(nlme)
library(multcomp)
library(stringr)

##NOTE: first set working directory to folder containing Native Diversity data##

#DATA IMPORT AND PREPATION FOR ANPP DATA##
##Read in raw data (Native Diversity ANPP and shelter key) as table data frames## 
NatDivANPP <- read.csv("ClimVar_NativeDiv_ANPP_20150424.csv") %>%
  tbl_df() %>%
    select(Plot..1.16., Subplot..L..XC..F.G.B..C..D., Species, Weight..g.) %>%
    na.omit() 
    names(NatDivANPP) = c("plot", "subplot", "species", "weight")

##Convert subplot and species from factor to character and back to factor to correct number of factor levels in each (too many initally)#
NatDivANPP <- mutate(NatDivANPP, subplot = as.character(subplot), species = as.character(species))
NatDivANPP <- mutate(NatDivANPP, subplot = as.factor(subplot), species = as.factor(species))

key <- read.csv("Shelter_key.csv")


##Add columns to shelter key frame to indicate rainfall treatment as binary data, where 1 = receives rain and 0 = rain excluded##
key <- mutate(key, spring = ifelse(treatment == "fallDry" | treatment =="controlRain", 1,0), 
              fall = ifelse(treatment == "fallDry" | treatment == "consistentDry", 0,1)) %>%
    mutate(spring=as.factor(spring), fall = as.factor(fall)) %>%
    tbl_df()


##Merge shelter key and density ANPP data so in one complete table, rename fields##
togNDanpp <- merge(NatDivANPP, key) %>%
    mutate(spec_treat=paste(treatment, species, sep="_")) %>%
    mutate(spec_treat = as.factor(spec_treat))
names(togNDanpp)=c("plot", "subplot", "species", "weight", "rain_treat", "shelterBlock", "shelter", "spring", "fall", "spec_treat")


##DATA ANALYSIS FOR NATIVE DIVERSITY ANPP##

#Visualization of Native Diversity ANPP data##
pdf("Native_ANPP_rainfallresponse.pdf")
ggplot(togNDanpp, aes(x=treatment, y=weight)) + geom_boxplot() + facet_wrap(~species) + labs(title="Native ANPP by treatment")
dev.off()

#Attempt at ANOVA on Native Diversity ANPP, does not yield anything significant#
NDmodel <- lme(weight~ rain_treat*species, random=~1 |shelterBlock, data=togNDanpp,
              contrasts=list(rain_treat=contr.treatment, species=contr.treatment))
anova(NDmodel)
summary(glht(NDmodel, linfct=mcp(rain_treat="Tukey", species="Tukey")))

#ANOVA on Native Diversity ANPP, using created treatment_species interaction variable#
NDmodel2 <- lme(weight~ spec_treat, random=~1 |shelterBlock, data=togNDanpp,
                contrasts=list(spec_treat=contr.treatment))
summary(glht(NDmodel2, linfct=mcp(spec_treat="Tukey")))
