##R script for spring 2015 Native Diversity data analysis##
library(tidyr)
library(dplyr)
library(ggplot2)
library(nlme)
library(multcomp)
library(stringr)

##NOTE: first set working directory to folder containing Native Diversity data##

#DATA IMPORT AND PREPATION FOR NATIVE COVER DATA##
##Read in raw data (Native Diversity Cover) as table data frames## 
NatDivCov <- read.csv("ClimVar_NativeDiv_Cover_20150423.csv", skip=3) %>%
    tbl_df()

#extract notes from native cover data and remove "X" from plot number#
natcovnotes<-NatDivCov[1,] %>%
    gather(plot, notes, X1:X16) %>%
    select(-Plot) %>%
    mutate(plot=extract_numeric(plot))

#extract vegetation data
natdivveg<-NatDivCov[-c(1,2,3,4,5,6,23),] %>% 
    gather(plot, cover, X1:X16, na.rm=T) %>%
    mutate(plot=extract_numeric(plot))

#import keys
vegkey<-read.csv("~/Dropbox/ClimVar/DATA/Plant_composition_data/Cover/Cover_EnteredData/ClimVar_cover_specieskey.csv")
shelterkey<-read.csv("~/Dropbox/ClimVar/DATA/Plant_composition_data/Cover/Cover_EnteredData/Shelter_key.csv")

#merge data with keys
natdivveg0<-merge(natdivveg, vegkey)
NDvegtidy<-merge(natdivveg0, shelterkey) %>%
    tbl_df() %>%
    filter(cover!="") %>%
    mutate(cover=as.numeric(cover))
unique(NDvegtidy$cover)


##DATA ANALYSIS FOR NATIVE DIVERSITY COVER##

#Visualization of Native Diversity Cover data##
#simple histogram to see distribution#
hist(NDvegtidy$cover)

#native status x treatment#
ggplot(NDvegtidy, aes(x=status, y=cover)) + geom_boxplot() + facet_wrap(~treatment) + labs(title="Native subplot cover by native status and treatment")

#functional group x treatment#
ggplot(NDvegtidy, aes(x=func, y=cover)) + geom_boxplot() + facet_wrap(~treatment) + labs(title="Native subplot cover by functional group and treatment")

#by species, across all treatments#
ggplot(NDvegtidy, aes(x=genus, y=cover, fill=factor(status))) + geom_boxplot() + labs(title="Native subplot cover by genus")

#to see cover across the spatial gradient of plots west to east (since plots towards the west side seemed drier)#
ggplot(NDvegtidy, aes(x=status, y=cover, fill=factor(plot))) + geom_boxplot() + facet_wrap(~plot) + labs(title="Native subplot cover by species")