#the purpose of this script is to clean up soil extracts data 
library(tidyverse)
library(stringr)
setwd("~/Dropbox/ClimVar/DATA/")

## Import soil extracts
soil <- read.csv("~/Dropbox/ClimVar/DATA/Soil Extracts/EXTRACTS_out_HEADERS.csv", na.strings="99999")
soil$Treatment<-as.factor(soil$Treatment)

#extract notes & remove notes from df
soil_notes <- soil[,14:15] 
soil<-soil %>% dplyr::select(-(12:19)) %>% 
  dplyr::rename(plot=Plot) %>% dplyr::rename(subplot = Treatment)%>%
  mutate(subplot = recode(subplot, "1"= "XC", "2"="C", "3"="L", "4"="G", "5"="F", "6"="B"))

## Import import shelter key
shelterkey <- read.csv("~/Dropbox/ClimVar/DATA/Plant_composition_data/Cover/Cover_EnteredData/Shelter_key.csv")

#merge data with shelter key
soil <- merge(soil, shelterkey) %>%
  tbl_df() 

write.csv(soil, "~/Dropbox/ClimVar/DATA/Soil Extracts/Soil_CleanedData/ClimVar_soil_2015.csv")






