#the purpose of this script is to extract and clean up cover data for litter, bare and rock
#this information is left out of the clean veg cover data
library(tidyverse)
library(stringr)

#########################
### 2015 COVER DATA #####
#########################

## Import import shelter key
shelterkey <- read.csv("~/Dropbox/ClimVar/DATA/Plant_composition_data/Cover/Cover_EnteredData/Shelter_key.csv")

### Read in 2015 cover data 
covdat_2015_0 <- read.csv("~/Dropbox/ClimVar/DATA/Plant_composition_data/Cover/Cover_EnteredData/ClimVar_Cover_20150508.csv", skip=3) %>%
  tbl_df()

covdat_compost_2015 <- read.csv("~/Dropbox/ClimVar/DATA/Plant_composition_data/Cover/Cover_EnteredData/ClimVar_CompostCover_20150611.csv", skip = 3) %>%
  tbl_df()

covdat_2015 <- left_join(covdat_2015_0, covdat_compost_2015) %>%
  tbl_df()

#extract notes
covdat_notes_2015 <- covdat_2015[1,] %>%
  gather(plot, notes, X1B:X16C) %>%
  dplyr::select(-Plot)

#extract litter data
litter_2015 <- covdat_2015[-c(1:2,8:87),] %>%
  gather(plotname, cover, X1B:X16C, na.rm=T) %>%
  mutate(plotname=substr(plotname, 2,5)) %>%
  mutate(plot=extract_numeric(plotname)) %>%
  mutate(subplot=substr(plotname,2,4), subplot=ifelse(plot > 9, substr(plotname, 3,4), subplot))

litter_2015$Plot<-as.factor(litter_2015$Plot)
levels(litter_2015$Plot)<-c("Percent_bare","Percent_forb","Percent_grass", "Percent_litter", "Litter_depth_cm" )
levels(litter_2015$Plot)

#merge data with shelter key
litter_2015 <- merge(litter_2015, shelterkey) %>%
  tbl_df() %>%
  filter(cover!="") %>%
  mutate(cover=as.numeric(cover),
         year = 2015)

sort(unique(litter_2015$Plot))
rm(covveg_2015, covdat_notes_2015, covdat_2015, covdat_2015_0, covdat_compost_2015)

#########################
### 2016 COVER DATA #####
#########################

covdat_2016 <- read.csv("~/Dropbox/ClimVar/DATA/Plant_composition_data/Cover/Cover_EnteredData/ClimVar_Cover_20160426.csv", skip=3) %>%
  tbl_df()

#extract notes
covdat_notes_2016 <- covdat_2016[1,] %>%
  gather(plot, notes) %>%
  filter(plot != "Plot") %>%
  tbl_df()

#extract litter data
litter_2016 <- covdat_2016[-c(1,7:74),] %>% 
  gather(plotname, cover, X1B:X16XC, na.rm=T) %>%
  mutate(plotname=substr(plotname, 2,5)) %>%
  mutate(plot=extract_numeric(plotname)) %>%
  mutate(subplot=substr(plotname,2,4), subplot=ifelse(plot > 9, substr(plotname, 3,4), subplot)) %>%
  mutate(Plot = as.character(Plot))

# merge data with shelter key
litter_2016 <- merge(litter_2016, shelterkey, all.x = T) %>%
  tbl_df() %>%
  filter(cover!="") %>%
  mutate(cover=as.numeric(cover)) %>%
  mutate(year = 2016)

sort(unique((litter_2016$Plot))) 

rm(covveg_2016, covdat_notes_2016, covdat_2016)

#########################
### 2017 COVER DATA #####
#########################

# import data
covdat_2017 <- read.csv("~/Dropbox/ClimVar/DATA/Plant_composition_data/Cover/Cover_EnteredData/ClimVar_Cover_20170526.csv") %>%
  tbl_df()

#no notes for 2017
#extract litter data
litter_2017 <- covdat_2017[-c(1,7:70),] %>% 
  gather(plotname, cover, X1F:X16F, na.rm=T) %>%
  mutate(plotname=substr(plotname, 2,5)) %>%
  mutate(plot=extract_numeric(plotname)) %>%
  mutate(subplot=substr(plotname,2,4), subplot=ifelse(plot > 9, substr(plotname, 3,4), subplot)) %>%
  mutate(Plot = as.character(Plot))

# check that all variables were retained
sort(unique(litter_2017$Plot))

# merge data with shelter key
litter_2017 <- merge(litter_2017, shelterkey, all.x = T) %>%
  tbl_df() %>%
  filter(cover!="") %>%
  mutate(cover=as.numeric(cover)) %>%
  mutate(year = 2017)


rm(covdat_2017)


##################################
### COMBINE YEARS, ADD ZEROS #####
##################################

#remove the legume treatment, only 1 year of data
dat.cover <- rbind(litter_2015, litter_2016, litter_2017) %>%
  tbl_df() %>%
  dplyr::select(-plotname) %>%
  dplyr::rename(type=Plot) %>%
  filter(subplot!="L")

vegkey <- dat.cover %>%
  dplyr::select(type) %>%
  unique()

# check for consistent names across years
sort(unique(dat.cover$type))

# add in 0s for absent species
zerokey <- expand.grid(unique(dat.cover$type), unique(dat.cover$plot), unique(dat.cover$subplot), unique(dat.cover$year))

names(zerokey) =c("type", "plot", "subplot", "year")

fullmat1 <- merge(zerokey, vegkey, all.x = T) 
fullmat2 <- merge(fullmat1, shelterkey, all.x = T)  


dat.cover_with0 <- merge(fullmat2, dat.cover, all.x = T) %>%
  mutate(cover = ifelse(is.na(cover), 0, cover)) %>%
  tbl_df() 
  

rm(fullmat1, fullmat2, vegkey, shelterkey, zerokey, veg_2015, veg_2016, veg_2017)

write.csv(dat.cover_with0, "~/Dropbox/ClimVar/DATA/Plant_composition_data/Cover/Cover_CleanedData/ClimVar_litter-cover.csv")






