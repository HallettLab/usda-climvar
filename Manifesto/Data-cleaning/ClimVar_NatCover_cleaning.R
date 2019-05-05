##R script for spring 2015 Native Diversity data analysis##
library(tidyr)
library(dplyr)
library(ggplot2)
library(nlme)
library(multcomp)
library(stringr)

##NOTE: first set working directory to folder containing Native Diversity data##

shelterkey <- read.csv("~/Dropbox/ClimVar/DATA/Plant_composition_data/Cover/Cover_EnteredData/Shelter_key.csv")

#DATA IMPORT AND PREPATION FOR NATIVE COVER DATA ##

################ 2015 #################

##Read in raw data (Native Diversity Cover) as table data frames## 
NatDivCov_2015 <- read.csv("ClimVar_NativeDiv_Cover_20150423.csv", skip=3) %>%
  tbl_df()

#extract notes from native cover data and remove "X" from plot number#
natcovnotes_2015 <- NatDivCov_2015[1,] %>%
  gather(plot, notes, X1:X16) %>%
  select(-Plot) %>%
  mutate(plot=extract_numeric(plot))

#extract vegetation data
natdivveg_2015 <- NatDivCov_2015[-c(1,2,3,4,5,6,23),] %>% 
  gather(plot, cover, X1:X16, na.rm=T) %>%
  mutate(plot=extract_numeric(plot))  %>%
  mutate(Plot = as.character(Plot),
         Plot = ifelse(Plot == "Lupinus microcarpus densiflorus ", "Lupinus microcarpus densiflorus", Plot))
         

vegkey_2015 <- read.csv("~/Dropbox/ClimVar/DATA/Plant_composition_data/Cover/Cover_EnteredData/ClimVar_cover_specieskey.csv")

#merge data with keys
natdivveg0_2015 <-left_join(natdivveg_2015, vegkey_2015) 

NDvegtidy_2015 <- merge(natdivveg0_2015, shelterkey) %>%
  tbl_df() %>%
  filter(cover!="") %>%
  mutate(cover=as.numeric(cover))  %>%
  select(-Plot) %>%
  mutate(year = 2015)


unique(NDvegtidy_2015$cover)



################ 2016 #################

##Read in raw data (Native Diversity Cover) as table data frames## 
NatDivCov_2016 <- read.csv("NativeDiversityCover_20160504.csv", skip=3) %>%
  tbl_df()

#extract notes from native cover data and remove "X" from plot number#
natcovnotes_2016 <- NatDivCov_2016[1,] %>%
  gather(plot, notes, X1:X16) %>%
  select(-Plot) %>%
  mutate(plot=extract_numeric(plot))

#extract vegetation data
natdivveg_2016 <- NatDivCov_2016[-c(1,2,3,4,5,6,7),-2] %>% 
  gather(plot, cover, X1:X16, na.rm=T) %>%
  mutate(plot=extract_numeric(plot)) %>%
  mutate(Plot = as.character(Plot),
         Plot = ifelse(Plot == "Lupinus microcarpus densiflorus ", "Lupinus microcarpus densiflorus", Plot),
         Plot = ifelse(Plot == "Anagallis arvensis", "Anagalis arvensis", Plot),
         Plot = ifelse(Plot == "Avena barbata", "Ave_bar", Plot),
         Plot = ifelse(Plot == "Avena fatua", "Ave_fat", Plot),
         Plot = ifelse(Plot == "Avenua fatua", "Ave_fat", Plot),
        Plot = ifelse(Plot == "Bromus hordeaceus", "Bro_hor", Plot),
         Plot = ifelse(Plot == "Lolium multiflorum", "Lol_mul", Plot),
         Plot = ifelse(Plot == "Taeniatherum caput-medusae", "Tae_med", Plot),
         Plot = ifelse(Plot == "Juncus bufonious", "Juncus bufonius", Plot),
         Plot = ifelse(Plot == "Unknown forb seedling", "Unknown forbs seedlings", Plot),
         Plot = ifelse(Plot == "Sherardia arvense", "Sherardia arvensis", Plot))
         
         


vegkey_2016 <- read.csv("~/Dropbox/ClimVar/DATA/Plant_composition_data/Cover/Cover_EnteredData/ClimVar_cover_specieskey_2016.csv")

#merge data with keys
natdivveg0_2016 <-left_join(natdivveg_2016, vegkey_2016)

# check all species retained
sort(unique(tbl_df(natdivveg0_2016 %>%
                     filter(is.na(species_name)))$Plot))


NDvegtidy_2016 <- merge(natdivveg0_2016, shelterkey) %>%
  tbl_df() %>%
  filter(cover!="") %>%
  mutate(cover=as.numeric(cover)) %>%
  select(-Plot, -Plot2, -ToDelete) %>%
  mutate(year = 2016)

unique(NDvegtidy_2016$cover)


###################################
### PUT 2015 AND 2016 TOGETHER ####
###################################

dat.nativecover <- rbind(NDvegtidy_2015, NDvegtidy_2016) 

vegkey_all <- dat.nativecover %>%
  select(species_name, species, genus, func, status) %>%
  unique()

zerokey <- expand.grid(unique(dat.nativecover$species_name), unique(dat.nativecover$plot), unique(dat.nativecover$year))
names(zerokey) =c("species_name", "plot", "year")

fullmat1 <- left_join(zerokey, vegkey_all)
fullmat2 <- left_join(fullmat1, shelterkey)


dat.cover_with0 <- merge(fullmat2, dat.nativecover, all.x = T) %>%
  mutate(cover = ifelse(is.na(cover), 0, cover)) %>%
  tbl_df() %>%
  mutate(subplot = "Native") %>%
  select(year, plot, subplot, treatment, shelterBlock, shelter, species_name, genus, species, func, status, cover)

write.csv(dat.cover_with0, "ClimVar_MasterCover_Native_1516.csv")
