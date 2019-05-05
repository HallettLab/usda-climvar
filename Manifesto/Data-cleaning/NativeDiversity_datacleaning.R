# clean and compile native diversity datasets
# author(s): LMH, CTW
# initated: May 2019

# script purpose:


# notes:
## this script retains pieces of older code LMH wrote for cleaning the native diversity cover datasets


# -- SETUP -----
rm(list=ls()) # clean environment
library(dplyr)
options(stringsAsFactors = F)
na_vals = c(" ", "", NA, "NA")

# set path to ClimVar plant data folder (main level)
datpath <- "~/Dropbox/ClimVar/DATA/Plant_composition_data/"
# list files in Native_Diversity entered data folder
datfiles <- list.files(paste0(datpath,"Native_Diversity/Native_Diversity_EnteredData/"), full.names = T)

# read in data:
# 2015 natdiv cover
cov2015 <- read.csv(datfiles[grep("cover.*2015", datfiles, ignore.case = T)], na.strings = na_vals, strip.white = T)
# 2015 natdiv anpp
anpp2015 <- read.csv(datfiles[grep("anpp.*2015", datfiles, ignore.case = T)], na.strings = na_vals, strip.white = T)
# 2016 natdiv cover
cov2016 <- read.csv(datfiles[grep("cover.*2016", datfiles, ignore.case = T)], na.strings = na_vals, strip.white = T)
# drought treatment lookup
shelterkey <- read.csv(paste0(datpath, "Shelter_key.csv"))

# read in and compile species key from ClimVar cover folder
spplist <- data.frame()
sppfiles <- list.files(paste0(datpath, "Cover/Cover_EnteredData/"), full.names = T)
# keep only spp list files (remove ClimVar cover dataset files)
sppfiles <- sppfiles[grep("specieskey", sppfiles)]
for(i in 1:length(sppfiles)){
  if(i == 1){
  temp_dat <- read.csv(sppfiles[i])
  colkeep <- colnames(temp_dat)
  }else{
    temp_dat <- read.csv(sppfiles[i])
    temp_dat <- temp_dat[colkeep]
  }
  #rbind to main dataset
  spplist <- rbind(spplist, temp_dat)
  rm(temp_dat)
  if(i == length(sppfiles)){
    # clean up final spplist
    spplist <- unique(spplist)
  }
}


# -- PREP COVER DATA -----
# 2015
#extract notes from native cover data and remove "X" from plot number#
natcovnotes_2015 <- NatDivCov_2015[4,] %>%
  gather(plot, notes, 1:16) %>%
  select(-Plot) %>%
  mutate(plot=extract_numeric(plot))

#extract vegetation data
natdivveg_2015 <- NatDivCov_2015[-c(1,2,3,4,5,6,23),] %>% 
  gather(plot, cover, X1:X16, na.rm=T) %>%
  mutate(plot=extract_numeric(plot))  %>%
  mutate(Plot = as.character(Plot),
         Plot = ifelse(Plot == "Lupinus microcarpus densiflorus ", "Lupinus microcarpus densiflorus", Plot))

