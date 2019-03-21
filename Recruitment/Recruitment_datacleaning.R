# ClimVar
# 2017 Recruitment data cleaning
# Script purpose:
## 1. Read in raw recruitment experiment data
## 2. Scale stem counts to per meter density (average subsamples when applicable)
## 3. Join rainfall treatment data
## 4. Write cleaned data to .csv in Cleaned_data DropBox folder

# Notes: 
## -- working directory set to DropBox/ClimVar/DATA/Plant_composition_data
## -- germination whole plot size = 25cm x 25cm, subsample size = 5cm x 5cm
## -- rule: defer to whole plot numbers when available (even when subsample data present)

# Script setup
library(tidyverse)
theme_set(theme_bw())

#########################
## 1. Read in raw data ##
#########################

recruit_dat <- read.csv("ComExpt17/ComExpt_EnteredData/recruitment_dat_spring2017.csv",
                        na.strings = c("", "NA", " "),
                        stringsAsFactors = FALSE)
str(recruit_dat)
# Issues:
## 1. NAs in sample.date are April samples (sample date in April unknown)
###  -- Make date a categorical?: April/May
## 2. Can remove: recorder, miniplot, notes, scaling columns, density columns
###  -- Set scaling and calculate density in cleaning script

treatment <- read.csv("Shelter_key.csv",
                      stringsAsFactors = FALSE) 
str(treatment) 

###################
## 2. Clean data ##
###################

recruit_clean <- recruit_dat %>%
  # remove single fecundity observation for BRNI, manually adjust in code later as its the only case
  # remove notes, will add disturbance data back in later
  dplyr::select(-recorder, -miniplot, -individuals.fecund, -heads.fecund, -wp.fecund, -notes) %>%
  mutate(sample_month = ifelse(species %in% c("TACA", "LOMU"), "May", "April")) %>%
  gather(metric, value, individuals.sub1:wp.sub2) %>%
  mutate(sample = ifelse(grepl("sub1", metric)==TRUE, 1,2),
         metric = gsub(".sub1|.sub2","", metric)) %>%
  filter(!is.na(value)) %>%
  spread(metric, value) %>%
  dplyr::select(-sample.date) %>%
  group_by(plot, species, sample_month) %>%
  mutate(double_sample = length(heads))
  # create rules for averaging stem counts and heads by plot and species
  # where whole plot censused, prioritize that number over 5x5cm subsample
  group_by(plot, species, sample_month) %>%
  mutate(stems = ifelse(wp==1, individuals,
                              mean(individuals, na.rm=TRUE)))
  #merge(treatment)

# create keyword search for disturbance (e.g. gopher, rain, browsing, background weeds)
keywords <- c("rain", "disturb", "infest", "invade", "gopher", "munch")

