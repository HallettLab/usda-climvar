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
  dplyr::select(-recorder, -miniplot, -notes, -scaling.sub1, -scaling.sub2, -plot.density.sub1, -plot.density.sub2) %>%
  mutate(sample_month = ifelse(species %in% c("TACA", "LOMU"), "May", "April")) %>%
  gather(metric, value, individuals.sub1:wp.sub2) %>%
  mutate(sample = ifelse(grepl("sub1", metric)==TRUE, 1,2),
         metric = gsub(".sub1|.sub2","", metric)) %>%
  filter(!is.na(value)) %>%
  spread(metric, value) %>%
  mutate(wp = ifelse(is.na(wp)==TRUE, 0, wp)) %>%
  dplyr::select(-sample.date) %>%
  merge(treatment)
