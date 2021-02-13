# prep USDA Compost TRHI biomass-fecund data for USDA ClimVar competition allometry
# author(s): ctw, jan 2021

# script purpose:
## read in USDA Compost competition TRHI data (collected in field by NS Spring 2020, processed @ Eugene)
## read in USDA Compost treatment key (CTW copied over to )
## subset to control plots (? need to screen and consider)..
# > screen for any differences between wet, control, and drought plots. only select TRHI from unammended/unfert plots.
# > 4 blocks x 1 nutrient control x 3 ppt trts = potentially 12 specimens, unless there are diffs by ppt
# > background TRHI is also included
## prep what's needed to run through allometric-seed_datacleaning.R with other species dat (apply similar colnames as other datasets)
## write out prepped dataset to ClimVar Dropbox/Competition/Data/EnteredData/outside_project_data_related/


# notes:
## TRHI clipped in April 2017 in USDA ClimVar Competition, too early for seed development
## to infill, borrowing data from USDA Compost competition experiment.
## For allometric relationships in USDA ClimVar, individuals of each species were selected on size gradient (small to large), 10 from "dry" plots (consistent dry or fall dry drought) and 10 from "wet" plots (control or spring dry)
## > CTW pulled individuals from subplots where they were growing as the background species (plenty to choose from), typically from low density plots, but selected haphazardly prioritizing healthy/representative plants to capture size gradient
## USDA Compost experiment is not set up exactly like that and no specimens were harvested specifically for allometry.. TRHI read into this script are the actual experiment data

# AM says 999 in phyto = a background specimen



# -- SETUP -----
library(tidyverse)
theme_set(theme_test())
options(stringsAsFactors = F)
na_vals = c("", " ", NA, "NA")
datpath <- "~/Dropbox/ClimVar/Competition/Data/" # path for CTW + LMH, can add for other users as needed

# list related compost files 
compostfiles <- list.files(paste0(datpath, "Competition_EnteredData/outside_project_data_related/"),full.names = T)

# read in compost dat
trhi <- read.csv(compostfiles[grep("winter2021", compostfiles)], na.strings = na_vals)
# trtment key
compost_key <- read.csv(compostfiles[grep("TreatmentK", compostfiles)], na.strings = na_vals)
# compost layout
layout <- read.csv(compostfiles[grep("layout", compostfiles)], na.strings = na_vals)


# review import
str(trhi)
summary(trhi) # one phyto entry has 5 indivs (stems).. that shouldn't be the case (AS planted 3 max).. if it's not a plot kept in subset, moot
str(compost_key)
str(layout)
summary(layout) # no empty rows, good.



# -- DATA PREP ----
# join trt info to trhi dat, drop fert plots, screen for diffs btwn ppt trts (cuidado con phyto v. background!)
trhi_sub <- left_join(trhi, compost_key) %>%
  # join background spp
  left_join(distinct(layout[c("plot", "subplot", "background")]))

# quick viz to see who has samples v. not and screen for obvi diffs before drop amended plots
# who haz?
trhi_sub %>%
  mutate(present = as.numeric(!is.na(num_stems))) %>%
  ggplot(aes(factor(block),factor(present))) +
  geom_point(aes(col = phyto == 999), position = position_jitter(width = 0.2, height = 0.1), alpha = 0.5) +
  scale_color_discrete(name = "Background") +
  facet_grid(ppt_trt~nut_trt)
# hm.. 
# > TRHI generally did not survive (at least not to harvest time in spring) in control ppt plots 
# > some survival/maturation in nutrient control plots, but survival rate best in compost which.. those specimens will likely be bigger than in other nut_trt plots)
# > there are no background specimen survivors in control nut_trt plots.. meaning we need to use phyto dat..

summary(is.na(trhi_sub$num_stems)) # only 29 specimens max
with(trhi_sub, sapply(split(num_stems, nut_trt), function(x) summary(is.na(x)))) # there are only 10 specimens from nutrient control plots.. 
# > I guess we take what we can get! specimens clipped for ClimVar allo likely were in intraspec comp.. not quite the same as phytos in interspec comp, but.. better than using plant treated with compost or fertilizer
# > note too: where more than one individual was clipped per phyto pos, it's all combined here (e.g. 3 stems, combined tot seeds for all 3. ClimVar allo specs were measured individually

# how many stems repped in trts?
with(subset(trhi_sub, !is.na(num_stems)), sapply(split(num_stems, nut_trt), sum)) # okie.. 72 stems for nutrient control
# how many stems in nutrient control per drought trt?
with(subset(trhi_sub, !is.na(num_stems) & nut_trt == "N"), sapply(split(num_stems, ppt_trt), sum))

# check diffs in biomass and seed fecund by ppt and nut trt (expect to see diff by nut trt, esp since compost trt retains soil moisture..)



  
