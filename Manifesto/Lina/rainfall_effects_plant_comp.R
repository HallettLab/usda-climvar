library(tidyverse)
library(ggplot2)
library(ggpubr)

#Rainfall treatment effects on plant composition
setwd("./Dropbox/ClimVar/DATA/Plant_composition_data")
veg_all <- read.csv("Cover/Cover_CleanedData/ClimVar_species-cover.csv", header = TRUE)

#Filter 2015 data
veg_2015 <- veg_all %>% 
  filter(year == 2015) %>%
  filter(subplot != "XC")

