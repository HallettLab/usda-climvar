# clean and compile SLA and leaf water content data from 2015 and 2016
# author(s): ctw

# notes: 
# this cleaning may have already been done by LMH and paired with LL or Brad trait data, but not clear in R scripts on GitHub
# CTW making this so it's clear that it's done
# leaf wet weights not collected on AVBA or ERBO in 2015 
# 2016 leaves not collected from any shelters; collected in front of ClimVar experiment or on grassy hillslope by garage/lab

# > ctw also collected plant height data for 2016 species, but that dataset lives in the traits subfolder



# -- SETUP ----
library(dplyr)
options(stringsAsFactors = F)
na_vals <- c(" ", "", NA, "NA")

# set path to climvar plant data
datpath <- "~/Dropbox/ClimVar/DATA/Plant_composition_data/"
# list files
datfiles <- list.files(paste0(datpath, "SLA/SLA_EnteredData"), full.names = T)

# read in data
# 2015 SLA + leaf trait
sla2015 <- read.csv(datfiles[grep("2015", datfiles)], na.strings = na_vals, strip.white = T)
# 2016 SLA + leaf trait
sla2016 <- read.csv(datfiles[grep("2016", datfiles)], na.strings = na_vals, strip.white = T)

# shelter treatment
treatment <- read.csv(paste0(datpath,"Shelter_Key.csv"))


# -- REVIEW DATA -----
#2015
str(sla2015) # wet weight not collected
summary(sla2015)
# review unique vals in select cols
sapply(sla2015[c("Date.collected","Subplot", "species")], unique)
#Subplot: what is S? control soil?; short species codes

#2016
str(sla2016) # wet weight collected
summary(sla2016)
sapply(sla2016[c("date_collected","species", "weighing_notes")], unique)
# all leaves collected outside experiment, no plot or subplot data; full species names

# columns to keep:
# species, rep, from where, date collected, wet wgt, dry wgt, leaf area, sla, compile notes?


# -- CLEAN DATASETS -----
clean2015 <- sla2015
#lower case all colnames except SLA
colnames(clean2015)[!colnames(clean2015) == "SLA"] <- casefold(colnames(clean2015)[!colnames(clean2015) == "SLA"]) 
clean2015 <- clean2015 %>%
  dplyr::select(date.collected:species, dry.weight:SLA) %>%
  # change names to same as 2016 names
  rename(date_collected = date.collected,
         leaf_area_cm2 = leaf_area,
         number = speciesno,
         code4 = species) %>%
  # recalculate dry weight to be sure correct
  mutate(dry_wgt_g = dry.weight - tin.weight,
         #remove letters from rep
         number = as.numeric(gsub("[a-z]", "", number, ignore.case = T)),
         # convert date to Date class
         date_collected = as.Date(date_collected, format = "%m/%d/%y"),
         # add wet weight for pairing with 2016 data
         wet_wgt_g = NA,
         # recode "S" subplot to XS
         subplot = gsub("S", "XS", subplot),
         # recode AV to AVBA, ERO to ERBO (only 2 species exist)
         code4 = ifelse(code4 == "AV", "AVBA", "ERBO"),
         species = ifelse(code4 == "AVBA", "Avena barbata", "Erodium botrys")) %>%
  # select final columns
  dplyr::select(plot, subplot, species, code4, number, date_collected, 
                wet_wgt_g, dry_wgt_g, leaf_area_cm2, SLA, species) %>%
  arrange(code4, plot, subplot, number)
  
# clean 2016
clean2016 <- sla2016 %>%
  # add cols in 2015 sla data that aren't in 2016 data
  mutate(plot = NA, 
         subplot = NA,
         # correct species misspellings
         species = gsub("Convulv", "Convolv", species),
         species = gsub("pulchra", "pulcher", species),
         code4 = casefold(paste0(substr(species, 1,2), # first two letters of genus
                                # first two letters of epithet (name after spaces)
                                substr(gsub("^[A-Z].+ ","",species),1,2)), 
                          upper = T),
         # convert date to Date class
         date_collected = as.Date(date_collected, format = "%m/%d/%y")) %>%
  rename(leaf_area_cm2 = area_cm2) %>%
  arrange(code4, number)

# reorder cols to match cleaned 2015 dataset 
clean2016 <- clean2016[colnames(clean2015)]


# -- COMPILE ALL DATASETS -----
# combine 2015 and 2016 datasets, add shelter info
all.tog <- rbind(clean2015, clean2016) %>%
  left_join(treatment) %>%
  mutate(plot = ifelse(is.na(plot), "outside experiment", plot)) %>%
  dplyr::select(species:number, plot, treatment:shelter, subplot:SLA)


# -- WRITE OUT -----
write.csv(all.tog, paste0(datpath, "SLA/SLA_CleanedData/SLA_leafwgts_2015-2016_cleaned.csv"), row.names = F)
