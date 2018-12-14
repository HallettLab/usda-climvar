library(tidyverse)
library(stringr)


#########################
### 2015 COVER DATA #####
#########################

## Import import keys
vegkey_2015 <- read.csv("~/Dropbox/ClimVar/DATA/Plant_composition_data/Cover/Cover_EnteredData/ClimVar_cover_specieskey_2015.csv")
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

#extract vegetation data
covveg_2015 <- covdat_2015[-c(1:8),] %>% 
  gather(plotname, cover, X1B:X16C, na.rm=T) %>%
  mutate(plotname=substr(plotname, 2,5)) %>%
  mutate(plot=extract_numeric(plotname)) %>%
  mutate(subplot=substr(plotname,2,4), subplot=ifelse(plot > 9, substr(plotname, 3,4), subplot))

#merge data with keys
veg0_2015 <- merge(covveg_2015, vegkey_2015)
veg_2015 <- merge(veg0_2015, shelterkey) %>%
  tbl_df() %>%
  filter(cover!="") %>%
  mutate(cover=as.numeric(cover),
         year = 2015) %>%
  mutate(species_name = as.character(species_name),
         species = as.character(species), 
         status = as.character(status),
         species_name = ifelse(species_name == "Galium aparine", "Galium parisiense", 
                               ifelse(species_name == "Torilis nodosa", "Torilis arvensis", species_name)),
         species = ifelse(species == "aparine", "parisiense", 
                          ifelse(species == "nodosa", "arvensis",species)),
         status = ifelse(species_name == "Galium parisiense", "non-native", status))

sort(unique(veg_2015$species_name))
rm(vegkey_2015, veg0_2015, covveg_2015, covdat_notes_2015, covdat_2015, covdat_2015_0, covdat_compost_2015)


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

#extract vegetation data
covveg_2016 <- covdat_2016[-c(1:6),] %>% 
  gather(plotname, cover, X1B:X16XC, na.rm=T) %>%
  mutate(plotname=substr(plotname, 2,5)) %>%
  mutate(plot=extract_numeric(plotname)) %>%
  mutate(subplot=substr(plotname,2,4), subplot=ifelse(plot > 9, substr(plotname, 3,4), subplot)) %>%
  mutate(Plot = as.character(Plot))

#import keys
vegkey_2016 <- read.csv("~/Dropbox/ClimVar/DATA/Plant_composition_data/Cover/Cover_EnteredData/ClimVar_cover_specieskey_2016.csv") %>%
  tbl_df() %>%
  mutate(Plot = as.character(Plot), Plot2 = as.character(Plot2)) %>%
  mutate(Plot = ifelse(Plot2 != "", Plot2, Plot)) %>%
  dplyr::select(-Plot2, -ToDelete) %>%
  mutate(species_name = as.character(species_name),
         species = as.character(species))

# merge data with species key
veg0_2016 <- merge(covveg_2016, vegkey_2016, all.x=T)


# check that all species were retained
length(intersect(vegkey_2016$Plot, veg0_2016$Plot))
length(unique(veg0_2016$Plot))

# merge data with shelter key
veg_2016 <- merge(veg0_2016, shelterkey, all.x = T) %>%
  tbl_df() %>%
  filter(cover!="") %>%
  mutate(cover=as.numeric(cover)) %>%
  mutate(year = 2016)

sort(unique((veg_2016$species_name))) 

rm(vegkey_2016, veg0_2016, covveg_2016, covdat_notes_2016, covdat_2016)



#########################
### 2017 COVER DATA #####
#########################

# import data
covdat_2017 <- read.csv("~/Dropbox/ClimVar/DATA/Plant_composition_data/Cover/Cover_EnteredData/ClimVar_Cover_20170526.csv") %>%
  tbl_df() %>%
  gather(plotname, cover, X1F:X16F, na.rm=T) %>%
  mutate(plotname=substr(plotname, 2,5)) %>%
  mutate(plot=extract_numeric(plotname)) %>%
  mutate(subplot=substr(plotname,2,4), subplot=ifelse(plot > 9, substr(plotname, 3,4), subplot)) %>%
  mutate(Plot = as.character(Plot)) %>%
  filter(!Plot%in%c("Percent_rock", "Percent_bare", "Percent_litter", "Litter_depth_cm",
                    "Percent_grass", "Percent_forb"))

# import keys
vegkey_2017 <- read.csv("~/Dropbox/ClimVar/DATA/Plant_composition_data/Cover/Cover_EnteredData/ClimVar_cover_specieskey_2017.csv") %>%
  tbl_df() %>%
  mutate(Plot = as.character(Plot), Plot2 = as.character(Plot2)) %>%
  mutate(Plot = ifelse(Plot2 != "", Plot2, Plot)) %>%
  dplyr::select(-Plot2, -ToDelete) %>%
  mutate(species_name = as.character(species_name))

# merge data with species key
veg0_2017 <- merge(covdat_2017, vegkey_2017, all.x=T)


# check that all species were retained
length(intersect(vegkey_2017$Plot, veg0_2017$Plot))
length(unique(veg0_2017$Plot))


# merge data with shelter key
veg_2017 <- merge(veg0_2017, shelterkey, all.x = T) %>%
  tbl_df() %>%
  filter(cover!="") %>%
  mutate(cover=as.numeric(cover)) %>%
  mutate(year = 2017)


rm(vegkey_2017, veg0_2017, covdat_2017)


##################################
### COMBINE YEARS, ADD ZEROS #####
##################################

dat.cover <- rbind(veg_2015, veg_2016, veg_2017) %>%
  tbl_df() %>%
  dplyr::select(-Plot)

vegkey <- dat.cover %>%
  dplyr::select(species_name, genus, species, func, status) %>%
  unique()

# check for consistent names across years
sort(unique(dat.cover$species_name))

# add in 0s for absent species
zerokey <- expand.grid(unique(dat.cover$species_name), unique(dat.cover$plot), unique(dat.cover$subplot), unique(dat.cover$year))

names(zerokey) =c("species_name", "plot", "subplot", "year")

 fullmat1 <- merge(zerokey, vegkey, all.x = T) 
 fullmat2 <- merge(fullmat1, shelterkey, all.x = T)  


dat.cover_with0 <- merge(fullmat2, dat.cover, all.x = T) %>%
  dplyr::select(-plotname) %>%
  mutate(cover = ifelse(is.na(cover), 0, cover)) %>%
  tbl_df()  %>%
  mutate(func2 = as.character(func),
         func2 = ifelse(genus == "Trifolium" | genus == "Vicia", "Nfixer", func2),
         func2 = as.factor(func2)) %>%
  filter(subplot != "L")

rm(fullmat1, fullmat2, vegkey, shelterkey, zerokey, veg_2015, veg_2016, veg_2017)

write.csv(dat.cover_with0, "~/Dropbox/ClimVar/DATA/Plant_composition_data/Cover/Cover_CleanedData/ClimVar_species-cover.csv")




