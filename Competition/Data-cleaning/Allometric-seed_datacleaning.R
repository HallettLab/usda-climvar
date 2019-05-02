# calculate biomass-fecundity allometric relationship by competiton species
# authors: LMH, CTW
# initate: Jan 2019

# script purpose:
# derive allometric relationship of biomass to seed production (fecundity) by species by fall drought treatment
# write out to competition cleaned data 

# notes: 
# specimens collected by fall dry and fall wet. check to see if treatment makes a difference on fecundity.
# more than 20 specimens collected for some species (e.g. april and late-season samples for forbs)


# -- SETUP -----
rm(list = ls()) # clean environment
library(readxl)
library(tidyverse)
library(broom)
options(stringsAsFactors = F)
theme_set(theme_bw())
na_vals = c("", " ", NA, "NA")
datpath <- "~/Dropbox/ClimVar/Competition/Data/"

# iterate through each tab in competition specimens 2017 workbook, read in and compile data.
# all tabs have same headers, each species x dry/wet on its own tab
# get spreadsheets
specsheets <- excel_sheets(paste0(datpath, "Competition_EnteredData/Competition_specimens_spring2017.xlsx"))
specsheets # first tab is metadata, tabs 2 through end are data tabs

# get colname order from metadata spreadsheet for ordering master compiled data frame (below)
specmeta <- read_excel(paste0(datpath, "Competition_EnteredData/Competition_specimens_spring2017.xlsx"), 
                       sheet = specsheets[1],
                       na = na_vals)
# clean up read in metadata table
# preserve variable descriptions and order
specmeta <- specmeta[23:nrow(specmeta),]
# set row one as colnames then remove
colnames(specmeta) <- specmeta[1,]; specmeta <- specmeta[-1,]

# read in cleaned biomass2 to join with allometric seed rates
comp.all <- read.csv(paste0(datpath, "/Competition_CleanedData/ClimVar_Comp_combined-biomass-2.csv"),
                     na.strings = na_vals, strip.white = T) 


# -- READ IN AND COMPILE ALL SPECIMEN TABS -----
# initiate df for storing specimen data
specdat <- data.frame()
for(i in specsheets[-1]){ # don't read in metadata tab
  print(paste("Reading in",i))
  tempdat <- read_excel(paste0(datpath,"Competition_EnteredData/Competition_specimens_spring2017.xlsx"), 
                        sheet = i, na = na_vals, trim_ws = T)
  #print(paste("Colnames for", i, "are:"))
  #print(colnames(tempdat))
  if(ncol(specdat) == 0){
    specdat <- rbind(specdat, tempdat)
  }else{
    #compare colnames in master vs read in dataset
    master_names <- colnames(specdat); tempdat_names <- colnames(tempdat)
    missing_master <- master_names[!master_names %in% tempdat_names]
    missing_tempdat <- tempdat_names[!tempdat_names %in% master_names]
    if(length(missing_master)>0){
      print(paste("Master colnames missing in", i, "being added to bind:"))
      print(missing_master)
      addtoTemp <- as.data.frame(matrix(ncol=length(missing_master), nrow = nrow(tempdat)))
      colnames(addtoTemp) <- missing_master
      tempdat <- cbind(tempdat, addtoTemp)
    }
    if(length(missing_tempdat)>0){
      print(paste(i, "colnames missing in master specimen data frame being added to bind:"))
      print(missing_tempdat)
      addtoMaster <- as.data.frame(matrix(ncol=length(missing_tempdat), nrow = nrow(specdat)))
      colnames(addtoMaster) <- missing_tempdat
      specdat <- cbind(specdat, addtoMaster)
    }
    # alphabetize temporarily to rbind; reorder cols when at last dataset read in
    tempdat <- tempdat[sort(colnames(tempdat))]
    specdat <- specdat[sort(colnames(specdat))]
    specdat <- rbind(specdat, tempdat)
  }
  if(i == specsheets[length(specsheets)]){
    # reorder columns
    specdat <- specdat[specmeta$Variable]
    print("All specimens imported and compiled in master data frame!")
    # clean up environment
    rm(addtoMaster, addtoTemp, tempdat, master_names, tempdat_names, missing_master, missing_tempdat, i)
  }
}


# -- SENSITIVITY CHECKS -----
# dry vs wet (does it matter? or can pool all specimens)
ggplot(specdat, aes(wgt_g, seeds, col = trt)) +
  geom_point(alpha = 0.5) +
  geom_smooth() +
  facet_grid(trt~species, scales = "free")

# april vs may samples for forbs (LACA)
subset(specdat, species == "LACA") %>%
  mutate(clip_month = ifelse(lubridate::month(clip_date)==4, "April", "May"),
         clip_month = ifelse(is.na(clip_month), "May", clip_month)) %>%
  ggplot(aes(wgt_g, seeds, col = clip_month)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se=F) +
  facet_wrap(~species, scales = "free") #hm.. april looks like better time to clip? maybe plants pooping out by may
# potentially should exclude may samples?

# check how slope changes if exclude may
# everything
summary(lm(seeds~wgt_g, data = subset(specdat, species == "LACA")))
# april only
summary(lm(seeds~wgt_g, data = subset(specdat, species == "LACA" & clip_date == as.Date("2017-04-19"))))


# -- SIMPLIFY (PER LMH) AND JOIN WITH CLEANED COMPETITION BIOMASS -----
# > note: ctw cleaned up code here significantly since all datasets compiled in for loop
# > LMH previously simplified each into separate datasets bc read in as separate list items
# put it together!
allo.tog <- specdat %>%
  select(species, trt, specimen, seeds, wgt_g) %>%
  filter(!is.na(wgt_g))

# graph it all together!
ggplot(allo.tog, aes(x=wgt_g, y=seeds, color = trt)) + geom_point() + geom_smooth(method = "lm", se =F) + 
  facet_wrap(~species, scales = "free") + geom_smooth(data = allo.tog, aes(x=wgt_g, y=seeds), method = "lm", se = F, color = "black")
# treatment doesn't seem to have any effect, can group all together in models

# long form way, don't care...
l <- lm(seeds~wgt_g, data = subset(allo.tog, species == "LACA"))
lacaout <- tidy(l) %>%
  mutate(species = "LACA")

l <- lm(seeds~wgt_g, data = subset(allo.tog, species == "AVFA"))
avfaout <- tidy(l) %>%
  mutate(species = "AVFA")

l <- lm(seeds~wgt_g, data = subset(allo.tog, species == "BRHO"))
brhoout <- tidy(l) %>%
  mutate(species = "BRHO")

l <- lm(seeds~wgt_g, data = subset(allo.tog, species == "VUMY"))
vumyout <- tidy(l) %>%
  mutate(species = "VUMY")

allo.out <- rbind(lacaout, avfaout, brhoout, vumyout)

allo.out_tomerge <- allo.out %>%
  select(species, term, estimate) %>%
  spread(term, estimate) %>%
  # keeps species codes to link to species lookup table (and traits potentially)
  #mutate(species = recode(species, LACA = "Lasthenia", AVFA = "Avena", BRHO = "Bromus", VUMY = "Vulpia")) %>%
  select(phyto = species,
         intercept = `(Intercept)`, 
         slope = wgt_g)

# join allometric seeding rates with biomass
comp.all2 <- left_join(comp.all, allo.out_tomerge) %>%
  mutate(seedsOut = intercept + slope*ind_weight_g)

write_csv(comp.all2, paste0(datpath, "Competition_CleanedData/ClimVar_Comp_fecundity.csv"))
