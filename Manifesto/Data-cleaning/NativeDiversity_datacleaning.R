# clean and compile native diversity datasets
# author(s): LMH, CTW
# initated: May 2019

# script purpose:
# read in native diveristy cover and anpp data
# build a native diversity species lookup table using climvar spp list and usda plants database
# clean and compile native diversity cover make long-form, tidy
# clean anpp separately
# write out both cleaned datasets and spp list to Native Diversity Cleaned Data subfolder


# notes:
## this script retains pieces of older code LMH wrote for cleaning the native diversity cover datasets
## requires package `request` to scrape usda plants database



# -- SETUP -----
rm(list=ls()) # clean environment
library(dplyr)
library(request) # to fetch usda plants database data
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
trtkey <- read.csv(paste0(datpath, "Shelter_key.csv"))

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


# -- FIX TYPOS IN SPP LIST -----
fixlist <- list(c("Anagalis", "Anagallis"), #bad spelling, correct spelling
                c("Cynosaurus", "Cynosurus"),
                c("Fillago", "Filago"),
                c("Vupia sp.", "Vulpia sp."))
# itate through fix list and replace values in species name and genus columns
for(i in 1:length(fixlist)){
  spplist[c("species_name", "genus")] <- sapply(spplist[c("species_name", "genus")], function(x)
    # replace bad spelling with correct spelling
    gsub(fixlist[[i]][1], fixlist[[i]][2], x))
}


# -- TIDY AND TRANSPOSE ENTERED COVER DATA -----
# loop through cover data, transpose, and append to master cover dataset
# initiate master data frame
cover_master <- data.frame()
nonspp_master <- data.frame()

# loop iterates through each cover dataset and adds new spp to master spp set
for(i in c("cov2015", "cov2016")){
  # read in dataset
  vegdat <- get(i)
  print(paste("Transposing and tidying", i, "dataset"))
  # remove blank rows
  vegdat <- vegdat[!is.na(vegdat[,1]),]
  # standardize colnames
  colnames(vegdat) <- c("col1", paste0("p", 1:(ncol(vegdat)-1)))
  # id where plot and cover data start
  plotpos <- grep("plot", vegdat[,1], ignore.case = T)
  
  
  # -- CREATE PLOT-DATE DATA TABLE -----
  # pull out recorder, sample date, and notes into separate data frame
  # notes row is always 1 below plot row in native cover datasets
  notes <- vegdat[1:(plotpos+1),]  
  notes <- data.frame(t(notes)) #transpose, change matrix class to data frame
  colnames(notes) <- casefold(gsub(":", "", notes[1,])) # set row 1 as column names, make lower case and remove colon
  notes <- notes[-1,] # remove row 1
  rownames(notes) <- seq(1,nrow(notes),1) # rename rownames in numeric sequence
  
  # change date from character to date format and add year column
  notes$date <- as.Date(notes$date, format = "%m/%d/%y")
  notes$yr <- as.numeric(substr(notes$date, 1,4))
  
  # join plot treatment info
  notes <- mutate(notes, plot = as.numeric(plot)) %>%
    left_join(trtkey, by = "plot") %>% #left_join preserves order of sampling, merge alphabetizes plots
    #reorder cols
    dplyr::select(plot, treatment:shelter, yr, date, recorder, notes)
  
  
  # -- PREP ABUNDANCE DATA FOR TIDYING AND TRANSPOSING -----
  # remove notes df rows from vegdat
  vegdat2 <- vegdat[c(plotpos:nrow(vegdat)),] %>%
    # remove notes row, and native/other headers
    subset(!grepl("notes| measurement", col1, ignore.case = TRUE))
  
  # id row where species abundance data starts
  spppos <- grep("Achil", vegdat2[,1])
  # remove any species rows that don't have any entries for abundance value
  allNAs <- apply(vegdat2[, 2:ncol(vegdat2)], 1, function(x) all(is.na(x)))
  # set rows up to litter depth as FALSE to preserve them (in case 0s not entered and rock or bare never encountered)
  allNAs[1:spppos] <- FALSE 
  # allNAs <- allNAs[!grepl("NA",names(allNAs))] # remove NAs created by ignoring rows 1 through litter depth
  vegdat2 <- vegdat2[!allNAs,] 
  
  # pull out non-species measurements (e.g. percent litter, litter depth)
  nonsppdat <- vegdat2[1:(spppos-1),]
  # clean up names (to match cov2016 format)
  nonsppdat[,1] <- gsub("%", "Percent_", nonsppdat[,1])
  nonsppdat[,1] <- gsub(" [/].*", "", nonsppdat[,1])
  nonsppdat[,1] <- gsub("non-native", "exotic", nonsppdat[,1])
  nonsppdat$col1[grep("depth", nonsppdat$col1)] <- "Litter_depth_cm"
  nonsppdat <- data.frame(t(nonsppdat))
  colnames(nonsppdat) <- casefold(nonsppdat[1,]) # set row 1 as column names, make lower case
  nonsppdat <- nonsppdat[-1,] # remove row 1
  rownames(nonsppdat) <- 1:nrow(nonsppdat) # rename rownames in numeric sequence
  nonsppdat <- mutate_all(nonsppdat, as.numeric) %>% # all vars numeric
    subset(!is.na(plot)) %>%
    gather(met, val, colnames(.)[2]:ncol(.)) %>%
    # add year
    mutate(yr = as.numeric(gsub("[a-z]+", "",i, ignore.case = T))) %>%
    dplyr::select(plot, yr, met, val) #reorder cols
  
  # create species-only dataframe
  vegdat3 <- vegdat2[c(1,spppos:nrow(vegdat2)),] 
  # clean up spp names based on climvar spp list
  vegdat3$col1 <- sapply(vegdat3$col1,function(x) ifelse(x %in% spplist$Plot, # if name is in the spp list
                                                         # replace with spp list species_name values
                                                         spplist$species_name[spplist$Plot == x],
                                                         # else do nothing
                                                         x))
  
  vegdat3 <- data.frame(t(vegdat3)) # transpose, change matrix to data frame
  colnames(vegdat3) <- vegdat3[1,] # set row 1 as colnames
  colnames(vegdat3)[1] <- casefold(colnames(vegdat3)[1])
  vegdat3 <- vegdat3[-1,] # remove row 1
  # renumber row names sequentially
  rownames(vegdat3) <- 1:nrow(vegdat3)
  # gather all species in single column
  vegdat3 <- gather(vegdat3, species, cover, 2:ncol(vegdat3)) %>%
    # remove row if species wasn't in that plot
    subset(!is.na(cover)) %>%
    # make plot and cover numeric classes
    mutate_at(vars(plot, cover), as.numeric) %>%
    # order by plot then species A-Z
    arrange(plot, species)
  
  # join notes and cover data
  clean_vegdat <- full_join(notes, vegdat3, by = "plot")
  
  
  # -- APPEND LONG FORM TO MASTER LONG FORM -----
  cover_master <- rbind(cover_master, clean_vegdat)
  nonspp_master <- rbind(nonspp_master, nonsppdat)
  
  # clean up if cover 2016
  if(i == "cov2016"){
    nonspp_master <- spread(nonspp_master, met, val)
    # grab nonspp columns names for sorting master data frame
    nonsppnames <- colnames(nonspp_master)
    # remove plot and yr
    nonsppnames <- nonsppnames[!nonsppnames %in% c("plot", "yr")]
    
    # put it all together
    cover_master <- left_join(cover_master, nonspp_master, by = c("plot", "yr"))
    
    # clean up
    rm(vegdat, vegdat2, vegdat3, notes, nonsppdat, clean_vegdat)
  }
}


# -- CREATE NATIVE COVER SPP LIST -----
# compile all unique native diversity species, starting with A. millefolium (first species on datasheet)
natdivspp <- c(cov2015$Sheet[grep("Achil", cov2015$Sheet):nrow(cov2015)],
               cov2016$Page[grep("Achil", cov2016$Page):nrow(cov2016)])

natdivspp <- sort(unique(natdivspp))
natdivspp
# are all native diversity cover species in the climvar cover species list?
summary(natdivspp %in% spplist$Plot)
# which native diversity species are not in the climvar cover species list?
natdivspp[!natdivspp %in% spplist$Plot]




