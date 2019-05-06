# clean and compile native diversity datasets
# author(s): LMH, CTW
# initated: May 2019

# script purpose:
# read in native diveristy cover and anpp data
# build a native diversity species lookup table using climvar spp list and usda plants database
# clean and compile native diversity cover make long-form, tidy
# clean anpp separately
# write out both cleaned datasets and spp list to ClimVar Dropbox Native_Diversity_CleanedData subfolder


# notes:
## 2016 natdiv anpp entered with 2016 anpp dataset (live in anpp folder)
## >> 2016 may anpp datasheets nowhere to be found on climvar dropbox! need to be rescanned 
## this script retains pieces of older code LMH wrote for cleaning the native diversity cover datasets
## requires package `request` to scrape usda plants database

# requirements to run script:
#devtools::install_github("sckott/request")



# -- SETUP -----
rm(list=ls()) # clean environment
library(dplyr)
library(tidyr)
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
# 2016 anpp (contains natdiv anpp)
anpp2016 <- read.csv(paste0(datpath, "ANPP/ANPP_EnteredData/ClimVar_ANPP_20160516.csv"), na.strings = na_vals)
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
  # move general survey note to its own row if the 2015 dataset
  if(i == "cov2015"){
    gennotes <- notes[grep("general survey", notes$notes, ignore.case = T),]
    gennotes$plot <- NA
    # extract general note from the first note
    gennotes$notes <- regmatches(gennotes$notes, regexpr("General .* %", gennotes$notes))
    # remove general note from plot 1 note
    notes$notes <- gsub("[(]G.*%[)] ","",notes$notes)
    # put together, general note in first row
    notes <- rbind(gennotes, notes)
  }
  
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
  # manual species corrections
  # change any Vulpia to Vulpia bromoides (VUBR is common grass in USDA compost project and was ID'd in nat div 2016)
  vegdat3$col1[grepl("Vulpia", vegdat3$col1)] <- "Vulpia bromoides"
  # correct JUBU spelling
  vegdat3$col1 <- gsub("bufonious", "bufonius", vegdat3$col1)
  # correct SHAR spelling
  vegdat3$col1 <- gsub("arvense", "arvensis", vegdat3$col1) # there's only ANAR and SHAR, same epithet for both
  
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
    cover_master <- left_join(cover_master, nonspp_master, by = c("plot", "yr")) %>%
      # rearrange colnames
      dplyr::select(colnames(notes), 
                    percent_native, percent_exotic, percent_rock, percent_bare,
                    percent_litter, litter_depth_cm, species, cover)
    
    # clean up
    rm(vegdat, vegdat2, vegdat3, notes, nonsppdat, clean_vegdat, nonspp_master, gennotes)
  }
}



# -- CREATE NATIVE COVER SPP LIST -----
# grab all identified species for scraping usda plants database
natdivspp <- sort(unique(cover_master$species)) %>% na.omit()
spplist_master <- data.frame(species = natdivspp, 
                             genus = gsub(" .*", "", natdivspp),
                             epithet = gsub("^[a-z]+ ", "", natdivspp, ignore.case = T),
                             unknown = ifelse(grepl("Unk", natdivspp), 1,0)) %>%
  mutate(code4 = casefold(paste0(substr(genus,1,2), substr(epithet,1,2)), upper = T))


# -- APPEND USDA PLANTS DATA -----
# specify vars desired from usda plants database (there are 134)
usda_plantvars <- c("Symbol","Accepted_Symbol_x","Scientific_Name_x","Common_Name","State_and_Province",
                    "Category","Family","Family_Common_Name","Duration","Growth_Habit","Native_Status")

spplist_master <- cbind(spplist_master, data.frame(matrix(nrow = nrow(spplist_master), ncol=length(usda_plantvars))))
colnames(spplist_master)[which(colnames(spplist_master) == "X1"):ncol(spplist_master)] <- usda_plantvars
# run usda plants api query to scrape species info
# NOTE!!: this will throw an error ["Client error: (400) Bad Request"] if a species is spelled incorrectly in the cover data (a good QA check)
# loop will issue warnings about cbind command providing more variables to replace than there are in 
for(p in spplist_master$species[spplist_master$unknown == 0]){
  print(paste("Pulling USDA Plants data for",p))
  temp_genus <- spplist_master$genus[spplist_master$species == p]
  # take only first word in epithet (e.g. ignore subspecies or variety)
  temp_epithet <- gsub(" [a-z]+$", "", spplist_master$epithet[spplist_master$species == p]) 
  
  # grab usda plants data
  if(grepl("ssp.", p)){
    temp_susbp <- gsub("^[A-Z].+ ssp. ", "", p)
    temp_epithet <- gsub(" .*", "", temp_epithet)
    templist <- api("https://plantsdb.xyz") %>%
      api_path(search) %>%
      api_query_(Genus = eval(temp_genus), Species = eval(temp_epithet), Subspecies = eval(temp_susbp))
  }
  # special case for medusahead (any hyphenated species epithet is a special case, search function bonks with hyphen)
  if(grepl("Taen", p)){
    templist <- api("https://plantsdb.xyz") %>%
      api_path(search) %>%
      api_query(Genus = Taeniatherum, Species = `caput-medusae`)
  }
  if(!grepl(" ssp[.]|-", p)){
    templist <- api("https://plantsdb.xyz") %>%
      api_path(search) %>%
      api_query_(Genus = eval(temp_genus), Species = eval(temp_epithet))
  }
  # isolate desired cols
  temp_df <- templist$data[1,colnames(templist$data) %in% usda_plantvars]
  # rematch to updated name if accepted symbol doesn't match symbol
  if(temp_df$Symbol != temp_df$Accepted_Symbol_x){
    templist2 <- api("https://plantsdb.xyz") %>%
      api_path(search) %>%
      api_query_(Symbol = eval(temp_df$Accepted_Symbol_x))
    update_df <- templist2$data[1,colnames(templist2$data) %in% usda_plantvars]
    temp_df[,which(colnames(temp_df)=="Common_Name"):ncol(temp_df)] <- update_df[,which(colnames(update_df)=="Common_Name"):ncol(update_df)]
  }
  # cleanup empty cells
  temp_df[temp_df==""] <- NA
  # append to master data frame
  #usdaplants_df <- rbind(usdaplants_df, temp_df)
  # add to spplist_master
  spplist_master[spplist_master$species == p,usda_plantvars] <- as.data.frame(temp_df)
  
  # clean up if last species
  if(p == max(spplist_master$species[spplist_master$unknown == 0])){
    # finish by adding fxnl_grp and simplified nativity col
    spplist_master$fxnl_grp[grepl("Gram", spplist_master$Growth_Habit)] <- "Grass"
    spplist_master$fxnl_grp[grepl("Forb", spplist_master$Growth_Habit, ignore.case = T)] <- "Forb"
    spplist_master$fxnl_grp[spplist_master$Family == "Fabaceae"] <- "N-fixer"
    
    spplist_master$nativity[grepl("L48 .I.",spplist_master$Native_Status)] <- "Exotic"
    spplist_master$nativity[grepl("L48 .N.",spplist_master$Native_Status)] <- "Native"
    
    # add functional group, lifespan to unknown forb
    spplist_master$fxnl_grp[spplist_master$code4 == "UNFO"] <- "Forb"
    
    # clean up environment
    rm(temp_df, templist, templist2, update_df, temp_epithet, temp_genus,p)
  }
}


# combine spplist with cleaned cover data and write out
cover_master2 <- left_join(cover_master,spplist_master) %>%
  #drop unknown column 
  dplyr::select(-c(unknown, State_and_Province, Native_Status)) %>%
  # for general note rows, specify "All" in plot col
  mutate(plot = ifelse(is.na(plot), "All", plot))




# -- CLEAN ANPP DATASET ----
# apply 2016 colnames to 2015 dataset (minus date sorted and sorter)
colnames(anpp2015) <- colnames(anpp2016)[!grepl("sorte", colnames(anpp2016), ignore.case = T)]
# remove blank rows
anpp2015 <- anpp2015[!is.na(anpp2015$Group),]
#subset 2016 anpp to native subplots only
anpp2016 <- anpp2016 %>%
  subset(Subplot == "Native") %>%
  # drop sorter columns
  dplyr::select(-c(Sorter, Date_sorted))

# rbind both datasets
master_anpp <- rbind(anpp2015, anpp2016) %>%
  rename_all(casefold) %>%
  rename(wgt_g = weight_g) %>%
  # convert date collected to Date class and add year column
  mutate(date_collected = as.Date(date_collected, format = "%m/%d/%y"),
         yr = as.numeric(substr(date_collected,1,4)),
         notes = ifelse(!is.na(notes_on_bag) & !is.na(notes_from_sorting),
                        paste(notes_on_bag, notes_from_sorting, sep = "; "),
                        ifelse(is.na(notes_on_bag), notes_from_sorting, notes_on_bag)),
         # convert any NAs in wgt_g to 0 (no native present)
         wgt_g = ifelse(is.na(wgt_g), 0, wgt_g),
         # change non-native to exotic
         group = gsub("Non-native", "Exotic", group)) %>%
  # add drought treatment data
  left_join(trtkey) %>%
  # reorder cols and drop page, line, and date weighed
  dplyr::select(plot, treatment:shelter, subplot, yr, date_collected, recorder, group, wgt_g, notes)


# -- FINISHING -----
# write out cleaned cover
write.csv(cover_master2, paste0(datpath, "Native_Diversity/Native_Diversity_CleanedData/NativeDiversity_Cover_2015-2016_cleaned.csv"), row.names = F)
# write out cleaned anpp
write.csv(master_anpp, paste0(datpath, "Native_Diversity/Native_Diversity_CleanedData/NativeDiversity_ANPP_2015-2016_cleaned.csv"), row.names = F)
# write out cleaned native diversity species list
write.csv(spplist_master, paste0(datpath, "Native_Diversity/Native_Diversity_CleanedData/NativeDiversity_SpeciesKey.csv"), row.names = F)
