# generate competition and germination species lookup table
# author: ctw (caitlin.t.white@colorado.edu)
# created: april 2019

# script purpose:
# compile species lookup table for competition and germination experiments 
# append USDA plants database data using Scott Chamberlain's request package
# write out each set to main level data folders of Competition and Recruitment experiments on ClimVar dropbox
# can be used for programmatically converting codes or names in datasets, in standardized way, or crosswalking traits datasets

# requirements to run script:
#devtools::install_github("sckott/request")


# -- SETUP ----
# clear environment
rm(list=ls())
# load libraries needed to use USDA plants api
library(readxl)
library(request) # to access USDA plants api
library(tibble) # request reads in data to tibble
library(dplyr) # request functions rely on pipes  
options(stringsAsFactors = F)
na_vals <- c(" ", "", NA, "NA")

# read in species list seeded in 2017 Competition and Germination (aka Recruitment) experiments
spp <- read_excel("~/Dropbox/ClimVar/Competition/Setup/CompGerm_seed_order_prep/SeedWeighing_Guide.xlsx",
           sheet = 1, na = na_vals, trim_ws = T)


# -- PREP DATA TO RUN THROUGH USDA PLANTS API -----
# trim to spp actually seeded
spp <- unique(spp[!is.na(spp$Species),c(1,3)])
colnames(spp) <- casefold(colnames(spp))
spp$species #remove: ERBO, "yellow star", "lupinus succulentis" (misspelled), VUMI
spp <- spp[!grepl("succulentis|yellow|botrys|michros", spp$species),]
#capitalize first letter in name
spp$species <- with(spp, paste0(casefold(substr(species,1,1), upper = T), substr(species, 2, nchar(species))))
spp$genus <- gsub(" .*", "", spp$species)
spp$epithet <- gsub("^[A-Z][a-z]+ ", "", spp$species)
# sp code used by LMH and CTW (upper case)
spp$code4 <- casefold(paste0(substr(spp$genus,1,2), substr(spp$epithet,1,2)), upper = T)
# sp code used by JL (to link with trait data) (lower case)
spp$code6 <- casefold(paste0(substr(spp$genus,1,3), substr(spp$epithet,1,3)))
# preserve experiment separately for writing out
destination <- spp[,1:2]
# remove experiment and duplicate species (seeded in both experiments)
spp <- spp[,!grepl("exp", colnames(spp))]; spp <- spp[!duplicated(spp$species),]


# -- RUN THROUGH USDA API, APPEND USDA PLANTS DATA -----
# specify vars desired from usda plants database (there are 134)
usda_plantvars <- c("Symbol","Accepted_Symbol_x","Scientific_Name_x","Common_Name","State_and_Province",
                    "Category","Family","Family_Common_Name","Duration","Growth_Habit","Native_Status")

spp2 <- cbind(spp, data.frame(matrix(nrow = nrow(spp), ncol=length(usda_plantvars))))
colnames(spp2)[which(colnames(spp2) == "X1"):ncol(spp2)] <- usda_plantvars
# run usda plants api query to scrape species info
# NOTE!!: this will throw an error ["Client error: (400) Bad Request"] if a species is spelled incorrectly in the cover data (a good QA check)
# loop will issue warnings about cbind command providing more variables to replace than there are in 
for(p in spp2$species){
  print(paste("Pulling USDA Plants data for",p))
  temp_genus <- spp2$genus[spp2$species == p]
  temp_epithet <- spp2$epithet[spp2$species == p] 
  
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
  # add to spp
  spp2[spp2$species == p,usda_plantvars] <- as.data.frame(temp_df)
  
  # finish by adding fxnl_grp and simplified nativity col
  if(p == spp2$species[length(spp2$species)]){
    spp2$fxnl_grp[grepl("Gram", spp2$Growth_Habit)] <- "Grass"
    spp2$fxnl_grp[grepl("Forb", spp2$Growth_Habit, ignore.case = T)] <- "Forb"
    spp2$fxnl_grp[spp2$Family == "Fabaceae"] <- "N-fixer"
    
    spp2$nativity[grepl("L48 .I.",spp2$Native_Status)] <- "Exotic"
    spp2$nativity[grepl("L48 .N.",spp2$Native_Status)] <- "Native"
    
    # alphabetize species
    spp2 <- arrange(spp2, species) # for some reason base order and sort won't work here.. must use dpylr function
    print("fin!")
  }
}

# review and write out to competition and germination data folders, dependent on which experiment seeded in
View(spp2)
write.csv(subset(spp2, species %in% destination$species[grepl("germ", destination$experiment)]),
          "~/Dropbox/ClimVar/Recruitment/Data/Recruitment_SpeciesKey.csv", row.names = F)
write.csv(subset(spp2, species %in% destination$species[grepl("comp", destination$experiment)]),
          "~/Dropbox/ClimVar/Competition/Data/Competition_SpeciesKey.csv", row.names = F)
