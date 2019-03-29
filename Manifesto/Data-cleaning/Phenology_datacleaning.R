# Combine and clean all phenology datasets for usda-climvar project
# created: mar 2019
# contact: caitlin.t.white@colorado.edu

# script purpose:
# bind all phenology data and associated photo filenames from 2015, 2016, and 2017 usda-climvar seasons
# standardize column names, subplot values
## > notes: 2017 season did not collect phenology data as in 2015-2016 but did take plot-level photos 3x during season
## > phenology surveyed in oct 2016 for native/DD and legume plots only in prep for competition and recruitment installment, include those data here
## > phenology photos as of 2019 mar 26 are in CTW's personal google drive (not enough dropbox space) and shared with LMH, will move to suding lab google drive later..

# manual corrections made to entered data files as a result of combing through:
# 1) date changed in 2015 data from 5/15/15 to 5/5/15
# 2) percent brown and percent green switched for plot 4B on 4/28/15 (should be 20% green, 78% brown)
# 3) data entry error for plot 1C 4/21/15 corrected from 75% bare to 15% bare
# 4) inserted lines in 2015 datasets for general plot-level notes (e.g. plot 1, subplot == NA)



# -- SETUP -----
rm(list = ls()) # clean environment
# load needed libraries
library(tidyverse)
library(readxl)
library(lubridate)
# modify default settings
options(stringsAsFactors = F)
theme_set(theme_bw())
na_vals <- c("" , " ", "NA", NA)

# set path to data in Climvar Dropbox
datpath <- "~/DropBox/ClimVar/DATA/"

#read in data
# for-loop to read in all pheno dats
pheno_dats <- list.files(paste0(datpath, "Plant_composition_data/Phenology/Phenology_EnteredData/"))
# remove csv files and Spring2017 workbook (no phenology data, only photos)
pheno_files <- pheno_dats[!grepl("csv|Spring2017", pheno_dats)]
#initiate df for all phenology data
pheno <- data.frame()
for(i in pheno_files){
  temp_dat <- read_excel(paste0(datpath, "Plant_composition_data/Phenology/Phenology_EnteredData/",i),
                       sheet = 1, na = na_vals, trim_ws = T)  
  #reorder columns and give standard name for binding
  colnames(temp_dat)[grepl("date", colnames(temp_dat), ignore.case = T)] <- "Date"
  colnames(temp_dat)[grepl("Plot", colnames(temp_dat), fixed = T)] <- "Plot"
  colnames(temp_dat)[grepl("Subplot", colnames(temp_dat), fixed = T)] <- "Subplot"
  colnames(temp_dat) <- gsub("[.]| ", "_", colnames(temp_dat))
  colnames(temp_dat) <- casefold(colnames(temp_dat))
  temp_dat <- temp_dat[c("date", "plot", "subplot", "percent_green", "percent_brown", "percent_bare", "notes")]
  # rbind to master phenology df
  pheno <- rbind(pheno, temp_dat)
}


# for-loop to read in all pheno photo keys
# remove csv files from pheno_dats and Spring 2017 xlsx (will add in later)
pheno_photo_dats <- pheno_dats[grepl("xls", pheno_dats) & !grepl("Spring2017", pheno_dats)]
#initiate df for all phenology data
pheno_photos <- data.frame()
for(i in pheno_photo_dats){
  temp_dat <- read_excel(paste0(datpath, "Plant_composition_data/Phenology/Phenology_EnteredData/",i),
                         sheet = 2, na = na_vals)  
  #reorder columns and give standard name for binding
  colnames(temp_dat) <- casefold(colnames(temp_dat))
  #temp_dat <- temp_dat[c("date", "plot", "subplot", ", "notes")]
  # rbind to master phenology df
  pheno_photos <- rbind(pheno_photos, temp_dat)
}
# add in Spring 2017 photo data (no phenology collected, but shelter photos taken 3 times during season)
## there are three photo keys in the spring 2017 workbook
for(i in 1:3){
  temp_dat <- read_excel(paste0(datpath, "Plant_composition_data/Phenology/Phenology_EnteredData/",pheno_dats[grepl("Spring2017", pheno_dats)]),
                         sheet = i, na = na_vals)  
  #reorder columns and give standard name for binding
  colnames(temp_dat) <- casefold(colnames(temp_dat))
  temp_dat <- temp_dat[!colnames(temp_dat)=="notes"]
  # rbind to master phenology df
  pheno_photos <- rbind(pheno_photos, temp_dat)
}



# -- PHENOLOGY CLEANING -----
# remove any rows where all cells are na
pheno <- subset(pheno, apply(pheno, 1, function(x) sum(is.na(x)))<ncol(pheno))

# check structure
str(pheno)
glimpse(pheno)
# what are the unique values?
lapply(pheno, function(x) sort(unique(x)))
# trim white space to standardize values before replacement
pheno <- mutate_all(pheno, trimws)
# change any "T" (trace) value to 0.5
# (25x25) to (25*25)/(100*200) (1mx2m subplot)
# < or > adjusted by -+ 2
pheno[c("percent_green", "percent_brown", "percent_bare")] <- sapply(pheno[c("percent_green", "percent_brown", "percent_bare")], function(x) ifelse(x == "T", 0.1, x)) 
pheno[c("percent_green", "percent_brown", "percent_bare")] <- sapply(pheno[c("percent_green", "percent_brown", "percent_bare")], function(x) ifelse(x == "(25x25)x2", (((25^2)/(100*200))*2)*100, x)) 
pheno[c("percent_green", "percent_brown", "percent_bare")] <- sapply(pheno[c("percent_green", "percent_brown", "percent_bare")], function(x) ifelse(grepl("25x25", x, ignore.case = T), ((25^2)/(100*200))*100, x)) 
pheno[c("percent_green", "percent_brown", "percent_bare")] <- sapply(pheno[c("percent_green", "percent_brown", "percent_bare")], 
                                                                     # if less than 1, assign 0.5
                                                                     function(x) ifelse(grepl("<1", x), 0.5,
                                                                                        #otherwise, if less than, subtract 2
                                                                                        ifelse(grepl("<",x), readr::parse_number(x)-2,
                                                                                               #if greater than, add 2; otherwise leave as is
                                                                                               ifelse(grepl(">", x), readr::parse_number(x)+2,x)))) 
# other clean up..
## standardize subplots (1st letter only, XC for control),
pheno$subplot <- sapply(pheno$subplot, function(x) ifelse(x == "Control", "XC", 
                                                          ifelse(x != "XC", substr(x,1,1), x)))
##coerce plot and percentages to numeric
pheno[c("plot", "percent_green", "percent_brown", "percent_bare")] <- sapply(pheno[c("plot", "percent_green", "percent_brown", "percent_bare")], as.numeric) 
## remove time from date, and coerce to date values
pheno$date <- as.Date(gsub(" [0-9].*", "", pheno$date), format = "%Y-%m-%d")
#review
glimpse(pheno) #fine
summary(pheno)
#arrange by date, plot, subplot
pheno <- arrange(pheno, date, plot, subplot)
#what are the plots with NAs (should just be general data collection notes)
pheno[is.na(pheno$plot),] #general data collection notes only, set plot and subplot as "All"
pheno$plot[is.na(pheno$plot)] <- "All"
pheno$subplot[is.na(pheno$subplot)] <- "All" # applies to system-wide and plot-wide notes


# visual check all looks fine
# note: intially plotted with outlier for pct bare (about 75% mid season), checked data sheets and had typed in "75" instead of "15"
# corrected typo in xlsx and csv on 3/25/19
pheno %>%
  gather(var, value, percent_green:percent_bare) %>%
  ggplot(aes(as.factor(date), value, group = date)) +
  geom_boxplot() +
  geom_point(alpha = 0.5) +
  facet_grid(var~lubridate::year(date), scales = "free_x") # looks fine

# --old code for correction to plot 1 C, won't run with corrected data
# #which row has 75% bare?
# pheno[pheno$percent_bare > 60,] #.. plot 1 C
# ggplot(subset(pheno, plot == 1 & subplot == "Compost"), aes(as.factor(date), percent_bare)) +
#   geom_point() # looks like data entry error.. reviewing data sheet, should be 15%. correcting raw entered data xlsx and csv..

# check by plot
# note: 4B decreased then increased in percent brown when first ran this ggplot code
# looking further.. seems like the percent green and percent brown values were switched. corrected raw data on 3/25/19
# ggplot looks better after correction
pheno %>%
  subset(lubridate::year(date) < 2016 & !is.na(subplot)) %>%
  gather(var, value, percent_green:percent_bare) %>%
  ggplot(aes(date, value, col = subplot)) +
  geom_path() +
  geom_point(alpha = 0.5) +
  facet_grid(plot~var, scales = "free_x") # okay (with corrected data)

pheno %>%
  #2016 only and remove rows with shelter-level notes
  subset(lubridate::year(date) == 2016 & !is.na(subplot)) %>%
  # remove october visit
  subset(month(`date`) < 10) %>%
  gather(var, value, percent_green:percent_bare) %>%
  ggplot(aes(date, value, col = subplot)) +
  geom_path() +
  geom_point(alpha = 0.5) +
  facet_grid(plot~var, scales = "free_x") # looks okay



# -- PHENOLOGY PHOTO CLEANING -----
# inspect data structure, values
glimpse(pheno_photos)
summary(pheno_photos)
lapply(pheno_photos[c("plot", "subplot", "datecaptured")], function(x) sort(unique(x)))

# standardize subplots
## change any NA in subplot to "shelter" (picture was for entire shelter)
pheno_photos$subplot[is.na(pheno_photos$subplot)] <- "All"
## abbreviate subplots Legume and Native/DD
pheno_photos$subplot <- with(pheno_photos, ifelse(subplot == "Legume", "L",
                                                  ifelse(grepl("Nati", subplot), "N/DD", subplot)))
## check
unique(pheno_photos$subplot) #good

# clean up date
## rename datecapture to date
colnames(pheno_photos)[grepl("date", colnames(pheno_photos))] <- "date"
## change date in photo key dataset from POSIXct to Date, remove timestamp
pheno_photos$date <- as.Date(gsub(" [0-9].*", "", pheno_photos$date), format = "%Y-%m-%d")
## check
unique(pheno_photos$date) #good

# change plot to character for joining with phenology dataset
pheno_photos$plot <- as.character(pheno_photos$plot)


# -- FINISHING -----
# join photo key data
pheno_master <- full_join(pheno, pheno_photos)
# final checks..
summary(pheno_master) 
lapply(select(pheno_master, date:percent_bare), function(x) sort(unique(x)))
# add a note to 2017 rows to indicate only photos collected, no data recorded
pheno_master$notes[year(pheno_master$date) == 2017] <- "No phenology data collected for Spring 2017 season, only photos of shelter/ambient plots"

# write out to cleaned phenology subfolder
write.csv(pheno_master, 
          file = paste0(datpath, "Plant_composition_data/Phenology/Phenology_CleanedData/ClimVar_Phenology_clean.csv"),
          row.names = F)
