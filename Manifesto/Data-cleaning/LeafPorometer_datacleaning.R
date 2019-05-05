# clean and compile 2015-2016 climvar porometer data
# author: ctw (caitlin.t.white@colorado.edu)
# initiated: may 2019

# script purpose:
# read in porometer raw datasets
# clean and compile both, append treatment data
# write out to climvar dropbox: Leaf_Poromter/Leaf_Porometer_CleanedData

# notes:
## porometer readings taken for A. barbata and E. botrys only in 2015 and 2016
## in 2015, taken from neutral subplots (control-plants, control-soil or both seeded)
## in 2016, taken from grass seeded, forb seeded, both seeded, and control-plants



# -- SETUP ----
library(readxl)
library(tidyverse)
library(lubridate)
options(stringsAsFactors = F)
theme_set(theme_bw())
na_vals = c("", " ", NA, "NA")

#set pathway to plant data on climvar dropbox
datpath <- "~/Dropbox/ClimVar/DATA/Plant_composition_data/"
# list files in porometer raw data
porfiles <- list.files(paste0(datpath, "Leaf_Porometer/Leaf_Porometer_RawData/"), full.names = T)
# keep raw data files only
porfiles <- porfiles[grep(".xls", porfiles)]

# read in data
# porometer
raw2015 <- read_excel(porfiles[1], sheet = 1, na = na_vals, trim_ws = T)
readme2015 <- read_excel(porfiles[1], sheet = 2, na = na_vals, trim_ws = T)
raw2016 <- read_excel(porfiles[2], sheet = 1, na = na_vals, trim_ws = T)
# drought treatment lookup
treatment <- read.csv(paste0(datpath,"Shelter_key.csv"), na.strings = na_vals)


# -- PREP AND COMPILE POROMETER DATA ----
# check readme
data.frame(readme2015)
# 5 reps per leaf per plot, sampled AVBA and ERBO, only sampled in Controls (XC, XS) or Both (B) subplots
glimpse(raw2015) # read in conductance and temp units to first row, date time correct
glimpse(raw2016) # read in conductance and temp units to first row, date time correct
# sample naming conventions consistent in both datasets

# rbind both, skipping first column, then adjust formatting
porclean <- rbind(raw2015[-1,], raw2016[-1,]) %>%
  # convert conductance and temp to numeric
  mutate_at(vars(Conductance, Temperature), as.double) %>%
  # lower case colnames, replace space with underscore
  rename_all(funs(casefold(gsub(" ", "_",.)))) %>%
  # split date from time
  mutate(date = parse_date(measurement_time, format = "%Y-%m-%d %H:%M:%S"),
         time = parse_time(measurement_time, format = "%Y-%m-%d %H:%M:%S"),
         year = year(date),
         # split sample ID in to plot, subplot, species and rep
         plot = parse_number(substr(sample_id, 1, 3)),
         subplot = str_extract(sample_id, "XC|XS|X|B|F|G|S"),
         # correct X in 2016 to XC
         subplot = ifelse(subplot == "X", "XC", subplot),
         code4 = str_extract(sample_id, "AV|ER"),
         # change sp abbreviation to 4-letter code
         code4 = ifelse(code4 == "AV", "AVBA", "ERBO"),
         # add full species name
         species = ifelse(code4 == "AVBA", "Avena barbata", "Erodium botrys"),
         rep = parse_number(gsub("P[0-9]+", "", sample_id))) %>%
  # join drought treatment data
  left_join(treatment) %>%
  #select final columns
  dplyr::select(plot, treatment:shelter, subplot, year, date:time, species, code4, rep, temperature, conductance) %>%
  # add units to conductance and temp colnames
  rename(cond_mmol_m2s = conductance,
         temp_c = temperature)

# check that all looks okay
str(porclean)
summary(porclean)
boxplot(porclean$cond_mmol_m2s)
boxplot(porclean$temp_c)

ggplot(porclean, aes(treatment, cond_mmol_m2s)) +
  geom_boxplot() +
  geom_point(aes(col = temp_c)) +
  facet_grid(year~code4) # looks okay, AV in 2015 lower vals than 2016 but similar trend


# -- WRITE OUT -----
write_csv(porclean, paste0(datpath, "Leaf_Porometer/Leaf_Porometer_CleanedData/LeafPorometer_2015-2016_clean.csv"))
