library(tidyverse)

## Set working directory
setwd("~/Dropbox/ClimVar/DATA/Plant_composition_data")

shelterkey <- read.csv("ANPP/ANPP_EnteredData/Shelter_key.csv")

#########################
##### 2015 ANPP DATA ####
#########################

aprdat_2015 <- read.csv("ANPP/ANPP_EnteredData/ClimVar_ANPP_20150322.csv") %>%
  tbl_df() %>%
  filter(!is.na(subplot)) %>%
  filter(!is.na(plot)) %>%
  mutate(func="Forb", func=ifelse(species=="AVBA" | species=="AVFA" | species == "BRDI" | species=="BRHO" | species=="HOMA" |
                                    species=="LOMU" | species== "TACA" | species == "BRADIS" | species == "CYDA" | 
                                    species == "VUMY" | species == "VUBR" | species == "VUSP" | species == "BRST" |species == "HOMU" | species == "BRMA" |
                                    species == "BRMI" | species == "GRASS" |
                                    species =="JUNCUS SP." | species == "BRODIA" | species == "JUBU", "Grass", func), 
         func=ifelse(species=="TRHI" | species == "TRSP" | species=="TRSPP" | 
                       species == "VISA" | species == "TRDU" | species == "TRGL" | species == "TRSU", "N-fixer",func),
         func = ifelse(species == "N/A", "All", func)) %>%
  mutate(weight_g=weight_g*16) %>%
  mutate(date = Date..collected., 
         harvest = "First",
         year = 2015) %>%
  select(year, date, harvest, plot, subplot, species, func, weight_g)


maydat_2015 <- read.csv("ANPP/ANPP_EnteredData/ClimVar_ANPP_20150419.csv") %>%
  tbl_df() %>%
  mutate(func="Forb", func=ifelse(species=="AVBA" | species=="AVFA" | species=="BRHO" | species=="HOMA" |
                                    species=="LOMU" | species== "TACA" | species == "BRADIS" | species == "CYDA" | 
                                    species == "VUMY" | species == "VUBR" | species == "BRST" |species == "HOMU" | species == "BRMI" |
                                    species =="JUNCUS SP." | species == "BRODIA" | species == "BRMA", "Grass", func), 
         func=ifelse(species=="TRHI" | species == "TRSP" | species=="TRSPP" | 
                       species == "VISA" | species == "TRDU" | species == "TRGL" | species == "TRSU", "N-fixer",func),
         func = ifelse(species == "N/A", "All", func)) %>%
  #func=ifelse(species=="CESO", "YellowStar", func), 
  #func=ifelse(species=="TACA", "Medusa", func)) %>%
  mutate(weight_g=weight_g*16) %>%
  mutate(date = Date..collected.,
         harvest = "Second",
         year = 2015) %>%
  select(year, date, harvest, plot, subplot, species, func, weight_g)

ANPP_2015_byspp <- merge(rbind(aprdat_2015, maydat_2015), shelterkey) %>%
 tbl_df() %>%
  mutate(focal=0, focal=ifelse(subplot=="G" & func=="Grass", 1, focal),
         focal=ifelse(subplot=="F" & func=="Forb", 1, focal),
         focal=ifelse(subplot=="B" & (func=="Forb" | func=="Grass"), 1, focal))

ANPP_2015_byfunc <- ANPP_2015_byspp %>%
  group_by(plot, year, date, harvest, subplot, func,treatment, shelterBlock, shelter) %>%
  summarize(weight_g_m = sum(weight_g)) %>%
  filter(subplot != "L", subplot != "C") %>%
  spread(func, weight_g_m, fill = 0) %>%
  gather(func, weight_g_m, 9:11) %>%
  tbl_df()

ANPP_2015_withC <- ANPP_2015_byspp %>%
  group_by(plot, year, date, harvest, subplot, treatment, shelterBlock, shelter) %>%
  summarize(weight_g_m = sum(weight_g)) %>%
  tbl_df()



#########################
##### 2016 ANPP DATA ####
#########################

aprdat_2016 <- read.csv("ANPP/ANPP_EnteredData/ClimVar_ANPP_20160416.csv") %>%
  tbl_df() %>%
  mutate(weight_g = ifelse(is.na(Weight_g), 0, Weight_g)) %>%
  mutate(weight_g = weight_g*16, 
         subplot = Subplot,
         func = Group,
         plot = Plot,
         date = Date_collected) %>%
  filter(!is.na(subplot)) %>%
  mutate(subplot = as.character(subplot),
         subplot = ifelse(subplot == "C", "XC", subplot)) %>%
  mutate(focal=0, focal=ifelse(subplot=="G" & func=="Grass", 1, focal),
         focal=ifelse(subplot=="F" & func=="Forb", 1, focal),
         focal=ifelse(subplot=="B" & (func=="Forb" | func=="Grass"), 1, focal)) %>%
  select(date, subplot, func, plot, weight_g, focal) %>%
  group_by(date, subplot, func, plot, focal) %>%
  summarize(weight_g = sum(weight_g)) %>%
  tbl_df() %>%
  mutate(year = 2016,
         harvest = "First")


maydat_2016 <- read.csv("ANPP/ANPP_EnteredData/ClimVar_ANPP_20160516.csv") %>%
  tbl_df() %>%
  mutate(weight_g = ifelse(is.na(Weight_g), 0, Weight_g)) %>%
  mutate(weight_g = weight_g*16, 
         subplot = Subplot,
         func = Group,
         plot = Plot,
         date = Date_collected) %>%
  mutate(focal=0, focal=ifelse(subplot=="G" & func=="Grass", 1, focal),
         focal=ifelse(subplot=="F" & func=="Forb", 1, focal),
         focal=ifelse(subplot=="B" & (func=="Forb" | func=="Grass"), 1, focal)) %>%
  select(date, subplot, func, plot, weight_g, focal) %>%
  mutate(year = 2016, harvest = "Second")


ANPP_2016_byfunc <- merge(rbind(aprdat_2016, maydat_2016), shelterkey) %>%
  tbl_df() %>%
  filter(subplot != "Native" & subplot != "C") %>%
  mutate(weight_g_m = weight_g) %>%
  select(-focal, -weight_g)



#########################
##### 2017 ANPP DATA ####
#########################

aprdat_2017 <- read.csv("ANPP/ANPP_EnteredData/ClimVar_ANPP_20170415.csv") %>%
  tbl_df() %>%
  mutate(weight_g = ifelse(is.na(wgt.g), 0, wgt.g)) %>%
  mutate(weight_g = weight_g*16, 
         func = group,
         date = date.clipped) %>%
  filter(!is.na(subplot),
         subplot != "C",
         subplot != "extra") %>%
  mutate(subplot = as.character(subplot),
         group = as.character(group)) %>%
  mutate(group = ifelse(group == "F", "Forb", group),
         group = ifelse(group == "N", "N-fixer", group),
         group = ifelse(group == "G", "Grass", group)) %>%
  mutate(focal= 0, focal=ifelse(subplot == "G" & func == "Grass", 1, focal),
         focal= ifelse(subplot=="F" & func=="Forb", 1, focal),
         focal= ifelse(subplot=="B" & (func=="Forb" | func=="Grass"), 1, focal)) %>%
  select(date, subplot, func, plot, weight_g, focal) %>%
  group_by(date, subplot, func, plot, focal) %>%
  summarize(weight_g = sum(weight_g)) %>%
  tbl_df() %>%
  mutate(year = 2017,
         harvest = "First")


maydat_2017 <- read.csv("ANPP/ANPP_EnteredData/ClimVar_ANPP_20170512.csv") %>%
  tbl_df() %>%
  tbl_df() %>%
  mutate(weight_g = ifelse(is.na(wgt.g), 0, wgt.g)) %>%
  mutate(weight_g = weight_g*16, 
         func = group,
         date = date.clipped) %>%
  mutate(subplot = as.character(subplot),
         group = as.character(group)) %>%
  mutate(group = ifelse(group == "F", "Forb", group),
         group = ifelse(group == "N", "N-fixer", group),
         group = ifelse(group == "G", "Grass", group)) %>%
  mutate(focal=0, focal=ifelse(subplot=="G" & func=="Grass", 1, focal),
         focal=ifelse(subplot=="F" & func=="Forb", 1, focal),
         focal=ifelse(subplot=="B" & (func=="Forb" | func=="Grass"), 1, focal)) %>%
  select(date, subplot, func, plot, weight_g, focal) %>%
  mutate(year = 2017, harvest = "Second")


ANPP_2017_byfunc <- merge(rbind(aprdat_2017, maydat_2017), shelterkey) %>%
  tbl_df() %>%
  filter(subplot != "Native" & subplot != "C") %>%
  mutate(weight_g_m = weight_g) %>%
  select(-focal, -weight_g)



###################################################################
#### PUT THEM ALL TOGETHER: BY FUNCTIONAL GROUP IN B, G, F, XC ####
###################################################################

ANPP_byfunc <- rbind(ANPP_2015_byfunc, ANPP_2016_byfunc, ANPP_2017_byfunc) %>%
  tbl_df() %>%
  mutate(subplot = as.character(subplot)) %>%
  group_by(year, subplot, plot, treatment) %>%
  mutate(totweight = sum(weight_g_m), propweight = weight_g_m/totweight) %>%
  tbl_df() %>%
  mutate(treatment=ordered(treatment, levels = c(consistentDry="consistentDry", springDry="springDry", fallDry="fallDry", controlRain="controlRain"))) 

write.csv(ANPP_byfunc, "ANPP/ANPP_CleanedData/ClimVar_ANPP-func-groups.csv", row.names = F)

############################################################
#### PUT THEM ALL TOGETHER: MAY BIOMASS; ALL SUBPLOTS ####
############################################################

maydat_2015_all <- ANPP_2015_withC %>%
  tbl_df() %>%
  filter(harvest != "First") %>%
  select(-harvest)


maydat_2016_all <- left_join(shelterkey, maydat_2016) %>%
  group_by(plot, year, date, harvest, subplot, treatment, shelterBlock, shelter) %>%
  summarize(weight_g_m = sum(weight_g)) %>%
  tbl_df()  %>%
  filter(harvest != "First") %>%
  select(-harvest)


maydat_2017_all <- left_join(shelterkey, maydat_2017) %>%
  group_by(plot, year, date, harvest, subplot, treatment, shelterBlock, shelter) %>%
  summarize(weight_g_m = sum(weight_g)) %>%
  tbl_df()  %>%
  filter(harvest != "First") %>%
  select(-harvest)


May_ANPP <- rbind(maydat_2015_all, maydat_2016_all, maydat_2017_all) %>%
  tbl_df() %>%
  filter(subplot != "L" & subplot != "Native")

write.csv(May_ANPP, "ANPP/ANPP_CleanedData/ClimVar_ANPP-peak.csv", row.names = F)
