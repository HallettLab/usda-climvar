ANPP <- read.csv("~/Dropbox/ClimVar/DATA/Plant_composition_data/ANPP/ANPP_CleanedData/ClimVar_ANPP-peak.csv", header = TRUE)
BNPP <- read.csv("~/Dropbox/ClimVar/DATA/Plant_composition_data/BNPP/BNPP_CleanedData/BNPP_MayHarvest_2015.csv", header = TRUE)
CWM <- read.csv("~/Dropbox/ClimVar/DATA/Plant_composition_data/Traits/Traits_ProcessedData-GH/ClimVar_trait-diversity-GH.csv", header = TRUE)

###How does total biomass relate to CWM traits?

##Add BNPP to ANPP for total biomass
library(tidyverse)
#Replace entry "20-Oct" to "10-20" in depth column in the BNPP dataset
BNPP$depth <- as.character(BNPP$depth) #set depth as character
BNPP <- BNPP %>%
  mutate(depth = replace(depth, depth == "20-Oct", "10-20")) %>%
  mutate(treatment_code = fct_recode(treatment, consDry = "consistentDry", control = "controlRain")) #make a column with shorter treatment names

#Calculate the aggregated BNPP in three levels of soil depth
BNPP1 <- BNPP %>%
  group_by(plot, subplot, treatment, treatment_code, shelterBlock, shelter, fall, spring) %>% #group data
  summarise(agg_BNPP = sum(bmass_g_m2))#sum BNPP 

#Filter 2015 ANPP 
ANPP1 <- ANPP %>%
  filter(year == 2015) %>%
  filter(subplot %in% c("B", "F", "G"))

#Join dataframes
Joined <- ANPP1 %>%
  left_join(BNPP1, by = c("plot", "subplot", "treatment", "shelterBlock", "shelter")) %>%
  left_join(CWM, by = c("plot", "subplot", "treatment", "shelterBlock", "shelter", "year", "spring", "fall")) %>%
  select(1,4,8,12,21:31) %>%
  mutate(total = weight_g_m + agg_BNPP)

#Standardize joined data
library(vegan)
stand_Joined_num <- decostand(Joined[,3:16], "standardize")
stand_Joined <- cbind(Joined[,1:2], stand_Joined_num)

##Backwards stepwise regression
model0 <- lm(total ~ CWM.Ht + CWM.SRLF , stand_Joined)
summary(model0)
AIC(model0)
