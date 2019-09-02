# Combine the clean background and phytometer data with predicted seeds
# authors: LMH, CTW
# initiated: Oct 2018 (will be modified as develop analyses)

# script purpose:
# compile all in one dataset, by competition plot, by plot, by treatment:
# phytometer, phytometer per plant weight, phytometer per plant flowers (LACA only)
# phytometer seeds projected (based on allometric rates)
# background competitor, density treatment, background anpp, background flowers, background density, background seeds added 

# notes:
# we don't need background projected fecundity right? or do we want that?
# as of 2019-05-02, still need to reconcile may forb flower counts vs. april for TRHI and ESCA
# also LACA on anpp samples vs LACA in stem counts



# -- SETUP -----
rm(list = ls()) # clean environment
library(tidyverse)
options(stringsAsFactors = F)
theme_set(theme_bw())
na_vals = c(" ", "", NA, "NA")

# set pathway to competition data on dropbox
datpath <- "~/Dropbox/ClimVar/Competition/Data/"

# read in cleaned data:
# cleaned background competitor density and anpp
background <- read.csv(paste0(datpath,"Competition_CleanedData/Competition_background_clean.csv"),
                       na.strings = na_vals, strip.white = T)
# cleaned phytometer stem count and anpp with predictions
phytodat <- read.csv(paste0(datpath,"Competition_CleanedData/Competition_phytometers_predicted.csv"),
                     na.strings = na_vals, strip.white = T)
# allometric seed rates (don't think needed, but just in case)
allodat <- read.csv(paste0(datpath,"Competition_CleanedData/Competition_allometric_clean.csv"),
                    na.strings = na_vals, strip.white = T)
# species lookup
spplist <-read.csv(paste0(datpath,"Competition_SpeciesKey.csv"),
                   na.strings = na_vals, strip.white = T)

# -- COMBINE PHYTO AND BACKGROUND DATASETS -----
# look at colnames and structure
glimpse(background)
glimpse(phytodat)
glimpse(allodat)
glimpse(spplist)

nrow(phytodat) #1056
# combine background and phytometers (there should only be 1056 when combine)
phytocomp <- left_join(phytodat, background) 
summary(phytocomp)
# which are the NAs in background density? (should only be control plots)
unique(phytocomp$backgroundspp[is.na(phytocomp$insitu_plot_bdensity)]) # Control (good)
# what are the plots with no background anpp?
unique(phytocomp$backgroundspp[is.na(phytocomp$plot_banpp)]) # Control (good), ESCA and LACA
# what are the ESCA and LACA plots?
background[is.na(background$plot_banpp),] # these are the two obs where stems recorded but anpp samples missing
# 2 x 5 phytometers = 10 plots, so that's the 10 extra NA rows in the combined dataset background anpp
phytocomp[(is.na(phytocomp$plot_banpp) & phytocomp$backgroundspp != "Control"),] #yup


# -- ADD IN ALLOMETRIC DATA, PROJECT FECUNDITY -----
phytocomp2 <- left_join(phytocomp, allodat[c("phytometer", "intercept", "slope")], by = "phytometer")
#reorganize columns
phytocomp2 <- phytocomp2[,c(colnames(background), colnames(phytodat)[4:11], "intercept", "slope", colnames(phytodat)[16:ncol(phytodat)])] %>%
  rename(background = species,
         bdensity= backgrounddensity,
         bcode4 = backgroundspp) %>%
  left_join(spplist[c("code4", "species")], by = c("phytometer" = "code4")) %>%
  rename(pcode4 = phytometer,
         phytometer = species) %>%
  dplyr::select(plot, falltreatment, treatment:shelter, 
                # background species columns
                background, bcode4:insitu_bdisturbed, seedsAdded:bdensity_flag, 
                # phytometer and predicted seeds out columns
                phytometer, pcode4:uprPI.95)

# clean up regression cols (4 signif digits)
combined <- phytocomp2 
combined[,grep("intercept",colnames(combined)):ncol(combined)] <- sapply(combined[,grep("intercept",colnames(combined)):ncol(combined)], function(x) round(x,4)) 

# -- FINISHING -----
# write out
write.csv(combined, paste0(datpath, "Competition_CleanedData/Competition_combined_clean.csv"), row.names = F)



##### OLD CODE #####
# # focus on background weights
# background.phyto <- background %>%
#   select(-density, -ind_flower, -seedsAdded, - tot_weight_g) %>%
#   tbl_df()
# 
# # focus on background densities
# background.density = background %>%
#   tbl_df() %>%
#   select(-ind_weight_g, -phyto) 
# 
# # clean up phytometer format
# phyto <- phyto.bmass %>%
#   tbl_df() %>%
#   mutate(ind_weight_g = ind.weight.g) %>%
#   select(-ind.weight.g,  -shelter, -disturbed) %>%
#   mutate(backgroundspp = as.character(backgroundspp),
#         backgrounddensity = as.character(backgrounddensity),
#          phyto = as.character(phyto)) %>%
#   mutate(backgroundspp = ifelse(is.na(backgroundspp), phyto, backgroundspp),
#          backgrounddensity = ifelse(is.na(backgrounddensity), "none", backgrounddensity))
# 
# # # make a file of disturbed plots
# disturbed <- phyto.bmass %>%
#   select(plot, backgroundspp, backgrounddensity, disturbed) %>%
#   unique()
# # disturbed <- phyto.bmass %>%
# #   select(plot, backgroundspp, backgrounddensity, chomped:rain_disturbed) %>%
# #   gather(category, disturbance, chomped:rain_disturbed) %>%
# #   mutate(disturbance = ifelse(disturbance != "", category, NA)) %>%
# #   select(-category) %>%
# #   filter(!is.na(disturbance))
# 
# # put it all together
# comp.dat0 <- left_join(rbind(phyto, background.phyto), disturbed) 
# comp.dat <- left_join(comp.dat0, background.density) %>%
#   mutate(treatment=ordered(treatment, levels = c(consistentDry="consistentDry", springDry="springDry", fallDry="fallDry", controlRain="controlRain"))) %>%
#   mutate(backgrounddensity=ordered(backgrounddensity, levels = c(none="none", low="low", high="high"))) %>%
#   mutate(competitor_density = density) %>%
#   select(-density)

# write.csv(comp.dat, "~/Dropbox/ClimVar/Competition/Data/Competition_CleanedData/ClimVar_Comp_combined-biomass-2.csv", row.names = F) %>%
#   tbl_df()


# # join allometric seeding rates with biomass
# comp.all2 <- left_join(comp.all, allo.out_tomerge) %>%
#   mutate(seedsOut = intercept + slope*ind_weight_g)
# 
# write_csv(comp.all2, paste0(datpath, "Competition_CleanedData/ClimVar_Comp_fecundity.csv"))