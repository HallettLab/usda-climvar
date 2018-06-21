require(gdata)
library(tidyverse)


# Read in data
background_stems <-read.xls("~/Dropbox/ClimVar/Competition/Data/Competition_EnteredData/Competition_background_spring2017.xlsx", sheet=2, header=T, na.strings="N/A") %>%
  tbl_df()

# Read in data
background_bio <-read.xls("~/Dropbox/ClimVar/Competition/Data/Competition_EnteredData/Competition_background_spring2017.xlsx", sheet=5, header=T, na.strings="#N/A!") %>%
  tbl_df() 

# Clean up biomass data
background_bio2 <- background_bio %>%
  mutate(dry_wgt_g = as.character(dry_wgt_g),
         dry_wgt_g = ifelse(dry_wgt_g == "NA", NA, dry_wgt_g),
         dry_wgt_g = as.numeric(dry_wgt_g)) %>%
  mutate(flower_bcount = flowers,
         stems_clipped = stems,
         harvest = ifelse(clip_date == "April 2017", "April", "May")) %>%
  dplyr::select(plot, background, plants_clipped, area_clipped_cm2, 
                stems_clipped, dry_wgt_g, flower_bcount, harvest) %>%
  separate(background, c("backgroundspp", "backgrounddensity"), sep = "_") %>%
  mutate(backgroundspp = recode(backgroundspp, AVFA = "Avena", BRHO = "Bromus", LACA = "Lasthenia",
                                ESCA = "Eschscholzia", TRHI = "Trifolium", VUMY = "Vulpia"),
         backgrounddensity = recode(backgrounddensity, LO = "low", HI = "high")) 
  


# Note: 1 low VUMY missing, think just was no material there, VUMY has stems for all others except 1 hi, 3 low, 3 high
# Bromus has no stem counts 
# tocheck <- background_bio2 %>%
#   filter(!is.na(area_clipped_cm2))


# Clean up stem data
background_stems2 <- background_stems %>%
  tbl_df() %>%
  mutate(area_counted_cm2 = sample.area.cm2,
         stems_counted = stems.5cm,
         flowers_counted = flowers,
        stalks_counted = stalks) %>%
  separate(background, c("backgroundspp", "backgrounddensity"), sep = "_") %>%
  dplyr::select(plot, backgroundspp, backgrounddensity, stems_counted, area_counted_cm2, flowers_counted, stalks_counted)

# put it all together
tog <- left_join(background_bio2, background_stems2) %>%
  mutate(stems_counted1 = ifelse(is.na(stems_counted), stems_clipped, stems_counted),
         stems_clipped1 = ifelse(is.na(stems_clipped), plants_clipped, stems_clipped),
         stems_clipped1 = ifelse(is.na(stems_clipped1) & area_counted_cm2 == 25, 
                                       stems_counted, stems_clipped1)) %>%
  mutate(ind_weight_g = dry_wgt_g/stems_clipped1) %>%
  ## density should be multiplied by 2500 for the full plot; but because the
  ## seeds were concentrated in the center i'm estimating that it was in fact a smaller area (1/4 that size)
  mutate(density = stems_counted1/area_counted_cm2*625) %>%
  mutate(flowers = ifelse(is.na(flowers_counted), flower_bcount, flowers_counted),
         ind_flower = flowers/stems_counted1) %>%
  select(plot, backgroundspp, backgrounddensity, harvest, ind_weight_g, density, ind_flower)


# Read in keys  
shelter.key <- read.csv("~/Dropbox/ClimVar/Competition/Data/Competition_EnteredData/Shelter_key.csv") %>%
  tbl_df()

seed.key <- read.csv("~/Dropbox/ClimVar/Competition/Data/Competition_EnteredData/CompExpt_seedingKey.csv") %>%
  tbl_df()

# Join it all
background_clean <-left_join(shelter.key, tog) %>%
  mutate(falltreatment = "wet",
         falltreatment = ifelse(treatment == "fallDry" | treatment == "consistentDry", "dry", falltreatment))

background_clean <-left_join(background_clean, seed.key) %>%
  mutate(perPersist = density/seedsAdded)

# Deal with the duplicate measures through averaging
background <- background_clean %>%
  mutate(phyto = backgroundspp) %>%
  select(plot, treatment:backgroundspp, backgrounddensity, phyto, falltreatment, ind_weight_g, density, seedsAdded) %>%
  group_by(plot, treatment, backgroundspp, backgrounddensity, phyto, falltreatment, shelterBlock) %>%
  summarize(ind_weight_g = mean(ind_weight_g), density = mean(density), seedsAdded = mean(seedsAdded))

write.csv(background, "~/Dropbox/ClimVar/Competition/Data/Competition_CleanedData/ClimVar_Comp_background-biomass.csv", row.names = F) %>%
  tbl_df()

