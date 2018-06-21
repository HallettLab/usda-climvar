# Read in data
background_stems <-read.xls("~/Dropbox/ClimVar/Competition/Data/Competition_EnteredData/Competition_background_spring2017.xlsx", sheet=2, header=T, na.strings="#N/A!")


# Convert for area
background_stems2 <- background_stems %>%
  tbl_df() %>%
  ## density should be multiplied by 2500 for the full plot; but because the
  ## seeds were concentrated in the center i'm estimating that it was in fact a smaller area (1/4 that size)
  mutate(density = stems.5cm/sample.area.cm2*625) %>%
  separate(background, c("backgroundspp", "backgrounddensity"), sep = "_") %>%
  dplyr::select(plot, backgroundspp, backgrounddensity, density)

# Read in keys  
shelter.key <- read.csv("~/Dropbox/ClimVar/Competition/Data/Competition_EnteredData/Shelter_key.csv") %>%
  tbl_df()

seed.key <- read.csv("~/Dropbox/ClimVar/Competition/Data/Competition_EnteredData/CompExpt_seedingKey.csv") %>%
  tbl_df()

# Join it all
back_stems_clean <-left_join(shelter.key, background_stems2) %>%
  mutate(falltreatment = "wet",
         falltreatment = ifelse(treatment == "fallDry" | treatment == "consistentDry", "dry", falltreatment))

back_stems_clean <-left_join(back_stems_clean, seed.key) %>%
  mutate(perPersist = density/seedsAdded)

## Not exporting yet because there is something weird about Vulpia - a lot of NAs


# Quick visual
ggplot((back_stems_clean), aes(x=backgrounddensity, y = perPersist, fill=treatment)) + geom_boxplot() + 
  facet_wrap( ~backgroundspp, scales = "free") #+ scale_y_log10()


ggplot((back_stems_clean), aes(x=backgrounddensity, y = perPersist, fill=falltreatment)) + geom_boxplot() + 
  facet_wrap( ~backgroundspp, scales = "free") #+ scale_y_log10()

