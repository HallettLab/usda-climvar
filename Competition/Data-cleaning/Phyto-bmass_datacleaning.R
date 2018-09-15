require(gdata)
library(tidyverse)


# Read in data
phyto.dat <-read.xls("~/Dropbox/ClimVar/Competition/Data/Competition_EnteredData/Competition_phytometers_spring2017.xlsx", sheet=3, header=T, na.strings="#N/A!") %>%
  tbl_df() 

# Get individual weights
phyto.dat2 <- phyto.dat %>%
  mutate(drywgt.g = as.character(dry_wgt_g),
         drywgt.g = ifelse(drywgt.g == "NA", NA, drywgt.g),
         drywgt.g = as.numeric(drywgt.g),
         no.plants = as.character(stems),
         no.plants = ifelse(no.plants == "NA", NA, no.plants),
         no.plants = as.numeric(no.plants)) %>%
  mutate(ind.weight.g = drywgt.g/no.plants,
         phyto = phytometer) %>%
  mutate(disturbed = as.character(disturbed),
    disturbed = ifelse(disturbed ==" 17 stems\\", "", disturbed),
    disturbed = as.numeric(as.character(disturbed))) %>%
  dplyr::select(plot, background, phyto, ind.weight.g, disturbed) %>%
  separate(background, c("backgroundspp", "backgrounddensity"), sep = "_") %>%
  filter(!is.na(plot)) %>%
  mutate(backgroundspp = ifelse(backgroundspp == "x Control x", NA, backgroundspp))
  
# Read in keys  
shelter.key <- read.csv("~/Dropbox/ClimVar/Competition/Data/Competition_EnteredData/Shelter_key.csv") %>%
  tbl_df()

# Put it all together
phyto.bmass <-left_join(phyto.dat2, shelter.key) %>%
  mutate(falltreatment = "wet",
         falltreatment = ifelse(treatment == "fallDry" | treatment == "consistentDry", "dry", falltreatment)) %>%
  mutate(backgroundspp = recode(backgroundspp, AVFA = "Avena", BRHO = "Bromus", LACA = "Lasthenia",
                                ESCA = "Eschscholzia", TRHI = "Trifolium", VUMY = "Vulpia"),
         phyto = recode(phyto, AVFA = "Avena", BRHO = "Bromus", LACA = "Lasthenia",
                                ESCA = "Eschscholzia", TRHI = "Trifolium", VUMY = "Vulpia"),
         backgrounddensity = recode(backgrounddensity, LO = "low", HI = "high")) 


write.csv(phyto.bmass, "~/Dropbox/ClimVar/Competition/Data/Competition_CleanedData/ClimVar_Comp_phytometer-biomass-2.csv", row.names = F) %>%
  tbl_df()

