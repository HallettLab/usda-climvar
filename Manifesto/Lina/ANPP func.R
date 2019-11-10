ANPP_fun <- read.csv("./Dropbox/ClimVar/DATA/Plant_composition_data/ANPP/ANPP_CleanedData/ClimVar_ANPP-func-groups.csv", header = TRUE)
library(tidyverse)

# facultative mech: are grasses growing more in mixed plots than grass-only 
ANPP_grass <- ANPP_fun %>%
  filter(func%in% c("Grass","G")) %>%
  filter(year == 2015)
ggplot(ANPP_grass, aes(x = subplot, y = weight_g_m)) +
  geom_boxplot() +
  facet_grid(~harvest) +
  labs(y = "Grass ANPP weight (g/m2)")
