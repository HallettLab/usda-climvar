## Combine the background and phytometer data
# Either source Phyto-bmass_datacleaning and Background-stem-biomass_datacleaning
library(gdata)
library(tidyverse)

# Or read it in:
background <- read.csv("~/Dropbox/ClimVar/Competition/Data/Competition_CleanedData/ClimVar_Comp_background-biomass-2.csv")
phyto.bmass <- read.csv("~/Dropbox/ClimVar/Competition/Data/Competition_CleanedData/ClimVar_Comp_phytometer-biomass-2.csv")

# focus on background weights
background.phyto <- background %>%
  select(-density, -ind_flower, -seedsAdded, - tot_weight_g) %>%
  tbl_df()

# focus on background densities
background.density = background %>%
  tbl_df() %>%
  select(-ind_weight_g, -phyto) 

# clean up phytometer format
phyto <- phyto.bmass %>%
  tbl_df() %>%
  mutate(ind_weight_g = ind.weight.g) %>%
  select(-ind.weight.g,  -shelter, -disturbed) %>%
  mutate(backgroundspp = as.character(backgroundspp),
        backgrounddensity = as.character(backgrounddensity),
         phyto = as.character(phyto)) %>%
  mutate(backgroundspp = ifelse(is.na(backgroundspp), phyto, backgroundspp),
         backgrounddensity = ifelse(is.na(backgrounddensity), "none", backgrounddensity))

# # make a file of disturbed plots
disturbed <- phyto.bmass %>%
  select(plot, backgroundspp, backgrounddensity, disturbed) %>%
  unique()
# disturbed <- phyto.bmass %>%
#   select(plot, backgroundspp, backgrounddensity, chomped:rain_disturbed) %>%
#   gather(category, disturbance, chomped:rain_disturbed) %>%
#   mutate(disturbance = ifelse(disturbance != "", category, NA)) %>%
#   select(-category) %>%
#   filter(!is.na(disturbance))

# put it all together
comp.dat0 <- left_join(rbind(phyto, background.phyto), disturbed) 
comp.dat <- left_join(comp.dat0, background.density) %>%
  mutate(treatment=ordered(treatment, levels = c(consistentDry="consistentDry", springDry="springDry", fallDry="fallDry", controlRain="controlRain"))) %>%
  mutate(backgrounddensity=ordered(backgrounddensity, levels = c(none="none", low="low", high="high"))) %>%
  mutate(competitor_density = density) %>%
  select(-density)

# write.csv(comp.dat, "~/Dropbox/ClimVar/Competition/Data/Competition_CleanedData/ClimVar_Comp_combined-biomass-2.csv", row.names = F) %>%
#   tbl_df()

## QUICK VISUALS

# take the average across blocks
comp.dat2 <- comp.dat %>%
  group_by(backgroundspp, backgrounddensity, phyto, falltreatment) %>%
  # filter(!is.na(disturbance)) %>%
  summarize(mean_weight = mean(ind_weight_g, na.rm = T))


# some graphs
ggplot(background.density, aes(x=falltreatment, y=density)) + geom_boxplot() + facet_grid(backgrounddensity~backgroundspp)

ggplot(subset(comp.dat2, backgrounddensity != "none"), aes(x=backgrounddensity, y = mean_weight, color = backgroundspp, group = backgroundspp)) + geom_point() +
  geom_line() + 
  facet_grid(phyto~falltreatment, scales = "free")

ggplot(comp.dat2, aes(x=backgrounddensity, y = mean_weight, color = falltreatment, group = falltreatment)) + geom_point() +
  geom_line() + 
  facet_grid(phyto~backgroundspp, scales = "free")
