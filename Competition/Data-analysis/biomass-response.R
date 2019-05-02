## Some preliminary scripts analyzing phytometer growth under different competition scenarios

library(tidyverse)
library(readr)

background <- read_csv("~/Dropbox/ClimVar/Competition/Data/Competition_CleanedData/ClimVar_Comp_background-biomass-2.csv")
comp.dat <- read_csv("~/Dropbox/ClimVar/Competition/Data/Competition_CleanedData/ClimVar_Comp_combined-biomass.csv")



# take the average across blocks
comp.dat2 <- comp.dat %>%
  group_by(backgroundspp, backgrounddensity, phyto, falltreatment) %>%
  # filter(!is.na(disturbance)) %>%
  summarize(mean_weight = mean(ind_weight_g, na.rm = T))


# some graphs
ggplot(background, aes(x=falltreatment, y=density)) + geom_boxplot() + facet_grid(backgrounddensity~backgroundspp)

ggplot(comp.dat2, aes(x=backgrounddensity, y = mean_weight, color = falltreatment, group = falltreatment)) + geom_point() +
  geom_line() + 
  facet_grid(phyto~backgroundspp, scales = "free")

ggplot(subset(comp.dat2, backgrounddensity != "none"), aes(x=falltreatment, y = mean_weight, color = backgroundspp, group = backgroundspp)) + geom_point() +
  geom_line() + 
  facet_grid(phyto~backgrounddensity, scales = "free")
ggsave("~/Dropbox/ClimVar/Competition/Figures/Exploratory/all-competition.pdf", width = 8, height = 10)


grass.dat <- comp.dat2 %>%
  filter(backgroundspp%in%c("Avena", "Bromus", "Vulpia"),
         phyto%in%c("Avena", "Bromus", "Vulpia")) 

ggplot(subset(grass.dat, backgrounddensity != "none"), aes(x=falltreatment, y = mean_weight, color = backgroundspp, group = backgroundspp)) + geom_point() +
  geom_line() + 
  facet_grid(phyto~backgrounddensity, scales = "free")
ggsave("~/Dropbox/ClimVar/Competition/Figures/Exploratory/grass-competition.pdf", width = 8, height = 10)



# -- VISUALIZE TRENDS -----
## what I want to have is 
# take the average across blocks
comp.dat2 <- ungroup(comp.all2) %>%
  group_by(backgroundspp, backgrounddensity, phyto, falltreatment, seedsAdded) %>%
  # filter(!is.na(disturbance)) %>%
  summarize(mean_weight = mean(ind_weight_g, na.rm = T),
            mean_seed = mean(seedsOut, na.rm = T)) %>%
  as.data.frame()


# some graphs of density and biomass
ggplot(comp.dat2, aes(x=falltreatment, y=density)) + 
  geom_boxplot() + facet_grid(backgrounddensity~backgroundspp)

ggplot(subset(comp.dat2, backgrounddensity != "none"), 
       aes(x=backgrounddensity, y = mean_weight, color = backgroundspp, group = backgroundspp)) + geom_point() +
  geom_line() + 
  facet_grid(phyto~falltreatment, scales = "free")

ggplot(comp.dat2, aes(x=backgrounddensity, y = mean_weight, color = falltreatment, group = falltreatment)) + geom_point() +
  geom_line() + 
  facet_grid(phyto~backgroundspp, scales = "free")


# some graphs of seeds
ggplot(comp.dat2, aes(x=falltreatment, y=mean_seed/seedsAdded)) + geom_boxplot() + 
  facet_grid(backgrounddensity~backgroundspp, scales = "free") 


ggplot(subset(comp.dat2, backgrounddensity != "none"), aes(x=backgrounddensity, y = mean_seed, color = backgroundspp, group = backgroundspp)) + geom_point() +
  geom_line() + 
  facet_grid(phyto~falltreatment, scales = "free")

ggplot(comp.dat2, aes(x=backgrounddensity, y = mean_seed, color = falltreatment, group = falltreatment)) + geom_point() +
  geom_line() + 
  facet_grid(phyto~backgroundspp, scales = "free")

# seems like:
# Avena is always the superior competitor
# Brome vs Vulpia is weather dependent: Brome can increase when rare under wet, Vuplia under dry
# Lasthenia is always the weakest competitor

