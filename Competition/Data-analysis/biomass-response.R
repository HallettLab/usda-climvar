# preliminary analysis of phytometer growth and fecundity under different competition scenarios
# authors: LMH, CTW
# initiated: Oct 2018 (modified over time)

# script purpose:
# ...

# notes:
# ...



# -- SETUP -----
rm(list = ls()) # clean environment
library(tidyverse)
options(stringsAsFactors = F)
theme_set(theme_bw())
na_vals = c(" ","", NA, "NA")

# set pathway to climvar dropbox competition data folder
datpath <- "~/Dropbox/ClimVar/Competition/Data/"

# read in data
# cleaned, combined competition dataset (cleaned background, phytometers, predicted vals, plot treatments)
comp.dat <- read.csv(paste0(datpath, "Competition_CleanedData/Competition_combined_clean.csv"),
                     na.strings = na_vals, strip.white = T)
# spp list
spplist <- read.csv(paste0(datpath, "Competition_SpeciesKey.csv"),
                    na.strings = na_vals, strip.white = T)


# -- VISUALS -----
# take the average across blocks
comp.dat2 <- comp.dat %>%
  group_by(background, bcode4, bdensity, phytometer, pcode4, falltreatment) %>%
  # filter(!is.na(disturbance)) %>%
  summarize(phyto_mean_weight = mean(p.ind.wgt.g, na.rm = T))


# some graphs
ggplot(comp.dat2, aes(x=falltreatment, y=bdensity)) + geom_boxplot() + facet_grid(backgrounddensity~backgroundspp)

ggplot(comp.dat2, aes(x=bdensity, y = phyto_mean_weight, color = falltreatment, group = falltreatment)) + geom_point() +
  geom_line() + 
  facet_grid(pcode4~bcode4, scales = "free")

ggplot(subset(comp.dat2, !is.na(bdensity)), aes(x=falltreatment, y = phyto_mean_weight, color = background, group = background)) + 
  geom_point() +
  geom_line() + 
  facet_grid(phytometer~bdensity, scales = "free")
ggsave("~/Dropbox/ClimVar/Competition/Figures/Exploratory/all-competition_May2019.pdf", width = 8, height = 10)


comp.dat2 %>%
  filter(bcode4 %in% spplist$code4[spplist$fxnl_grp == "Grass"],
         pcode4 %in% spplist$code4[spplist$fxnl_grp == "Grass"],
         !is.na(bdensity)) %>%
  ggplot(aes(x=falltreatment, y = phyto_mean_weight, color = background, group = background)) + geom_point() +
  geom_line() + 
  facet_grid(phytometer~bdensity, scales = "free")
ggsave("~/Dropbox/ClimVar/Competition/Figures/Exploratory/grass-competition_May2019.pdf", width = 8, height = 10)



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

