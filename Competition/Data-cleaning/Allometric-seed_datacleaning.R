library(gdata)
library(pbapply)
library(tidyverse)
library(broom)
library(readr)

## batch read in the different tabs of the seed weight xlsx file
source("https://gist.github.com/schaunwheeler/5825002/raw/3526a15b032c06392740e20b6c9a179add2cee49/xlsxToR.r")

alldat <- xlsxToR("~/Dropbox/ClimVar/Competition/Data/Competition_EnteredData/Competition_specimens_spring2017.xlsx", header = TRUE)
alldat2 <- alldat[-1]

# unlist (for unnecessary ease)
for (i in seq(alldat2))
  assign(paste0("df", i), alldat2[[i]]) 

# format each (this is clunky but whatevs)
laca_d <- df1 %>%
  tbl_df() %>%
  select(species, trt, seeds, wgt_g)  %>%
  filter(!is.na(wgt_g)) %>%
  mutate(wgt_g = as.numeric(as.character(wgt_g)))

laca_w <- df2 %>%
  select(species, trt, seeds, wgt_g)  %>%
  filter(!is.na(wgt_g)) %>%
  mutate(wgt_g = as.numeric(as.character(wgt_g)))

ggplot(laca_d, aes(wgt_g, seeds)) + geom_point() + geom_point(data = laca_w, aes(wgt_g, seeds), color = "blue")


avfa_d <- df3 %>%
  tbl_df() %>%
  select(species, trt, seeds, wgt_g)  %>%
  filter(!is.na(wgt_g)) %>%
  mutate(wgt_g = as.numeric(as.character(wgt_g)))

avfa_w <- df4 %>%
  select(species, trt, seeds, wgt_g)  %>%
  filter(!is.na(wgt_g)) %>%
  mutate(wgt_g = as.numeric(as.character(wgt_g)))

ggplot(avfa_d, aes(wgt_g, seeds)) + geom_point() + geom_point(data = avfa_w, aes(wgt_g, seeds), color = "blue")

brho_d <- df5 %>%
  tbl_df() %>%
  select(species, trt, seeds, wgt_g)  %>%
  filter(!is.na(wgt_g)) %>%
  mutate(wgt_g = as.numeric(as.character(wgt_g)))

brho_w <- df6 %>%
  select(species, trt, seeds, wgt_g)  %>%
  filter(!is.na(wgt_g)) %>%
  mutate(wgt_g = as.numeric(as.character(wgt_g)))

ggplot(brho_d, aes(wgt_g, seeds)) + geom_point() + geom_point(data = brho_w, aes(wgt_g, seeds), color = "blue")


vumy_d <- df7 %>%
  tbl_df() %>%
  select(species, trt, seeds, wgt_g)  %>%
  filter(!is.na(wgt_g)) %>%
  mutate(wgt_g = as.numeric(as.character(wgt_g)))

vumy_w <- df8 %>%
  select(species, trt, seeds, wgt_g)  %>%
  filter(!is.na(wgt_g)) %>%
  mutate(wgt_g = as.numeric(as.character(wgt_g)))

ggplot(vumy_d, aes(wgt_g, seeds)) + geom_point() + geom_point(data = vumy_w, aes(wgt_g, seeds), color = "blue")

# put it together!
allo.tog <- rbind(laca_d, laca_w, avfa_d, avfa_w, brho_d, brho_w, vumy_d, vumy_w)

# graph it all together!
ggplot(allo.tog, aes(x=wgt_g, y=seeds, color = trt)) + geom_point() + geom_smooth(method = "lm", se =F) + 
  facet_wrap(~species, scales = "free") + geom_smooth(data = allo.tog, aes(x=wgt_g, y=seeds), method = "lm", se = F, color = "black")

# long form way, don't care...
l <- lm(seeds~wgt_g, data = subset(allo.tog, species == "LACA"))
lacaout <- tidy(l) %>%
  mutate(species = "LACA")


l <- lm(seeds~wgt_g, data = subset(allo.tog, species == "AVFA"))
avfaout <- tidy(l) %>%
  mutate(species = "AVFA")

l <- lm(seeds~wgt_g, data = subset(allo.tog, species == "BRHO"))
brhoout <- tidy(l) %>%
  mutate(species = "BRHO")

l <- lm(seeds~wgt_g, data = subset(allo.tog, species == "VUMY"))
vumyout <- tidy(l) %>%
  mutate(species = "VUMY")

allo.out <- rbind(lacaout, avfaout, brhoout, vumyout)

allo.out_tomerge <- allo.out %>%
  select(species, term, estimate) %>%
  spread(term, estimate) %>%
  mutate(species = recode(species, LACA = "Lasthenia", AVFA = "Avena", BRHO = "Bromus", VUMY = "Vulpia")) %>%
  select(phyto = species,
         intercept = `(Intercept)`, 
         slope = wgt_g)


comp.all <- read_csv( "~/Dropbox/ClimVar/Competition/Data/Competition_CleanedData/ClimVar_Comp_combined-biomass-2.csv") 

comp.all2 <- left_join(comp.all, allo.out_tomerge) %>%
  mutate(seedsOut = intercept + slope*ind_weight_g)

write_csv(comp.all2, "~/Dropbox/ClimVar/Competition/Data/Competition_CleanedData/ClimVar_Comp_fecundity.csv")

## what I want to have is 

# take the average across blocks
comp.dat2 <- comp.all2 %>%
  tbl_df() %>%
  group_by(backgroundspp, backgrounddensity, phyto, falltreatment, seedsAdded) %>%
  # filter(!is.na(disturbance)) %>%
  summarize(mean_weight = mean(ind_weight_g, na.rm = T),
            mean_seed = mean(seedsOut, na.rm = T))


# some graphs of density and biomass
ggplot(background.density, aes(x=falltreatment, y=density)) + geom_boxplot() + facet_grid(backgrounddensity~backgroundspp)

ggplot(subset(comp.dat2, backgrounddensity != "none"), aes(x=backgrounddensity, y = mean_weight, color = backgroundspp, group = backgroundspp)) + geom_point() +
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

