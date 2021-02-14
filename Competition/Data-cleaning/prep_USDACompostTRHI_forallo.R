# prep USDA Compost TRHI biomass-fecund data for USDA ClimVar competition allometry
# author(s): ctw, jan 2021

# script purpose:
## read in USDA Compost competition TRHI data (collected in field by NS Spring 2020, processed @ Eugene)
## read in USDA Compost treatment key (CTW copied over to )
## subset to control plots (? need to screen and consider)..
# > screen for any differences between wet, control, and drought plots. only select TRHI from unammended/unfert plots.
# > 4 blocks x 1 nutrient control x 3 ppt trts = potentially 12 specimens, unless there are diffs by ppt
# > background TRHI is also included
## prep what's needed to run through allometric-seed_datacleaning.R with other species dat (apply similar colnames as other datasets)
## write out prepped dataset to ClimVar Dropbox/Competition/Data/EnteredData/outside_project_data_related/


# notes:
## TRHI clipped in April 2017 in USDA ClimVar Competition, too early for seed development
## to infill, borrowing data from USDA Compost competition experiment.
## For allometric relationships in USDA ClimVar, individuals of each species were selected on size gradient (small to large), 10 from "dry" plots (consistent dry or fall dry drought) and 10 from "wet" plots (control or spring dry)
## > CTW pulled individuals from subplots where they were growing as the background species (plenty to choose from), typically from low density plots, but selected haphazardly prioritizing healthy/representative plants to capture size gradient
## USDA Compost experiment is not set up exactly like that and no specimens were harvested specifically for allometry.. TRHI read into this script are the actual experiment data

# AM says 999 in phyto = a background specimen



# -- SETUP -----
rm(list = ls())
library(tidyverse)
theme_set(theme_test())
options(stringsAsFactors = F)
na_vals = c("", " ", NA, "NA")
datpath <- "~/Dropbox/ClimVar/Competition/Data/" # path for CTW + LMH, can add for other users as needed

# list related compost files 
compostfiles <- list.files(paste0(datpath, "Competition_EnteredData/outside_project_data_related/"),full.names = T)

# read in compost dat
# trhi data per stem (individual)
trhi_ind <- read.csv(compostfiles[grep("_individ", compostfiles)], na.strings = na_vals)
# trhi data per sampling point (aggregated data, at level of phytometer or background)
trhi_agg <- read.csv(compostfiles[grep("_agg", compostfiles)], na.strings = na_vals)
# trtment key
compost_key <- read.csv(compostfiles[grep("TreatmentK", compostfiles)], na.strings = na_vals)
# compost layout
layout <- read.csv(compostfiles[grep("layout", compostfiles)], na.strings = na_vals)


# review import
## individuals
str(trhi_ind)
summary(trhi_ind)
## aggregated data
str(trhi_agg)
## treatment and design info
str(compost_key)
str(layout)
summary(layout) # no empty rows, good.


# note TRHI in plot 10, pos 6 misentered with subplot = 10;
# only two subplots in plot 10 with TRHI @ 6 = subplot 5 and 6
# subplot 5 is present in dataset so changing subplot 10 to 6
# > note there are 5 stems entered for this phyto pos.. there should only be 3 max and NS wrote 3 on datasheet so not sure where extra 2 came from
# > keeping 5 stems as is bc no way to correct since plant biomass data pooled
trhi_ind$subplot[trhi_ind$plot == 10 & trhi_ind$phyto == 6] # to confirm only subplots 5 and 10 present
layout[layout$plot == 10 & layout$phyto == "TRHI", c("subplot", "phytonum")] # confirm pos 6 in subplots 5 and 6
trhi_ind$subplot[trhi_ind$plot == 10 & trhi_ind$phyto == 6 & trhi_ind$subplot == 10] <- 6
trhi_agg$subplot[trhi_agg$plot == 10 & trhi_agg$phyto == 6 & trhi_agg$subplot == 10] <- 6


# -- SCREEN DATA -----
# join trt info to trhi dat, drop fert plots, screen for diffs btwn ppt trts (cuidado con phyto v. background!)

# 1. individual data screen -----
trhi_ind <- left_join(trhi_ind, compost_key) %>%
  # join background spp
  left_join(distinct(layout[c("plot", "subplot", "background")]))
# does everyone have a background val?
summary(is.na(trhi_ind$background)) # good
# do all 999s match up with background == TRHI?
distinct(subset(trhi_ind, phyto == 999 | background == "TRHI", c("phyto", "background"))) # good

# quick viz to see who has samples v. not and screen for obvi diffs before drop amended plots
# who haz?
trhi_ind %>%
  mutate(present = as.numeric(!is.na(stem))) %>%
  ggplot(aes(factor(block),factor(present))) +
  geom_point(aes(col = phyto == 999), position = position_jitter(width = 0.2, height = 0.1), alpha = 0.5) +
  scale_color_discrete(name = "Background") +
  facet_grid(ppt_trt~nut_trt)
# looks like enough background stems in control nutrient.. altho no stems from block 1 but that's ok. we didn't sample allo specs from all blocks in ClimVar 

summary(!is.na(trhi_ind$stem)) # no survivorship or missing between 29 phyto positions and background
with(trhi_ind, sapply(split(stem, fulltrt), function(x) summary(!is.na(x)))) 
# > there are 87 specimens from nutrient control plots


# how many stems from background or control phytos?
with(subset(trhi_ind, !is.na(stem) & background %in% c("TRHI", "Control")), sapply(split(stem, fulltrt), sum))
# how many stems in nutrient control per drought trt?
with(subset(trhi_ind, !is.na(stem) & nut_trt == "N" & background %in% c("TRHI", "Control")), sapply(split(stem, paste(background, ppt_trt)), sum))
# there is enough in background to match what we did for climvar, don't need to use control plot phytos
# > do still check for diffs between biomass and seed production

# > note.. we don't have per stem plant biomass.. poop. may need to use aggregated data after all..
# > maybe compare individual stem to seed mass relationship and then agg stem to agg seed mass relationship to get a sense of how representative agg biomass:seed count will be of indiv biomass:seed count? 

# sum indiv data by stem (flower and seed count split per stem)
# > replace "no" in notes with 0 in appropriate column first
trhi_indstems <- trhi_ind %>%
  replace_na(list(stem = 0, flower = 0, num_seeds = 0, seed_mass_g = 0)) %>%
  #group_by(names(trhi_ind)[!names(trhi_ind) %in% c("notes", "flower", "num_seeds", "seed_mass_g")]) %>%
  group_by(id, plot, subplot, background, phyto, stem) %>%
  summarise(flowers = sum(flower),
            seeds = sum(num_seeds),
            seed_mass_g = sum(seed_mass_g)) %>%
  ungroup() %>%
  # rejoin plot treatment info
  left_join(compost_key)

# out of curiosity, what is the mean per-stem seed count and seed mass?
with(trhi_indstems, sapply(split(seeds, fulltrt), mean)) # productive in fertilizer
with(trhi_indstems, sapply(split(seed_mass_g, fulltrt), mean)) # interesting.. moisture influence suggested

ggplot(subset(trhi_indstems, stem >0), aes(block, seeds)) +
  stat_summary() +
  facet_grid(ppt_trt~nut_trt)

ggplot(subset(trhi_indstems, stem >0), aes(block, seed_mass_g)) +
  stat_summary() +
  facet_grid(ppt_trt~nut_trt)

ggplot(subset(trhi_indstems, stem >0), aes(ppt_trt, seed_mass_g, col = nut_trt)) +
  stat_summary(position = position_jitter(width = 0.2, height = 0)) # will be interesting for AS

# check breakdown of individual seeds to mass
ggplot(subset(trhi_indstems), aes(seeds, seed_mass_g)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_grid(ppt_trt~nut_trt)

# nutrient control only for heads up on any diffs between ppts treats
ggplot(subset(trhi_indstems), aes(seeds, seed_mass_g, col = ppt_trt)) +
  geom_point() +
  geom_smooth(aes(fill = ppt_trt), method = "lm", formula = "y ~ 0 + x") +
  facet_wrap(~nut_trt) # interesting..

# > def don't use samples from C or F for climvar allo relationship
# try to only use XC and D in control nutrient plot if possible (if there are enough pooled samples)



# 2. aggregate data screen ------
# add treatment and design info to agg data
trhi_agg <- left_join(trhi_agg, compost_key) %>%
  # join background spp
  left_join(distinct(layout[c("plot", "subplot", "background")]))
# does everyone have a background val?
summary(is.na(trhi_agg$background)) # good

# pull nutrient control samples
trhi_sub <- subset(trhi_agg, nut_trt == "N")
  
# check for diffs in TRHI background vs. TRHI phyto v comp vs. TRHI phyto in control
mutate(trhi_sub, growgroup = ifelse(grepl("Con|TRHI", background), background, "Competitor")) %>%
  ggplot(aes(biomass_tot, seeds_tot, col = growgroup)) +
  geom_point() +
  geom_smooth(aes(fill = growgroup), method = "lm", formula = "y ~ 0 + x") +
  facet_wrap(~ppt_trt)
# > I don't see much of a difference between seed count ~ biomass by growing conditions (whether sample is TRHI background spec, TRHI phyto in comp or TRHI phyto in control)
# think there are enough data points in D and in XC to drop watered treatment (doing that aligns better with climvar methods)
# to be sure, check counts
with(subset(trhi_sub, !is.na(num_stems)), sapply(split(num_stems, ppt_trt), length)) # yup. we just took 10 per "dry" (fall drought/consistent drought) and "wet" (ambient, spring drought)





# -- PREP DATA FOR ALLO SCRIPTS -----
# create dataset structured and with headers to match other allometric spec data read in
# cols to keep: species, ppt_trt, seeds_tot, biomass_tot, plot, background (just so it's there in case anyone curious about provenance)
# > keep flower count too just in case
trhi_sub <- subset(trhi_sub, ppt_trt != "W") %>% # drop watered plots
  # create col for provenance of data
  mutate(project = "USDA Compost", harvest_date = "May 2020") %>%
  select(project, harvest_date, plot, nut_trt, ppt_trt, background, id, species, id,  num_stems, num_flowers, seeds_tot, biomass_tot) %>%
  # drop empty
  subset(!is.na(num_stems) & !is.na(biomass_tot) ) # two plots that have stems present noted are missing biomass and seed dat
 
# .. should we take mean so output reps mass and seed count per stem like the rest?.... hm 
# does it change anything?
# as is
with(subset(trhi_sub, ppt_trt == "XC"), summary(lm(seeds_tot ~ 0 + biomass_tot)))
with(subset(trhi_sub, ppt_trt == "D"), summary(lm(seeds_tot ~ 0 + biomass_tot)))
with(trhi_sub, summary(lm(seeds_tot ~ 0 + biomass_tot*ppt_trt))) # 11 additional seeds per additional gram in ambient moisture compared to drought

# if things were at per stem level..
trhi_sub$perstem_flowers <- with(trhi_sub, round(num_flowers/num_stems, 3))
trhi_sub$perstem_seeds <- with(trhi_sub, round(seeds_tot/num_stems, 3))
trhi_sub$perstem_mass <- with(trhi_sub, round(biomass_tot/num_stems, 3))

# derived per stem
with(subset(trhi_sub, ppt_trt == "XC"), summary(lm(perstem_seeds ~ 0 + perstem_mass))) # as is = +107, per stem = +102
with(subset(trhi_sub, ppt_trt == "D"), summary(lm(perstem_seeds ~ 0 + perstem_mass)))  # as is = +109, per stem = +112
with(trhi_sub, summary(lm(perstem_seeds ~ 0 + perstem_mass*ppt_trt))) # as is = +11 seeds in control, per stem = +10 .. close enough
# > slight differences

# write out everything so it's all available, that was go just modify in allo script if we decide to use switch which numbers we use
write_csv(trhi_sub, paste0(datpath, "Competition_EnteredData/outside_project_data_related/USDACompost_Competition_TRHIspecimens_spring2020_prepped.csv"))
