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
                     na.strings = na_vals, strip.white = T) %>%
  tbl_df()

# and bring back allodat
allodat <- read.csv(paste0(datpath,"Competition_CleanedData/Competition_allometric_clean.csv"),
                    na.strings = na_vals, strip.white = T) %>%
  rename(bcode4 = phytometer)

phyto.dat <- comp.dat %>%
  select(plot:bdensity, insitu_bdisturbed, phytometer, pcode4, insitu_pstems, p_totwgt_seedfit ) %>%
  filter(!is.na(p_totwgt_seedfit)) %>%
  mutate(p_totwgt_seedfit = ifelse(insitu_pstems == 0, 0, p_totwgt_seedfit)) %>%
  select(-phytometer, -bcode4) %>%
  rename(stemsIn = insitu_pstems, seedsOut = p_totwgt_seedfit, disturbed = insitu_bdisturbed,
         block = shelterBlock) 

back.dat0 <- comp.dat %>%
  filter(!is.na(background)) %>%
  select(plot:bdensity,insitu_bdisturbed, b.ind.wgt.g, seedsAdded, insitu_plot_bdensity) %>%
  unique() 

back.dat <- left_join(back.dat0, allodat) %>%
  mutate(p_totwgt_seedfit = intercept + b.ind.wgt.g*insitu_plot_bdensity*slope) %>%
  mutate(p_totwgt_seedfit = ifelse(insitu_plot_bdensity == 0, 0, p_totwgt_seedfit)) %>%
  select(plot:bdensity,insitu_bdisturbed,seedsAdded, p_totwgt_seedfit, insitu_plot_bdensity) %>%
  rename(seedsIn = seedsAdded, seedsOut = p_totwgt_seedfit, pcode4 = bcode4,
         stemsIn = insitu_plot_bdensity, disturbed = insitu_bdisturbed,
         block = shelterBlock) %>%
  filter(!is.na(seedsOut)) %>%
  select(plot, falltreatment, treatment, block, shelter, background, bdensity, disturbed, pcode4, stemsIn, seedsOut)


## check that always the right number of seeds added 
ggplot(back.dat, aes(x= seedsIn, y = stemsIn)) + geom_point() + geom_abline(intercept = 0, slope = 1)

# calculate stems in
stems.in <- rbind(phyto.dat, back.dat) %>%
  filter(pcode4 != "TRHI") %>%
  select(-seedsOut) %>%
  mutate(varnew = paste(pcode4, "stemsIn", sep = "_")) %>%
  select(-pcode4) %>%
  spread(varnew, stemsIn, fill = 0) %>%
  select(plot, falltreatment, treatment, block, disturbed:VUMY_stemsIn)

# calculate seeds out 
seeds.out <-  rbind(phyto.dat, back.dat) %>%
  filter(pcode4 != "TRHI") %>%
  select(plot, falltreatment, treatment, block, seedsOut, pcode4, background, bdensity) %>%
  rename(species = pcode4)

stemsin.seedsout <- left_join(seeds.out, stems.in)

rm(allodat, back.dat, back.dat0, comp.dat, phyto.dat, seeds.out, stems.in)

