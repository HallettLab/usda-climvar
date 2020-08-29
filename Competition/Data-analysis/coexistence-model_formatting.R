# -- SETUP -----
rm(list = ls()) # clean environment
library(tidyverse)
options(stringsAsFactors = F)
theme_set(theme_bw())
na_vals = c(" ","", NA, "NA")

# set pathway to climvar dropbox competition data folder

# specify dropbox pathway (varies by user -- tweak when share with Caitlin)
if(file.exists("~/Dropbox/Shared/ClimVar/Competition/Data/")){
  # LGS
  datpath <- "~/Dropbox/Shared/ClimVar/Competition/Data/"
  # LMH
}else{
  datpath <- "~/Dropbox/ClimVar/Competition/Data/"
}

### FECUNDITY ### 
# read in data
# cleaned, combined competition dataset (cleaned background, phytometers, predicted vals, plot treatments)
comp.dat <- read.csv(paste0(datpath, "Competition_CleanedData/Competition_combined_clean.csv"),
                     na.strings = na_vals, strip.white = T) %>%
  tbl_df()

# and bring back allodat
allodat <- read.csv(paste0(datpath,"Competition_CleanedData/Competition_allometric_clean.csv"),
                    na.strings = na_vals, strip.white = T)

# isolate the phytometer data
phyto.dat <- comp.dat %>%
  select(plot:bdensity, insitu_bdisturbed, phytometer, pcode4, insitu_pstems, p_seedfit, pfit_source) %>%
  filter(!is.na(p_seedfit)) %>%
  mutate(p_seedfit = ifelse(insitu_pstems == 0, 0, p_seedfit)) %>%
  select(-phytometer, -bcode4, -shelter) %>%
  rename(stemsIn = insitu_pstems, seedsOut = p_seedfit, disturbed = insitu_bdisturbed,
         block = shelterBlock) %>%
  mutate(seedsIn = ifelse(pcode4 == "AVFA" | pcode4 == "TRHI", 10, 12),
         seedsIn = ifelse(pcode4 == "LACA" | pcode4 == "ESCA", 15, seedsIn),
         seedsIn = ifelse(stemsIn > seedsIn, stemsIn, seedsIn)) %>%
  # code to subset seedfit based on individual pwgt or total pwgt
  subset(grepl("ind", pfit_source)) %>% # uncomment this line if want based on indidivual phyto wgt
  #subset(grepl("tot", pfit_source)) %>% # uncomment this line if want based on total phyto wgt
  select(-pfit_source)

# isolate the background data
back.dat0 <- comp.dat %>%
  filter(!is.na(background)) %>%
  select(plot:bdensity,insitu_bdisturbed, b.ind.wgt.g, seedsAdded, insitu_plot_bdensity) %>%
  unique() 

back.dat <- left_join(back.dat0, allodat, by = c("bcode4" = "species")) %>%
  mutate(p_seedfit = intercept + b.ind.wgt.g*insitu_plot_bdensity*slope) %>%
  mutate(p_seedfit = ifelse(insitu_plot_bdensity == 0, 0, p_seedfit)) %>%
  select(plot:bdensity,insitu_bdisturbed,seedsAdded, p_seedfit, insitu_plot_bdensity) %>%
  rename(seedsIn = seedsAdded, seedsOut = p_seedfit, pcode4 = bcode4,
         stemsIn = insitu_plot_bdensity, disturbed = insitu_bdisturbed,
         block = shelterBlock) %>%
  filter(!is.na(seedsOut)) %>%
  select(plot, falltreatment, treatment, block, background, bdensity, disturbed, pcode4, stemsIn, seedsOut, seedsIn)


## check that always the right number of seeds added 
ggplot(back.dat, aes(x= seedsIn, y = stemsIn, col = pcode4, shape = bdensity)) + geom_point() + geom_abline(intercept = 0, slope = 1)

# calculate stems in
stems.in <- rbind(phyto.dat, back.dat) %>%
  filter(pcode4 != "TRHI") %>%
  select(-seedsOut, -seedsIn) %>%
  mutate(varnew = paste(pcode4, "stemsIn", sep = "_")) %>%
  select(-pcode4) %>%
  spread(varnew, stemsIn, fill = 0) %>%
  select(plot, falltreatment, treatment, block, disturbed:VUMY_stemsIn, background, bdensity)


# calculate seeds in
seeds.added <- rbind(phyto.dat, back.dat) %>%
  filter(pcode4 != "TRHI") %>%
  select( -stemsIn, -seedsOut) %>%
  mutate(varnew = paste(pcode4, "seedsIn", sep = "_")) %>%
  select(-pcode4) %>%
  spread(varnew, seedsIn, fill = 0) %>%
  select(plot, falltreatment, treatment, block, disturbed:VUMY_seedsIn, background, bdensity)


# calculate seeds out 
seeds.out <-  rbind(phyto.dat, back.dat) %>%
  filter(pcode4 != "TRHI") %>%
  select(plot, falltreatment, treatment, block, seedsOut, pcode4, background, bdensity) %>%
  rename(species = pcode4)

# stemsin.seedsout for model
stemsin.seedsout <- left_join(seeds.out, stems.in)


# seedsin.seedsout for model
seedsin.seedsout <- left_join(seeds.out, seeds.added)





# If we want to break up by recruitment windows ---------------------------


### RECRUITMENT ###
recruit.dat <- read.csv(paste0(datpath, "Competition_CleanedData/ClimVar_Comp_background-recruitment.csv"),
                     na.strings = na_vals, strip.white = T) %>%
  tbl_df() %>%
  rename(block = shelterBlock,
         genus = backgroundspp, 
         bdensity = backgrounddensity,
         stemsRecruit = density) %>%
  select(plot, falltreatment, treatment, block, genus, bdensity, stemsRecruit) 



tomerge <- back.dat %>%
  select(background) %>%
  unique() %>%
  mutate(background2 = background) %>%
  separate(background2, c("genus", "species"))

back.dat <- left_join(back.dat, tomerge)

### SURVIVAL ###
together <- left_join(recruit.dat, back.dat) %>%
  filter(!is.na(background), background != "Trifolium hirtum") %>%
  mutate(perRecruit = stemsRecruit/seedsIn)



ggplot(together, aes(x= seedsIn, y = stemsRecruit, color = falltreatment)) + geom_point() +
  facet_wrap(~background,  scales = "free") + 
  geom_abline(intercept = 0, slope = 1)

ggplot(together, aes(x= genus, y = perRecruit, color = falltreatment)) + geom_boxplot() +
  facet_wrap(~bdensity)

ggplot(together, aes(x= stemsRecruit, y = stemsIn, color = falltreatment)) + geom_point() +
  facet_wrap(~background,  scales = "free") + 
  geom_abline(intercept = 0, slope = 1)

# remove unnecessary
rm(allodat, back.dat, back.dat0, comp.dat, phyto.dat, seeds.out, stems.in, seeds.added, tomerge, recruit.dat, together)

  