# clean and compile background competitor biomass and plot density
# authors: LMH, CTW

# script purpose:


# notes:
# sometimes stems subsampled but flowers counted in whole plot.. need to discuss how to treat that if use flower data
# avena sampled in 10x10cm when subsampled (because stems thicker), all else in 5x5cm


# -- SETUP ----
rm(list=ls()) # clean environment
require(readxl)
library(tidyverse)
options(stringsAsFactors = F)
na_vals <- c(" ", "", NA, "NA")

# set competition data pathway
datpath <- "~/Dropbox/ClimVar/Competition/Data/"
# pathway to entered daata
comppath <- paste0(datpath, "Competition_EnteredData/Competition_background_spring2017.xlsx")

# Read in data
# background competitor metatadata
background_meta <- read_excel(comppath, sheet = "metadata", na = na_vals, skip = 21)
# background competitor stem counts
background_stems <-read_excel(comppath, sheet=2, na = na_vals)
# background competitor biomass
background_bio <-read_excel(comppath, sheet=3, na = na_vals)
# treatment
shelter.key <- read.csv("~/Dropbox/ClimVar/Competition/Data/Shelter_key.csv", na.strings = na_vals)
# seeding lookup table
seed.key <- read.csv("~/Dropbox/ClimVar/Competition/Data/Competition_EnteredData/ClimVar_CompGerm_seedingKey.csv", na.strings = na_vals)


# -- SENSITIVITY CHECKS -----
# look at dates collected by species
lapply(split(background_bio$clip_date, background_bio$background), unique)
# all forbs have late season clips in addition to main april clip, look at difference in wgts by month

subset(background_bio, grepl("TRH|ESC|LAC", background) & dry_wgt_g >0 & !is.na(dry_wgt_g)) %>%
  mutate(clip_month = ifelse(grepl("Apr", clip_date), "April", "May"),
         background = gsub("_.*", "", background)) %>%
  ggplot(aes(stems, dry_wgt_g, col = clip_month)) +
  geom_jitter(alpha = 0.6, width = 0.1) +
  facet_wrap(~background)

subset(background_bio, grepl("TRH|ESC|LAC", background) & dry_wgt_g >0 & !is.na(dry_wgt_g)) %>%
  mutate(clip_month = ifelse(grepl("Apr", clip_date), "April", "May"),
         background = gsub("_.*", "", background)) %>%
  ggplot(aes(clip_month, dry_wgt_g, group = clip_month)) +
  geom_boxplot() +
  geom_jitter(pch=1, width = 0.1) +
  facet_wrap(~background)

subset(background_stems, grepl("TRI|ESC|LAs", background, ignore.case = T) & stems >0 & !is.na(stems)) %>%
  mutate(sample_month = ifelse(grepl("Apr", sample_date), "April", "May"),
         background = gsub("_.*", "", background)) %>%
  ggplot(aes(sample_month, stems, group = sample_month)) +
  geom_boxplot() +
  geom_jitter(pch=1, width = 0.1) +
  facet_wrap(~background) # stems only sampled in situ in april

# > not much of a difference for anpp stems late season, but there is a difference for ESCA and TRI ANPP in may..
# > cleaned phytometer dataset only keeps april samples.. so i guess keep only april here too?
# > can acknowledge TRHI and ESCA developed more in May, maybe include a supplemental figure if it makes sense


## FLOWERS OBSERVED BY SAMPLE DATE ##
subset(background_bio, grepl("TRH|ESC|LAC", background) & !is.na(flowers)) %>%
  mutate(clip_month = ifelse(grepl("Apr", clip_date), "April", "May"),
         background = gsub("_.*", "", background)) %>%
  ggplot(aes(clip_month, flowers, group = clip_month)) +
  geom_boxplot() +
  geom_jitter(pch=1, width = 0.1) +
  facet_wrap(~background) 
#hm.. seeds not counted for ESCA or TRHI (no allometric relationship, so i guess it doesn't matter)
# no difference for LACA (only 1 may data point for LACA anyway)

# how many LACA obs have broken stalks?
summary(is.na(background_stems$broken_LACA_stalks[grepl("LA", background_stems$background, ignore.case = T)]))
length(grep("LA", background_stems$background, ignore.case = T)) # 7 of 32 LACA have broken stems..

# > come back to this.. if use flower data for anything, talk to LMH on how to treat
# > also, after looking at data for ESCA and TRI, just keep flower data for LACA since sporadic for ESCA and TRHI (also flower data on those more from May... unless we want to project flowers based on stems?)

# plot flowers by stem count
subset(background_bio,!is.na(flowers)) %>%
  mutate(clip_month = ifelse(grepl("Apr", clip_date), "April", "May"),
         background = gsub("_.*", "", background)) %>%
  ggplot(aes(stems, flowers, col = clip_month)) +
  geom_jitter(alpha = 0.6, width = 0.2, height = 0) +
  facet_wrap(~background)

# plot flowers by whether counted in whole plot or subsample
subset(background_stems,!is.na(flowers_wp)) %>%
  mutate(background = gsub("_.*", "", background)) %>%
  ggplot(aes(as.factor(flowers_wp), flowers)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.6, width = 0.2, height = 0, col = "grey50") +
  facet_wrap(~background) # if you project subsampled flowers out to the plot level, that's going to be quite a bit more..
# compare projected flowers
subset(background_stems,!is.na(flowers_wp)) %>%
  mutate(background = gsub("_.*", "", background),
         project_flowers = ifelse(flowers_wp == 0, flowers * 75, flowers)) %>% # multiple 3/4s of plot, being conservative
  ggplot(aes(as.factor(flowers_wp), project_flowers)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.6, width = 0.2, height = 0, col = "grey50") +
  facet_wrap(~background) # that seems pretty unrealistic..

# try scaling up to half plot, split by low and high density to see if more reasonable?
subset(background_stems, grepl("La", background) & !is.na(flowers_wp)) %>%
  mutate(project_flowers = ifelse(flowers_wp == 0, flowers * 50, flowers)) %>%
  ggplot(aes(as.factor(flowers_wp), project_flowers, col = broken_LACA_stalks)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.6, width = 0.2, height = 0) +
  facet_wrap(~background) # idk.. low density seems more realistic i guess.. high still seems like a lot, but a lot of plants seeded so maybe?
# there are broken LACA stalks in both whole plot and subsample counts, doesn't seem to be a huge factor
# > bigger consideration is whether appropriate to assume can scale up flowers directly to plot level..
# > maybe try flowers per stem and scale that?

subset(background_stems, grepl("La", background) & !is.na(stems) & stems > 0) %>%
  mutate(flowers = ifelse(is.na(flowers), 0, flowers),
         flor_per_stem = flowers/stems) %>%
  ggplot(aes(background, flor_per_stem)) +
  geom_boxplot() +
  geom_jitter(aes(col = as.factor(flowers_wp)), width = 0.2, height = 0) # looking better





# -- CHECK AREA SEEDED -----
# nicolai wrote in the notes he left a 10cm buffer on each side when seeding (after CTW left SFREC)
# when CTW and Nicolai were seeding together, left about a 5-cm marging all on sides
# this would be..
area_seeded <- (50*50) - (2*(40*10)+ 2*(20*10))
area_seeded
# CTW recorded percent cover when clipping AV and VUMY to get a sense of coverage
# look at coverage
subset(background_stems, grepl("Ave|Bro|Vul", background)) %>%
  mutate(backgrounddensity = gsub(".*_", "", background),
         backgroundspp = gsub("_.*", "", background)) %>%
  ggplot(aes(backgrounddensity, pct_cover_wp)) +
  geom_point() +
  facet_wrap(~backgroundspp)

# converted to area cm2 (comparing to 1300cm^2)
subset(background_stems, grepl("Ave|Bro|Vul", background)) %>%
  mutate(backgrounddensity = gsub(".*_", "", background),
         backgroundspp = gsub("_.*", "", background)) %>%
  ggplot(aes(backgrounddensity, (pct_cover_wp/100)*2500)) +
  geom_boxplot() +
  geom_hline(yintercept = area_seeded, col = "red") +
  geom_jitter(pch = 1, width = 0.2, height = 0) +
  ggtitle("Background species area covered vs. estimated area seeded (red line)") +
  facet_wrap(~backgroundspp)
# > 2019-05-01: sent LMH an email with this figure asking what she thinks we should do..
# > will use 1300cm^2 seeded for now so can press on with code

# avena sampled in 10x10cm when subsampled, all else in 5x5cm 


# -- PREP DATASETS ----- 
# decisions made for prepping data (for now, as of 2019-05-01):
# area seeded = 1300cm^2 (10cm buffers on all sides)
# scale flowers up directly
# exclude late-season data for forbs (since phytometer data only considers april counts)

# Clean up biomass data
background_bio2 <- background_bio

  mutate(flowers_clipped = ifelse(flowers_wp == 1, flowers, flowers*(area_seeded/25)),
         stems_clipped = ifelse(stems_wp),
         harvest = ifelse(clip_date == "April 2017", "April", "May")) %>%
  mutate(backgroundspp = gsub("_.*", "", background),
         backgrounddensity = ifelse(grepl("LO", background), "low", "high")) %>%
  dplyr::select(plot, backgroundspp, backgrounddensity, area_clipped_cm2, 
              stems_clipped, dry_wgt_g, flower_bcount, harvest)


# Note: 1 low VUMY missing, think just was no material there, VUMY has stems for all others except 1 hi, 3 low, 3 high
# Bromus has no stem counts 
# tocheck <- background_bio2 %>%
#   filter(!is.na(area_clipped_cm2))


# Clean up stem data
background_stems2 <- background_stems %>%
  mutate(stems_scale = ifelse(stems_wp == 0, area_seeded/stems_area_cm2, 1),
         flor_scale = ifelse(flowers_wp == 0, area_seeded/flowers_area_cm2, 1),
         insitu_bstems = stems * stems_scale,
         insitu_bflowers = flowers * flor_scale,
         backgroundspp = gsub("-.*", "", background),
         backgrounddensity = gsub(".*_", "", background)) %>%
  dplyr::select(plot, backgroundspp, backgrounddensity, insitu_bstems, insitu_bflowers, disturbed)
# one bromus stem count is missing, but percent cover and ANPP area present.. could be infilled by estimating stems ~ ANPP + pct_cover for BRHO

BRHO_bdat <- subset(background_stems, grepl("Bro", background)) %>%
  mutate(type = ifelse(grepl("lo", background), "BRHO_LO", "BRHO_HI")) %>%
  left_join(background_bio, by = c("plot", "comp_plot", "type" = "background")) %>%
  left_join(shelter.key)

summary(lm(stems.y ~ dry_wgt_g + pct_cover_wp, data = BRHO_bdat))
summary(lm(stems.y ~ pct_cover_wp, data = BRHO_bdat)) # stronger with just pct_cover alone
# 5.4376+(0.264*.15) = 5.4772.. in the biomass data, 3 BRHO entries have weights of 0.41g (1 is the missing BRHO, and the other two have 6 stems)
# one of those 0.41g entries was a fall dry with 10% cover and low density, the missing brho is consistent dry with 15% cover and is low density..
# maybe just infill 6 stems as well for the missing BRHO? so we don't lost data point, an explain we infilled through best estimation available
background_stems2$insitu_bstems[is.na(background_stems2$insitu_bstems)] <- 6 * (area_seeded/25)
# yield 312 stems in the comp plot; the 0.41 fall dry plot has 6 stems in a 5x5cm subsample (and 7 in 5x5cm subsample 2 [see plot notes])
# BUT whole plot only had 41 stems (Nikolai counted all) because coverage was patchy.. so 312 could be a big overestimation for BRHO LO on plot 14

# check corr between brho hi and brho lo
brhohi <- subset(BRHO_bdat, type == "BRHO_HI") %>% dplyr::select(plot, stems.x) %>% rename(hi_stems = stems.x)
brholo <- subset(BRHO_bdat, type == "BRHO_LO") %>% dplyr::select(plot, stems.x) %>% rename(lo_stems = stems.x)
test <- merge(brhohi, brholo)
cor.test(~ test$lo_stems + test$hi_stems)
plot(test$hi_stems ~ test$lo_stems) #hm.. will keep at 6


# -- JOIN STEMS AND BIOMASS -----

tog <- left_join(background_bio2, background_stems2) %>%
   mutate(ind_weight_g = dry_wgt_g/stems_clipped) %>%
#   ## density should be multiplied by 2500 for the full plot; but because the
#   ## seeds were concentrated in the center i'm estimating that it was in fact a smaller area (1/4 that size)
#   mutate(density = stems_counted/area_counted_cm2*625) %>%
#
#   # mutate(stems_counted1 = ifelse(is.na(stems_counted), stems_clipped, stems_counted),
#   # stems_clipped1 = ifelse(is.na(stems_clipped), plants_clipped, stems_clipped),
#   # stems_clipped1 = ifelse(is.na(stems_clipped1) & area_counted_cm2 == 25,
#   #                               stems_counted, stems_clipped1)) %>%
#   #  mutate(ind_weight_g = dry_wgt_g/stems_clipped1) %>%
#   ## density should be multiplied by 2500 for the full plot; but because the
#   ## seeds were concentrated in the center i'm estimating that it was in fact a smaller area (1/4 that size)
#   ## might make this an if/else if the plot is subsampled or completely sampled...
#   #  mutate(density = stems_counted1/area_counted_cm2*625) %>%
#   mutate(tot_weight_g = ind_weight_g*density) %>%
#   mutate(flowers = ifelse(is.na(flowers_counted), flower_bcount, flowers_counted),
#          ind_flower = flowers/stems_counted) %>%
#   select(plot, backgroundspp, backgrounddensity, harvest, ind_weight_g, density, tot_weight_g, ind_flower)




# Join it all
background_clean <-left_join(shelter.key, tog) %>%
  mutate(falltreatment = "wet",
         falltreatment = ifelse(treatment == "fallDry" | treatment == "consistentDry", "dry", falltreatment))

background_clean <-left_join(background_clean, seed.key) %>%
  mutate(perPersist = density/seedsAdded)


  


# -- FINISHING -----
write.csv(background, "~/Dropbox/ClimVar/Competition/Data/Competition_CleanedData/ClimVar_Comp_background-biomass-2.csv", row.names = F) %>%
  tbl_df()



####### OLD CODE (LMH WROTE EARLIER BUT NOT USING) ###############
# used some but not all of this above (e.g. not scaling LMH used)
# put it all together
# tog <- left_join(background_bio2, background_stems2) %>%
#   mutate(stems_clipped = as.numeric(as.character(stems_clipped)),
#          stems_counted = as.numeric(as.character(stems_counted)),
#          area_clipped_cm2 = as.numeric(as.character(area_clipped_cm2))) %>%
#   mutate(ind_weight_g = dry_wgt_g/stems_clipped) %>%
#   ## density should be multiplied by 2500 for the full plot; but because the
#   ## seeds were concentrated in the center i'm estimating that it was in fact a smaller area (1/4 that size)
#   mutate(density = stems_counted/area_counted_cm2*625) %>%
#   
#   # mutate(stems_counted1 = ifelse(is.na(stems_counted), stems_clipped, stems_counted),
#   # stems_clipped1 = ifelse(is.na(stems_clipped), plants_clipped, stems_clipped),
#   # stems_clipped1 = ifelse(is.na(stems_clipped1) & area_counted_cm2 == 25, 
#   #                               stems_counted, stems_clipped1)) %>%
#   #  mutate(ind_weight_g = dry_wgt_g/stems_clipped1) %>%
#   ## density should be multiplied by 2500 for the full plot; but because the
#   ## seeds were concentrated in the center i'm estimating that it was in fact a smaller area (1/4 that size)
#   ## might make this an if/else if the plot is subsampled or completely sampled...
#   #  mutate(density = stems_counted1/area_counted_cm2*625) %>%
#   mutate(tot_weight_g = ind_weight_g*density) %>%
#   mutate(flowers = ifelse(is.na(flowers_counted), flower_bcount, flowers_counted),
#          ind_flower = flowers/stems_counted) %>%
#   select(plot, backgroundspp, backgrounddensity, harvest, ind_weight_g, density, tot_weight_g, ind_flower) 



# no need to average because not keeping late-season samples for forbs in cleaned dataset
# # Deal with the duplicate measures through averaging
# background <- background_clean %>%
#   mutate(phyto = backgroundspp) %>%
#   select(plot, treatment:backgroundspp, backgrounddensity, phyto, falltreatment, ind_weight_g, tot_weight_g, density, seedsAdded, ind_flower) %>%
#   group_by(plot, treatment, backgroundspp, backgrounddensity, phyto, falltreatment, shelterBlock) %>%
#   # summarize(ind_weight_g = mean(ind_weight_g), density = mean(density), seedsAdded = mean(seedsAdded)) %>%
#   group_by(plot, treatment, backgroundspp, backgrounddensity, phyto, falltreatment, shelterBlock) %>%
#   mutate(maxval = max(ind_weight_g),
#          ismax = ifelse(ind_weight_g == maxval, 1, 0)) %>%
#   filter(ismax == 1 | is.na(ismax)) %>%
#   summarize(ind_weight_g = mean(ind_weight_g),
#             density = mean(density),
#             seedsAdded = mean(seedsAdded),
#             tot_weight_g = mean(tot_weight_g),
#             ind_flower = mean(ind_flower))