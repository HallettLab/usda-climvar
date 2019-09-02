# clean and compile background competitor biomass and plot density
# authors: LMH, CTW

# script purpose:
# read in background anpp and stems datasets, with background metadata and shelter key
# clean ANPP: calculate individual weight
# clean background density: project stem counts to plot level (using estimate actual area seeded)
# combine both datasets and project total plot background ANPP based on individual weight and projected plot stem counts
# join shelter treatment dataset
# write combined dataset to Competition_CleanedData folder

# notes:
# sometimes stems subsampled but flowers counted in whole plot.. need to discuss how to treat that if use flower data
# avena sampled in 10x10cm when subsampled (because stems thicker), all else in 5x5cm


# -- SETUP ----
rm(list=ls()) # clean environment
require(readxl)
library(tidyverse)
options(stringsAsFactors = F)
theme_set(theme_bw())
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
  geom_jitter(alpha = 0.6, width = 0.1, height = 0) +
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
# when CTW and Nicolai were seeding together, left about a 5-cm margin all on sides
# this would be..
area_seeded <- (50*50) - (2*(50*10)+ 2*(30*10))
area_seeded
# CTW recorded percent cover when clipping AV and VUMY to get a sense of coverage
# look at coverage
subset(background_stems, grepl("Ave|Bro|Vul", background)) %>%
  mutate(backgrounddensity = gsub(".*_", "", background),
         backgroundspp = gsub("_.*", "", background)) %>%
  ggplot(aes(backgrounddensity, pct_cover_wp)) +
  geom_point() +
  facet_wrap(~backgroundspp)

# converted to area cm2 (comparing to 900cm^2)
subset(background_stems, grepl("Ave|Bro|Vul", background)) %>%
  mutate(backgrounddensity = gsub(".*_", "", background),
         backgroundspp = gsub("_.*", "", background)) %>%
  ggplot(aes(backgrounddensity, (pct_cover_wp/100)*2500)) +
  geom_boxplot() +
  geom_hline(yintercept = area_seeded, col = "red") +
  geom_jitter(pch = 1, width = 0.2, height = 0) +
  labs(x = "Background density",
       y = "Plot area covered (cm^2)") +
  ggtitle("Background species area covered vs. estimated area seeded (red line)") +
  facet_wrap(~backgroundspp)
# > 2019-05-01: sent LMH an email with this figure asking what she thinks we should do..
# > will use 900cm^2 seeded for now so can press on with code

# avena sampled in 10x10cm when subsampled, all else in 5x5cm 


# -- QUALITY CHECK DATASETS/TROUBLESHOOT NAs ----- 
# decisions made for prepping data (for now, as of 2019-05-01):
# area seeded = 900cm^2 (10cm buffers on all sides)
# scale flowers up directly
# exclude late-season data for forbs (since phytometer data only considers april counts)

# look at repeat samples once more..
background_bio %>%
  mutate(ID = paste(plot, background)) %>%
  subset(ID %in% ID[sample2 == 1]) %>% 
  View() # only ESCA, TRHI, LACA, none had flowers in April, but most in May 

# check for NAs in stem and ANPP counts
# anpp dataset
summary(background_bio) # 3 NAs in background bio
background_bio[is.na(background_bio$stems),c("plot", "background", "stems", "dry_wgt_g", "disturbed", "QA_flag", "QA_notes")]
background_stems[background_stems$plot %in% c(3,10,14) & 
                   grepl("Esc.*hi|Las.*lo|Bro.*lo", background_stems$background), ]
# stem count dataset
summary(background_stems)
background_stems[is.na(background_stems$stems),] # only the same bromus_low from above

# infill missing brho_low stem count in anpp dataset
# one bromus stem count is missing, but percent cover and ANPP area present.. could be infilled by estimating stems ~ ANPP + pct_cover for BRHO
BRHO_bdat <- subset(background_stems, grepl("Bro", background)) %>%
  mutate(type = ifelse(grepl("lo", background), "BRHO_LO", "BRHO_HI")) %>%
  left_join(background_bio, by = c("plot", "comp_plot", "type" = "background")) %>%
  left_join(shelter.key)
# simple models
summary(lm(stems.y ~ dry_wgt_g + pct_cover_wp, data = BRHO_bdat))
summary(lm(stems.y ~ pct_cover_wp, data = BRHO_bdat)) # stronger with just pct_cover alone
# 5.43765+(0.26468*.15) = 5.477352.. in the biomass data, 3 BRHO entries have weights of 0.41g (1 is the missing BRHO, and the other two have 6 stems)
# one of those 0.41g entries was a fall dry with 10% cover and low density, the missing brho is consistent dry with 15% cover and is low density..
# maybe just infill 6 stems as well for the missing BRHO? so we don't lost data point, an explain we infilled through best estimation available

# check corr between brho hi and brho lo
brhohi <- subset(BRHO_bdat, type == "BRHO_HI") %>% dplyr::select(plot, stems.x, stems_wp) %>% rename(hi_stems = stems.x, hi_wp = stems_wp)
brholo <- subset(BRHO_bdat, type == "BRHO_LO") %>% dplyr::select(plot, stems.x, stems_wp) %>% rename(lo_stems = stems.x, lo_wp = stems_wp)
test <- merge(brhohi, brholo) %>%
  # keep only plots where stems counted at same scale
  subset(lo_wp == hi_wp)
summary(lm(test$lo_stems ~ test$hi_stems))
plot(test$lo_stems ~ test$hi_stems, main = "Plot-level background BRHO stems in HI\nvs. BRHO stems in LO") #hm.. will keep at 6

# infill p14 bromus_low stem count and anpp stem count (are same for BRHO and VUMY) with 6
background_bio$stems[is.na(background_bio$stems) & background_bio$plot == 14] <- 6
background_stems$stems[is.na(background_stems$stems)] <- 6
# fill in area clipped/counted too (so script below projects plot density [won't if there are NAs])
background_stems$stems_area_cm2[is.na(background_stems$stems_area_cm2) & grepl("Br.*lo",background_stems$background)] <- 25
background_stems$stems_wp[is.na(background_stems$stems_wp) & grepl("Br.*lo",background_stems$background)] <- 0


# -- PREP DATASETS ----
# Clean up biomass data
background_bio2 <- background_bio %>%
  # remove 2nd sample
  subset(sample2 != 1) %>%
  rename(bflowers_anpp = flowers,
         bstems_anpp = stems,
         disturbed_banpp = disturbed) %>%
  mutate(b.ind.wgt.g = ifelse(bstems_anpp == 0, 0, round(dry_wgt_g/bstems_anpp,4)), # round to 4 places since box scale only weighed out to 4
         backgroundspp = gsub("_.*", "", background),
         backgrounddensity = ifelse(grepl("LO", background), "low", "high")) %>%
  dplyr::select(plot, backgroundspp, backgrounddensity, 
              bstems_anpp, dry_wgt_g, b.ind.wgt.g, bflowers_anpp, area_clipped_cm2, disturbed_banpp)


# Clean up stem data
background_stems2 <- background_stems %>%
  mutate(stems_scale = ifelse(stems_wp == 0, area_seeded/stems_area_cm2, 1),
         flor_scale = ifelse(flowers_wp == 0, area_seeded/flowers_area_cm2, 1),
         insitu_plot_bdensity = stems * stems_scale,
         insitu_plot_bflowers = flowers * flor_scale,
         backgroundspp = gsub("_.*", "", background),
         backgrounddensity = gsub(".*_", "", background)) %>%
  rename(insitu_bdisturbed = disturbed) %>%
  dplyr::select(plot, backgroundspp, backgrounddensity, insitu_plot_bdensity, insitu_plot_bflowers, insitu_bdisturbed)

# recode backgroundspp with code4
code4 <- c("AVFA", "BRHO", "ESCA", "LACA", "TRHI", "VUMY")
genera <- sort(unique(background_stems2$backgroundspp))
for(i in 1:length(genera)){
  background_stems2$backgroundspp[background_stems2$backgroundspp == genera[i]] <- code4[i]
}
# check
print(data.frame(orig = background_stems$background, newspp = background_stems2$backgroundspp, newdens =background_stems2$backgrounddensity))
# looks good


# -- JOIN STEMS AND BIOMASS -----
# join individual anpp and projected plot stem density
tog <- left_join(background_bio2, background_stems2) %>%
  # project plot ANPP
  # for VUMY and BRHO, whether multiply ANPP * area_seeded/area_clipped or multiply individual dry weight * project plot density, you get the same number
  # > bc area clipped was area stem counted (all other background spp, we only clipped 3 individuals for ANPP)
  mutate(plot_banpp = round(b.ind.wgt.g * insitu_plot_bdensity,4)) %>%
  dplyr::select(plot, backgroundspp, backgrounddensity, b.ind.wgt.g, plot_banpp, disturbed_banpp, insitu_plot_bdensity, insitu_plot_bflowers, insitu_bdisturbed)
  
# do the projected data seem reasonable?
summary(tog)
ggplot(tog, aes(insitu_plot_bdensity, plot_banpp)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm") +
  facet_wrap(backgroundspp ~ backgrounddensity, scales = "free")
# it seems reasonable based on what CTW saw in the field:
# low density plots typically had heartier background plants (i.e. higher ANPP plants) bc they weren't in competition with themselves
# for example, the BRHO high density plot having a negative trend with higher stem density makes sense in contrast with BRHO low

# look at flowers
ggplot(subset(tog, backgroundspp %in% c("TRHI", "LACA", "ESCA")), aes(insitu_plot_bdensity, insitu_plot_bflowers)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm") +
  facet_wrap(backgroundspp ~ backgrounddensity, scales = "free", nrow = 3)
# may need to revisit incorporating may flower data for forbs.. but maybe not bc only have seed data for LACA anyway
# but could consider incorporating flowers counted in ANPP weighing
# plus also dealing with broken LACA stems..
# > 2019-05-02: go with this for now..



# -- JOIN TREATMENT AND SEEDING DATA -----
# prep seeding key for joining
bseeding <- subset(seed.key, grepl("comp", experiment, ignore.case = T)) %>%
  dplyr::select(species:treatment, seeds_per_plot) %>%
  rename(seedsAdded = seeds_per_plot)
  
# join treatment data and seeding data
background_clean <- left_join(shelter.key, tog) %>%
  mutate(falltreatment = ifelse(grepl("fallDry|consistentDry", treatment), "dry", "wet")) %>% 
  #join seeding key data as final QA check
  left_join(bseeding, by = c("backgroundspp" = "code4", "backgrounddensity" = "treatment")) %>%
  # calculate percent seeded that survived
  mutate(perPersist = round((insitu_plot_bdensity/seedsAdded)*100, 2),
         # flag any projected densities that are unrealistic
         bdensity_flag = ifelse(insitu_plot_bdensity > seedsAdded, 1, NA))

summary(background_clean)
summary(background_clean$bdensity_flag == 1) #10 observations exceed realistic values based on numbers seeded
# what are these flagged density observations?
ggplot(subset(background_clean, perPersist > 100), aes(backgroundspp, perPersist, col = backgrounddensity)) +
  geom_jitter(width = 0.25, height = 0, alpha = 0.8, size = 2) +
  ggtitle("Percent persistence in background seeded\nthat exceeds individuals seeded (i.e. unrealistic)")
# mostly AV and BRHO (not counted by CTW), 1 TRHI.. all in low density plots..

# perhaps using percent plot covered could help with better density estimates? not sure what to do about TRHI

# what does everything else look like?
ggplot(data = subset(background_clean, perPersist < 100), aes(backgrounddensity, perPersist)) +
  #geom_point(data = subset(background_clean, perPersist > 100), aes(backgrounddensity, perPersist), col = "red") +
  geom_boxplot() +
  geom_jitter(width = 0.25, height = 0, pch = 1) +
  geom_text(data = unique(background_clean[c("backgrounddensity", "backgroundspp", "seedsAdded")]), 
            aes(x = backgrounddensity, y = 80, label = paste0("+",unique(seedsAdded))),nudge_x = -0.2, col = "dodgerblue3") +
  labs(x = "Background density",
       y = "Percent persisted (plot density/seeds added)",
       title = "Background competitor persistence (# seeds added in blue text)",
       subtitle = paste("1 unrealistic AVFA low persistence value:",background_clean$perPersist[background_clean$perPersist>100], "(not shown)")) +
  facet_wrap(~backgroundspp)
# LACA and VUMY seem too low on survivorship (esp LACA.. but those seeds are tiny)


# -- FINISHING -----
# write out to competition cleaned data folder
write.csv(background_clean, "~/Dropbox/ClimVar/Competition/Data/Competition_CleanedData/Competition_background_clean.csv", row.names = F)










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