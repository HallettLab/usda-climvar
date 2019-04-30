# clean phytometer biomass dataset
# authors: LMH, CTW
# created: Jan 2019

# script purpose:
# calculate average individual phytometer weights per species per competition plot
# > necesitates making data QA decisions..
# (TO DO [later, extra]: create addition cleaned up dataset for just AVFA that includes culm and pct greenness data)

# script steps:
## read in phytometer biomass and shelter key
## check phytometer biomass sensitivity to position and, for forbs, later season clip date (i.e. sample2 == 1)
## calculate individual plant wgts from (sample drymass / sample stem count)
## join with shelter key and write out to competition cleaned data folder


# -- SETUP ----
rm(list=ls()) # clean environment
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
options(stringsAsFactors = F)
theme_set(theme_bw())
na_vals <- c(" ", "", NA, "NA")

# set path to main level competition data folder
datpath <- "~/Dropbox/ClimVar/Competition/Data/"

# read in data
# phytometer metadata (for variable definition reference)
phyto.meta <- read_excel(paste0(datpath, "Competition_EnteredData/Competition_phytometers_spring2017.xlsx"), 
                         sheet="metadata", na = na_vals, trim_ws = T, skip = 20)
# phytometer biomass
phyto.dat <- read_excel(paste0(datpath, "Competition_EnteredData/Competition_phytometers_spring2017.xlsx"), 
                        sheet="phyto_anpp", na = na_vals, trim_ws = T)
# shelter key
shelter.key <- read.csv(paste0(datpath,"Shelter_key.csv"), na.strings = na_vals, strip.white = T)


# -- SENSITIVITY CHECKS ON RAW DATA -----
## PHYTOMETER POSITION ##
# noticeable differences in wgt by sp per competitor by position?
# don't expect any differences but just in case...
for(i in sort(unique(phyto.dat$phytometer))){
  p <- subset(phyto.dat, phytometer == i & sample2 == 0 & is.na(disturbed)) %>%
    ggplot(aes(as.factor(position), dry_wgt_g)) +
    geom_boxplot(alpha = 0.6) +
    ggtitle(paste(i, "biomass sensitivity to phytometer position")) +
    facet_wrap(~background)
  print(p)
}
# plot all together
subset(phyto.dat, sample2 == 0 & is.na(disturbed)) %>%
  ggplot(aes(as.factor(position), dry_wgt_g, col = phytometer)) +
  stat_summary(fun.y = mean, geom = "point", alpha = 0.6, pch =1) +
  ggtitle(paste("Phytometer position sensitivity: mean biomass per species per position")) +
  facet_wrap(~background)

# plot phytometer weights against all background competitors and densities in same panel 
subset(phyto.dat, sample2 == 0 & is.na(disturbed)) %>%
  ggplot(aes(as.factor(position), dry_wgt_g, group = background)) +
  stat_summary(fun.y = mean, geom = "point", alpha = 0.6) +
  ggtitle(paste("Phytometer position sensitivity: mean biomass per species per position")) +
  facet_wrap(~phytometer, scales = "free_y")

# > no visually noticeable influences in position. good!


## LATE SEASON CLIP FOR FORBS ##
# e.g. did we clip ESCA or TRHI too early? how to adjust if so?
unique(phyto.dat$phytometer[phyto.dat$sample2==1]) #ESCA, TRHI and VUMY (ignore VUMY, all way clipped in May)
data.frame(phyto.dat[phyto.dat$phytometer == "VUMY" & phyto.dat$sample2 == 1,])
# looked at data, definitely ignore this row, rep 2 of VUMY at the same position in same plot, and significantly more stems.. anomaly

#plot ESCA and TRHI 
subset(phyto.dat, phytometer %in% c("ESCA", "TRHI")) %>%
  mutate(clip_date = ifelse(!grepl("pril", clip_date), "May", "April")) %>%
  ggplot(aes(clip_date, dry_wgt_g)) +
  geom_boxplot() +
  #geom_point(alpha = 0.5) +
  facet_wrap(phytometer ~ background, scales = "free_y")
# esca seems to have developed more in may (i.e. heavier) *IF* it competed (e.g. april-may difference not the same in AVFA_HI vs AVFA_LO)
# control may esca also heavier .. ESCA is a later season forb tho
# trhi more consistently pooped out by may, except for in low density inferior grasses (VUMY ANd BRHO... but VUMY AND BRHO are also earlier season grasses than AVFA, and TRHI can hold on longer)

# > best we can do is acknowledge esca could have been clipped in may (since mid-later season forb) if results prove funky


## LOOK FOR OUTLIERS IN BIOMASS OR STEM COUNTS ##
#biomass
ggplot(subset(phyto.dat, sample2 == 0), aes(background, dry_wgt_g)) +
  geom_point(alpha = 0.5) +
  theme(axis.text.x = element_text(angle = 90)) +
  facet_wrap(~ phytometer, scales = "free_y")
# BRHO (laca_hi), ESCA (vumy_hi), and LACA (trhi_lo) all have 1 outlier..
# check ranges
with(subset(phyto.dat, sample2 == 0), lapply(split(dry_wgt_g, phytometer), function(x) range(x, na.rm=T)))
# range for values greater than 0
with(subset(phyto.dat, sample2 == 0 & dry_wgt_g>0), lapply(split(dry_wgt_g, phytometer), function(x) range(x, na.rm=T)))

#stem counts
ggplot(subset(phyto.dat, sample2 == 0), aes(background, stems)) +
  geom_point(alpha = 0.5) +
  theme(axis.text.x = element_text(angle = 90)) +
  facet_wrap(~ phytometer, scales = "free_y")
# BRHO (laca_hi), ESCA (vumy_hi), and LACA (trhi_lo) all have 1 outlier..
# check ranges
with(subset(phyto.dat, sample2 == 0), lapply(split(stems, phytometer), function(x) range(x, na.rm=T)))
# range for values greater than 0
with(subset(phyto.dat, sample2 == 0 & stems>0), lapply(split(stems, phytometer), function(x) range(x, na.rm=T)))
# not as bad
# check boxplots
ggplot(subset(phyto.dat, sample2 == 0), aes(background, stems)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90)) +
  facet_wrap(~ phytometer, scales = "free_y") #meh.. one AVFA outlier in BRHO_LO.. AV competitive.. metadata notes up to roughly 10 AV seeded, so outlier number is possible


# > general conclusions about data quality:
## > exclude sample2
## > no worries about clip date or phyto position
## > check how outliers in ESCA, TRHI, and, LACA biomass play out in individual weights


# -- CLEAN UP AND AGGREGATE DATA -----
# Get individual weights
phyto.dat2 <- as.data.frame(subset(phyto.dat, sample2 == 0)) # excude 2nd rep samples
phyto.dat2$ind.weight.g <- with(phyto.dat2, round(dry_wgt_g/stems,6)) # round to 6 decimal places (altho our analytic scale could only measure out to 4..)
# check NaNs generated (is it all from 0 biomass wgt? [i.e. nothing grew])
phyto.dat2$dry_wgt_g[is.nan(phyto.dat2$ind.weight.g)] # yes
# change NaNs to 0
phyto.dat2$ind.weight.g[is.nan(phyto.dat2$ind.weight.g)] <- 0
# split background competition and background seeding density
phyto.dat2$backgroundspp <- substr(phyto.dat2$background,1,4)
phyto.dat2$backgrounddensity <- gsub("[A-Z]{4}_", "", phyto.dat2$background)
# clean up background control names
phyto.dat2$backgroundspp[grepl("Co", phyto.dat2$background)] <- "Control"
phyto.dat2$backgrounddensity[grepl("Co", phyto.dat2$background)] <- NA
# change NA in disturbed to 0 (so true binary var)
phyto.dat2$disturbed[is.na(phyto.dat2$disturbed)] <- 0

# select final cols for clean dataset
phyto.dat3 <- phyto.dat2 %>%
  dplyr::select(plot, backgroundspp, backgrounddensity, phytometer, dry_wgt_g, stems, ind.weight.g, disturbed)
# rename cols to lmh's preferred names
colnames(phyto.dat3)[colnames(phyto.dat3) %in% c("phytometer", "dry_wgt_g", "stems")] <- c("phyto", "drywgt.g", "no.plants")


# -- FINISHING ----
# join with shelter key and clean up
phyto.bmass <-left_join(phyto.dat3, shelter.key, by = "plot") %>%
  mutate(backgrounddensity = recode(backgrounddensity, LO = "low", HI = "high"),
         falltreatment = ifelse(treatment %in% c("fallDry", "consistentDry"), "dry", "wet")) #springDry received ambient rainfall in fall
# note: LMH recoded species codes, but CTW keeping as they are in case want to pair with species descriptive info or traits
  # mutate(backgroundspp = recode(backgroundspp, AVFA = "Avena", BRHO = "Bromus", LACA = "Lasthenia",
  #                               ESCA = "Eschscholzia", TRHI = "Trifolium", VUMY = "Vulpia"),
  #        phyto = recode(phyto, AVFA = "Avena", BRHO = "Bromus", LACA = "Lasthenia",
  #                       ESCA = "Eschscholzia", TRHI = "Trifolium", VUMY = "Vulpia"),


# write out to dropbox competition cleaned data folder
#write.csv(phyto.bmass, paste0(datpath, "Competition_CleanedData/ClimVar_Comp_phytometer-biomass-2.csv"), row.names = F)
# give more informative name
write.csv(phyto.bmass, paste0(datpath, "Competition_CleanedData/ClimVar_Comp_phytometer_biomass_clean.csv"), row.names = F)

