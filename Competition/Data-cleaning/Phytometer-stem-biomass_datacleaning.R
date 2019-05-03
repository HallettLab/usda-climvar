# clean and compile phytometer biomass and in situ stem count
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
## join stems counted in field and shelter key data
## write out to competition cleaned data folder

# note:
# should probably keep notes from anpp and stems datasets so can make decisions pre-analysis about how to treat browsed samples and confirm missing samples 
# but CTW too lazy to write that in right now, just want to get all scripts working together first



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
# phytometer in situ stem count
phyto.stems <- read_excel(paste0(datpath, "Competition_EnteredData/Competition_phytometers_spring2017.xlsx"), 
                        sheet="phyto_stems", na = na_vals, trim_ws = T)
# shelter key
shelter.key <- read.csv(paste0(datpath,"Shelter_key.csv"), na.strings = na_vals, strip.white = T)


# -- SENSITIVITY CHECKS ON RAW DATA -----
## PHYTOMETER POSITION ##
# noticeable differences in wgt by sp per competitor by position?
# don't expect any differences but just in case...
for(i in sort(unique(phyto.dat$phytometer))){
  p <- subset(phyto.dat, phytometer == i & sample2 == 0 & is.na(disturbed)) %>%
    ggplot(aes(as.factor(position), dry_wgt_g)) +
    geom_point(alpha = 0.6) +
    ggtitle(paste(i, "biomass sensitivity to phytometer position")) +
    facet_wrap(~background, scales = "free_y")
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

# check stem count by position
for(i in sort(unique(phyto.stems$phytometer))){
  p <- subset(phyto.stems, phytometer == i & is.na(disturbed)) %>%
    ggplot(aes(as.factor(position), stems)) +
    geom_point(alpha = 0.6) +
    ggtitle(paste(i, "stem count sensitivity to phytometer position")) +
    facet_wrap(~background, scales = "free_y")
  print(p)
}
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
phyto.dat2$p.ind.wgt.g <- with(phyto.dat2, round(dry_wgt_g/stems,6)) # round to 6 decimal places (altho our analytic scale could only measure out to 4..)
# check NaNs generated (is it all from 0 biomass wgt? [i.e. nothing grew])
phyto.dat2$dry_wgt_g[is.nan(phyto.dat2$p.ind.wgt.g)] # yes
# change NaNs to 0
phyto.dat2$p.ind.wgt.g[is.nan(phyto.dat2$p.ind.wgt.g)] <- 0
# split background competition and background seeding density
phyto.dat2$backgroundspp <- substr(phyto.dat2$background,1,4)
phyto.dat2$backgrounddensity <- gsub("[A-Z]{4}_", "", phyto.dat2$background)
# clean up background control names
phyto.dat2$backgroundspp[grepl("Co", phyto.dat2$background)] <- "Control"
phyto.dat2$backgrounddensity[grepl("Co", phyto.dat2$background)] <- NA
# change NA in disturbed to 0 (so true binary var)
phyto.dat2$disturbed[is.na(phyto.dat2$disturbed)] <- 0

# prep stem dataset for merging with biomass
phyto.stems2 <- dplyr::select(phyto.stems, plot, background:stems, disturbed) %>%
  rename(insitu_pstems = stems,
         insitu_pdisturbed = disturbed) %>%
  mutate(insitu_pdisturbed = ifelse(is.na(insitu_pdisturbed), 0, insitu_pdisturbed),
         # split background species from density treatment
         backgroundspp = gsub("_.*", "", background),
         backgrounddensity = ifelse(grepl("lo", background), "LO", 
                                    ifelse(grepl("hi", background), "HI", NA)),
         # change genera to 4-letter code
         backgroundspp = recode(backgroundspp, Avena = "AVFA", Bromus = "BRHO", Lasthenia = "LACA",
                                Eschscholzia = "ESCA", Vulpia = "VUMY", Trifolium = "TRHI"),
         backgroundspp = ifelse(grepl("Cont", backgroundspp), "Control", backgroundspp),
         # change genera to 4-letter code in phytometer col
         phytometer = recode(phytometer, Avena = "AVFA", Bromus = "BRHO", Lasthenia = "LACA",
                             Eschscholzia = "ESCA", Vulpia = "VUMY", Trifolium = "TRHI"))


# combine biomass and in situ stem count datasets, select final cols for clean dataset
phyto.dat3 <- phyto.dat2 %>%
  dplyr::select(plot, backgroundspp, backgrounddensity, phytometer, dry_wgt_g, stems, p.ind.wgt.g, disturbed) %>%
  # rename stem column in ANPP dataset to not confused with in situ stem column
  rename(pANPP_stems = stems,
         pANPP_disturbed = disturbed) %>%
  # join stem counts 
  left_join(phyto.stems2[!colnames(phyto.stems2) %in% c("background", "position")]) %>%
  mutate(p_totwgt = round(p.ind.wgt.g*insitu_pstems,4))

# logic check: do any ANPP stem counts exceed field stem counts? (should not)
summary(phyto.dat3$pANPP_stems > phyto.dat3$insitu_pstems)
# what are the NAs?
View(phyto.dat3[is.na(phyto.dat3$pANPP_stems > phyto.dat3$insitu_pstems),])
# these are true NAs for ANPP:
# 1) VUMY was planted in place of BRHO
# 2) ESCA wasn't clipped, but stems counted
# 3) TRHI was missing (didn't see it in plot photo, no sample for it in lab.. but since stem count exists, must have not been clipped)

# rename cols to lmh's preferred names -- not yet
# colnames(phyto.dat3)[colnames(phyto.dat3) %in% c("phytometer", "dry_wgt_g", "ANPP_stems", "field_stems")] <- c("phyto", "drywgt.g", "ANPP.no.plants", "field.no.stems")


# -- FINISHING ----
# join with shelter key and clean up
phyto.bmass <-left_join(phyto.dat3, shelter.key, by = "plot") %>%
  mutate(backgrounddensity = recode(backgrounddensity, LO = "low", HI = "high"),
         falltreatment = ifelse(treatment %in% c("fallDry", "consistentDry"), "dry", "wet")) #springDry received ambient rainfall in fall

# write out to dropbox competition cleaned data folder
#write.csv(phyto.bmass, paste0(datpath, "Competition_CleanedData/ClimVar_Comp_phytometer-biomass-2.csv"), row.names = F)
# give more informative name
write.csv(phyto.bmass, paste0(datpath, "Competition_CleanedData/Competition_phytometers_clean.csv"), row.names = F)

