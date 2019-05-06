#Caitlin White
#EBIO/INSTAAR
#September 2016

#Data exploration of SFREC USDA shelter project native diversity data from Yr 1 and Yr2

### ------> NOTE <------ ####
# CTW cleaned up this script May 2019 so it runs (i.e. minimal update)..
# incorporated code from other scripts to consolidate scripts
# needs more work..
# *****************************



#Load necessary libraries
library(tidyverse)
options(stringsAsFactors = F)
theme_set(theme_bw())
na_vals <- c(" ", "", NA, "NA")

datpath <- "~/Dropbox/Plant_composition_data/DATA/Native_Diversity/"

#read in data
natdat <- read.csv(paste0(datpath,"Native_Diversity_CleanedData/NativeDiversity_Cover_2015-2016_cleaned.csv"), 
                   na.strings = na_vals, strip.white = T)
str(natdat)


## -- Data prep
# parking piece of older code, fix later:
##Add columns to shelter key frame to indicate rainfall treatment as binary data, where 1 = receives rain and 0 = rain excluded##
# key <- mutate(key, spring = ifelse(treatment == "fallDry" | treatment =="controlRain", 1,0), 
#               fall = ifelse(treatment == "fallDry" | treatment == "consistentDry", 0,1)) %>%
#   mutate(spring=as.factor(spring), fall = as.factor(fall)) %>%
#   tbl_df()



# treatment as factor (per LMH)
natdat <- mutate(natdat, treatment=ordered(treatment, levels = c(consistentDry="consistentDry", springDry="springDry", fallDry="fallDry", controlRain="controlRain")))

#add lifehistory variable (exported unique species and found almost all (natives included are annual), specify exceptions otherwise assign "annual" to all)
#natdat$lifehistory

## -- LMH visuals -----
# Looks like Clarkia, Eschsholzia and Achillea are the main drivers of "natives do better in wet falls" pattern
# Might want to reconsider Mimulus given it's low recruitment
ggplot(subset(natdat, nativity=="Native"), aes(x=treatment, y=cover, fill=as.factor(yr))) + geom_boxplot() + facet_wrap(~species, scales = "free")
#ggplot(subset(dat.nat, nativity!="" & nativity=="Exotic"), aes(x=treatment, y=cover, fill=as.factor(yr))) + geom_boxplot() + facet_wrap(~species, scales = "free")


## -- Exploratory data visuals
# 1) what are general trends by year, by nativity, by func, by lifeform?
ggplot(subset(natdat, !is.na(cover) & !is.na(nativity))) + geom_boxplot(aes(nativity, cover)) + facet_grid(yr ~ treatment)
ggplot(subset(natdat, !is.na(cover) & !is.na(nativity))) + geom_boxplot(aes(fxnl_grp, cover)) + facet_grid(yr ~ treatment)
ggplot(subset(natdat, !is.na(cover) & !is.na(nativity))) + geom_boxplot(aes(Duration, cover)) + facet_grid(yr ~ treatment)
ggplot(subset(natdat, !is.na(cover) & !is.na(nativity))) + geom_boxplot(aes(nativity, cover, fill = Duration)) + facet_grid(yr ~ treatment)

ggplot(subset(natdat, nativity=="Native"), aes(plot,cover)) + geom_point(aes(color=nativity)) + 
  stat_smooth() +
  theme_bw() +
  facet_wrap(~ treatment) 

ggplot(subset(natdat, fxnl_grp=="forb"), aes(plot,cover, fill=nativity)) + 
  #geom_point() + 
  #stat_smooth() +
  stat_summary(fun.y=sum, geom="bar", alpha=0.5) +
  theme_bw() +
  facet_grid(treatment ~ year) 

ggplot(natdat, aes(plot,cover, fill=Duration)) + 
  #geom_jitter() + 
  #stat_smooth() +
  stat_summary(fun.y=sum, geom="bar", alpha=0.5) +
  theme_bw() +
  facet_grid(treatment ~ year) 


##-- Summarize data on subplot level, by func, Duration, and nativity
natdat.subplot <- natdat %>%
  subset(!is.na(cover) & !is.na(nativity) & !is.na(fxnl_grp)) %>%
  group_by(yr, plot, treatment, shelterBlock, shelter, fxnl_grp, nativity, Duration) %>%
  summarise(cover=sum(cover))

ggplot(natdat.subplot, aes(plot,cover, group=nativity, fill=nativity)) + 
  stat_summary(fun.y=sum, geom="bar", alpha=0.5) +
  #geom_jitter(aes(color=nativity), size=5, width=1) + 
  #stat_smooth() +
  theme_bw() +
  facet_grid(fxnl_grp ~ yr) 


#create dataset for funcal means for plotting horizontal lines
fxnl.sum <- natdat.subplot %>%
  group_by(yr, plot, treatment, shelter, fxnl_grp) %>%
  summarise(sum.cover=sum(cover), nobs=length(cover))

fxnl.means <- fxnl.sum %>%
  group_by(yr, treatment, shelter, fxnl_grp) %>%
  summarise(mean.cover=mean(sum.cover), se.cover= sd(sum.cover)/sqrt(length(sum.cover)), nobs=length(sum.cover))

ggplot(natdat.subplot, aes(plot,cover, group=fxnl_grp, fill=fxnl_grp)) + 
  stat_summary(fun.y=sum, geom="bar", alpha=0.5) +
  #geom_jitter(aes(color=nativity), size=5, width=1) + 
  #stat_smooth() +
  theme_bw() +
  facet_grid(treatment ~ yr) +
  geom_hline(aes(yintercept=mean.cover, color=fxnl_grp), fxnl.means)

#create dataset for native nativity means for plotting horizontal lines
nativity.sum <- natdat.subplot %>%
  group_by(yr, plot, treatment, shelter, nativity) %>%
  summarise(sum.cover=sum(cover), nobs=length(cover))

nativity.means <- nativity.sum %>%
  group_by(yr, treatment, shelter, nativity) %>%
  summarise(mean.cover=mean(sum.cover), se.cover= sd(sum.cover)/sqrt(length(sum.cover)), nobs=length(sum.cover))

ggplot(natdat.subplot, aes(plot,cover, group=nativity, fill=nativity)) + 
  stat_summary(fun.y=sum, geom="bar", alpha=0.5) +
  #geom_jitter(aes(color=nativity), size=5, width=1) + 
  #stat_smooth() +
  theme_bw() +
  facet_grid(treatment ~ yr) +
  geom_hline(aes(yintercept=mean.cover, color=nativity), nativity.means)

#create dataset for percent cover, summed by plot, nativity, and func
natdat.nozero <- natdat[natdat$cover != 0 & !is.na(natdat$cover),]

natdat.stat.fxn <- natdat.nozero %>%
  subset(!is.na(nativity) & !is.na(fxnl_grp)) %>%
  group_by(yr, plot, treatment, shelterBlock, shelter, fxnl_grp, nativity) %>%
  summarise(cover=sum(cover), nobs=length(nativity))

natdat.stat.fxn <- unite_(natdat.stat.fxn, "group", c("nativity", "fxnl_grp"), sep=" ", remove=F)
natdat.stat.fxn$group <- as.factor(natdat.stat.fxn$group)

ggplot(natdat.stat.fxn) + geom_boxplot(aes(nativity, cover, fill=group)) + 
  facet_grid(yr ~ treatment) + 
  labs(y="Cover (%)", x=NULL) +
  theme(axis.title.y=element_text(face="bold", size=20),
        axis.text.y=element_text(size=16),
        axis.text.x = element_blank(),
        strip.text=element_text(face="bold", size=16),
        legend.text=element_text(face="bold", size=16),
        legend.direction="horizontal",
        legend.position="bottom",
        legend.title=element_blank())

ggplot(natdat.stat.fxn) + geom_boxplot(aes(nativity, nobs, fill=group)) + 
  facet_grid(yr ~ treatment) + 
  labs(y="Count unique species", x=NULL) +
  scale_y_continuous(limits = c(0,10), breaks=c(2,4,6,8,10), expand=c(0,0)) +
  theme(axis.title.y=element_text(face="bold", size=20),
        axis.text.y=element_text(size=16),
        axis.text.x = element_blank(),
        strip.text=element_text(face="bold", size=16),
        legend.text=element_text(face="bold", size=16),
        legend.direction="horizontal",
        legend.position="bottom",
        legend.title=element_blank())

#what species drive patterns each year and each treatment? for natives and Exotics?
ggplot(subset(natdat, nativity == "Native" & !is.na(nativity) & !is.na(cover))) + geom_boxplot(aes(nativity, cover, fill=species)) + 
  facet_grid(yr ~ treatment) + 
  labs(y="Cover (%)", x=NULL) +
  theme_bw() + 
  theme(axis.title.y=element_text(face="bold", size=20),
        axis.text.y=element_text(size=16),
        axis.text.x = element_blank(),
        strip.text=element_text(face="bold", size=16),
        legend.text=element_text(face="bold", size=16),
        legend.direction="horizontal",
        legend.position="bottom",
        legend.title=element_blank())

ggplot(subset(natdat, nativity == "Native" & !is.na(nativity) & !is.na(cover)), aes(species, cover, fill=species)) + 
  #geom_point() + 
  stat_summary(fun.y="mean", geom="bar") +
  stat_summary(fun.y = mean,
               fun.ymin = function(x) mean(x) - sd(x), 
               fun.ymax = function(x) mean(x) + sd(x), 
               geom = "pointrange",
               color="black") +
  facet_grid(yr ~ treatment) + 
  labs(y="Cover (%)", x=NULL) +
  theme_bw() + 
  theme(axis.title.y=element_text(face="bold", size=20),
        axis.text.y=element_text(size=16),
        axis.text.x = element_blank(),
        strip.text=element_text(face="bold", size=16),
        legend.text=element_text(face="bold", size=16),
        #legend.direction="horizontal",
        #legend.position="bottom",
        legend.title=element_blank())

ggplot(subset(natdat, nativity == "Exotic" & !is.na(nativity) & !is.na(cover)), aes(species, cover, fill=species)) + 
  #geom_point() + 
  stat_summary(fun.y="mean", geom="bar") +
  stat_summary(fun.y = mean,
               fun.ymin = function(x) mean(x) - sd(x), 
               fun.ymax = function(x) mean(x) + sd(x), 
               geom = "pointrange",
               color="black") +
  facet_grid(yr ~ treatment) + 
  labs(y="Cover (%)", x=NULL) +
  theme_bw() + 
  theme(axis.title.y=element_text(face="bold", size=20),
        axis.text.y=element_text(size=16),
        axis.text.x = element_blank(),
        strip.text=element_text(face="bold", size=16),
        legend.text=element_text(face="bold", size=16),
        #legend.direction="horizontal",
        #legend.position="bottom",
        legend.title=element_blank())

ggplot(data=subset(natdat, nativity == "Exotic" & fxnl_grp == "Grass" & !is.na(cover)), aes(species, cover, fill=species)) + 
  #geom_boxplot() +
  #geom_point() + 
  stat_summary(fun.y="mean", geom="bar") +
  stat_summary(fun.y = mean,
               fun.ymin = function(x) mean(x) - sd(x), 
               fun.ymax = function(x) mean(x) + sd(x), 
               geom = "pointrange",
               color="black") +
  facet_grid(yr ~ treatment) + 
  labs(y="Cover (%)", x=NULL) +
  theme_bw() + 
  theme(axis.title.y=element_text(face="bold", size=20),
        axis.text.y=element_text(size=16),
        axis.text.x = element_blank(),
        strip.text=element_text(face="bold", size=16),
        legend.text=element_text(face="bold", size=16),
        #legend.direction="horizontal",
        #legend.position="bottom",
        legend.title=element_blank())

ggplot(data=subset(natdat, nativity == "Exotic" & fxnl_grp == "Forb"), aes(species, cover, fill=species)) + 
  #geom_boxplot() + 
  stat_summary(fun.y="mean", geom="bar") +
  stat_summary(fun.y = mean,
               fun.ymin = function(x) mean(x) - sd(x), 
               fun.ymax = function(x) mean(x) + sd(x), 
               geom = "pointrange",
               color="black") +
  facet_grid(yr ~ treatment) + 
  labs(y="Cover (%)", x=NULL) +
  theme_bw() + 
  theme(axis.title.y=element_text(face="bold", size=20),
        axis.text.y=element_text(size=16),
        axis.text.x = element_blank(),
        strip.text=element_text(face="bold", size=16),
        legend.text=element_text(face="bold", size=16),
        #legend.direction="horizontal",
        #legend.position="bottom",
        legend.title=element_blank())

#Spp driving patterns:
#Natives: CLAM, ESCA (Minor roles: COHE, GICA, NEME, JUBU, LUMI)
#Exotic grasses: AVBA, LOMU, TACA
#Exotic forbs: ERBO, TRSU, CESO (Minor roles: TRHI, CEGL, SIGA, VIA, HYGL, HYRA)

##plot again focusing on spp driving patterns
natdat.drivers <- subset(natdat, epithet == "barbata" | epithet == "multiflorum" | epithet == "caput-medusae" |
                           epithet == "botrys" | epithet == "subterraneum" | epithet == "solstitialis" | epithet == "sativa" |
                           epithet == "amoena" | epithet == "californica" | epithet == "heterophylla" | epithet == "bufonius")
natdat.drivers <- unite_(natdat.drivers, "group", c("nativity", "fxnl_grp"), sep=" ", remove=F)        
natdat.drivers$species <- factor(natdat.drivers$species, levels=unique(natdat.drivers$species[order(natdat.drivers$group)], ordered=T))

ggplot(natdat.drivers, aes(species, cover, fill=species)) + 
  #geom_boxplot() + 
  scale_fill_manual(values = c("indianred1", "hotpink1", "lightcoral", "lightsalmon",
                               "darkgreen", "darkseagreen", "forestgreen", "lightgreen",
                               "royalblue1", "skyblue1", "steelblue4")) +
  stat_summary(fun.y="mean", geom="bar") +
  stat_summary(fun.y = mean,
               fun.ymin = function(x) mean(x) - (sd(x)/sqrt(length(x))), 
               fun.ymax = function(x) mean(x) + (sd(x)/sqrt(length(x))), 
               geom = "pointrange",
               color="black") +
  facet_grid(yr ~ treatment) + 
  labs(y="Mean cover (%)", x=NULL) +
  theme_bw() + 
  theme(axis.title.y=element_text(face="bold", size=20),
        axis.text.y=element_text(size=16),
        axis.text.x = element_blank(),
        strip.text=element_text(face="bold", size=16),
        legend.text=element_text(face="bold", size=14),
        legend.direction="horizontal",
        legend.position="bottom",
        legend.title=element_blank())
