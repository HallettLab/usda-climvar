#Caitlin White
#CU Boulder, EBIO & INSTAAR
#October 2016

#Native Diversity ANPP preliminary exploration

### ------> NOTE <------ ####
# CTW cleaned up this script May 2019 so it runs (i.e. minimal update)..
# incorporated code from other scripts to consolidate scripts
# needs more work..
# *****************************


#load needed libraries
library(dplyr)
library(ggplot2)
options(stringsAsFactors = F)
theme_set(theme_bw())
na_vals <- c(" ", "", "NA", NA)

#set pathway to cover data
datpath <- "~/Dropbox/ClimVar/DATA/Plant_composition_data/Native_Diversity/"

#read in data
#Native Diversity ANPP
Native.ANPP <- read.csv(paste0(datpath,"Native_Diversity_CleanedData/NativeDiversity_ANPP_2015-2016_cleaned.csv"), 
                        na.strings = na_vals, strip.white = T)
#Nat Diversity cover
Native.cover <- read.csv(paste0(datpath,"Native_Diversity_CleanedData/NativeDiversity_Cover_2015-2016_cleaned.csv"), 
                         na.strings = na_vals, strip.white = T)


# -- EXPLORE DATA -----
# parking piece of older code, fix later:
##Add columns to shelter key frame to indicate rainfall treatment as binary data, where 1 = receives rain and 0 = rain excluded##
# key <- mutate(key, spring = ifelse(treatment == "fallDry" | treatment =="controlRain", 1,0), 
#               fall = ifelse(treatment == "fallDry" | treatment == "consistentDry", 0,1)) %>%
#   mutate(spring=as.factor(spring), fall = as.factor(fall)) %>%
#   tbl_df()


#Summarize cover data by plot, nativity, and year to add in to ANPP data
Native.cover.grouped <- subset(Native.cover, plot != "All") %>%
  mutate(plot = as.numeric(plot)) %>%
  group_by(yr, plot, treatment, shelterBlock, shelter, nativity) %>%
  summarise(cover = sum(cover)) %>%
  ungroup()
  
Native.all <- full_join(Native.cover.grouped, Native.ANPP, by = c("plot", "treatment", "shelterBlock", "shelter", "yr", "nativity" = "group"))
#change "NA" to 0 values for wgt_g (so every cover point will have pair)
Native.all$wgt_g[is.na(Native.all$wgt_g)] <- 0
Native.all$nativity[is.na(Native.all$nativity)] <- "Native"
Native.all$cover[is.na(Native.all$cover) & grepl("ESCA", Native.all$notes)] <- 0.5
Native.all$cover[is.na(Native.all$cover)] <- 0

#Visualize data
ggplot(Native.all) + geom_boxplot(aes(nativity, wgt_g, fill=nativity)) + 
  facet_grid(yr ~ treatment) + 
  labs(y="Weight (g)", x=NULL) +
  theme_bw() + 
  coord_flip() +
  theme(axis.title.x=element_text(face="bold", size=20),
        axis.text.x=element_text(size=16),
        axis.text.y = element_blank(),
        strip.text=element_text(face="bold", size=16),
        legend.text=element_text(face="bold", size=16),
        # legend.direction="horizontal",
        # legend.position= c(0.8, 0),
        legend.title=element_blank())


ggplot(Native.all, aes(cover, wgt_g, color=nativity)) + geom_point() +
  stat_smooth(method="lm", se=F) +
  facet_grid(yr ~ treatment)

ggplot(Native.all, aes(nativity, wgt_g)) + 
  geom_boxplot() +
  geom_jitter(alpha = 0.8, width = 0.2, height = 0) +
  facet_grid(yr ~ treatment)

#can't look at just native vs. non-native forb ANPP because don't have that data (only grouped non-native species)


#ANOVA tests for ANPP and cover grouped native vs. non-native
#Just consider 2015 data, 2016 has issues due to grassy, wetter year and lack of seeding/weeding Fall 2015
cover.anova <- aov(cover ~ nativity + treatment + shelterBlock + nativity*treatment, data = subset(Native.all, yr == "2015"))
summary(cover.anova)

anpp.anova <- aov(wgt_g ~ nativity + treatment + nativity*treatment, data = subset(Native.all, yr == "2015"))
summary(anpp.anova)

lm.cover <- lm(cover ~ nativity + treatment + nativity*treatment, data = subset(Native.all, yr == "2015"))
summary(lm.cover)

lm.anpp <- lm(wgt_g ~ nativity + treatment + nativity*treatment, data = subset(Native.all, yr == "2015"))
summary(lm.anpp)

fxnl.cover.anova <- aov(cover ~ nativity + fxnl_grp + treatment + nativity*treatment + fxnl_grp*treatment + nativity*fxnl_grp,
                        data = subset(Native.cover, yr == "2015"))
summary(fxnl.cover.anova)

lm.fxnl.cover <- lm(cover ~ nativity + fxnl_grp + treatment + nativity*treatment + fxnl_grp*treatment + nativity*fxnl_grp,
                    data = subset(Native.cover, yr == "2015"))
summary(lm.fxnl.cover)


# parking older code from 2015, need to update..
##DATA ANALYSIS FOR NATIVE DIVERSITY ANPP##

#Visualization of Native Diversity ANPP data##
# pdf("Native_ANPP_rainfallresponse.pdf")
# ggplot(togNDanpp, aes(x=treatment, y=weight)) + geom_boxplot() + facet_wrap(~species) + labs(title="Native ANPP by treatment")
# dev.off()
# 
# #Attempt at ANOVA on Native Diversity ANPP, does not yield anything significant#
# NDmodel <- lme(weight~ rain_treat*species, random=~1 |shelterBlock, data=togNDanpp,
#                contrasts=list(rain_treat=contr.treatment, species=contr.treatment))
# anova(NDmodel)
# summary(glht(NDmodel, linfct=mcp(rain_treat="Tukey", species="Tukey")))
# 
# #ANOVA on Native Diversity ANPP, using created treatment_species interaction variable#
# NDmodel2 <- lme(weight~ spec_treat, random=~1 |shelterBlock, data=togNDanpp,
#                 contrasts=list(spec_treat=contr.treatment))
# summary(glht(NDmodel2, linfct=mcp(spec_treat="Tukey")))
