# Clean and compile climvar 2017 recruitment data
# author(s): CTW (caitlin.t.white@colorado.edu)
# intiated: May 2019 (code may change over course of analysis)

# Script purpose:
## 1. Read in raw recruitment experiment data
## 2. Scale stem counts to miniplot scale
## 3. Join rainfall treatment data and seeding key
## 4. Calculate percent recruitment
## 5. Write cleaned data .csv to Recruitment_CleanedData subfolder in ClimVar/Recruitment dropbox

# Notes: 
## -- germination whole plot size = 25cm x 25cm, subsample size = 5cm x 5cm
## -- rule: defer to whole plot numbers when available (even when subsample data present)

# > 2019-05-04, CTW: there are some projected plot stem counts and percent recruitment values that are unrealistic
# > 14 values in total, most B. nigra. Cleaning script includes code for 2 QA figures and writes out figs to figures folder on dropbox
# > sent figs to Julie and Lauren for feedback...



# -- SETUP -----
rm(list = ls())
library(readxl)
library(tidyverse)
options(stringsAsFactors = F)
theme_set(theme_bw())
na_vals = c(" ", "", NA, "NA")

# set path to recruitment data folder
datpath <- "~/Dropbox/ClimVar/Recruitment/Data/"

# read in data
# metadata (variable description)
recruit_meta <- read_excel(paste0(datpath, "Recruitment_EnteredData/recruitment_stems_spring2017.xlsx"),
                           sheet = "metadata",
                           na = na_vals, skip = 22)
# stem counts
recruit_dat <- read_excel(paste0(datpath, "Recruitment_EnteredData/recruitment_stems_spring2017.xlsx"),
                          sheet = "recruit_stems",
                          na = na_vals)
#plot notes
plot_notes <- read_excel(paste0(datpath, "Recruitment_EnteredData/recruitment_stems_spring2017.xlsx"),
                         sheet = "recruit_plot_notes",
                         na = na_vals)
# seeding key
seedkey <- read.csv(paste0(datpath, "Recruitment_EnteredData/ClimVar_CompGerm_seedingKey.csv"),
                    na.strings = na_vals, strip.white = T)
# species lookup table
spplist <- read.csv(paste0(datpath, "Recruitment_SpeciesKey.csv"),
                    na.strings = na_vals, strip.white = T)
# drought treatment key
treatment <- read.csv(paste0(datpath,"Shelter_key.csv")) 


# -- REVIEW RECRUITMENT DATA -----
str(recruit_dat)
summary(recruit_dat)
# what the NAs in stems?
data.frame(recruit_dat[is.na(recruit_dat$stems),]) # heads counted at different scale than stems
unique(recruit_dat$sample_date) # april date generic, may dates sprecific.. got coerced into character in data import
# what's in the stems_fecund columns
unique(recruit_dat$stems_fecund)
recruit_dat[recruit_dat$stems_fecund %in% 5:6,] # only 2 observations, BRNI; can drop this column from cleaned dataset
# review variables
View(recruit_meta)
# review plot notes
plot_notes # qualitative.. don't need to keep in cleaned dataset



# -- CLEAN DATA -----
recruit_clean <- recruit_dat %>%
  # convert sample dates to April/May (search by numeric string in case any NAs)
  mutate(sample_date = ifelse(grepl("[0-9]{5}", sample_date), "May 2017", sample_date)) %>%
  dplyr::select(-c(page,recorder, stems_fecund)) %>%
  # rename spcode to pair with spplist (later on)
  rename(code4 = spcode,
         # rename miniplot_wp to stem count wp (will put head count scale in separate column later)
         stemct_wp = miniplot_wp) %>%
  arrange(plot, code4, sample_date, sample2) %>%
  # create plot-species ID for splitting dataset
  mutate(ID = paste0(plot, code4),
         # create head count scale column so can keep all data in same row (e.g. if heads and stems counted at different scales)
         headct_wp = ifelse(is.na(heads), NA, stemct_wp)) %>% # if no heads counted, head scale is NA, else equal to stem count scale (for now, repeat samples need treatment)
  #move ID to first col, head count scale column after stem count scale columns
  dplyr::select(ID, plot:stemct_wp, headct_wp, sample2:QA_notes)

# remove obs that have multiple samples (i.e sample2 = 1) and deal with those separately, then add back in
repeats <- recruit_clean %>%
  subset(ID %in% recruit_clean$ID[duplicated(recruit_clean$ID)])

# remove repeat samples from main clean dataset (add back in when cleaned)
recruit_clean2 <- anti_join(recruit_clean, repeats)  

# proceed to cleaning repeated samples...
## rules: 
## if taken in same month, prioritize whole plot counts over subsamples (5x5cm)
## if two subsamples in same month, use average
## prioritize april over may samples (if repeat in may -- most were just checks on april counts)
## > can put in notes that saw x number of plants in may
## > could be that plants were miss-ID'd in April (e.g. grasses not yet developed) and therefore over counter
## > but could also be that number of plants were there in April, but some had died by May (altho CTW would probably have seen the dead plant)
## > either way, may re-checks not consistently done for all miniplots so can't apply to all april values, hence stick with april values
## if heads and stems counted at different scales, keep both rows (need to scale up differently)

View(repeats)  
# may require manual cleaning (want to append notes, varies by case.. there are 17 pairs)
# separate each sample event into its own data frame
repeats1 <- subset(repeats, sample2 == 0)
repeats2 <- subset(repeats, sample2 == 1)
# start with first pair..
#1 1ELGL: append may notes to april obseration
repeats_clean <- repeats1[1,]; repeats_clean$field_notes <- gsub("2", "2 in miniplot", paste(repeats2$field_notes[1], "in May")) 

#2 "3BRCA": add note about may count
repeats_clean[2,] <- repeats1[2,]; repeats_clean$field_notes[2] <- paste0(repeats1$field_notes[2],"; CTW counts ", repeats2$stems[2], " in miniplot in May")

#3 "3BRNI": 2 april subsamples, average stems and heads and add note
repeats_clean[3,] <- repeats1[3,]; repeats_clean$stems[3] <- mean(c(repeats1$stems[3], repeats2$stems[3])); repeats_clean$heads[3] <- mean(c(repeats1$heads[3], repeats2$heads[3]))
repeats_clean$field_notes[3] <- paste0("average of 2 april subsamples (stems, heads): ", 
                                       # april subsample 1
                                       str_flatten(repeats1[3, c("stems", "heads")], collapse = ", "), "; ",
                                       # april subsample 2
                                       str_flatten(repeats2[3, c("stems", "heads")], collapse = ", "))

#4 "3CLAM": average two subsamples from april
repeats_clean[4,] <- repeats1[4,]; repeats_clean$stems[4] <- mean(c(repeats1$stems[4], repeats2$stems[4])) # no head counts
repeats_clean$field_notes[4] <- gsub("^patchy", paste0("patchy; average of 2 april subsamples (stems): ", 
                                                       repeats1$stems[4], ", ", repeats2$stems[4]), repeats1$field_notes[4])

#5 "3MECA": append note about may count
repeats_clean[5,] <- repeats1[5,]; repeats_clean$field_notes[5] <- paste0(repeats1$field_notes[5], "; May sampling: CTW counts ", repeats2$stems[5], " stems. ", repeats2$field_notes[5])

#6 "3TACA": defer to whole plot count; append subsample note
repeats_clean[6,] <- repeats1[6,]; repeats_clean$field_notes[6] <- paste(repeats1$field_notes[6], repeats2$field_notes[6], sep = ". ")

#7 "4MECA": subsample note already in primary observation notes, just remove second sample
repeats_clean[7,] <- repeats1[7,]

#8 "4TACA": stems counted in 5x5cm, heads were whole plot; assign head number to primary row and change headct_wp to 1
repeats_clean[8,] <- repeats1[8,]
repeats_clean$heads[8] <- repeats2$heads[8]; repeats_clean$headct_wp[8] <- repeats2$headct_wp[8]

#9 "5ELGL": append may count and may field notes to april observation
repeats_clean[9,] <- repeats1[9,]; repeats_clean$field_notes[9] <- paste("May sampling: CTW also counts", repeats2$stems[9], "stems.", repeats2$field_notes[9] )

#10 "5STPU": append may count to april field notes
repeats_clean[10,] <- repeats1[10,]; repeats_clean$field_notes[10] <- paste("May sampling: CTW also counts", repeats2$stems[10], "stems.")

#11 "5TRHI": append may count and may field notes to april observation
repeats_clean[11,] <- repeats1[11,]; repeats_clean$field_notes[11] <- paste("May sampling: CTW counts", repeats2$stems[11], "stems.", repeats2$field_notes[11] )

#12 "6BRNI": average may subsamples, append note about averaging
repeats_clean[12,] <- repeats1[12,]; repeats_clean$stems[12] <- mean(c(repeats1$stems[12], repeats2$stems[12])); repeats_clean$heads[12] <- mean(c(repeats1$heads[12], repeats2$heads[12]))
repeats_clean$field_notes[12] <- paste(paste0("Average of 2 april subsamples (stems, heads): ", 
                                              # april subsample 1
                                              str_flatten(repeats1[12, c("stems", "heads")], collapse = ", "), "; ",
                                              # april subsample 2
                                              str_flatten(repeats2[12, c("stems", "heads")], collapse = ", ")), 
                                       # original field note after averaging note
                                       repeats1$field_notes[12], sep = ". ")

#13 "7TACA": Append subsample counts to field notes
repeats_clean[13,] <- repeats1[13,]
repeats_clean$field_notes[13] <- paste(repeats1$field_notes[13], 
                                       paste0("5x5cm patch subsample (stems, heads): ", 
                                              # stems
                                              str_flatten(repeats2[13, c("stems", "heads")], collapse = ", "))) 

#14 "9CLAM": defer to whole plot count, append subsample count to field notes (ahead of field notes)
repeats_clean[14,] <- repeats1[14,]; repeats_clean$field_notes[14] <- paste(repeats2$field_notes[14], repeats1$field_notes[14], sep = ". ")
# clean up note
repeats_clean$field_notes[14] <- gsub("[(?].* entry", "", repeats_clean$field_notes[14]) # take out comment to LMH
repeats_clean$field_notes[14] <- gsub("pilla/", "pillar/", repeats_clean$field_notes[14]) # fix typo

#15 "11CLAM": 2 april subsamples, average and add note to field notes
repeats_clean[15,] <- repeats1[15,]; repeats_clean$stems[15] <- mean(c(repeats1$stems[15], repeats2$stems[15])) # no heads counted
repeats_clean$field_notes[15] <- paste0("average of 2 april subsamples (stems only): ", 
                                        # april subsample 1
                                        str_flatten(c(repeats1$stems[15],repeats2$stems[15]), collapse = ", "))

#16 "12BRNI": defer to whole plot stem count, append subsample count to field notes, 
# calculate heads based on flowers/stem in subsample and add 30 to account for outlier, round up to nearest whole number
repeats_clean[16,] <- repeats1[16,]
# heads in whole plot = (stems whole plot minus 1) * (subsample 2 heads/stem) +30 flower outlier, all rounded to nearest whole number
repeats_clean$heads[16] <- round((repeats1$stems[16]-1)*(repeats2$heads[16]/repeats2$stems[16]) + 30, 0)
repeats_clean$headct_wp[16] <- 1 # calculated head count reps whole plot
# add note
repeats_clean$field_notes[16] <- paste0(repeats1$field_notes[16],". April subsample (stems, heads): ", repeats2$stems[16], ", ", repeats2$heads[16],
                                        ". Whole plot heads calculated by multiplying whole plot stems (minus 1) by heads per stem from subsample, adding 30 for outlier, and rounding to nearest whole number.")

#17 "15TRIN": stems counted in 5x5cm, heads were whole plot; assign head number to primary row and change headct_wp to 1
repeats_clean[17,] <- repeats1[17,]
repeats_clean$heads[17] <- repeats2$heads[17]; repeats_clean$headct_wp[17] <- repeats2$headct_wp[17]
repeats_clean$field_notes[17] <- paste(repeats1$stems[17], "stems counted in subsample,", 
                                       repeats2$heads[17], "heads counted in whole miniplot.")

# check that all miniplot are accounted for
if(sum(!repeats1$ID %in% unique(repeats_clean$ID))>0){
  print("Stop! Repeat sample miniplots missing from repeats_clean data frame:")
  print(repeats1$ID[!repeats1$ID %in% unique(repeats_clean$ID)])
} else{
  # rejoin repeat samples with main dataset
  recruit_clean3 <- rbind(recruit_clean2, repeats_clean) %>%
    arrange(plot, miniplot) %>%
    dplyr::select(-c(miniplot, sample2)) #remove miniplot and sample2 cols (no more repeat samples in dataset)
  # final check all obs accounted for
  if(length(recruit_clean3$ID) == length(unique(recruit_clean$ID))){
    # clean up
    rm(repeats1, repeats2, recruit_clean, recruit_clean2, repeats, repeats_clean)
    print("Cleaned up repeat observations joined with main dataset. Proceed to scaling up stem and head count with recruits_clean3!")
  }else{
    print("Miniplots missing:")
    print(!unique(recruit_clean$ID) %in% recruit_clean$ID)
  }
}



# -- SCALE UP STEM AND HEAD COUNTS -----
# assign plot area for scaling up (miniplot = 25x25cm)
# assume 1cm buffer on all sides
plot_area <- (25*25) - (2*(1*25) + 2*(1*23)) # this should correspond to whatever we think the area seeded actually was
plot_scale <- plot_area/(5*5) # subsample = 5x5cm

recruit_clean3$plot_stems <- with(recruit_clean3, ifelse(stemct_wp == 0, stems*plot_scale, stems))
recruit_clean3$plot_heads <- with(recruit_clean3, ifelse(headct_wp == 0, heads*plot_scale, heads))
# manual corrections to scaling based on field notes..
# review notes
recruit_clean3$field_notes[!is.na(recruit_clean3$field_notes)]
# 2 plots need manual corrections to scaling: 6LOMU, 7LOMU
# only 10x20cm area covered in whole plot
data.frame(recruit_clean3[grepl("10x20", recruit_clean3$field_notes),])
fix <- grep("10x20", recruit_clean3$field_notes)
recruit_clean3$plot_stems[fix] <- recruit_clean3$stems[fix] *((10*20)/(5*5))
recruit_clean3$plot_heads[fix] <- recruit_clean3$heads[fix] *((10*20)/(5*5))
# confirm correct
data.frame(recruit_clean3[grepl("10x20", recruit_clean3$field_notes),]) #yes

# review QA notes (okay to drop?)
recruit_clean3[!is.na(recruit_clean3$QA_notes),11:12] # LMH was supposed to revew all.. assume done by now. okay to drop


# -- ADD SEEDING, CALCULATE RECRUITMENT -----
recruit_clean4 <- recruit_clean3 %>%
  dplyr::select(-c(ID, QA_flag, QA_notes)) %>%
  left_join(seedkey[seedkey$treatment == "recruitment", c("code4", "seeds_per_plot")]) %>%
  rename(seedsAdded = seeds_per_plot) %>%
  mutate(pctRecruit = round(100 * (plot_stems/seedsAdded),2),
         flag_recruit = ifelse(pctRecruit>100, 1, 0),
         disturbed = ifelse(is.na(disturbed), 0, 1))

# look at distribution of percent recruitment by species
ggplot(recruit_clean4, aes(code4, pctRecruit)) +
  geom_hline(aes(yintercept = 100), col = "red") +
  annotate(geom = "text", x = "LUSU", y = 105, vjust = 0, label = "100% recruited", color = "red") +
  geom_boxplot() +
  geom_point(data = subset(recruit_clean4, flag_recruit == 1), aes(code4, pctRecruit), col = "red", alpha = 0.8) +
  labs(x = "Species",
       y = "Percent recruited\n(#indiviuals/#seeded)*100",
       title = paste("Percent recruitment QA:", 
                     length(recruit_clean4$flag_recruit[recruit_clean4$flag_recruit == 1]),  
                     "values exceed 100%"),
       subtitle = paste(length(recruit_clean4$flag_recruit[recruit_clean4$flag_recruit == 1 & recruit_clean4$code4 == "BRNI"]),
                        "of 16 BRNI observations flagged; flagged LOMU in plot invaded by outside LOMU"))
# it doesn't look that bad.. a few species have outliers, BRNI looks overcounted on nearly half
# LOMU outlier is due to outside LOMU invading plot, CTW could have overcounted TACA on some plots (e.g. one noted as being dense)
# CLAM was pretty dense on some plot, not sure what to do about high TRIN..
# write out figure to Recruitment dropbox folder
ggsave(paste0(datpath,"Recruitment_Figures/QA_pctRecruit_boxplots_byspecies.pdf"),
       height = 4, width = 6, scale = 1.1)


# number of flags by species
sapply(split(recruit_clean4$flag_recruit, recruit_clean4$code4), function(x) sum(x==1))


# -- JOIN DROUGHT TREATMENT AND SPPLIST DATA -----
# compile final dataset
#final set
recruit_combined <- recruit_clean4 %>%
  left_join(treatment) %>%
  left_join(spplist) %>%
  mutate(falltreatment = ifelse(grepl("Rain|spring", treatment), "wet", "dry"))

# is there any patterns to percent recruitment flags by drought treatment? (e.g. wet plots tend to have overcounts)
ggplot(subset(recruit_combined, flag_recruit == 1), aes(treatment, pctRecruit, col = code4)) +
  geom_jitter(aes(shape = as.factor(disturbed)), size = 2, height = 0, width = 0.2, alpha = 0.8) +
  labs(x = "Drought treatment",
       y = "Percent recruited\n(#indiviuals/#seeded)*100",
       title = "Most overcounts in fall wet plots (except Brassica nigra)",
       subtitle = "High LOMU = plot invaded by outside LOMU") +
  scale_color_discrete(name = "Species") +
  scale_shape_discrete(name = "Disturbed?")

# sort of .. wetter plots (not dry in fall) tend to have overcounts 
# BRNI overcounted in all treatment types, all other species in wetter plots (except 1 TRIN plot in fall dry)
# .. control plots were super invaded generally
# write out figure to Recruitment dropbox figures folder
ggsave(paste0(datpath,"Recruitment_Figures/QA_flagged_pctRecruit_values.pdf"),
       height = 4, width = 6)


# -- FINISHING ----
# re-arrange columns and write out
names(recruit_combined)
recruit_combined <- recruit_combined %>%
  dplyr::select(plot, falltreatment, treatment:shelter, code4, code6, stems, stemct_wp, plot_stems, 
               heads, headct_wp, plot_heads, seedsAdded:flag_recruit, sample_date:disturbed, species:nativity) %>%
  #drop a few USDA Plants database cols
  dplyr::select(-c(State_and_Province, Native_Status))

#write out
write.csv(recruit_combined, paste0(datpath, "Recruitment_CleanedData/Recruitment_cleaned.csv"), row.names = F)
