# generate seeding key (lookup table) for competition and recruitment experiments
# mar 2019
# contact: caitlin.t.white@colorado.edu

#script purpose:
# create lookup table for # seeds added to competition background low and high density treatments + recruitment miniplots in 2017 experiment.
# competition low density background and recruitment plots were seeded at rate of 4g per m^2
# competition high density background seeded at rate of 16g/m^2
# germination in ambient seeded at 4g/m^2 (which is why some germination experiment spp have seeding calculations at 1g (for 0.5x0.5m^2 GA plots))

# final lookup will contain:
# 1) species name
# 2) species 4-letter code
# 3) species 6-letter code (to link with JL data)
# 3) pertinent experiment
# 4) pertinent treatment (hi/low for comp exp; recruitment or GA (germination in ambient) for germ exp)
# 5) prescribed seed mass rate
# 6) plot area
# 7) weight seeded per plot (adjusted for viability and seed attachments on AVFA and TACA)
# 8) mean seed weight per species
# 9) seeds added (estimated per seeds added)

# write out final seeding lookup table to recruitment and competition entered data folders; this script gets pushed to github repo



# -- SETUP -----
library(readxl)
library(tidyverse)
options(stringsAsFactors = F)
theme_set(theme_bw())
na_vals <- c("", " ", NA, "NA")

#set absolute path to dropbox folder to retrieve and write out data
## --> can vary by user, enter code for your path below <--
dropbox_path <- "../../Dropbox/ClimVar/" #ctw's path

#read in excel data:
#seeding instructions from LMH to CTW
# >note: there are multiple rows (weights) per species in the germination experiment.. which row is correct?
seed_prep <- read_excel(paste0(dropbox_path, "Competition/Setup/CompGerm_seed_order_prep/SeedWeighing_Guide.xlsx"),
                        sheet = 1, na = na_vals, trim_ws = T)
#ctw per seed weights for AV, ERO, and TACA
seed_wgt <- read_excel(paste0(dropbox_path, "Competition/Setup/CompGerm_seed_order_prep/SeedWeighing_Guide.xlsx"),
                       sheet = 2, na = na_vals, trim_ws = T)
#kellen ryan seed weights **error in clarkia amoena
larson_seedwgt <- read_excel(paste0(dropbox_path, "Recruitment/Data/Recruitment_SeedWeights/CA_dryseedmass.xlsx"),
                                    sheet = 1, na = na_vals, trim_ws = T)

  
# -- DATA PREP -----
str(seed_prep)
#remove unneeded cols
seed_prep <- seed_prep[!colnames(seed_prep) %in% c("X__1", "X__2", "X__3", "X__4", "X__5")]
#review note for AVFA, ERO, and TACA:
seed_prep$Experiment[48]

# retain cols needed in ryan data
str(larson_seedwgt)
# check perseedwt
summary((larson_seedwgt$drymass_g/larson_seedwgt$no_seeds) - larson_seedwgt$perseedwt) #should be all 0s.. yes!
# check reps
sort(unique(larson_seedwgt$rep))

# check variation per species
ggplot(larson_seedwgt, aes(species, perseedwt)) +
  geom_point(size = 2, alpha = 0.5) +
  stat_summary(col = "orchid", alpha = 0.5) +
  ggtitle("Recruitment: raw weights, per species, multiple reps each (black), with mean and se (purple)") +
  coord_flip() +
  theme(plot.title = element_text(size = 12))
# variability all fairly small, AVFA varies a little more and ERO (but ERBO not used in any experiment)

larson_seedmeans <- larson_seedwgt %>%
  #correct spelling for clarkia amoena
  mutate(species = ifelse (species == "claarm", "claamo", species)) %>%
  #keeping all reps? since methods varied for number of seeds weighed.. more info the better
  group_by(species) %>%
  summarise(mean_perseed_wgt_g = mean(perseedwt),
            se_perseed_wgt_g = sd(perseedwt)/sqrt(length(perseedwt)))

# per seed weight for AVFA, ERO and TACA (field-collected seeds)
str(seed_wgt)
field_seeds <- seed_wgt[,1:3] %>%
  mutate(rep = readr::parse_number(Species),
         Species = gsub(" |[0-9]", "", Species))
colnames(field_seeds)[1:3] <- c("species", "gross_wgt_g", "seed_only_wgt_g")
field_seeds <- field_seeds[c("species", "rep", "gross_wgt_g", "seed_only_wgt_g")]

# look at variability
field_seeds %>%
  gather(var, val, gross_wgt_g:seed_only_wgt_g) %>%
  ggplot(aes(0, val)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.5, width = 0.1) +
  stat_summary(col = "orchid") + #by default plots mean with 1 se
  labs(y = "Weight (g)") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  facet_grid(species ~ var, scales = "free_y")

# ctw calculated seed masses for AVFA, ERBO and TACA.. close to Julie's, will keep Julie's weights since same methods and intrumentation with other species seed weights
# ctw's weights are slightly greater than Julie's
field_seed_means <- field_seeds %>%
  group_by(species) %>%
  summarise(mean_perseed_wgt_g = mean(gross_wgt_g),
            se_perseed_wgt_g = sd(gross_wgt_g)/sqrt(length(gross_wgt_g)))



# -- COMPILE -----  
# 1) prep germination seeding rates
germ_seeding <- subset(seed_prep, grepl("germ", Experiment)) %>%
  # correct treatment for GA plots
  mutate(Treatment = ifelse(conversion == 0.25, "GA (germination in ambient)", "recruitment")) %>%
  # remove spp not used in experiment
  subset(!grepl("yello|microstach|botrys|lentis", Species)) %>%
  dplyr::select(Experiment:'g/m2', finalwgt) %>%
  mutate(plot_area = ifelse(Treatment == "recruitment", 0.25*0.25, 0.5*0.5),
         code6 = paste0(substr(Species,1,3), substr(gsub("[a-z]+ ", "", Species), 1,3))) %>%
  left_join(larson_seedmeans[c("species", "mean_perseed_wgt_g")], by = c("code6" = "species"))
# clean up colnames
colnames(germ_seeding) <- casefold(colnames(germ_seeding))
colnames(germ_seeding)[grepl("m2" ,colnames(germ_seeding))] <- "seed_rate_g_m2"
colnames(germ_seeding)[grepl("final" ,colnames(germ_seeding))] <- "final_wgt_seeded_g"

# 2) prep competition seeding rates
comp_seeding <- subset(seed_prep, grepl("com", Experiment)) %>%
  # remove spp not used in experiment
  subset(!grepl("yello|microstach|botrys|lentis", Species)) %>%
  dplyr::select(Experiment:'g/m2', finalwgt) %>%
  mutate(Treatment = gsub("_.+", "", Treatment), #remove alt designation, since actually used
           plot_area = 0.5*0.5,
         code6 = paste0(substr(Species,1,3), substr(gsub("[a-z]+ ", "", Species), 1,3))) %>%
  left_join(larson_seedmeans[c("species", "mean_perseed_wgt_g")], by = c("code6" = "species"))
# clean up colnames
colnames(comp_seeding) <- casefold(colnames(comp_seeding))
colnames(comp_seeding)[grepl("m2" ,colnames(comp_seeding))] <- "seed_rate_g_m2"
colnames(comp_seeding)[grepl("final" ,colnames(comp_seeding))] <- "final_wgt_seeded_g"

# 3) combine germ and comp tables
seed_key <- rbind(comp_seeding, germ_seeding) %>%
  mutate(seeds_per_plot = round(final_wgt_seeded_g/mean_perseed_wgt_g))
#capitalize first letter
seed_key$species <- with(seed_key, paste0(casefold(substr(species, 1,1), upper = T),
                  substr(species, 2, nchar(species))))
seed_key$code4 <- with(seed_key, casefold(paste0(substr(species, 1,2), substr(gsub("[A-Z][a-z]+ ", "", species),1,2)), upper = T))
#reorder columns
seed_key <- seed_key %>%
  dplyr::select(experiment:species,code4, code6, treatment:seed_rate_g_m2, plot_area, final_wgt_seeded_g, mean_perseed_wgt_g, seeds_per_plot)
# add m2 to plot_area name for unit info
colnames(seed_key)[colnames(seed_key) == "plot_area"] <- "plot_area_m2"

#write out
# to recruitment folder
write.csv(seed_key, paste0(dropbox_path, "Recruitment/Data/Recruitment_EnteredData/ClimVar_CompGerm_seedingKey.csv"), row.names = F)
# to competition folder
write.csv(seed_key, paste0(dropbox_path, "Competition/Data/Competition_EnteredData/ClimVar_CompGerm_seedingKey.csv"), row.names = F)
