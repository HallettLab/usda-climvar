ANPP <- read.csv("./Dropbox/ClimVar/DATA/Plant_composition_data/ANPP/ANPP_CleanedData/ClimVar_ANPP-peak.csv", header = TRUE)
BNPP <- read.csv("./Dropbox/ClimVar/DATA/Plant_composition_data/BNPP/BNPP_CleanedData/BNPP_MayHarvest_2015.csv", header = TRUE)
aCWM <- read.csv("./Dropbox/ClimVar/DATA/Plant_composition_data/BNPP/BNPP_CleanedData/above_traits.csv", header = TRUE)
bCWM <- read.csv("./Dropbox/ClimVar/DATA/Plant_composition_data/BNPP/BNPP_CleanedData/Belowground_CWM_traits.csv", header = TRUE)
env <- read.csv("./Dropbox/ClimVar/DATA/Plant_composition_data/BNPP/BNPP_CleanedData/env_2015.csv", header = TRUE)

###How does total biomass relate to CWM traits?

##Add BNPP to ANPP for total biomass
library(tidyverse)
#Replace entry "20-Oct" to "10-20" in depth column in the BNPP dataset
BNPP$depth <- as.character(BNPP$depth) #set depth as character
BNPP <- BNPP %>%
  mutate(depth = replace(depth, depth == "20-Oct", "10-20")) %>%
  mutate(treatment_code = fct_recode(treatment, consDry = "consistentDry", control = "controlRain")) #make a column with shorter treatment names

#Calculate the aggregated BNPP in three levels of soil depth
BNPP1 <- BNPP %>%
  group_by(plot, subplot, treatment, shelterBlock) %>% #group data
  summarise(agg_BNPP = sum(bmass_g_m2)) %>% #sum BNPP
  select(plot, subplot, treatment, shelterBlock, agg_BNPP)

#Filter 2015 ANPP 
ANPP1 <- ANPP %>%
  filter(year == 2015) %>%
  filter(subplot %in% c("B", "F", "G")) %>%
  select(plot, subplot, treatment, shelterBlock, weight_g_m)

#Select aboveground CWM
above_CWM <- aCWM %>%
  select(shelterBlock, subplot, treatment, CWM.Ht, CWM.LDMC, CWM.SLA)

#Select belowground CWM
below_CWM <- bCWM %>%
  select(shelterBlock, subplot, treatment, nbsp, CWM.Dens, CWM.DiamC, CWM.SRLC, CWM.SRLF, CWM.PropF)

#Join dataframes
Joined <- ANPP1 %>%
  left_join(BNPP1, by = c("plot", "subplot", "treatment", "shelterBlock")) %>%
  left_join(above_CWM, by = c("subplot", "treatment", "shelterBlock")) %>%
  left_join(below_CWM, by = c("subplot", "treatment", "shelterBlock")) %>%
  mutate(total = weight_g_m + agg_BNPP)

#Standardize joined data
library(vegan)
stand_Joined_num <- decostand(Joined[,5:16], "standardize")
stand_Joined <- cbind(Joined[,1:4], stand_Joined_num)

#OR log transform joined data
#log_Joined_num <- log(Joined[, 5:16])
#log_Joined <- cbind(Joined[,1:4], log_Joined_num)

#Subset data by subplot
stand_both <- stand_Joined %>%
  filter(subplot == "B")
stand_forb <- stand_Joined %>%
  filter(subplot == "F") 
stand_grass <- stand_Joined %>%
  filter(subplot == "G")

##Backwards stepwise regression
library(MASS)
model0 <- lm(total ~ CWM.Ht + CWM.LDMC + CWM.SLA + CWM.Dens + CWM.DiamC + CWM.SRLC + CWM.SRLF + CWM.PropF, stand_both)
step_both <- stepAIC(model0, direction = "backward", trace = FALSE)
step_both$anova
model1 <- lm(total ~ CWM.Ht + CWM.LDMC + CWM.SLA + CWM.DiamC + CWM.SRLC + CWM.SRLF, stand_both)
summary(model1)

model2 <- lm(total ~ CWM.Ht + CWM.LDMC + CWM.SLA + CWM.Dens + CWM.DiamC + CWM.SRLC + CWM.SRLF + CWM.PropF, stand_forb)
step_forb <- stepAIC(model2, direction = "backward", trace = FALSE)
step_forb$anova
model3 <- lm(total ~ CWM.LDMC + CWM.SRLC + CWM.SRLF, stand_forb)
summary(model3)

model4 <- lm(total ~ CWM.Ht + CWM.LDMC + CWM.SLA + CWM.Dens + CWM.DiamC + CWM.SRLC + CWM.SRLF + CWM.PropF, stand_grass)
step_grass<- stepAIC(model4, direction = "backward", trace = FALSE)
step_grass$anova
model5 <- lm(total ~ CWM.Dens + CWM.SRLF + CWM.PropF, stand_grass)
summary(model5)

model6 <- lm(total ~CWM.Ht + CWM.LDMC + CWM.SLA + CWM.Dens + CWM.DiamC + CWM.SRLC + CWM.SRLF + CWM.PropF, stand_Joined)
step_total <- stepAIC(model6, direction = "backward", trace = FALSE)
step_total$anova
model7 <- lm(total ~ CWM.Ht + CWM.DiamC + CWM.SRLF, stand_Joined)
summary(model7)

###How does the environment affect BNPP?
#Select env of interest
env_data <- env %>%
  dplyr::select(plot, subplot, treatment, shelterBlock, SOC, SIC, MBC, MBIC, NH4, NO3.NO2, avg_sm, Litter_depth_cm, Percent_bare, Percent_litter)

#Join Biomass and env_data
Joined_env <- ANPP1 %>%
  left_join(BNPP1, by = c("plot", "subplot", "treatment", "shelterBlock")) %>%
  left_join(env_data, by = c("plot", "subplot", "treatment", "shelterBlock")) %>%
  mutate(total = weight_g_m + agg_BNPP)

#Log transform data
log_Joined_env_num <- log(Joined_env[,5:17]) ##There are NaN/Inf
log_Joined_env <- cbind(Joined_env[,1:4], log_Joined_env_num)

#Instead standardize data
library(vegan)
stand_env_num <- decostand(Joined_env[,5:17], "standardize")
stand_env <- cbind(Joined_env[,1:4], stand_env_num)

#Subset data by subplot
stand_both_env <- stand_env %>%
  filter(subplot == "B")
stand_forb_env <- stand_env %>%
  filter(subplot == "F") %>%
  filter(!is.na(avg_sm))
stand_grass_env <- stand_env %>%
  filter(subplot == "G")

#Backwards stepwise regression
env0 <- lm(agg_BNPP ~ SOC + SIC + MBC + MBIC + NH4 + NO3.NO2 + avg_sm + Litter_depth_cm + Percent_bare + Percent_litter, stand_both_env)
step_both_env <- stepAIC(env0, direction = "backward", trace = FALSE)
step_both_env$anova
env1 <- lm(agg_BNPP ~ Percent_bare + Percent_litter, stand_both_env)
summary(env1)

env2 <- lm(agg_BNPP ~ SOC + SIC + MBC + MBIC + NH4 + NO3.NO2 + avg_sm + Litter_depth_cm + Percent_bare + Percent_litter, stand_forb_env)
step_forb_env <- stepAIC(env2, direction = "backward", trace = FALSE)
step_forb_env$anova
env3 <- lm(agg_BNPP ~ MBC + NH4 + NO3.NO2 + Litter_depth_cm, stand_forb_env)
summary(env3)

env4 <- lm(agg_BNPP ~ SOC + SIC + MBC + MBIC + NH4 + NO3.NO2 + avg_sm + Litter_depth_cm + Percent_bare + Percent_litter, stand_grass_env)
step_grass_env <- stepAIC(env2, direction = "backward", trace = FALSE)
step_grass_env$anova
env5 <- lm(agg_BNPP ~ MBC + NH4 + NO3.NO2 + Litter_depth_cm, stand_grass_env)
summary(env5)

###How does the environment affect total biomass?
env6 <- lm(total ~ SOC + SIC + MBC + MBIC + NH4 + NO3.NO2 + avg_sm + Litter_depth_cm + Percent_bare + Percent_litter, stand_both_env)
step_both_total <- stepAIC(env6, direction = "backward", trace = FALSE)
step_both_total$anova
env7 <- lm(total ~ MBC + MBIC + NH4 + avg_sm + Litter_depth_cm + Percent_bare, stand_both_env)
summary(env7)

env8 <- lm(total ~ SOC + SIC + MBC + MBIC + NH4 + NO3.NO2 + avg_sm + Litter_depth_cm + Percent_bare + Percent_litter, stand_forb_env)
step_forb_total <- stepAIC(env8, direction = "backward", trace = FALSE)
step_forb_total$anova
env9 <- lm(total ~ MBIC + NO3.NO2 + avg_sm + Litter_depth_cm + Percent_bare, stand_forb_env)
summary(env9)

env10 <- lm(total ~ SOC + SIC + MBC + MBIC + NH4 + NO3.NO2 + avg_sm + Litter_depth_cm + Percent_bare + Percent_litter, stand_grass_env)
step_grass_total <- stepAIC(env10, direction = "backward", trace = FALSE)
step_grass_total$anova
env11 <- lm(total ~ MBC, stand_grass_env)
summary(env11)

###How does the environment affect aboveground biomass?
env12 <- lm(weight_g_m ~ SOC + SIC + MBC + MBIC + NH4 + NO3.NO2 + avg_sm + Litter_depth_cm + Percent_bare + Percent_litter, stand_grass_env)
step_grass_above <- stepAIC(env12, direction = "backward", trace = FALSE)
step_grass_above$anova
env13 <- lm(weight_g_m ~ 1, stand_grass_env)
summary(env13)

env14 <- lm(weight_g_m ~ SOC + SIC + MBC + MBIC + NH4 + NO3.NO2 + avg_sm + Litter_depth_cm + Percent_bare + Percent_litter, stand_forb_env)
step_forb_above <- stepAIC(env14, direction = "backward", trace = FALSE)
step_forb_above$anova
env15 <- lm(weight_g_m ~ NO3.NO2, stand_forb_env)
summary(env15)

env16 <- lm(weight_g_m ~ SOC + SIC + MBC + MBIC + NH4 + NO3.NO2 + avg_sm + Litter_depth_cm + Percent_bare + Percent_litter, stand_both_env)
step_both_above <- stepAIC(env16, direction = "backward", trace = FALSE)
step_both_above$anova
env17 <- lm(weight_g_m ~ 1, stand_both_env)
summary(env17)

#Diversity vs BNPP/Total Biomass plot
library(ggplot2)
ggplot(Joined, aes(x = nbsp, y = agg_BNPP, col = subplot)) +
  geom_point() +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE)+
  labs(x = "Species Richness", y = "BNPP (g) in depth 0-30 cm") +
  annotate("text", x= 6, y = 550, label = "R2 = 0.04", color = "#F8766D") +
  annotate("text", x= 6, y = 350, label = "R2 = 0.005", color = "#00BFC4") +
  annotate("text", x= 11, y = 200, label = "R2 = 0.02", color = "#00ba38") +
  theme_bw() +
  labs(col = "treatment")

ggplot(Joined, aes(x = nbsp, y = total, col = subplot)) +
  geom_point() +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE)+
  labs(x = "Species Richness", y = "Total biomass (g)") +
  theme_classic()
Both <- lm(agg_BNPP ~ nbsp, data = log_both)
summary(Both)
Forb <- lm(agg_BNPP ~ nbsp, data = log_forb)
summary(Forb)
Grass <- lm(agg_BNPP ~ nbsp, data = log_grass)
summary(Grass)

#ANPP and BNPP by subplot
library(ggpubr)
p1 <- ggplot(Joined, aes(x = subplot, y = weight_g_m, fill = subplot)) +
  geom_boxplot() +
  theme_bw() +
  annotate("text", x = 1, y = 1000, label = "a", size = 4) +
  annotate("text", x = 2, y = 1000, label = "b", size = 4) +
  annotate("text", x = 3, y = 1000, label = "ab", size = 4) +
  labs(y = "ANPP g/m2", x = "composition treatment") +
  scale_fill_discrete(name = "treatment", labels = c("Mixed", "Forb", "Grass")) +
  scale_x_discrete(labels = c("Mixed", "Forb", "Grass"))

p2 <- ggplot(Joined, aes(x = subplot, y = agg_BNPP, fill = subplot)) +
  geom_boxplot() +
  theme_bw() +
  annotate("text", x = 1, y = 800, label = "a", size = 4) +
  annotate("text", x = 2, y = 800, label = "b", size = 4) +
  annotate("text", x = 3, y = 800, label = "b", size = 4) +
  labs(y = "BNPP g/m2 depth 0-30 cm", x = "composition treatment") +
  scale_fill_discrete(name = "treatment", labels = c("Mixed", "Forb", "Grass")) +
  scale_x_discrete(labels = c("Mixed", "Forb", "Grass"))

ggarrange(p1, p2, ncol = 2, nrow = 1, 
          common.legend = TRUE, legend = "right",
          align = "v",labels = c("a)", "b)"))

fitA <- lm(weight_g_m~subplot, Joined)
TukeyHSD(aov(fitA))
fitB <- lm(agg_BNPP~subplot, Joined)
TukeyHSD(fitB)

#Biomass by rain treatment
p3 <- ggplot(Joined, aes(x = treatment, y = agg_BNPP, fill = treatment)) +
        geom_boxplot() +
        theme_bw() +
        labs(y = "BNPP g/m2 depth 0-30 cm", x = "precip treatment", fill = "treatment") +
        scale_fill_manual(values= c("#0070b8", "#b2c7e4", "#fccaaf", "#c85b23"))

p4 <- ggplot(Joined, aes(x = treatment, y = weight_g_m, fill = treatment)) +
        geom_boxplot() +
        theme_bw() +
        theme(legend.position = "none") +
        labs(y = "ANPP g/m2", x = "precip treatment", fill = "treatment") +
        scale_fill_manual(values= c("#0070b8", "#b2c7e4", "#fccaaf", "#c85b23"))

ggarrange(p3, p4, ncol = 2, nrow = 1, 
          common.legend = TRUE, legend = "right",
          align = "v",labels = c("a)", "b)"))

fit1 <- lm(weight_g_m~treatment*subplot, Joined)
anova(fit1)
fit2 <- lm(agg_BNPP~treatment*subplot, Joined)
anova(fit2)

#Regression ANPP ~ CWM.Ht
Ht_all <- lm(weight_g_m ~ CWM.Ht, Joined)
summary(Ht_all) #significant
Ht_both<- lm(weight_g_m ~ CWM.Ht, Joined%>%filter(subplot == "B")) 
summary(Ht_both) #significant
Ht_forb <- lm(weight_g_m ~ CWM.Ht, Joined%>%filter(subplot == "F"))
summary(Ht_forb) #not significant
Ht_grass <- lm(weight_g_m ~ CWM.Ht, Joined%>%filter(subplot == "G"))
summary(Ht_grass) #not significant

#Regression ANPP ~ CWM.LDMC
LDMC_all <- lm(weight_g_m ~ CWM.LDMC, Joined)
summary(LDMC_all) #significant
LDMC_both<- lm(weight_g_m ~ CWM.LDMC, Joined%>%filter(subplot == "B")) 
summary(LDMC_both) #not significant
LDMC_forb <- lm(weight_g_m ~ CWM.LDMC, Joined%>%filter(subplot == "F"))
summary(LDMC_forb) #not significant
LDMC_grass <- lm(weight_g_m ~ CWM.LDMC, Joined%>%filter(subplot == "G"))
summary(LDMC_grass) #not significant

#Regression ANPP ~ CWM.SLA
SLA_all <- lm(weight_g_m ~ CWM.SLA, Joined)
summary(SLA_all) #significant
SLA_both<- lm(weight_g_m ~ CWM.SLA, Joined%>%filter(subplot == "B")) 
summary(SLA_both) #not significant
SLA_forb <- lm(weight_g_m ~ CWM.SLA, Joined%>%filter(subplot == "F"))
summary(SLA_forb) #not significant
SLA_grass <- lm(weight_g_m ~ CWM.SLA, Joined%>%filter(subplot == "G"))
summary(SLA_grass) #not significant



