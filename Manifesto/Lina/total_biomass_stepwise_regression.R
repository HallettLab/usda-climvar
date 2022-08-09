ANPP <- read.csv("./Dropbox/ClimVar/DATA/Plant_composition_data/ANPP/ANPP_CleanedData/ClimVar_ANPP-peak.csv", header = TRUE)
BNPP <- read.csv("./Dropbox/ClimVar/DATA/Plant_composition_data/BNPP/BNPP_CleanedData/BNPP_MayHarvest_2015.csv", header = TRUE)
CWM <- read.csv("./Dropbox/ClimVar/DATA/Plant_composition_data/BNPP/BNPP_CleanedData/Belowground_CWM_traits.csv", header = TRUE)
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
  dplyr::select(plot, subplot, treatment, shelterBlock, agg_BNPP)

#Filter 2015 ANPP 
ANPP1 <- ANPP %>%
  filter(year == 2015) %>%
  filter(subplot %in% c("B", "F", "G")) %>%
  dplyr::select(plot, subplot, treatment, shelterBlock, weight_g_m)

#Standardize trait data
library(vegan)
stand_CWM_num <- decostand(CWM[,15:22], "standardize")
stand_CWM <- cbind(CWM[,1:5], stand_CWM_num)

#Join dataframes
stand_Joined <- ANPP1 %>%
  left_join(BNPP1, by = c("plot", "subplot", "treatment", "shelterBlock")) %>%
  left_join(stand_CWM, by = c("plot", "subplot", "treatment", "shelterBlock")) %>%
  mutate(total = weight_g_m + agg_BNPP)
Joined <- ANPP1 %>%
  left_join(BNPP1, by = c("plot", "subplot", "treatment", "shelterBlock")) %>%
  left_join(CWM, by = c("plot", "subplot", "treatment", "shelterBlock")) %>%
  mutate(total = weight_g_m + agg_BNPP)
##Backwards stepwise regression
library(MASS)
model0 <- lm(total ~CWM.Ht + CWM.LDMC + CWM.SLA + CWM.Dens + CWM.DiamC + CWM.SRLC + CWM.SRLF + CWM.PropF, stand_Joined)
step_total <- stepAIC(model0, direction = "backward", trace = FALSE)
step_total$anova
model1 <- lm(total ~ CWM.LDMC + CWM.DiamC + CWM.SRLC + CWM.SRLF + CWM.PropF, stand_Joined)
summary(model1)

model2 <- lm(weight_g_m ~CWM.Ht + CWM.LDMC + CWM.SLA , stand_Joined)
step_total <- stepAIC(model2, direction = "backward", trace = FALSE)
step_total$anova
model3 <- lm(weight_g_m ~CWM.Ht , stand_Joined)
summary(model3)

model4 <- lm(agg_BNPP ~CWM.Dens + CWM.DiamC + CWM.SRLC + CWM.SRLF + CWM.PropF , stand_Joined)
step_total <- stepAIC(model4, direction = "backward", trace = FALSE)
step_total$anova
model5 <- lm(weight_g_m ~CWM.DiamC, stand_Joined)
summary(model5)
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

#ANPP and BNPP by subplot Fig1
library(ggpubr)
p1 <- ggplot(Joined, aes(x = subplot, y = weight_g_m, fill = subplot)) +
  geom_boxplot() +
  theme_classic() +
  theme(legend.position = "none") +
  annotate("text", x = 1, y = 800, label = "a", size = 4) +
  annotate("text", x = 2, y = 800, label = "b", size = 4) +
  annotate("text", x = 3, y = 800, label = "ab", size = 4) +
  labs(y = bquote('ANPP'~(g/m^2)), x = "") +
  ylim(120, 800)+
  scale_fill_manual(name = "Treatment", labels = c("Mixed", "Forb", "Grass"),  values= c("#fc8d62", "#66c2a5", "#8da0cb")) +
  scale_x_discrete(labels = c("Mixed", "Forb", "Grass"))

p2 <- ggplot(Joined, aes(x = subplot, y = agg_BNPP, fill = subplot)) +
  geom_boxplot() +
  theme_classic() +
  theme(legend.position = "none") +
  annotate("text", x = 1, y = 800, label = "a", size = 4) +
  annotate("text", x = 2, y = 800, label = "b", size = 4) +
  annotate("text", x = 3, y = 800, label = "b", size = 4) +
  labs(y = bquote('BNPP'~(g/m^2)), x = "Composition treatment") +
  ylim(120, 800)+
  scale_fill_manual(name = "Treatment", labels = c("Mixed", "Forb", "Grass"),  values= c("#fc8d62", "#66c2a5", "#8da0cb")) +
  scale_x_discrete(labels = c("Mixed", "Forb", "Grass"))

pt_subplot <- ggplot(Joined, aes(x = subplot, y = total, fill = subplot)) +
  geom_boxplot() +
  theme_classic() +
  theme(legend.position = "none") +
  annotate("text", x = 1, y = 1600, label = "a", size = 4) +
  annotate("text", x = 2, y = 1600, label = "b", size = 4) +
  annotate("text", x = 3, y = 1600, label = "b", size = 4) +
  labs(y = bquote('Total Biomass'~(g/m^2)), x = "") +
  #ylim(120, 800)+
  scale_fill_manual(name = "Treatment", labels = c("Mixed", "Forb", "Grass"),  values= c("#fc8d62", "#66c2a5", "#8da0cb")) +
  scale_x_discrete(labels = c("Mixed", "Forb", "Grass"))

fitA <- aov(weight_g_m~subplot, Joined)
TukeyHSD(aov(fitA))
fitB <- aov(agg_BNPP~subplot, Joined)
TukeyHSD(fitB)
fitC <- aov(total~subplot, Joined)
TukeyHSD(fitC)

#Biomass by rain treatment 
Joined$treatment <- factor(Joined$treatment, levels = c("controlRain",  "springDry", "fallDry","consistentDry"))
p3 <- ggplot(Joined, aes(x = treatment, y = agg_BNPP)) +
        geom_boxplot() +
        geom_jitter(aes(x = treatment, y = agg_BNPP, color = subplot)) +
        theme_classic() +
        ylim(100, 700)+
        annotate("text", x = 2.5, y = 700, label = "NS", size = 4) +
        labs(y = bquote('BNPP'~(g/m^2)), x = "Rainfall treatment", fill = "treatment") +
        scale_color_manual( name = "Treatment", labels = c("Mixed", "Forb", "Grass"),  values= c("#fc8d62", "#66c2a5", "#8da0cb")) +
        scale_x_discrete(labels = c("Control\n", "Spring\n Dry\n", "Fall\n Dry\n", "Consistent\n Dry\n"))

p4 <- ggplot(Joined, aes(x = treatment, y = weight_g_m)) +
        geom_boxplot() +
        geom_jitter(aes(x = treatment, y = weight_g_m, color = subplot)) +
        theme_classic() +
        theme(legend.position = "none") +
        ylim(100, 700)+
        annotate("text", x = 2.5, y = 700, label = "NS", size = 4) +
        labs(y = bquote('ANPP'~(g/m^2)), x = "", fill = "treatment") +
        scale_color_manual( values= c("#fc8d62", "#66c2a5", "#8da0cb")) +
        #scale_fill_manual(name = "Treatment", labels = c("Control","Spring Dry", "Fall Dry", "Consistent Dry"), values= c("#0070b8", "#b2c7e4", "#fccaaf", "#c85b23")) +
        scale_x_discrete(labels = c("Control\n", "Spring\n Dry\n", "Fall\n Dry\n", "Consistent\n Dry\n"))

pt_rain <- ggplot(Joined, aes(x = treatment, y = total)) +
  geom_boxplot() +
  theme_classic() +
  geom_jitter(aes(x = treatment, y = total, color = subplot)) +
  theme(legend.position = "none") +
  ylim(200, 1500)+
  annotate("text", x = 2.5, y = 1500, label = "NS", size = 4) +
  labs(y = bquote('Total Biomass'~(g/m^2)), x = "", fill = "treatment") +
  scale_color_manual( values= c("#fc8d62", "#66c2a5", "#8da0cb")) +
  #scale_fill_manual(name = "Treatment", labels = c("Control","Spring Dry", "Fall Dry", "Consistent Dry"), values= c("#0070b8", "#b2c7e4", "#fccaaf", "#c85b23")) +
  scale_x_discrete(labels = c("Control\n", "Spring\n Dry\n", "Fall\n Dry\n", "Consistent\n Dry\n"))

fit1 <- aov(weight_g_m~treatment, Joined)
TukeyHSD(fit1)
fit2 <- aov(agg_BNPP~treatment, Joined)
TukeyHSD(fit2)
legend_1 <- get_legend(p1)
legend_2 <- get_legend(p3)
p1 <- p1 + theme(legend.position = "none")
p3 <- p3 + theme(legend.position = "none")
ggarrange(p1, p2, pt_subplot, 
          p4, p3, pt_rain, 
          ncol = 3, nrow = 2, 
          align = "v", common.legend = TRUE, legend = "right",
          labels = c("a)", "b)", "c)",  "d)", "e)", "f)"), 
          widths=c(1, 1 ,1))

###Relationships bw CWM of three aboveground traits and ANPP
p1a <- ggplot(data = Joined, aes(x = CWM.Ht, y = weight_g_m)) +
  geom_point(aes(color = as.factor(subplot))) +
  theme_classic() +
  xlab("CWM Height") +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE, color = "black")+
  labs(y= NULL) +
  theme(legend.position="none") +
  scale_color_manual(name = "treatment", labels = c("Mixed", "Forb", "Grass"), values = c("#fc8d62", "#66c2a5", "#8da0cb") )+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))

p2a <- ggplot(data = Joined, aes(x = CWM.SLA, y = weight_g_m)) +
  geom_point(aes(color = as.factor(subplot))) +
  theme_classic() +
  xlab("CWM Specific Leaf Area") +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE, color = "black")+
  labs(y= NULL) +
  theme(legend.position= "none")+
  scale_color_manual(name = "treatment", labels = c("Mixed", "Forb", "Grass"), values = c("#fc8d62", "#66c2a5", "#8da0cb") )+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))

p3a <- ggplot(data = Joined, aes(x = CWM.LDMC, y = weight_g_m)) +
  geom_point(aes(color = as.factor(subplot))) +
  theme_classic() +
  xlab("CWM Leaf Dry Matter Content") +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE, color = "black")+
  labs(y= NULL) +
  theme(legend.position = "none")+
  scale_color_manual(name = "treatment", labels = c("Mixed", "Forb", "Grass"), values = c("#fc8d62", "#66c2a5", "#8da0cb") )+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))


figurea <- ggarrange(p1a, p2a, p3a, 
                    ncol =3, nrow =1, common.legend = TRUE, legend = "bottom",
                    align = "v",labels = c("a)", "b)", "c)"))
annotate_figure(figurea, 
                left = text_grob("ANPP (g/m2)", rot = 90))

###Relationships bw CWM of five root traits and BNPP
p1 <- ggplot(data = Joined, aes(x = CWM.Dens, y = agg_BNPP)) +
  geom_point(aes(color = as.factor(subplot))) +
  theme_classic() +
  xlab("CWM Root Density") +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE, color = "black")+
  labs(y= NULL) +
  theme(legend.position="none") +
  scale_color_manual(name = "treatment", labels = c("Mixed", "Forb", "Grass"), values = c("#fc8d62", "#66c2a5", "#8da0cb") )+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))

p2 <- ggplot(data = Joined, aes(x = CWM.SRLF, y = agg_BNPP)) +
  geom_point(aes(color = as.factor(subplot))) +
  theme_classic() +
  xlab("CWM SRLF") +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE, color = "black")+
  labs(y= NULL) +
  theme(legend.position= "none")+
  scale_color_manual(name = "treatment", labels = c("Mixed", "Forb", "Grass"), values = c("#fc8d62", "#66c2a5", "#8da0cb") )+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))

p3 <- ggplot(data = Joined, aes(x = CWM.SRLC, y = agg_BNPP)) +
  geom_point(aes(color = as.factor(subplot))) +
  theme_classic() +
  xlab("CWM SRLC") +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE, color = "black")+
  labs(y= NULL) +
  theme(legend.position = "none")+
  scale_color_manual(name = "treatment", labels = c("Mixed", "Forb", "Grass"), values = c("#fc8d62", "#66c2a5", "#8da0cb") )+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))

p4 <- ggplot(data = Joined, aes(x = CWM.DiamC, y = agg_BNPP)) +
  geom_point(aes(color = as.factor(subplot))) +
  theme_classic() +
  xlab("CWM Diameter Coarse") +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE, color = "black")+
  labs(y= NULL) +
  theme(legend.position = "none")+
  scale_color_manual(name = "treatment", labels = c("Mixed", "Forb", "Grass"), values = c("#fc8d62", "#66c2a5", "#8da0cb") )+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))

p5 <- ggplot(data = Joined, aes(x = CWM.PropF, y = agg_BNPP)) +
  geom_point(aes(color = as.factor(subplot))) +
  theme_classic() +
  xlab("CWM Proportion of Fine") +
  labs(y= NULL) +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE, color = "black")+
  theme(legend.position = "none")+
  scale_color_manual(name = "treatment", labels = c("Mixed", "Forb", "Grass"), values = c("#fc8d62", "#66c2a5", "#8da0cb") )+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))

figure <- ggarrange(p1, p2, p3, p4, p5,
                    ncol =3, nrow =2, common.legend = TRUE, legend = "bottom",
                    align = "v",labels = c("a)", "b)", "c)", "d)", "e)"))
annotate_figure(figure, 
                left = text_grob("BNPP (g/m2) depth 0-30 cm", rot = 90))


#Regression ANPP ~ CWM.Ht Table S5
Ht_all <- lm(weight_g_m ~ CWM.Ht, Joined)
summary(Ht_all) #significant
Ht_both<- lm(weight_g_m ~ CWM.Ht, Joined%>%filter(subplot == "B")) 
summary(Ht_both) #not significant
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
summary(SLA_all) #not significant
SLA_both<- lm(weight_g_m ~ CWM.SLA, Joined%>%filter(subplot == "B")) 
summary(SLA_both) #not significant
SLA_forb <- lm(weight_g_m ~ CWM.SLA, Joined%>%filter(subplot == "F"))
summary(SLA_forb) #not significant
SLA_grass <- lm(weight_g_m ~ CWM.SLA, Joined%>%filter(subplot == "G"))
summary(SLA_grass) #not significant

#Regression BNPP ~ CWM.Dens
Dens_all <- lm(agg_BNPP ~ CWM.Dens, Joined)
summary(Dens_all) #not significant
Dens_both<- lm(agg_BNPP ~ CWM.Dens, Joined%>%filter(subplot == "B")) 
summary(Dens_both) #not significant
Dens_forb <- lm(agg_BNPP ~ CWM.Dens, Joined%>%filter(subplot == "F"))
summary(Dens_forb) #not significant
Dens_grass <- lm(agg_BNPP ~ CWM.Dens, Joined%>%filter(subplot == "G"))
summary(Dens_grass) #not significant

#Regression BNPP ~ CWM.SRLF
SRLF_all <- lm(agg_BNPP ~ CWM.SRLF, Joined)
summary(SRLF_all) #not significant
SRLF_both<- lm(agg_BNPP ~ CWM.SRLF, Joined%>%filter(subplot == "B")) 
summary(SRLF_both) #not significant
SRLF_forb <- lm(agg_BNPP ~ CWM.SRLF, Joined%>%filter(subplot == "F"))
summary(SRLF_forb) #not significant
SRLF_grass <- lm(agg_BNPP ~ CWM.SRLF, Joined%>%filter(subplot == "G"))
summary(SRLF_grass) #not significant

#Regression BNPP ~ CWM.SRLC
SRLC_all <- lm(agg_BNPP ~ CWM.SRLC, Joined)
summary(SRLC_all) #not significant
SRLC_both<- lm(agg_BNPP ~ CWM.SRLC, Joined%>%filter(subplot == "B")) 
summary(SRLC_both) #not significant
SRLC_forb <- lm(agg_BNPP ~ CWM.SRLC, Joined%>%filter(subplot == "F"))
summary(SRLC_forb) #not significant
SRLC_grass <- lm(agg_BNPP ~ CWM.SRLC, Joined%>%filter(subplot == "G"))
summary(SRLC_grass) #not significant

#Regression BNPP ~ CWM.DiamC
DiamC_all <- lm(agg_BNPP ~ CWM.DiamC, Joined)
summary(DiamC_all) #not significant
DiamC_both<- lm(agg_BNPP ~ CWM.DiamC, Joined%>%filter(subplot == "B")) 
summary(DiamC_both) #not significant
DiamC_forb <- lm(agg_BNPP ~ CWM.DiamC, Joined%>%filter(subplot == "F"))
summary(DiamC_forb) #not significant
DiamC_grass <- lm(agg_BNPP ~ CWM.DiamC, Joined%>%filter(subplot == "G"))
summary(DiamC_grass) #significant

#Regression BNPP ~ CWM.PropF
PropF_all <- lm(agg_BNPP ~ CWM.PropF, Joined)
summary(PropF_all) #not significant
PropF_both<- lm(agg_BNPP ~ CWM.PropF, Joined%>%filter(subplot == "B")) 
summary(PropF_both) #not significant
PropF_forb <- lm(agg_BNPP ~ CWM.PropF, Joined%>%filter(subplot == "F"))
summary(PropF_forb) #not significant
PropF_grass <- lm(agg_BNPP ~ CWM.PropF, Joined%>%filter(subplot == "G"))
summary(PropF_grass) #not significant

#Table S3 Mixed model of ppt and comp effects on biomass 
m_anpp <- lme(weight_g_m ~ treatment*subplot, random=~1|shelterBlock, Joined, na.action=na.exclude)
anova(m_anpp)
m_bnpp <-lme(agg_BNPP ~ treatment*subplot, random=~1|shelterBlock, Joined, na.action=na.exclude)
anova(m_bnpp)
m_total<- lme(total ~ treatment*subplot, random=~1|shelterBlock, Joined, na.action=na.exclude)
anova(m_total)

#Table S4 T-statisitcs and tukey of comp effect on biomass
summary(lm(weight_g_m~treatment, Joined))
summary(lm(agg_BNPP~treatment, Joined))
summary(lm(total~treatment, Joined))

