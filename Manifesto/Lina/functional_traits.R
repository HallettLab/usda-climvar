setwd("./Dropbox/ClimVar/DATA/Plant_composition_data/BNPP/BNPP_CleanedData")
BNPP0 <- read.csv("BNPP_MayHarvest_2015.csv", header = TRUE)
traits <- read.csv("traits_2015.csv", header = TRUE)
library(tidyverse)

#BNPP summary
BNPP0$depth <- as.character(BNPP0$depth) #set depth as character
BNPP <- BNPP0 %>%
  unique() %>% #Remove potential duplicates
  filter(!is.na(bmass_g_m2)) %>% #Check for any missing values
  mutate(depth = replace(depth, depth == "20-Oct", "10-20")) %>%
  mutate(treatment_code = fct_recode(treatment, consDry = "consistentDry", control = "controlRain")) %>% #make a column with shorter treatment names 
  group_by(plot, subplot, treatment, treatment_code, shelterBlock, shelter, fall, spring) %>% #group data
  summarise(agg_BNPP = sum(bmass_g_m2))#sum BNPP 

#join functional traits and BNPP summaries
BNPP_traits <- BNPP %>%
  right_join(traits, by = c("plot", "subplot", "shelterBlock", "treatment"))
BNPP_traits$AvDom <- as.factor(BNPP_traits$AvDom)

summary <- BNPP_traits %>%
  group_by(subplot, treatment) %>%
  summarize(mean_PropF = mean(CWM.PropF),
            mean_SRLF = mean(CWM.SRLF),
            mean_SRLC = mean(CWM.SRLC),
            mean_DiamC = mean(CWM.DiamC),
            mean_Dens = mean(CWM.Dens),
            mean_SLA = mean(CWM.SLA),
            mean_LDMC = mean(CWM.LDMC),
            mean_Ht = mean(CWM.Ht),
            mean_BNPP = mean(agg_BNPP))

#visualize BNPP vs Proportion of fine roots
p1 <- ggplot(BNPP_traits, aes(x = CWM.PropF, y = agg_BNPP, col = subplot, shape = treatment)) +
  geom_point() +
  theme_classic() +
  theme(legend.position = "none")

g1 <- ggplot(summary, aes(x = mean_PropF, y  = mean_BNPP, col = treatment, shape = subplot)) +
  geom_point() +
  theme_classic() #+
  theme(legend.position = "none")

#visualize BNPP vs Specific root length fine
p2 <- ggplot(BNPP_traits, aes(x = CWM.SRLF, y = agg_BNPP, col = subplot, shape = treatment)) +
  geom_point() +
  theme_classic() +
  theme(legend.position = "none")

g2 <- ggplot(summary, aes(x = mean_SRLF, y  = mean_BNPP, col = subplot, shape = treatment)) +
  geom_point() +
  theme_classic() +
  theme(legend.position = "none")

#visualize BNPP vs Specific root length coarse
p3 <- ggplot(BNPP_traits, aes(x = CWM.SRLC, y = agg_BNPP, col = subplot, shape = treatment)) +
  geom_point() +
  theme_classic() +
  theme(legend.position = "none")

g3 <- ggplot(summary, aes(x = mean_SRLC, y  = mean_BNPP, col = subplot, shape = treatment)) +
  geom_point() +
  theme_classic() +
  theme(legend.position = "none")

#visualize BNPP vs Coarse Root Diameter
p4 <- ggplot(BNPP_traits, aes(x = CWM.DiamC, y = agg_BNPP, col = subplot, shape = treatment)) +
  geom_point() +
  theme_classic() +
  theme(legend.position = "none")

g4 <- ggplot(summary, aes(x = mean_DiamC, y  = mean_BNPP, col = subplot, shape = treatment)) +
  geom_point() +
  theme_classic() +
  theme(legend.position = "none")

#visualize BNPP vs Root Density
p5 <- ggplot(BNPP_traits, aes(x = CWM.Dens, y = agg_BNPP, col = subplot, shape = treatment)) +
  geom_point() +
  theme_classic() +
  theme(legend.position = "none")

g5 <- ggplot(summary, aes(x = mean_Dens, y  = mean_BNPP, col = subplot, shape = treatment)) +
  geom_point() +
  theme_classic() +
  theme(legend.position = "none")

#visualize BNPP vs Specific Leaf Area
p6 <- ggplot(BNPP_traits, aes(x = CWM.SLA, y = agg_BNPP, col = subplot, shape = treatment)) +
  geom_point() +
  theme_classic() +
  theme(legend.position = "none")

g6 <- ggplot(summary, aes(x = mean_SLA, y  = mean_BNPP, col = subplot, shape = treatment)) +
  geom_point() +
  theme_classic() +
  theme(legend.position = "none")

#visualize BNPP vs Leaf Dry Matter Content
p7 <- ggplot(BNPP_traits, aes(x = CWM.LDMC, y = agg_BNPP, col = subplot, shape = treatment)) +
  geom_point() +
  theme_classic() +
  theme(legend.position = "none")

g7 <- ggplot(summary, aes(x = mean_LDMC, y  = mean_BNPP, col = subplot, shape = treatment)) +
  geom_point() +
  theme_classic() 

#visualize BNPP vs Height
p8 <- ggplot(BNPP_traits, aes(x = CWM.Ht, y = agg_BNPP, col = subplot, shape = treatment)) +
  geom_point() +
  theme_classic() +
  theme(legend.position = "none")

g8 <- ggplot(summary, aes(x = mean_Ht, y  = mean_BNPP, col = subplot, shape = treatment)) +
  geom_point() +
  theme_classic() 

library(ggpubr)
ggpubr::ggarrange(g1, g2, g3, g4, g5, g6, g7, g8, ncol = 4, nrow = 2,
             common.legend = TRUE,
             legend = "bottom")

#visualize BNPP vs avena cover 
ggplot(BNPP_traits, aes(x = AvCover, y = agg_BNPP, col = subplot, shape = treatment)) +
  geom_point() +
  theme_classic()  #no trend

#visualize BNPP vs avena dominance
ggplot(BNPP_traits, aes(x = AvDom, y = agg_BNPP, col = subplot)) +
  geom_boxplot() +
  theme_bw() +
  facet_wrap(~treatment) #see spring dry mixed plot 

#visualize BNPP vs percent bare 
ggplot(BNPP_traits, aes(x = Percent_bare, y = agg_BNPP, col = subplot, shape = treatment)) +
  geom_point() +
  theme_classic()

#visualize BNPP vs Litter depth
ggplot(BNPP_traits, aes(x = Litter_depth_cm, y = agg_BNPP, col = subplot, shape = treatment)) +
  geom_point() +
  theme_classic()

