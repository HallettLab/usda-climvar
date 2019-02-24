setwd("./Dropbox/ClimVar/DATA/Decagon data")
setwd("~/Dropbox/ClimVar/DATA/Decagon data")

library(tidyverse)
library(ggplot2)
library(lubridate)
library(gridExtra)

data <- read.csv("ClimVar_sm_2015.csv", header = TRUE)
data$month <- as.factor(data$month)
levels(data$month)
levels(data$treatment)

#visualize sm
ggplot(data = data, aes(x = month, y = sm, fill = treatment)) + geom_boxplot()

#summarize sm 
CV <- function(x){(sd(x)/mean(x))*100}
SE <- function(x){(sd(x)/sqrt(length(x)))}

sm_data <- data %>%
  filter(sm >0) %>%
  filter(subplot == "B"| subplot == "G" | subplot == "F") %>%
  mutate(date = ymd_hms(date)) %>%
  filter(date > '2014-11-13 02:00:00',
         date < '2015-05-10 18:00:00') %>%
  mutate(season = ifelse(month %in% c(11:12, 1), "fall", "spring")) 

summary_sm <- sm_data %>%
  group_by(treatment, plot, subplot, shelterBlock, season) %>%
  summarise(mean_sm = mean(sm), se_sm = SE(sm), sd_sm = sd(sm), cv_sm = CV(sm))

#plot soil moisture by season
ggplot(data = sm_data, aes(x = season, y = sm, col = treatment)) +
  geom_boxplot() +
  theme_classic()

#soil moisture affects BNPP?
setwd("./Dropbox/ClimVar/DATA/Plant_composition_data/BNPP/BNPP_CleanedData")
BNPP0 <- read.csv("BNPP_MayHarvest_2015.csv", header = TRUE)
BNPP0$depth <- as.character(BNPP0$depth) #set depth as character
#BNPP summary
BNPP <- BNPP0 %>%
  unique() %>% #Remove potential duplicates
  filter(!is.na(bmass_g_m2)) %>% #Check for any missing values
  mutate(depth = replace(depth, depth == "20-Oct", "10-20")) %>%
  mutate(treatment_code = fct_recode(treatment, consDry = "consistentDry", control = "controlRain")) %>% #make a column with shorter treatment names 
  group_by(plot, subplot, treatment, treatment_code, shelterBlock, shelter, fall, spring) %>% #group data
  summarise(agg_BNPP = sum(bmass_g_m2))#sum BNPP 
#join soil moisture and BNPP summaries
joined <- BNPP %>%
  right_join(summary_sm, by = c("plot", "subplot", "treatment", "shelterBlock"))
#subset joined data
sm_fall <- joined %>%
  filter(season == "fall")
sm_spring <- joined %>%
  filter(season == "spring")

#plot BNPP against mean soil moisture
ggplot(sm_fall, aes(x=mean_sm, y=agg_BNPP, col = subplot, shape = treatment)) +
  geom_point() +
  theme_classic() +
  labs(x = "Mean Soil Moisture", y = expression(paste("Aggregated BNPP (g/m"^2,") in soil depth 0-30 cm")))
ggplot(sm_spring, aes(x=mean_sm, y=agg_BNPP, col=subplot, shape = treatment)) +
  geom_point() +
  theme_classic() +
  labs(x = "Mean Soil Moisture", y = expression(paste("Aggregated BNPP (g/m"^2,") in soil depth 0-30 cm")))

#add linear regression lines
ggplot(sm_fall, aes(x=mean_sm, y=agg_BNPP, col = subplot)) +
    geom_point() +
    geom_smooth(method = lm , size = 1, se = FALSE, fullrange = TRUE) +
    theme_classic() +
    labs(x = "Mean Soil Moisture", y = expression(paste("Aggregated BNPP (g/m"^2,") in soil depth 0-30 cm")))
ggplot(sm_spring, aes(x=mean_sm, y=agg_BNPP, col=subplot)) +
  geom_point() +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = TRUE)+
  theme_classic() +
  labs(x = "Mean Soil Moisture", y = expression(paste("Aggregated BNPP (g/m"^2,") in soil depth 0-30 cm")))

#plot BNPP against cv soil moisture
ggplot(sm_fall, aes(x=cv_sm, y=agg_BNPP, col = subplot, shape = treatment)) +
  geom_point() +
  theme_classic() +
  labs(x = "CV Soil Moisture", y = expression(paste("Aggregated BNPP (g/m"^2,") in soil depth 0-30 cm")))
ggplot(sm_spring, aes(x=cv_sm, y=agg_BNPP, col=subplot, shape = treatment)) +
  geom_point() +
  theme_classic() +
  labs(x = "CV Soil Moisture", y = expression(paste("Aggregated BNPP (g/m"^2,") in soil depth 0-30 cm")))

#how do traits trade-off with soil moisture?
setwd("./Dropbox/ClimVar/DATA/Plant_composition_data/BNPP/BNPP_CleanedData")
traits <- read.csv("traits_2015.csv", header = TRUE)

sm_trait <- traits %>%
  right_join(summary_sm, by = c("plot", "subplot", "treatment", "shelterBlock")) %>%
  group_by(subplot, treatment, season) %>%
  summarize(mean_PropF = mean(CWM.PropF),
            mean_SRLF = mean(CWM.SRLF),
            mean_SRLC = mean(CWM.SRLC),
            mean_DiamC = mean(CWM.DiamC),
            mean_Dens = mean(CWM.Dens),
            mean_SLA = mean(CWM.SLA),
            mean_LDMC = mean(CWM.LDMC),
            mean_Ht = mean(CWM.Ht),
            mean_sm = mean(mean_sm))

Fall <- sm_trait %>%
  filter(season == "fall")
Spring <- sm_trait %>%
  filter(season == "spring")

f1 <- ggplot(Fall, aes(x=mean_sm, y = mean_PropF, col = subplot, shape = treatment)) +
  geom_point() +
  theme_classic() +
  labs(x = "Fall Mean Soil Moisture", y = "CWM.PropF") +
  theme(legend.position = "none")

s1 <- ggplot(Spring, aes(x=mean_sm, y = mean_PropF, col = subplot, shape = treatment)) +
  geom_point() +
  theme_classic() +
  labs(x = "Spring Mean Soil Moisture", y = "CWM.PropF") +
  theme(legend.position = "none")

f2 <- ggplot(Fall, aes(x=mean_sm, y = mean_SRLF, col = subplot, shape = treatment)) +
  geom_point() +
  theme_classic() +
  labs(x = "Fall Mean Soil Moisture", y = "CWM.SRLF") +
  theme(legend.position = "none")

s2 <- ggplot(Spring, aes(x=mean_sm, y = mean_SRLF, col = subplot, shape = treatment)) +
  geom_point() +
  theme_classic() +
  labs(x = "Spring Mean Soil Moisture", y = "CWM.SRLF") +
  theme(legend.position = "none")

f3 <- ggplot(Fall, aes(x=mean_sm, y = mean_SRLC, col = subplot, shape = treatment)) +
  geom_point() +
  theme_classic() +
  labs(x = "Fall Mean Soil Moisture", y = "CWM.SRLC") +
  theme(legend.position = "none")

s3 <- ggplot(Spring, aes(x=mean_sm, y = mean_SRLC, col = subplot, shape = treatment)) +
  geom_point() +
  theme_classic() +
  labs(x = "Spring Mean Soil Moisture", y = "CWM.SRLC") +
  theme(legend.position = "none")

f4 <- ggplot(Fall, aes(x=mean_sm, y = mean_DiamC, col = subplot, shape = treatment)) +
  geom_point() +
  theme_classic() +
  labs(x = "Fall Mean Soil Moisture", y = "CWM.DiamC") +
  theme(legend.position = "none")

s4 <- ggplot(Spring, aes(x=mean_sm, y = mean_DiamC, col = subplot, shape = treatment)) +
  geom_point() +
  theme_classic() +
  labs(x = "Spring Mean Soil Moisture", y = "CWM.DiamC") +
  theme(legend.position = "none")

f5 <- ggplot(Fall, aes(x=mean_sm, y = mean_Dens, col = subplot, shape = treatment)) +
  geom_point() +
  theme_classic() +
  labs(x = "Fall Mean Soil Moisture", y = "CWM.Dens") +
  theme(legend.position = "none")

s5 <- ggplot(Spring, aes(x=mean_sm, y = mean_Dens, col = subplot, shape = treatment)) +
  geom_point() +
  theme_classic() +
  labs(x = "Spring Mean Soil Moisture", y = "CWM.Dens") +
  theme(legend.position = "none")

f6 <- ggplot(Fall, aes(x=mean_sm, y = mean_SLA, col = subplot, shape = treatment)) +
  geom_point() +
  theme_classic() +
  labs(x = "Fall Mean Soil Moisture", y = "CWM.SLA") +
  theme(legend.position = "none")

s6 <- ggplot(Spring, aes(x=mean_sm, y = mean_SLA, col = subplot, shape = treatment)) +
  geom_point() +
  theme_classic() +
  labs(x = "Spring Mean Soil Moisture", y = "CWM.SLA") +
  theme(legend.position = "none")

f7 <- ggplot(Fall, aes(x=mean_sm, y = mean_LDMC, col = subplot, shape = treatment)) +
  geom_point() +
  theme_classic() +
  labs(x = "Fall Mean Soil Moisture", y = "CWM.LDMC") +
  theme(legend.position = "none")

s7 <- ggplot(Spring, aes(x=mean_sm, y = mean_LDMC, col = subplot, shape = treatment)) +
  geom_point() +
  theme_classic() +
  labs(x = "Spring Mean Soil Moisture", y = "CWM.LDMC") +
  theme(legend.position = "none")

f8 <- ggplot(Fall, aes(x=mean_sm, y = mean_Ht, col = subplot, shape = treatment)) +
  geom_point() +
  theme_classic() +
  labs(x = "Fall Mean Soil Moisture", y = "CWM.Ht") +
  theme(legend.position = "none")

s8 <- ggplot(Spring, aes(x=mean_sm, y = mean_Ht, col = subplot, shape = treatment)) +
  geom_point() +
  theme_classic() +
  labs(x = "Spring Mean Soil Moisture", y = "CWM.Ht") +
  theme(legend.position = "none")

library(ggpubr)
ggpubr::ggarrange(f1, f2, f3, f4, f5, f6, f7, f8, nrow=2, ncol= 4, common.legend = TRUE, legend = "bottom")
ggpubr::ggarrange(s1, s2, s3, s4, s5, s6, s7, s8, nrow=2, ncol= 4, common.legend = TRUE, legend = "bottom")
