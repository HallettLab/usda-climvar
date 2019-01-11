setwd("./Dropbox/ClimVar/DATA/Decagon data")
setwd("~/Dropbox/ClimVar/DATA/Decagon data")

library(tidyverse)
library(ggplot2)

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
  mutate(season = ifelse(month %in% 9:12, "fall", ifelse(month %in% 2:5, "spring", ifelse(month %in% 6:8, "summer", "fall")))) %>%
  group_by(season,treatment, subplot) %>%
  summarise(mean_sm = mean(sm), se_sm = SE(sm), sd_sm = sd(sm), cv_sm = CV(sm))

#plot soil moisture by season
ggplot(data = sm_data, aes(x = season, y = mean_sm, col = treatment, shape = subplot)) +
  geom_point() +
  #geom_errorbar(aes(ymin=mean_sm-se_sm, ymax=mean_sm+se_sm), width=.1) +
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
  group_by(subplot, treatment) %>%
  summarise(mean_BNPP = mean(agg_BNPP), se_BNPP = SE(agg_BNPP), sd_BNPP = sd(agg_BNPP), cv_BNPP = CV(agg_BNPP)) %>%
  right_join(sm_data, by = c("subplot", "treatment"))
#subset joined data
sm_fall <- joined %>%
  filter(season == "fall")
sm_spring <- joined %>%
  filter(season == "spring")
sm_summer <- joined %>%
  filter(season == "summer")
#plot BNPP against soil moisture
ggplot(sm_fall, aes(x=mean_sm, y= mean_BNPP, col = treatment, shape = subplot)) +
  geom_point() +
  theme_classic()
ggplot(sm_spring, aes(x=mean_sm, y=mean_BNPP, col=treatment, shape=subplot)) +
  geom_point() +
  theme_classic()
ggplot(sm_summer, aes(x=mean_sm, y=mean_BNPP, col=treatment, shape=subplot)) +
  geom_point() +
  theme_classic()

