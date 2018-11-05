#set working directory
setwd("C:/Users/Lina/Dropbox/ClimVar/DATA/Plant_composition_data/BNPP/BNPP_CleanedData")
#load the dataset
BNPP0 <- read.csv("BNPP_MayHarvest_2015.csv", header = TRUE)

#load package "tidyverse"
library("tidyverse")
#set data as a tibble
BNPP <- as.tibble(BNPP0)
BNPP
levels(BNPP$treatment)

#replace entry "20-Oct" to "10-20" in depth column
BNPP$depth <- as.character(BNPP$depth)
BNPP <- BNPP %>%
  mutate(depth = replace(depth, depth == "20-Oct", "10-20"))

#visualize data with a boxplot
library(ggplot2)
ggplot(BNPP, aes(x=treatment, y=bmass_g_m2, fill=depth, color=depth)) +
  geom_boxplot() +
  theme_classic()

####ANOVA biomass by treatment 
install.packages("multcomp")
library(multcomp)

#subset the data by depth
Top10 <- BNPP %>%
  filter(depth == "0-10")
Mid10 <- BNPP %>%
  filter(depth == "10-20")
Bottom10 <- BNPP %>%
  filter(depth == "20-30")

#name variables
Top10_trt <- Top10$treatment
Top10_biomass <- Top10$bmass_g_m2
Mid10_trt <- Mid10$treatment
Mid10_biomass <- Mid10$bmass_g_m2
Bottom10_trt <- Bottom10$treatment
Bottom10_biomass <- Bottom10$bmass_g_m2

#ANOVA
Top10ANOVA <- aov(Top10_biomass ~ Top10_trt)
summary(Top10ANOVA)
boxplot(Top10_biomass ~ Top10_trt)
summary(glht(Top10ANOVA, linfct = mcp(Top10_trt = "Tukey")))

Mid10ANOVA <- aov(Mid10_biomass ~ Mid10_trt)
summary(Mid10ANOVA)
boxplot(Mid10_biomass ~ Mid10_trt)

Bottom10ANOVA <- aov(Bottom10_biomass ~ Bottom10_trt)
summary(Bottom10ANOVA)
boxplot(Bottom10_biomass ~ Bottom10_trt)



