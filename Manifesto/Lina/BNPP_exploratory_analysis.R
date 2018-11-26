####DATA SET UP

#set working directory
setwd("C:/Users/Lina/Dropbox/ClimVar/DATA/Plant_composition_data/BNPP/BNPP_CleanedData")
#load the dataset
BNPP0 <- read.csv("BNPP_MayHarvest_2015.csv", header = TRUE)

setwd("C:/Users/Lina/Dropbox/ClimVar/DATA/Plant_composition_data/ANPP/ANPP_CleanedData")
ANPP0 <- read.csv("ClimVar_ANPP-peak.csv", header = TRUE)

#set data as a tibble

library("tidyverse")

BNPP <- as.tibble(BNPP0)
str(BNPP)

ANPP <- as.tibble(ANPP0)
str(ANPP)

#replace entry "20-Oct" to "10-20" in depth column
BNPP$depth <- as.character(BNPP$depth)
BNPP <- BNPP %>%
  mutate(depth = replace(depth, depth == "20-Oct", "10-20"))

#calculate total BNPP 
BNPP1 <- BNPP %>%
  group_by(plot, subplot, treatment, shelterBlock, shelter, fall, spring) %>%
  summarise(agg_BNPP = sum(bmass_g_m2))

#calculate root:shoot ratio
joined <- ANPP %>%
  filter(year == 2015) %>%
  filter(subplot == "F" | subplot == "B" | subplot == "G") %>%
  inner_join(BNPP1, by = c( "plot", "subplot", "treatment", "shelterBlock")) %>%
  mutate(root_shoot = agg_BNPP/weight_g_m)

####DATA VISUALIZATION

library(ggplot2)
library(ggpubr)

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

#plot biomass vs. treatment 
par(mfrow = c(1,3))
boxplot(Top10_biomass ~ Top10_trt, main = "Depth 0-10")
boxplot(Mid10_biomass ~ Mid10_trt, main = "Depth 10-20")
boxplot(Bottom10_biomass ~ Bottom10_trt, main = "Depth 20-30")

#plot biomass vs. treatment (group by depth)
ggplot(BNPP, aes(x=treatment, y=bmass_g_m2, fill=depth, color=depth)) +
  geom_boxplot() +
  theme_classic()

#plot biomass vs. treatment (group by functional group)
Top10_plot <- ggplot(Top10, aes(x=treatment, y=bmass_g_m2, fill=subplot, color=subplot)) +
  geom_boxplot() +
  theme_classic()

Mid10_plot <- ggplot(Mid10, aes(x=treatment, y=bmass_g_m2, fill=subplot, color=subplot)) +
  geom_boxplot() +
  theme_classic()

Bottom10_plot <- ggplot(Bottom10, aes(x=treatment, y=bmass_g_m2, fill=subplot, color=subplot)) +
  geom_boxplot() +
  theme_classic()

ggarrange(Top10_plot, Mid10_plot, Bottom10_plot, 
          labels = c( " Depth 0-10" , " Depth 10-20", "Depth 20-30"),
          ncol = 1,
          nrow = 3)

#plot aggregated biomass vs. treatment by functional group
total_biomass <- ggplot(BNPP1, aes(x=treatment, y=agg_BNPP, fill = subplot)) +
  geom_boxplot() +
  theme_classic()

#plot root:shoot ratio vs. treatment by functional group 
root_shoot_plot <- ggplot(joined, aes(x=treatment, y=root_shoot, fill = subplot)) +
  geom_boxplot() +
  theme_classic()

####ANALYSIS

library(nlme) 
library(multcomp)

###One way ANOVA biomass ~ treatment
Top10ANOVA <- aov(Top10_biomass ~ Top10_trt)
summary(Top10ANOVA)

Mid10ANOVA <- aov(Mid10_biomass ~ Mid10_trt)
summary(Mid10ANOVA)

Bottom10ANOVA <- aov(Bottom10_biomass ~ Bottom10_trt)
summary(Bottom10ANOVA)

###Randomized block ANOVA biomass ~ functional group * treatment (block as random)
#check assumptions for randomized block ANOVA

#Normality of the response variable at each level of the factor
BNPP.agg <- with(BNPP, aggregate(data.frame(bmass_g_m2), 
                                    by = list(A = treatment, B = subplot), mean))
boxplot(bmass_g_m2 ~ A, BNPP.agg) 
#Homogeneity of variance
with(BNPP.agg, plot(tapply(bmass_g_m2, A, var), 
                      tapply(bmass_g_m2, A, mean)))

#determine whether the design is balanced
!is.list(replications(bmass_g_m2 ~ treatment + subplot + depth, BNPP)) #returns FALSE
replications(bmass_g_m2 ~ treatment + subplot + depth, BNPP) #unbalanced 

#fit mixed model using lme
Top10.lme <- lme(data = Top10, Top10_biomass ~ Top10_trt * subplot, random = ~1| shelterBlock)
anova(Top10.lme)
summary(Top10.lme)
plot(resid(Top10.lme) ~fitted(Top10.lme))

Mid10.lme <- lme(data = Mid10, Mid10_biomass ~ Mid10_trt * subplot, random = ~1| shelterBlock)
anova(Mid10.lme)
summary(Mid10.lme)
plot(resid(Mid10.lme) ~fitted(Mid10.lme))

Bottom10.lme <- lme(data = Bottom10, Bottom10_biomass ~ Bottom10_trt * subplot, random = ~1| shelterBlock)
anova(Bottom10.lme)
summary(Bottom10.lme)
plot(resid(Bottom10.lme) ~fitted(Bottom10.lme))

total.lme <- lme(data = joined, agg_BNPP ~ treatment * subplot, random = ~1| shelterBlock)
anova(total.lme)
summary(total.lme)

#subset dataset
Both <- joined %>%
  filter(subplot == "B")
Grass <- joined %>%
  filter(subplot == "G")
Forb <- joined %>%
  filter(subplot == "F")

#