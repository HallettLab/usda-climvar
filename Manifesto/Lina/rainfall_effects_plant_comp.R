library(tidyverse)
library(ggplot2)
library(ggpubr)

#function for standard error
se <- function(x){
  sd(x)/sqrt(length(x))# this is a function for calculating standard error
} 

#Rainfall and comp treatment effects on plant composition
setwd("./Dropbox/ClimVar/DATA/Plant_composition_data")
veg_all <- read.csv("Cover/Cover_CleanedData/ClimVar_species-cover.csv", header = TRUE)

#Filter 2015 data
veg_2015 <- veg_all %>% 
  filter(year == 2015) %>%
  filter(subplot != "XC") %>%
  filter(subplot != "C") %>%
  group_by(shelterBlock, treatment, subplot, func2) %>%
  summarise(cover = sum(cover)/100)

#Summarize veg comp by treatments and functional groups
veg_func2 <- veg_2015 %>%
  group_by(treatment, subplot, func2) %>%
  summarise(mean = mean(cover), 
            se = se(cover))

#Plot them
veg_func2$treatment <- factor(veg_func2$treatment, levels =c("controlRain",  "springDry", "fallDry","consistentDry"))
veg_func2$subplot <- factor(veg_func2$subplot, levels = c("B", "F", "G"), labels = c("Mixed", "Forb", "Grass"))
ggplot(veg_func2, aes(x = treatment, y = mean, col = func2))+
  facet_wrap(~subplot)+
  geom_point(position = position_dodge(0.5), size = 2)+
  theme_bw() +
  labs(y = "Cover", x = "Rainfall Treatment", color = "Functional Group") +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0, position = position_dodge(0.5), size =1)+
  scale_color_manual( labels = c("Forb", "Grass", "Legume" ), values = c("seagreen3", "dodgerblue","mediumpurple"))+
  scale_x_discrete(labels = c("Control", "Spring", "Fall", "Consistent"))