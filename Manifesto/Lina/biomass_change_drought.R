#Goal: To understand biomass response to drought by composition 
library(tidyverse)
library(ggplot2)

#function for standard error
calcSE<-function(x){
  x <- x[!is.na(x)]
  sd(x)/sqrt(length(x))}

#Load ANPP1 and BNPP1 from "total_biomass_stepwise_regression.R"
#calculate relative change in biomass compared to control

avg_control_ANPP <- ANPP1 %>%
  filter(treatment == "controlRain") %>%
  group_by(subplot)%>%
  summarize(mean_controlRain = mean(weight_g_m))

avg_control_BNPP <- BNPP1 %>%
  filter(treatment == "controlRain") %>%
  group_by(subplot) %>%
  summarize(mean_controlRain = mean(agg_BNPP))

change_ANPP <- left_join(ANPP1, avg_control_ANPP) %>%
  group_by(subplot, treatment) %>%
  mutate(change_ANPP = (weight_g_m-mean_controlRain)/mean_controlRain) %>%
  filter(treatment != "controlRain") %>%
  summarize(mean_change_ANPP = mean(change_ANPP),
            se_change_ANPP = calcSE(change_ANPP))
  
change_BNPP <- left_join(BNPP1, avg_control_BNPP) %>%
  group_by(subplot, treatment) %>%
  mutate(change_BNPP = (agg_BNPP-mean_controlRain)/mean_controlRain) %>%
  filter(treatment != "controlRain") %>%
  summarize(mean_change_BNPP = mean(change_BNPP),
            se_change_BNPP = calcSE(change_BNPP))

ggplot(change_ANPP, aes(x = subplot, y = mean_change_ANPP))+
  geom_point()+
  facet_grid(~treatment)
