#Goal: To understand biomass response to drought by composition 
library(tidyverse)
library(ggplot2)
library(ggpubr)

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
  mutate(change_ANPP = (weight_g_m-mean_controlRain)/mean_controlRain) 

summary_change_ANPP <- change_ANPP %>%
  filter(treatment != "controlRain") %>%
  summarize(mean_change_ANPP = mean(change_ANPP),
            se_change_ANPP = calcSE(change_ANPP))
summary_change_ANPP$treatment <- factor(summary_change_ANPP$treatment, levels = c( "springDry", "fallDry","consistentDry"),
                                        labels = c("Spring\n Dry", "Fall\n Dry", "Consistent\n Dry"))
summary_change_ANPP$subplot <- factor(summary_change_ANPP$subplot, levels = c( "B", "F","G"),
                                        labels = c("Mixed", "Forb", "Grass"))

change_BNPP <- left_join(BNPP1, avg_control_BNPP) %>%
  group_by(subplot, treatment) %>%
  mutate(change_BNPP = (agg_BNPP-mean_controlRain)/mean_controlRain)

summary_change_BNPP <- change_BNPP %>%
  filter(treatment != "controlRain") %>%
  summarize(mean_change_BNPP = mean(change_BNPP),
            se_change_BNPP = calcSE(change_BNPP))
summary_change_BNPP$treatment <- factor(summary_change_BNPP$treatment, levels = c( "springDry", "fallDry","consistentDry"), 
                                        labels = c("Spring\n Dry", "Fall\n Dry", "Consistent\n Dry"))
summary_change_BNPP$subplot <- factor(summary_change_BNPP$subplot, levels = c( "B", "F","G"),
                                      labels = c("Mixed", "Forb", "Grass"))
    
#plot drought treatment effects on ANPP and BNPP 
Fig_drought_ANPP <- ggplot(summary_change_ANPP, aes(x = treatment, y = mean_change_ANPP, col = treatment))+
                          geom_point(size=4)+
                          labs(y=bquote("Relative change ANPP\n compared to control"), color = "Treatment")+
                          facet_grid(~subplot)+
                          geom_hline(aes(yintercept=0), color="grey") +
                          theme(axis.title.x=element_blank(),
                                text = element_text(size=14),
                                panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                panel.background = element_blank(),
                                axis.line = element_line(colour = "black"),
                                panel.border = element_rect(colour = "black", fill = NA, size = 1.2))+
                          geom_errorbar(aes(ymin = mean_change_ANPP-se_change_ANPP, ymax = mean_change_ANPP+se_change_ANPP), size = 0.5, 
                                        width = 0,position=position_dodge(0.9))+
                          scale_x_discrete(labels = c("Spring\n Dry", "Fall\n Dry", "Consistent\n Dry"))+
                          scale_color_manual(name = "Treatment", labels = c("Spring Dry", "Fall Dry","Consistent Dry" ), values= c( "#b2c7e4", "#fccaaf", "#c85b23"))

Fig_drought_BNPP <- ggplot(summary_change_BNPP, aes(x = treatment, y = mean_change_BNPP, col = treatment))+
                          geom_point(size=4)+
                          labs(y=bquote("Relative change BNPP\n compared to control"), color = "Treatment")+
                          facet_grid(~subplot)+
                          geom_hline(aes(yintercept=0), color="grey") +
                          theme(axis.title.x=element_blank(),
                                text = element_text(size=14),
                                panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                panel.background = element_blank(),
                                axis.line = element_line(colour = "black"),
                                panel.border = element_rect(colour = "black", fill = NA, size = 1.2))+
                          geom_errorbar(aes(ymin = mean_change_BNPP-se_change_BNPP, ymax = mean_change_BNPP+se_change_BNPP), size = 0.5, 
                                        width = 0,position=position_dodge(0.9))+
                          scale_x_discrete(labels = c("Spring\n Dry", "Fall\n Dry", "Consistent\n Dry"))+
                          scale_color_manual(name = "Treatment", labels = c("Spring Dry", "Fall Dry","Consistent Dry" ), values= c( "#b2c7e4", "#fccaaf", "#c85b23"))


ggarrange(Fig_drought_ANPP, Fig_drought_BNPP, nrow = 2, common.legend = TRUE, legend = "right")
