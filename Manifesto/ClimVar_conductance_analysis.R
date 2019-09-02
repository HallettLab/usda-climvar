# simple analysis of climvar 2015 + 2016 leaf poromter data
# created: jun 2015 (lmh); updated script may 2019 (ctw)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)

datpath <- "~/Dropbox/ClimVar/DATA/Plant_composition_data/Leaf_Porometer/"
dat<-read.csv(paste0(datpath, "Leaf_Porometer_CleanedData/LeafPorometer_2015-2016_clean.csv")) %>%
  tbl_df() 

boxplot(cond_mmol_m2s~treatment*species, data=dat)

myaov<-aov(cond_mmol_m2s~treatment*species + Error(shelterBlock/rep), data=dat)
summary(myaov)

myaov2<-aov(cond_mmol_m2s~treatment*species*year + Error(shelterBlock/rep), data=dat)
summary(myaov2)

dat_graph<-dat %>%
  group_by(year, treatment, species, plot) %>%
  summarize(cond=mean(cond_mmol_m2s)) %>%
  tbl_df()%>%
  group_by(year, treatment, species) %>%
  summarize(meancond=mean(cond), second=sd(cond)/sqrt(length(cond))) %>%
  tbl_df() %>%
  mutate(treatment=ordered(treatment, levels = c(consistentDry="consistentDry", springDry="springDry", fallDry="fallDry", controlRain="controlRain")))

cond_graph<-ggplot(dat_graph, aes(x=treatment, y=meancond, group=species, fill=species)) + 
  geom_bar(stat="identity", position="dodge") + theme_bw() + 
  geom_errorbar(aes(ymax = meancond+second, ymin = meancond-second), width=.25, position=position_dodge(width=0.9)) + 
  labs(x="Treatment", y=(expression("Stomatal conductance (mmol/(m"^{2}~ "*s))")), fill="Species") +
  theme(text = element_text(size=22),
        # move legend to white space in first panel to give more space to data
        legend.position = c(.12, .9),
        legend.background = element_rect(color = "black")) +
  facet_grid(~year)
  
tiff(paste0(datpath, "Leaf_Porometer_Figures/ClimVar_conductance_2015-2016.tiff"), width=900, height=600)
cond_graph
dev.off()
