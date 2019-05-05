library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)

mydat<-read.csv("ClimVar_porometer_2015.csv") %>%
  tbl_df()  %>%
  mutate(plot=extract_numeric(plot)) 

shelterkey<-read.csv("Shelter_key.csv")

dat<-merge(mydat, shelterkey) %>%
  tbl_df()

boxplot(cond~treatment*species, data=dat)

myaov<-aov(cond~treatment*species + Error(shelterBlock/replicate), data=dat)
summary(myaov)

dat_graph<-dat %>%
  group_by(treatment, species, plot) %>%
  summarize(cond=mean(cond)) %>%
  tbl_df()%>%
  group_by(treatment, species) %>%
  summarize(meancond=mean(cond), second=sd(cond)/sqrt(length(cond))) %>%
  tbl_df() %>%
  mutate(treatment=ordered(treatment, levels = c(consistentDry="consistentDry", springDry="springDry", fallDry="fallDry", controlRain="controlRain"))) %>%
  mutate(species2="Avena", species2=ifelse(species=="ER", "Erodium", species2))

cond_graph<-ggplot(dat_graph, aes(x=treatment, y=meancond, group=species2, fill=species2)) + 
  geom_bar(stat="identity", position="dodge") + theme_bw() + 
  geom_errorbar(aes(ymax = meancond+second, ymin = meancond-second), width=.25, position=position_dodge(width=0.9)) + 
  labs(x="Treatment", y=(expression("Stomatal conductance (mmol/(m"^{2}~ "*s))")), fill="Species") +
  theme(text = element_text(size=25))
  
tiff("ClimVar_conductance.tiff", width=900, height=600)
cond_graph
dev.off()
