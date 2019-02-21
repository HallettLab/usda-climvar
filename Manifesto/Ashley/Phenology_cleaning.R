library(tidyverse)

## Set working directory
setwd("~/Dropbox/ClimVar/DATA/Plant_composition_data")
shelterkey <- read.csv("ANPP/ANPP_EnteredData/Shelter_key.csv")

###############
#2015 Phenology
###############

pheno_2015_1 <- read.csv("Phenology/Phenology_EnteredData/ClimVar_Phenology_20150414.csv") %>% 
  dplyr::select(-Notes, -Recorder) %>% rename(plot=Plot)%>% rename(subplot=Subplot) %>% 
  add_column(time=1)
pheno_2015_1$Percent.Green <- as.numeric(gsub("[^[:alnum:]///' ]", "", pheno_2015_1$Percent.Green))
pheno_2015_1$Percent.Brown <- as.numeric(gsub("[^[:alnum:]///' ]", "", pheno_2015_1$Percent.Brown))
pheno_2015_1 <- left_join(shelterkey, pheno_2015_1)
pheno_2015_2 <- read.csv("Phenology/Phenology_EnteredData/ClimVar_Phenology_20150421.csv")%>% 
  dplyr::select(-Notes, -Recorder) %>% rename(plot=Plot)%>% rename(subplot=Subplot) %>%
  mutate(subplot = recode(subplot, "Control"= "XC", "Compost"="C", "Grass"="G", "Forb"="F", "Both"="B"))%>% 
  add_column(time=2)
pheno_2015_2$Percent.Green <- as.numeric(gsub("[^[:alnum:]///' ]", "", pheno_2015_2$Percent.Green))
pheno_2015_2$Percent.Brown <- as.numeric(gsub("[^[:alnum:]///' ]", "", pheno_2015_2$Percent.Brown))
pheno_2015_2 <- left_join(shelterkey, pheno_2015_2)
pheno_2015_3 <- read.csv("Phenology/Phenology_EnteredData/ClimVar_Phenology_20150428.csv")%>% 
  dplyr::select(-Notes, -Recorder) %>% rename(plot=Plot)%>% rename(subplot=Subplot) %>%
  rename(Date=Date.collected)%>%
  mutate(subplot = recode(subplot, "Control"= "XC", "Compost"="C", "Grass"="G", "Forb"="F", "Both"="B", "Both "="B"))%>% 
  add_column(time=3)
pheno_2015_3$Percent.Green <- as.numeric(gsub("[^[:alnum:]///' ]", "", pheno_2015_3$Percent.Green))
pheno_2015_3$Percent.Brown <- as.numeric(gsub("[^[:alnum:]///' ]", "", pheno_2015_3$Percent.Brown))
pheno_2015_3 <- left_join(shelterkey, pheno_2015_3)
pheno_2015_4 <- read.csv("Phenology/Phenology_EnteredData/ClimVar_Phenology_20150505.csv") %>% 
  dplyr::select(-Notes, -Recorder, -Page, -Line) %>% rename(plot=Plot..1.16.) %>% 
  rename(subplot=Subplot..L..XC..F.G.B..C..D.) %>%
  rename(Date=Date..collected.) %>% add_column(time=4)
pheno_2015_4 <- left_join(shelterkey, pheno_2015_4)

pheno_2015<-rbind(pheno_2015_1, pheno_2015_2, pheno_2015_3, pheno_2015_4)

pheno_2015_xc<-pheno_2015 %>% filter(subplot=="XC") %>%
  group_by(time, treatment) %>%
  summarize(meanPG=mean(Percent.Green), sePG=sd(Percent.Green)/sqrt(length(Percent.Green)))

ggplot(data=pheno_2015_xc, aes(x=time, y=meanPG, color=treatment))+
  geom_point(aes(cex=1.5))+
  geom_errorbar(aes(ymax = meanPG+sePG, ymin = meanPG-sePG), width=.25)+
  geom_line()+
  labs(x="Time", y="Percent Green") +
  theme(text = element_text(size=15))+
  theme_bw()+
  scale_color_manual(values = c( "sienna", "royalblue2","peachpuff","lightsteelblue3"), guide = guide_legend(title = "Treatment"))

###############
#2016 Phenology
###############
pheno_2016_1 <- read.csv("Phenology/Phenology_EnteredData/ClimVar_Phenology_20160427.csv") %>% 
  dplyr::select(-Notes, -Recorder) %>% rename(plot=Plot)%>% rename(subplot=Subplot) %>% 
    rename(Percent.Green=Percent_green) %>% rename(Percent.Brown=Percent_brown) %>% 
    rename(Percent.Bare=Percent_bare) %>% rename(Date=DateTime_start) %>%
    add_column(time=1) %>% filter(subplot!="NA")
pheno_2016_1 <- left_join(shelterkey, pheno_2016_1)
pheno_2016_2 <- read.csv("Phenology/Phenology_EnteredData/ClimVar_Phenology_20160502.csv") %>% 
  dplyr::select(-Notes, -Recorder) %>% rename(plot=Plot)%>% rename(subplot=Subplot) %>%
  rename(Percent.Green=Percent_green) %>% rename(Percent.Brown=Percent_brown) %>% 
  rename(Percent.Bare=Percent_bare) %>% rename(Date=DateTime_start) %>%
  mutate(subplot = recode(subplot, "Control"= "XC", "Compost"="C", "Grass"="G", "Forb"="F", "Both"="B"))%>% 
  add_column(time=2)%>% filter(subplot!="NA")
pheno_2016_2 <- left_join(shelterkey, pheno_2016_2)
pheno_2016_3 <- read.csv("Phenology/Phenology_EnteredData/ClimVar_Phenology_20160509.csv")%>% 
dplyr::select(-Notes, -Recorder) %>% rename(plot=Plot)%>% rename(subplot=Subplot) %>%
  rename(Percent.Green=Percent_green) %>% rename(Percent.Brown=Percent_brown) %>% 
  rename(Percent.Bare=Percent_bare) %>% rename(Date=DateTime_start) %>%
  add_column(time=3)%>% filter(subplot!="NA")
pheno_2016_3 <- left_join(shelterkey, pheno_2016_3)
pheno_2016_4 <- read.csv("Phenology/Phenology_EnteredData/ClimVar_Phenology_20160516.csv") %>% 
  dplyr::select(-Notes, -Recorder) %>% rename(plot=Plot)%>% rename(subplot=Subplot) %>%
  rename(Percent.Green=Percent_green) %>% rename(Percent.Brown=Percent_brown) %>% 
  rename(Percent.Bare=Percent_bare) %>% rename(Date=DateTime_start)  %>% add_column(time=4)%>% 
  mutate(Percent.Green = recode(Percent.Green, "T"= "0"))%>% 
  filter(subplot!="NA")
pheno_2016_4$Percent.Green <- as.numeric(gsub("[^[:alnum:]///' ]", "", pheno_2016_4$Percent.Green))
pheno_2016_4 <- left_join(shelterkey, pheno_2016_4)

pheno_2016<-rbind(pheno_2016_1, pheno_2016_2, pheno_2016_3, pheno_2016_4)

pheno_2016_xc<-pheno_2016 %>% filter(subplot=="XC") %>%
  group_by(time, treatment) %>%
  summarize(meanPG=mean(Percent.Green), sePG=sd(Percent.Green)/sqrt(length(Percent.Green)))

ggplot(data=pheno_2016_xc, aes(x=time, y=meanPG, color=treatment))+
  geom_point(aes(cex=1.5))+
  geom_errorbar(aes(ymax = meanPG+sePG, ymin = meanPG-sePG), width=.25)+
  geom_line()+
  labs(x="Time", y="Percent Green") +
  theme(text = element_text(size=15))+
  theme_bw()+
  scale_color_manual(values = c( "sienna", "royalblue2","peachpuff","lightsteelblue3"), guide = guide_legend(title = "Treatment"))

pheno_all<-rbind(pheno_2015,pheno_2016)

pheno_all_xc<-pheno_all %>% filter(subplot=="XC") %>%
  group_by(time, treatment) %>%
  summarize(meanPG=mean(Percent.Green), sePG=sd(Percent.Green)/sqrt(length(Percent.Green)))

ggplot(data=pheno_all_xc, aes(x=time, y=meanPG, color=treatment))+
  geom_point(aes(cex=1.5))+
  geom_errorbar(aes(ymax = meanPG+sePG, ymin = meanPG-sePG), width=.25)+
  geom_line()+
  labs(x="Time", y="Percent Green") +
  theme(text = element_text(size=15))+
  theme_bw()+
  scale_color_manual(values = c( "sienna", "royalblue2","peachpuff","lightsteelblue3"), guide = guide_legend(title = "Treatment"))


