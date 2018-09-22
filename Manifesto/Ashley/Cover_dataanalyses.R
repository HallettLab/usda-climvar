library(tidyverse)
library(nlme)
library(ggplot2)
library(dplyr)

setwd("~/Dropbox/ClimVar/DATA/Plant_composition_data")
cover<-read.csv("Cover/Cover_CleanedData/ClimVar_species-cover.csv")
levels(cover$plot)
levels(cover$species_name)
levels(cover$subplot)
levels(cover$year)
levels(cover$genus)
levels(cover$species)
levels(cover$func)
levels(cover$status)
levels(cover$treatment)
levels(cover$shelterBlock)
levels(cover$shelter)
levels(cover$func2)

#change plots, years, shelter to factors
cover[,'plot'] <- as.factor(as.character(cover[,'plot']))
cover[,'year'] <- as.factor(as.character(cover[,'year']))
cover[,'shelter'] <- as.factor(as.character(cover[,'shelter']))

#remove compost treatment
cover_noC<-filter(cover, subplot!='C')

####################################################
##rerun following block of code from Lauren's old code:
vegAv <- cover_noC %>%
  group_by(plot, subplot, treatment, shelterBlock) %>%
  mutate(totcover = sum(cover)) %>%
  mutate(AvCover = cover/totcover * 100) %>%
  filter(species_name == "Avena barbata") %>%
  dplyr::select(plot, subplot, AvCover) 

##SPECIES WITH MAX ABUNDANCE >5
vegabund <- cover_noC %>%
  group_by (plot, subplot, func, genus, treatment, shelterBlock) %>%
  summarize(cover=sum(cover)) %>%
  tbl_df() %>%
  group_by(genus) %>%
  mutate(maxcover=max(cover)) %>%
  tbl_df() %>%
  filter(maxcover>5)

ggplot(subset(vegabund, func=="grass" & (subplot=="XC" | subplot=="G" | subplot=="B")), aes(x=treatment, y=cover, fill=genus)) + 
  geom_bar(stat="identity", position="dodge") + 
  facet_grid(shelterBlock~subplot) +
  theme_bw() + 
  labs(x="Treatment", y="Percent cover", fill="Functional
       group") +
  theme(text = element_text(size=20))

ggplot(subset(vegabund, func=="forb" & (subplot=="XC" | subplot=="F" | subplot=="B")), aes(x=treatment, y=cover, fill=genus)) + 
  geom_bar(stat="identity", position="dodge") + 
  facet_grid(shelterBlock~subplot) +
  theme_bw() + 
  labs(x="Treatment", y="Percent cover", fill="Functional
       group") +
  theme(text = element_text(size=20))

###ALL GRASSES AND FORBS
gf<-cover_noC %>%
  group_by(plot, subplot, func, treatment, shelterBlock) %>%
  summarize(cover=sum(cover)) %>%
  tbl_df()

gfproportion <- gf %>%
  group_by(plot, subplot, treatment, shelterBlock) %>%
  mutate(totcover = sum(cover)) %>%
  mutate(percentForb = (cover/totcover)*100) 
  filter(func == "forb") %>%
  dplyr::select(-func)

gf_graphic <- gf %>%
  group_by(subplot, func, treatment) %>%
  summarize(meancover=mean(cover), secover=sd(cover)/sqrt(length(cover)))

ggplot(gf_graphic, aes(x=treatment, y=meancover, fill=func)) + 
  geom_bar(stat="identity", position="dodge") + 
  facet_wrap(~subplot) +
  theme_bw() + 
  geom_errorbar(aes(ymax = meancover+secover, ymin = meancover-secover), width=.25, position=position_dodge(width=0.9)) + 
  labs(x="Treatment", y="Percent cover", fill="Functional
       group") +
  theme(text = element_text(size=20))

ggplot((gfproportion), aes(x=treatment, y = percentForb)) + 
  geom_boxplot()

a<-lme(percentForb ~ treatment, random=~1|shelterBlock/subplot, data=gfproportion,
         contrasts=list(treatment=contr.treatment))

summary(a)        
summary(glht(a,linfct=mcp(treatment="Tukey")), alternative="Bonferonni")

gfproportion_graph <- gfproportion %>%
  group_by(treatment) %>%
  summarize(meanprop=mean(percentForb), seprop=sd(percentForb)/sqrt(length(percentForb)))

ggplot(gfproportion_graph, aes(x=treatment, y=meanprop)) + 
  geom_bar(stat="identity", position="dodge", fill = "gray") + 
  theme_bw() + 
  geom_errorbar(aes(ymax = meanprop+seprop, ymin = meanprop-seprop), width=.25, position=position_dodge(width=0.9)) + 
  labs(x="Treatment", y="Percent forbs") +
  theme(text = element_text(size=20))

ggplot(gfproportion, aes(x=percentForb, y=totcover, color = treatment)) + 
  geom_point() + geom_smooth(method = "lm")
############################################

#calculate proportion grass
gfproportionG <- gf %>%
  group_by(plot, subplot, treatment, shelterBlock) %>%
  mutate(totcover = sum(cover)) %>%
  mutate(percentGrass = (cover/totcover)*100) %>%
  filter(func == "grass") %>%
  dplyr::select(-func) %>%
  dplyr::select(-cover)

#calc prop forb
gfproportionF <- gf %>%
  group_by(plot, subplot, treatment, shelterBlock) %>%
  mutate(totcover = sum(cover)) %>%
  mutate(percentForb = (cover/totcover)*100) %>%
  filter(func == "forb") %>%
  dplyr::select(-func) %>%
  dplyr::select(-cover)

gf2<-cover_noC %>%
  group_by(plot, subplot, species_name, func, treatment, shelterBlock) %>%
  summarize(cover=sum(cover)) %>%
  tbl_df()

#proportion avena
vegAv2 <- gf2 %>%
  group_by(plot, subplot, treatment, shelterBlock) %>%
  mutate(totcover = sum(cover)) %>%
  mutate(AvCover = cover/totcover * 100) %>%
  filter(species_name == "Avena barbata") %>%
  dplyr::select(-func,-species_name,-cover)


# merge prop grass with forb prop and prop avena
gfprop_all<- merge(gfproportionG, gfproportionF, all.x = T) %>%
  tbl_df() 

gfprop_all<-merge(gfprop_all, vegAv2, all.x = T) %>%
  tbl_df() 

#create a stacked bar plot
#ggplot(gfprop_all, aes(fill=condition, y=value, x=specie)) +
  #geom_bar( stat="identity")

#does proportion of Avena drive community structure?
#is cover specific to site (XC and control rain)?
#does functional diversity/evenness/richness change with precipitation variability (drought vs. control) or seasonability (drought treatments)?
#how does functional diversity relate to coefficient of variation for soil moisture?
