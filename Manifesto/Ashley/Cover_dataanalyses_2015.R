library(tidyverse)
library(nlme)
library(ggplot2)
library(dplyr)
library(multcomp)
library(lsmeans)
library(MuMIn)
library(vegan)
library(RColorBrewer)

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

#remove compost treatment and select 2015
cover_noC_2015<-filter(cover, subplot!='C', year == '2015')

####################################################
##rerun following block of code from Lauren's old code:
vegAv_2015 <- cover_noC_2015 %>%
  group_by(plot, subplot, treatment, shelterBlock) %>%
  mutate(totcover = sum(cover)) %>%
  mutate(AvCover = cover/totcover * 100) %>%
  filter(species_name == "Avena barbata") %>%
  dplyr::select(plot, subplot, AvCover) 

##SPECIES WITH MAX ABUNDANCE >5
vegabund_2015 <- cover_noC_2015 %>%
  group_by (plot, subplot, func, genus, treatment, shelterBlock) %>%
  summarize(cover=sum(cover)) %>%
  tbl_df() %>%
  group_by(genus) %>%
  mutate(maxcover=max(cover)) %>%
  tbl_df() %>%
  filter(maxcover>5)

ggplot(subset(vegabund_2015, func=="grass" & (subplot=="XC" | subplot=="G" | subplot=="B")), aes(x=treatment, y=cover, fill=genus)) + 
  geom_bar(stat="identity", position="dodge") + 
  facet_grid(shelterBlock~subplot) +
  theme_bw() + 
  labs(x="Treatment", y="Percent cover", fill="Genus") +
  theme(text = element_text(size=20))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(subset(vegabund_2015, func=="forb" & (subplot=="XC" | subplot=="F" | subplot=="B")), aes(x=treatment, y=cover, fill=genus)) + 
  geom_bar(stat="identity", position="dodge") + 
  facet_grid(shelterBlock~subplot) +
  theme_bw() + 
  labs(x="Treatment", y="Percent cover", fill="Genus") +
  theme(text = element_text(size=20))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

###ALL GRASSES AND FORBS
gf_2015<-cover_noC_2015 %>%
  group_by(plot, subplot, func, treatment, shelterBlock, year) %>%
  summarize(cover=sum(cover)) %>%
  tbl_df()

gfproportion_2015 <- gf_2015 %>%
  group_by(plot, subplot, treatment, shelterBlock, year) %>%
  mutate(totcover = sum(cover)) %>%
  mutate(percentForb = (cover/totcover)*100) %>%
  filter(func == "forb") %>%
  dplyr::select(-func)

gf_graphic_2015 <- gf_2015 %>%
  group_by(subplot, func, treatment) %>%
  summarize(meancover=mean(cover), secover=sd(cover)/sqrt(length(cover)))

ggplot(gf_graphic_2015, aes(x=treatment, y=meancover, fill=func)) + 
  geom_bar(stat="identity", position="dodge") + 
  facet_wrap(~subplot) +
  theme_bw() + 
  geom_errorbar(aes(ymax = meancover+secover, ymin = meancover-secover), width=.25, position=position_dodge(width=0.9)) + 
  labs(x="Treatment", y="Percent cover", fill="Functional
       group") +
  theme(text = element_text(size=20))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

gfproportion_2015_noXC<-gfproportion_2015 %>% filter(subplot!="XC")

ggplot((gfproportion_2015_noXC), aes(x=treatment, y = percentForb, fill=subplot)) + 
  geom_boxplot()

a<-lme(percentForb ~ treatment, random=~1|shelterBlock/subplot, data=gfproportion_2015_noXC,
       contrasts=list(treatment=contr.treatment))

summary(a) 
anova(a)
r.squaredGLMM(a)
summary(glht(a,linfct=mcp(treatment="Tukey")), alternative="Bonferonni")


gfproportion_graph_2015_noXC <- gfproportion_2015 %>%
  group_by(treatment) %>% filter(subplot!="XC")%>%
  summarize(meanprop=mean(percentForb), seprop=sd(percentForb)/sqrt(length(percentForb)))

ggplot(gfproportion_graph_2015, aes(x=treatment, y=meanprop)) + 
  geom_bar(stat="identity", position="dodge", fill = "gray") + 
  theme_bw() + 
  geom_errorbar(aes(ymax = meanprop+seprop, ymin = meanprop-seprop), width=.25, position=position_dodge(width=0.9)) + 
  labs(x="Treatment", y="Percent forbs") +
  theme(text = element_text(size=20))

ggplot(gfproportion_2015_noXC, aes(x=percentForb, y=totcover, color = treatment))+
  geom_point() + geom_smooth(method = "lm")

ggplot(gfproportion_2015, aes(x=percentForb, y=totcover, color = treatment))+
  facet_grid(~subplot)+
  geom_point() + geom_smooth(method = "lm")
############################################

#calculate proportion grass for composition treatments
gfproportionG_2015 <- gf_2015 %>%
  group_by(plot, subplot, treatment, shelterBlock, year) %>%
  mutate(totcover = sum(cover)) %>%
  mutate(percentGrass = (cover/totcover)*100) %>%
  filter(func == "grass", subplot!= "XC") %>%
  dplyr::select(-func)

b<-lme(percentGrass ~ treatment, random=~1|shelterBlock/subplot, data=gfproportionG_2015,
       contrasts=list(treatment=contr.treatment))

summary(b) 
anova(b)
r.squaredGLMM(b) #1% explained by fixed effects
summary(glht(b,linfct=mcp(treatment="Tukey")), alternative="Bonferonni")

ggplot((gfproportionG_2015), aes(x=treatment, y = percentGrass, fill=subplot)) + 
  geom_boxplot()

#calc prop forb
gfproportionF_2015 <- gf_2015 %>%
  group_by(plot, subplot, treatment, shelterBlock, year) %>%
  mutate(totcover = sum(cover)) %>%
  mutate(percentForb = (cover/totcover)*100) %>%
  filter(func == "forb") %>%
  dplyr::select(-func) %>%
  dplyr::select(-cover)

gf2_2015<-cover_noC_2015 %>%
  group_by(plot, year, subplot, species_name, func, treatment, shelterBlock) %>%
  mutate(cover=sum(cover)) %>%
  tbl_df()


#proportion avena
vegAv2_2015 <- gf2_2015 %>%
  group_by(plot, subplot, treatment, shelterBlock, year) %>%
  mutate(totcover = sum(cover)) %>%
  mutate(AvCover = cover/totcover * 100) %>%
  filter(species_name == "Avena barbata") %>%
  dplyr::select(-func,-species_name,-cover, -genus, -species, -status, -func2, -X)

#proportion lolium
vegLol_2015 <- gf2_2015 %>%
  group_by(plot, subplot, treatment, shelterBlock, year)%>%
  mutate(totcover=sum(cover))%>%
  mutate(LolCover=cover/totcover*100) %>%
  filter(species_name == "Lolium multiflorum") %>%
  dplyr::select(-func,-species_name,-cover, -genus, -species, -status, -func2, -totcover, -X)

#proportion medusahead
vegTae_2015 <- gf2_2015 %>%
  group_by(plot, subplot, treatment, shelterBlock, year)%>%
  mutate(totcover=sum(cover))%>%
  mutate(TaeCover=cover/totcover*100) %>%
  filter(genus == "Taeniatherum") %>%
  dplyr::select(-func,-species_name,-cover, -genus, -species, -status, -func2, -totcover, -X)

# merge prop grass with forb prop and prop avena
gfprop_all_2015<- merge(gfproportionG_2015, gfproportionF_2015, all.x = T) %>%
  tbl_df() 

gfprop_all_2015<-merge(gfprop_all_2015, vegAv2_2015, all.x = T) %>%
  tbl_df() 

gfprop_all_2015 <-merge(gfprop_all_2015,vegLol_2015, all.x=T)%>%
  tbl_df() 

gfprop_all_2015 <-merge(gfprop_all_2015,vegTae_2015, all.x=T)%>% 
  tbl_df()

SE <- function(x){(sd(x)/sqrt(length(x)))}
#MeanCov_2015<-gfprop_all_2015 %>% group_by(treatment)%>%summarise(meanFCover=mean(percentForb))
#CovFse_2015<-aggregate(percentForb ~ treatment, data= gfprop_all_2015, FUN = SE)
#colnames(CovFse_2015)[colnames(CovFse_2015)=="percentForb"] <- "Fcov.SE"
#from ANNP_dataanalyses.R
#MeanCov2<- merge(MeanFunc, MeanCov)
#MeanCov2<-merge(MeanCov2,CovFse)

#ggplot(MeanCov2, aes(x=sm_cv, y=meanFCover, color=treatment))+
  #geom_errorbar(aes(ymin=meanFCover-Fcov.SE, ymax=meanFCover+Fcov.SE), width=1, position=pd)+
  #geom_point(position=pd)+
  #geom_smooth(method='lm')+
  #labs(x="Soil Moisture CV", y="Percent Forb Cover")

gf3_2015<-cover_noC_2015 %>%
  group_by(plot, year, subplot, treatment, shelterBlock) %>%
  mutate(totcover=sum(cover)) %>%
  tbl_df()


gf4_2015 <- gf3_2015 %>%
  group_by(year, subplot, treatment, shelterBlock, genus, species_name, func2) %>%
  summarize(percentCov = cover/totcover)%>%
  filter(subplot=="XC")%>%
  #filter(percentCov>0.10)%>%
  tbl_df()

gf5_2015 <- gf3_2015 %>%
  group_by(year, subplot, treatment, shelterBlock, genus, species_name, func2) %>%
  summarize(percentCov = cover/totcover)%>%
  #filter(percentCov>0.10)%>%
  tbl_df()

gf2_graph_2015 <- gf4_2015 %>%
  group_by(treatment,shelterBlock,genus,func2)%>%
  summarise(cover=mean(percentCov), secover=sd(percentCov)/sqrt(length(percentCov)))
  
gf3_graph_2015 <- gf4_2015 %>%
  group_by(treatment,genus, func2) %>%
  summarise(cover=mean(percentCov), secover=sd(percentCov)/sqrt(length(percentCov)))

gf5_graph_2015 <- gf5_2015 %>%
  group_by(treatment, subplot, genus, func2)%>%
  filter(subplot!="XC")%>%
  summarise(cover=mean(percentCov), secover=sd(percentCov)/sqrt(length(percentCov)))

gf6_graph_2015 <- gf5_2015 %>%
  group_by(treatment, subplot, genus, func2, shelterBlock)%>%
  filter(subplot!="XC")%>%
  summarise(cover=mean(percentCov), secover=sd(percentCov)/sqrt(length(percentCov)))

#create a stacked bar plot for XC
ggplot(gf3_graph_2015, aes(fill=genus, y=cover, x=treatment)) +
  geom_bar( stat="identity")

gf2_graph_2015 <- gf2_graph_2015 %>% arrange(func2, genus) %>% filter(cover>0.01, genus!="NA")
gf2_graph_2015$func2 <- factor(gf2_graph_2015$func2, c("forb","Nfixer","grass"))
gf2_graph_2015$genus<- factor(gf2_graph_2015$genus, c("Centaurea", "Convolvulus", "Erodium","Hypochaeris","Senecio", "Vicia","Avena","Bromus","Cynodon","Hordeum","Lolium", "Taeniatherum"))
ggplot(gf2_graph_2015, aes(fill=genus, y=cover, x=treatment, color=func2)) +
  #geom_bar( stat="identity")+
  theme_bw()+
  scale_fill_manual(values = c("orange", "orangered", "firebrick","indianred4","indianred", "green4", "lightblue", "skyblue2", "skyblue4", "dodgerblue3", "royalblue3","navy"))+
  geom_bar( stat="identity", position='stack')+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  facet_wrap(~shelterBlock)

ggplot(gf2_graph_2015, aes(fill=genus, y=cover, x=shelterBlock)) +
  geom_bar( stat="identity")+
  theme_bw()+
  scale_fill_manual(values = c("orange", "orangered", "firebrick","indianred4","indianred", "green4", "lightblue", "skyblue2", "skyblue4", "dodgerblue3", "royalblue3","navy"))+
  geom_bar( stat="identity", position='stack')+
  theme(axis.text.x = element_text(hjust = 0.5))
  

#now for control only using genera with > 1% cover for XC
gf3_graph_2015 <- gf3_graph_2015 %>% arrange(func2, genus) %>% filter(cover>0.01)
levels(gf3_graph_2015$genus)
gf3_graph_2015$func2 <- factor(gf3_graph_2015$func2, c("forb","Nfixer","grass"))
gf3_graph_2015$genus<- factor(gf3_graph_2015$genus, c("Centaurea", "Convolvulus", "Erodium","Hypochaeris","Senecio", "Vicia","Avena","Bromus","Cynodon","Hordeum","Lolium", "Taeniatherum"))
ggplot(gf3_graph_2015, aes(fill=genus, colour=func2,  y=cover, x=treatment:func2)) +
  theme_bw()+
  scale_fill_manual(values = c("orange", "orangered", "firebrick","indianred4","indianred", "green4", "lightblue", "skyblue2", "skyblue4", "dodgerblue3", "royalblue3","navy"))+
  geom_bar( stat="identity", position='stack')+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#lets look at cover for the composition treatments
gf5_graph_2015 <- gf5_graph_2015 %>% arrange(func2, genus) %>% filter(cover>0.01) %>% filter(genus!="NA")
levels(gf5_graph_2015$genus)
gf5_graph_2015$func2 <- factor(gf5_graph_2015$func2, c("forb","Nfixer","grass"))
gf5_graph_2015$genus<- factor(gf5_graph_2015$genus, c("Anagalis", "Centaurea","Cerastium", "Erodium","Hypochaeris", "Rumex","Sherardia","Silene", "Trifolium","Vicia","Avena","Bromus","Cynodon","Hordeum","Lolium", "Taeniatherum", "Vulpia"))
ggplot(gf5_graph_2015, aes(fill=genus, colour=func2,  y=cover, x=treatment)) +
  theme_bw()+
  facet_wrap(~subplot)+
  scale_fill_manual(values = c("goldenrod","orange", "darkorange2", "firebrick","indianred4","indianred", "saddlebrown", "palegreen", "green4", "lightblue", "skyblue2", "skyblue4", "dodgerblue3", "royalblue3","navy", "black"))+
  geom_bar( stat="identity", position='stack')+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

gf6_graph_2015 <- gf6_graph_2015 %>% arrange(func2, genus) %>% filter(cover>0.01) %>% filter(genus!="NA", genus!="Unknown")
levels(gf6_graph_2015$genus)
gf6_graph_2015$func2 <- factor(gf6_graph_2015$func2, c("forb","Nfixer","grass"))
gf6_graph_2015$genus<- factor(gf6_graph_2015$genus, c("Anagalis", "Carduus","Centaurea","Cerastium", "Erodium","Hypochaeris","Lactuca", "Lythrum", "Rumex","Sonchrum","Sherardia","Silene","Zeltnera", "Trifolium","Vicia","Avena","Brachypodium","Bromus","Cynodon","Hordeum","Lolium", "Taeniatherum", "Vulpia"))
ggplot(gf6_graph_2015, aes(fill=genus, colour=func2,  y=cover, x=treatment)) +
  theme_bw()+
  facet_wrap(~subplot*shelterBlock)+
  scale_fill_manual(values = c("goldenrod","orange", "orangered", "darkorange2", "firebrick","indianred4","indianred2","brown4", "indianred", "saddlebrown", "khaki3", "palegreen", "green4", "lightblue", "lightblue4", "skyblue2", "skyblue4", "dodgerblue3", "royalblue3","navy", "black"))+
  geom_bar( stat="identity", position='stack')+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


library(ggpubr)
#does proportion of Avena drive community structure?
#remove XC for comparison to Lina's data
gfprop_noXC_2015<-gfprop_all_2015%>% filter(subplot!="XC")
ggplot(gfprop_noXC_2015, aes(x=AvCover, y=totcover, color=treatment))+
  geom_point()+
  geom_smooth(method='lm')+
  labs(x="Percent Avena", y="Total Cover")

ggplot(gfprop_noXC_2015, aes(x=AvCover, y=percentGrass, color=treatment))+
  geom_point()+
  geom_smooth(method='lm' )+
  labs(x="Percent Avena", y="Percent Grass Cover")

ggplot(gfprop_noXC_2015, aes(x=AvCover, y=LolCover, color=treatment, shape=subplot))+
  geom_point()+
  #geom_smooth(method='lm' )+
  #facet_wrap(~treatment)+
  labs(x="Percent Avena", y="Percent Lolium")

gfprop_noF_2015<-gfprop_noXC_2015 %>% filter(subplot!="F")
gfprop_noF_season_2015 <- gfprop_noXC_2015 %>% filter(treatment=="fallDry","springDry")

ggscatterhist(gfprop_noF_2015, x = "AvCover", y = "LolCover",
              color = "treatment", size = 3, alpha = 0.6,
              margin.params = list(fill = "treatment", color = "black", size = 0.2))

ggscatterhist(gfprop_noF_2015, x = "AvCover", y = "TaeCover",
              color = "treatment", size = 3, alpha = 0.6,
              margin.params = list(fill = "treatment", color = "black", size = 0.2))

ggplot(gfprop_noXC_2015, aes(x=subplot, y=percentForb, color=treatment))+
  geom_boxplot()

gfprop_graph_2015<- gfprop_all_2015 %>% gather(func, percent, 9:13)%>%
  filter(subplot!="XC", func!="percentForb", func!="percentGrass")%>%
  tbl_df()

ggplot(gfprop_graph_2015, aes(x=treatment, y=percent, fill=func))+
  facet_wrap(~subplot)+
  geom_boxplot()

ggplot(gfprop_noXC_2015, aes(x=AvCover, y=totcover, color=treatment))+
  #facet_grid(~treatment)+
  geom_point()+
  geom_smooth(method='lm', se=FALSE)

ggplot(gfprop_noXC_2015, aes(x=TaeCover, y=totcover, color=treatment))+
  #facet_grid(~treatment)+
  geom_point()+
  geom_smooth(method='lm', se=FALSE)

gfprop_graph2_2015<-gfprop_all_2015 %>% gather(func, percent, 11:13) %>% filter(subplot=="B")
ggplot(gfprop_graph2_2015, aes(x=percent, y=totcover, group=func, color=func))+
  facet_grid(~treatment)+
  geom_point()+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method='lm', se=FALSE)

ggplot(gfprop_noXC_2015, aes(x=treatment, y=AvCover))+
  geom_boxplot()

#does treatment affect % cover for Lolium?
ggplot(gfprop_noXC_2015, aes(x=treatment, y=LolCover))+
  geom_boxplot()

g7<-lme(LolCover ~ treatment, random=~1|subplot/shelterBlock, gfprop_noXC_2015, na.action=na.exclude)
summary(g7)
anova(g7) 
r.squaredGLMM(g7) #12% of variation explained by fixed effects, 30% by whole model (spacial variation?)
qqnorm(residuals(g7))
qqline(residuals(g7))
shapiro.test(residuals(g7))
#not normally distributed, try log transformation
g8<-lme(log(LolCover+1) ~treatment, random=~1|subplot/shelterBlock, gfprop_noXC_2015, na.action=na.exclude)
summary(g8)
anova(g8)
r.squaredGLMM(g8) #14% of variation explained by fixed effects, 32% by whole model (spatial variation?)
qqnorm(residuals(g8))
qqline(residuals(g8))
shapiro.test(residuals(g8))
#normal
gLS8<-lsmeans(g8, ~treatment, adjust="Tukey")
contrast(gLS8, "pairwise") #control rain is greater than spring dry

#does treatment affect medusahead?
g5w<-lme(TaeCover ~ treatment, random=~1|subplot/shelterBlock, gfprop_noXC_2015, na.action=na.exclude)
summary(g5w)
anova(g5w)
r.squaredGLMM(g5w) #14% of variation explained by fixed effects, 25% by whole model (interannual variation?)
qqnorm(residuals(g5w))
qqline(residuals(g5w))
shapiro.test(residuals(g5w))
#not normally distributed, try log transformation
g6w<-lme(log(TaeCover+1) ~treatment, random=~1|shelterBlock/subplot, gfprop_noXC_2015, na.action=na.exclude)
summary(g6w)
anova(g6w) #main effect of treatment; p=0.0031
r.squaredGLMM(g6w) #19% of variation explained by fixed effects, 47% by whole model (spatial or subplot variation?)
qqnorm(residuals(g6w))
qqline(residuals(g6w))
shapiro.test(residuals(g6w))
#normal
gLS6w<-lsmeans(g6w, ~treatment, adjust="tukey")
contrast(gLS6w, "pairwise") #differences b/w control rain and fall dry; spring dry and fall dry
#more medusahead in fall dry
#treatment is signifiant, more medusahead in fall dry
ggplot(gfprop_noXC_2015, aes(x=treatment, y=TaeCover))+
  geom_boxplot()
#medusahead released from competition in fall dry?

#does treatment affect avena
av<-lme(AvCover ~ treatment, random=~1|shelterBlock/subplot, gfprop_noXC_2015, na.action=na.exclude)
summary(av)
anova(av)
r.squaredGLMM(av) #2% of variation explained by fixed effects, 26% by whole model (interannual variation?)
qqnorm(residuals(av))
qqline(residuals(av))
shapiro.test(residuals(av))
#not normally distributed, try log transformation
av2<-lme(log(AvCover+1) ~treatment, random=~1|shelterBlock/subplot, gfprop_noXC_2015, na.action=na.exclude)
summary(av2)
anova(av2)
r.squaredGLMM(av2) #2% of variation explained by fixed effects, 34% by whole model (interannual variation?)
qqnorm(residuals(av2))
qqline(residuals(av2))
shapiro.test(residuals(av2))
#normal
LSav<-lsmeans(av2, ~treatment, adjust="tukey")
contrast(LSav, "pairwise") #no differences for avena
ggplot(gfprop_noXC_2015, aes(x=treatment, y=AvCover))+
  geom_boxplot()

#does treatment affect total % cover ?
tot<-lme(totcover ~ treatment, random=~1|shelterBlock/subplot, gfprop_noXC_2015, na.action=na.exclude)
summary(tot)
anova(tot)
r.squaredGLMM(tot) #4% of variation explained by fixed effects, 37% by whole model (interannual variation?)
qqnorm(residuals(tot))
qqline(residuals(tot))
shapiro.test(residuals(tot))
#normal
LStot<-lsmeans(tot, ~treatment, adjust="tukey")
contrast(LStot, "pairwise") #no differences for avena
ggplot(gfprop_noXC_2015, aes(x=treatment, y=totcover))+
  geom_boxplot()

gfprop_noF_season_2015 <- gfprop_graph_2015 %>% filter(treatment!="consistentDry") %>% filter(subplot!="F", subplot!="XC")
gfprop_noF_season_graph_2015 <- gfprop_noF_season_2015 %>%
  group_by(func, treatment) %>%
  summarize(meancover=mean(percent), secover=sd(percent)/sqrt(length(percent)))
gfprop_noF_season_graph_2015<-gfprop_noF_season_graph_2015%>%filter(func!="percentForb")

ggplot(gfprop_noF_season_graph_2015, aes(x=treatment, y=meancover, color=func, group=func))+
  geom_line()+
  geom_point()+
  theme_bw()+
  geom_errorbar(aes(ymin=meancover-secover, ymax=meancover+secover), width=.2,position=position_dodge(0.05))+
  labs(x="Treatment", y="Mean Cover (%)")

#make plots to match Lina's BNPP plots
se <- function(x) sqrt(var(x)/length(x)) #create a function for SE
gf_prop_all_long <- gfprop_noXC_2015 %>% gather(func, cover, -plot, -subplot, -treatment,-shelterBlock, -year,-shelter)
#note: to include XC on these plots use object: gfprop_all_2015
summary_cover_shelter <- gf_prop_all_long %>%
  filter(treatment == "consistentDry" | treatment == "controlRain", func!="cover") %>%
  group_by(shelter, func) %>% #group by shelter and functional groups
  summarise(mean = mean(cover), #summarise by mean and SE
            SE = se(cover))

summary_cover_fall <- gf_prop_all_long %>%
  filter(treatment == "controlRain" | treatment == "fallDry", func!="cover") %>%
  group_by(treatment, func) %>% #group by fall rain treatment and functional groups
  summarise(mean = mean(cover), #summarise by mean and SE
            SE = se(cover))

summary_cover_spring <- gf_prop_all_long %>%
  filter(treatment == "controlRain" | treatment == "springDry", func!="cover") %>%
  group_by(treatment, func) %>% #gropu by spring rain treatment and functional groups
  summarise(mean = mean(cover), #summarise by mean and SE
            SE = se(cover))


#interaction plot of shelter treatment and functional groups
cover_shelter <- ggplot(summary_cover_shelter, aes(x = as.factor(shelter), y = mean, group = func, color = func)) +
  geom_line(size = 1.0) + #add lines
  geom_point(size = 2.0) + #add points
  geom_errorbar(aes(ymin = mean-SE, ymax = mean+SE), width = 0.1, size = 1.0) + #add error bars
  theme_classic() + #simple, classic theme 
  labs(x = "Shelter", y = expression(paste("% Cover"))) + #label the axis
  scale_color_manual(values = c("#999999","#E69F00","#56B4E9", "green", "Red", "purple")) + #specify color 
  theme(legend.position = "none") + #remove the legend 
  scale_y_continuous(limits = c(0, 150)) 
 # annotate("text", x= 1.5, y = 535, label = "Mixed", color = "#999999", angle = -40) +
 # annotate("text", x= 1.35, y = 450, label = "Grass", color = "#56B4E9", angle = 40) +
  #annotate("text", x= 1.5, y = 485, label = "Control", color = "Red", angle = -10) +
  #annotate("text", x= 1.5, y = 358, label = "Forb", color = "#E69F00", angle = 18)
cover_shelter

#interaction plot of fall rain treatment and functional groups
cover_fall <- 
  ggplot(summary_cover_fall, aes(x = treatment, y = mean, group = func, color = func)) +
  geom_line(size = 1.0) + #add lines
  geom_point(size = 2.0) + #add points
  geom_errorbar(aes(ymin = mean-SE, ymax = mean+SE), width = 0.1, size = 1.0) + #add error bars
  theme_classic() + #simple, classic theme 
  labs(x = "Fall Rain") + #label the axis
  scale_color_manual(values = c("#999999","#E69F00","#56B4E9", "green", "Red", "purple")) + #specify color
  theme(legend.position = "none") + #remove the legend
  scale_y_continuous(labels = NULL, name = NULL, limits = c(0,150))  #remove y-axis label
  #scale_x_discrete(limits=c(1,0)) + #change order of discrete x scale 
  #annotate("text", x= 1.5, y = 560, label = "Mixed", color = "#999999", angle = -10) +
  #annotate("text", x= 1.5, y = 445, label = "Grass", color = "#56B4E9", angle = 28) +
  #annotate("text", x= 1.5, y = 515, label = "Control", color = "Red", angle = -31) +
  #annotate("text", x= 1.5, y = 325, label = "Forb", color = "#E69F00", angle = -5)
cover_fall

#interaction plot of spring rain treatment and functional groups
cover_spring <- ggplot(summary_cover_spring, aes(x = treatment, y = mean, group = func, color = func)) +
  geom_line(size = 1.0) + #add lines
  geom_point(size = 2.0) + #add points
  geom_errorbar(aes(ymin = mean-SE, ymax = mean+SE), width = 0.1, size = 1.0) + #add error bars
  theme_classic() + #simple, classic theme 
  labs(x = "Spring Rain", y = expression(paste("Aggregated BNPP (g/m"^2,") in soil depth 0-30 cm"))) + #label the axis
  scale_color_manual(values = c("#999999","#E69F00","#56B4E9", "green", "Red", "purple"), labels = c("Mixed", "Forb dominant", "Grass dominant")) + #legend colors and labels 
  theme(legend.position = "none") + #remove the legend
  scale_y_continuous(labels = NULL, name = NULL, limits = c(0,150)) + #remove y-axis label
  annotate("text", x= 2, y = 130, label = "Avena", color = "#999999", angle = 0) +
  annotate("text", x= 2, y = 140, label = "Forb", color = "#56B4E9", angle = 0) +
  annotate("text", x= 2, y = 145, label = "Grass", color = "green", angle = 0) +
  annotate("text", x= 2, y = 120, label = "Medusahead", color = "Red", angle = 0) +
  annotate("text", x= 2, y = 150, label = "Total", color = "purple", angle = 0) +
  annotate("text", x= 2, y = 125, label = "Lolium", color = "#E69F00", angle = 0)
cover_spring

#compile interaction plots 
grid.arrange(cover_shelter, cover_fall, cover_spring, ncol = 3, widths = c(1.5,1.2,1.2))


#How does ANPP relate to total cover?
#recreate objects from ANPP analyses
May_ANPP<-read.csv("Plant_composition_data/ANPP/ANPP_CleanedData/ClimVar_ANPP-peak.csv")
#create subset with 2015 data only, remove compost subplot
May_ANPP_2015<-filter(May_ANPP, year=='2015', subplot !="C")
May_2015_XC<-filter(May_ANPP_2015, subplot=='XC')

May_2015 <- merge(May_ANPP_2015, gfprop_all_2015) %>% dplyr::select(-cover)

#create objects from community analysis
data<- cover %>% dplyr::select(-X, -species, -genus, -status, -func, -func2) %>% spread(species_name, cover)
data<-data %>% dplyr::select(-Unknown)
levels(cover$plot)
str(data)
levels(data$treatment)
levels(data$year)
levels(data$subplot)

data <- tibble::rowid_to_column(data, "subplot")
data2 <- filter(data, year =="2015", subplot!="C", subplot!="XC") %>% arrange(treatment, shelterBlock)

Treatment<-data2[,4]
Year<-data2[,3]
data2$ID <- seq.int(nrow(data2))
plotnames<-data2[,64]
cover.Bio<- data2 %>% dplyr::select(-c(1:6), -ID)
rownames(cover.Bio)<-plotnames
#check for empty rows
cover.Biodrop<-cover.Bio[rowSums(cover.Bio[, (1:57)]) ==0, ] #no empty rows, next step not needed
# remove empty rows cover.Biodrop<-cover.Bio[rowSums(cover.Bio[, (1:58)])  >0, ]

cover.rowsums <- rowSums(cover.Bio [1:57])

cover.relrow <- data.frame(cover.Bio /cover.rowsums)
cover.colmax<-sapply(cover.Bio ,max)
cover.relcolmax <- data.frame(sweep(cover.Bio ,2,cover.colmax,'/'))
cover.pa <- cover.Bio %>% mutate_each(funs(ifelse(.>0,1,0)), 1:57)
#calculate shannon-weiner
May_2015$H <- diversity(cover.Bio)
#calculate pielou's J
May_2015$J <- May_2015$H/log(specnumber(cover.Bio))

#make bray-curtis dissimilarity matrix
spp.bcd <- vegdist(cover.relrow)
#or
spp.pa.bcd<-vegdist(cover.pa, binary=T) #this calcs distances based on presence/absence, just in case

#run NMS ordination
spp.mds0 <-isoMDS(spp.bcd) #runs nms only once
spp.mds0  #by default 2 dimensions returned, stress is 18.995
ordiplot(spp.mds0)

#prefer to run multiple NMS ordinations
#with different starting configurations, and choose the best
#this function does that, and in addition does several other steps too
#including: 
#standardizing the data (though fuction call below turns this off with autotransform=F)
#calculating distance matrix (default bray-curtis)
#running NMDS with random starts
#rotation of axes to maximize variance of site scores on axis 1
#calculate species scores based on weighted averaging

help(metaMDS)
spp.mds<-metaMDS(cover.relrow, trace = FALSE, autotransform=T, trymax=100, k=4) #runs several with different starting configurations
#trace= TRUE will give output for step by step what its doing
#default is 2 dimensions, can put k=4 for 4 dimensions
spp.mds #solution converged after 30 tries, stress = 10
summary(spp.mds)

#plot results
stressplot(spp.mds, spp.bcd) #stressplot to show fit, fit is decent
ordiplot(spp.mds)
spscores1<-scores(spp.mds,display="sites",choices=1)
spscores2<-scores(spp.mds,display="sites",choices=2)
tplots<-data2[,2]
tplot_levels<-levels(tplots)
spscoresall<-data.frame(tplots,spscores1,spscores2)

#to color grasses and forbs labels:
colspec<- rep(c("plum1", "plum1", "palegreen", "palegreen", "palegreen", "palegreen", "palegreen", "palegreen", "palegreen", "palegreen", "palegreen", "plum1", "plum1", "plum1", "plum1","plum1", "palegreen", "palegreen", "plum1", "plum1", "plum1", "plum1", "plum1", "plum1", "plum1", "palegreen", "palegreen", "palegreen", "plum1", "plum1", "palegreen", "plum1", "plum1", "plum1", "palegreen", "plum1", "plum1", "plum1", "plum1", "plum1", "plum1", "plum1", "plum1", "plum1", "darkgreen", "plum1", "plum4", "plum4", "plum4", "plum4", "plum4", "plum4", "plum1", "plum4", "palegreen","palegreen", "plum1"))

#plots colored based on treatment
xc.plot <- ordiplot(spp.mds,choices=c(1,2), type = "none")   #Set up the plot
colsT2 <- rep(c("darkred","deepskyblue","goldenrod1", "Magenta"), each = 12) #color based on drought treatment
cols1 <- rep(c("darkred","deepskyblue","goldenrod1", "Magenta"))
shapes <- rep(c(15, 8, 17 ), each=1) #shapes on subplot
points(spscoresall$NMDS1,spscoresall$NMDS2,col=colsT2,pch=shapes) 
text(spp.mds, display = "species", cex=0.8, col=colspec) #label species
# add legend for treatment
legend("topright",legend=levels(Treatment), col=cols1, pch=19, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)
# add legend for year
legend("topleft",legend=levels(as.factor(as.character(data2$subplot))), col="black", pch=shapes, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)

permanova1 <- adonis(spp.bcd~data2$treatment*data2$subplot, perm=100, method="bray")
permanova1
#treatment does not significantly drive communities

#does diversity/evenness/richness change with precipitation variability (drought vs. control) or seasonability (drought treatments)?
Hm<-lme(H ~ treatment*subplot, random=~1|shelterBlock, May_2015, na.action=na.exclude)
summary(Hm)
anova(Hm)
r.squaredGLMM(Hm) #25% of variation explained by fixed effects, 40% by whole model (spatial variation?)
qqnorm(residuals(Hm))
qqline(residuals(Hm))
shapiro.test(residuals(Hm))
#normal
hLS<-lsmeans(Hm, ~subplot)
contrast(hLS, "pairwise") #forb is more diverse

ggplot(May_2015, aes(x=subplot, y=H))+
  #facet_wrap(~treatment)+
  geom_boxplot()
  #geom_smooth(method = "lm", se=FALSE)

Jm<-lme(J ~ treatment*subplot, random=~1|shelterBlock, May_2015, na.action=na.exclude)
summary(Jm)
anova(Jm)
r.squaredGLMM(Jm) #25% of variation explained by fixed effects, 40% by whole model (spatial variation?)
qqnorm(residuals(Jm))
qqline(residuals(Jm))
shapiro.test(residuals(Jm))
#normal
jLS<-lsmeans(Jm, ~subplot)
contrast(jLS, "pairwise") #forb is more even

ggplot(May_2015, aes(y=weight_g_m, x=percentForb, color = subplot))+
  geom_point()+
  geom_smooth(method="lm",se=F)
  #scale_colour_gradient(low="magenta", high="plum4")

#how does diversity differ by block?
ggplot(May_2015, aes(x=shelterBlock, y=H))+
  facet_wrap(~subplot)+
  geom_boxplot()

cover.Bio2<- data %>% filter(year=="2015", subplot=="XC") %>% dplyr::select(-c(1:6)) 
May_2015_XC$H<-diversity(cover.Bio2)
May_2015_XC$J <- May_2015_XC$H/log(specnumber(cover.Bio2))
Hm2<-lme(H ~ treatment, random=~1|shelterBlock, May_2015_XC, na.action=na.exclude)
summary(Hm2)
anova(Hm2)
r.squaredGLMM(Hm2) #25% of variation explained by fixed effects, 40% by whole model (spatial variation?)
qqnorm(residuals(Hm2))
qqline(residuals(Hm2))
shapiro.test(residuals(Hm2))
#normal
hLS2<-lsmeans(Hm2, ~treatment)
contrast(hLS2, "pairwise") #no differences in diversity by treatment for XC

Jm2<-lme(J ~ treatment, random=~1|shelterBlock, May_2015_XC, na.action=na.exclude)
summary(Jm2)
anova(Jm2)
r.squaredGLMM(Jm2) #25% of variation explained by fixed effects, 40% by whole model (spatial variation?)
qqnorm(residuals(Jm2))
qqline(residuals(Jm2))
shapiro.test(residuals(Jm2))
#not normal
jLS2<-lsmeans(Jm2, ~treatment)
contrast(jLS2, "pairwise") #forb is more even

###check trait differences
#does functional diversity change with precipitation variability (drought vs. control) or seasonability (drought treatments)?
traits<-read.csv("Traits/Traits_ProcessedData-GH/ClimVar_trait-diversity-GH.csv")
traits_2015<-traits %>% filter(subplot!="XC", subplot!="C", year=="2015")
May_2015<-merge(May_2015,traits_2015)
May_ANPP<-merge(May_ANPP,traits) 
May_ANPP <- May_ANPP %>%filter( subplot!="C")
May_all_XC <- May_ANPP %>% filter (subplot=="XC")

ggplot(May_ANPP, aes(y=RaoQ, x=forbCover, color = treatment, shape=subplot))+
  facet_wrap(~as.factor(as.character(year)))+
  geom_point()+
  geom_smooth(method="lm",se=F)
  #scale_colour_gradient(low="magenta", high="plum4")

gf_graph_all <- May_ANPP %>%
  group_by(subplot, treatment, year) %>%
  summarize(Grass=mean(grassCover), GrassSE=sd(grassCover)/sqrt(length(grassCover)))

gf_graph_forb <- May_ANPP %>%
  group_by(subplot, treatment, year) %>%
  summarize(Forb=mean(forbCover), ForbSE=sd(forbCover)/sqrt(length(forbCover)))

gf_graph_all <- merge(gf_graph_all, gf_graph_forb)
gf_graph_all <- gf_graph_all %>% gather(func, cover, -ForbSE, -GrassSE, -subplot, -treatment, -year) %>%
  gather(func2, SE, -func, -cover, -subplot, -treatment, -year) %>% dplyr::select(-func2)

ggplot(gf_graph_all, aes(x=treatment, y=cover, group = func, color=func)) + 
  geom_errorbar(aes(ymax = cover+SE, ymin = cover-SE), width=.25) + 
  geom_point() + 
  theme_bw() + 
  facet_wrap(~subplot*year, ncol=3)+
  labs(x="Treatment", y="Cover %") +
  theme(text = element_text(size=20))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#does functional diversity differ by subplot or treatment?
Fd<-lme(FDis ~ treatment*subplot, random=~1|shelterBlock, May_2015, na.action=na.exclude)
summary(Fd)
anova(Fd)
r.squaredGLMM(Fd) #42% of variation explained by fixed effects, 47% by whole model (spatial variation?)
qqnorm(residuals(Fd))
qqline(residuals(Fd))
shapiro.test(residuals(Fd))
#close to normal
Fd2<-lme(log(FDis+1) ~ treatment*subplot, random=~1|shelterBlock, May_2015, na.action=na.exclude)
summary(Fd2)
anova(Fd2)
r.squaredGLMM(Fd2) #39% of variation explained by fixed effects, 45% by whole model (spatial variation?)
qqnorm(residuals(Fd2))
qqline(residuals(Fd2))
shapiro.test(residuals(Fd2))
FdLS<-lsmeans(Fd2, ~subplot*treatment)
contrast(FdLS, "pairwise") #forb is more diverse

ggplot(May_2015, aes(x=subplot, y=FDis))+
  facet_wrap(~treatment)+
  geom_boxplot()

Rao<-lme(RaoQ ~ treatment*subplot, random=~1|shelterBlock, May_2015, na.action=na.exclude)
summary(Rao)
anova(Rao)
r.squaredGLMM(Rao) #51% of variation explained by fixed effects, 55% by whole model (spatial variation?)
qqnorm(residuals(Rao))
qqline(residuals(Rao))
shapiro.test(residuals(Rao))
#normal
rLS<-lsmeans(Rao, ~subplot)
contrast(rLS, "pairwise") #forb is more diverse

ggplot(May_2015, aes(x=subplot, y=RaoQ))+
  geom_boxplot()


#how does functional diversity differ by treatment for control plots?
May_2015_XC2<-May_all_XC %>% filter(year=="2015")
Rao2<-lme(RaoQ ~ treatment, random=~1|shelterBlock, May_2015_XC2, na.action=na.exclude)
summary(Rao2)
anova(Rao2)
r.squaredGLMM(Rao2) #17% of variation explained by fixed effects, 52% by whole model (spatial variation?)
qqnorm(residuals(Rao2))
qqline(residuals(Rao2))
shapiro.test(residuals(Rao2))
#normal
r2LS<-lsmeans(Rao, ~treatment)
contrast(r2LS, "pairwise") 

Fd2<-lme(FDis ~ treatment, random=~1|year/shelterBlock, May_2015_XC2, na.action=na.exclude)
summary(Fd2)
anova(Fd2)
r.squaredGLMM(Fd2) #14% of variation explained by fixed effects, 48% by whole model (spatial variation?)
qqnorm(residuals(Fd2))
qqline(residuals(Fd2))
shapiro.test(residuals(Fd2))
#normal
fd2LS<-lsmeans(Fd2, ~treatment)
contrast(fd2LS, "pairwise") #no differences

ggplot(May_all_XC, aes(x=shelterBlock, y=FDis))+
  facet_wrap(~year)+
  geom_boxplot()
#block A has higher functional dissimilarity
#if niche complementarity is in effect, we'd expect block A to have highest ecosystem function

#let's test if block A has highest FDis


ggplot(May_all_XC, aes(x=treatment, y=RaoQ, color=treatment))+
  facet_wrap(~year)+
  geom_boxplot()+
  geom_smooth(method="lm", se=F)


#does functional diversity stabilize the community?
#how does functional diversity relate to coefficient of variation for soil moisture?
