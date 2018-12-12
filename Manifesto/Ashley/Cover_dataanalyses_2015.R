library(tidyverse)
library(nlme)
library(ggplot2)
library(dplyr)
library(multcomp)
library(lsmeans)
library(MuMIn)

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

ggplot((gfproportion_2015), aes(x=treatment, y = percentForb)) + 
  geom_boxplot()

a<-lme(percentForb ~ treatment, random=~1|shelterBlock/subplot, data=gfproportion_2015,
       contrasts=list(treatment=contr.treatment))

summary(a)        
summary(glht(a,linfct=mcp(treatment="Tukey")), alternative="Bonferonni")

gfproportion_graph_2015 <- gfproportion_2015 %>%
  group_by(treatment) %>%
  summarize(meanprop=mean(percentForb), seprop=sd(percentForb)/sqrt(length(percentForb)))

ggplot(gfproportion_graph_2015, aes(x=treatment, y=meanprop)) + 
  geom_bar(stat="identity", position="dodge", fill = "gray") + 
  theme_bw() + 
  geom_errorbar(aes(ymax = meanprop+seprop, ymin = meanprop-seprop), width=.25, position=position_dodge(width=0.9)) + 
  labs(x="Treatment", y="Percent forbs") +
  theme(text = element_text(size=20))

ggplot(gfproportion_2015, aes(x=percentForb, y=totcover, color = treatment))+
  geom_point() + geom_smooth(method = "lm")

ggplot(gfproportion_2015, aes(x=percentForb, y=totcover, color = treatment))+
  facet_grid(~subplot)+
  geom_point() + geom_smooth(method = "lm")
############################################

#calculate proportion grass
gfproportionG_2015 <- gf_2015 %>%
  group_by(plot, subplot, treatment, shelterBlock, year) %>%
  mutate(totcover = sum(cover)) %>%
  mutate(percentGrass = (cover/totcover)*100) %>%
  filter(func == "grass") %>%
  dplyr::select(-func)

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

ggplot(gf2_graph_2015, aes(fill=genus, y=cover, x=shelterBlock, color=func2)) +
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
gf5_graph_2015$genus<- factor(gf5_graph_2015$genus, c("Anagalis", "Centaurea","Cerastium", "Convolvulus", "Erodium","Hypochaeris", "Rumex","Sherardia","Silene", "Trifolium","Vicia","Avena","Bromus","Cynodon","Hordeum","Lolium", "Taeniatherum", "Vulpia"))
ggplot(gf5_graph_2015, aes(fill=genus, colour=func2,  y=cover, x=treatment:func2)) +
  theme_bw()+
  facet_wrap(~subplot)+
  scale_fill_manual(values = c("goldenrod","orange", "darkorange2", "orangered", "firebrick","indianred4","indianred", "saddlebrown", "palegreen", "green4", "lightblue", "skyblue2", "skyblue4", "dodgerblue3", "royalblue3","navy", "black"))+
  geom_bar( stat="identity", position='stack')+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


library(ggpubr)
#does proportion of Avena drive community structure?
ggplot(gfprop_all_2015, aes(x=AvCover, y=totcover, color=treatment))+
  geom_point()+
  geom_smooth(method='lm')+
  labs(x="Percent Avena", y="Total Cover")

ggplot(gfprop_all_2015, aes(x=AvCover, y=percentGrass, color=treatment))+
  geom_point()+
  geom_smooth(method='lm' )+
  labs(x="Percent Avena", y="Percent Grass Cover")

ggplot(gfprop_all_2015, aes(x=AvCover, y=LolCover, color=treatment, shape=subplot))+
  geom_point()+
  #geom_smooth(method='lm' )+
  #facet_wrap(~treatment)+
  labs(x="Percent Avena", y="Percent Lolium")

gfprop_noF_2015<-gfprop_all_2015 %>% filter(subplot!="F")
gfprop_noF_season_2015 <- gfprop_all_2015 %>% filter(treatment=="fallDry","springDry")

ggscatterhist(gfprop_noF_2015, x = "AvCover", y = "LolCover",
              color = "treatment", size = 3, alpha = 0.6,
              margin.params = list(fill = "treatment", color = "black", size = 0.2))

ggscatterhist(gfprop_noF, x = "AvCover", y = "TaeCover",
              color = "treatment", size = 3, alpha = 0.6,
              margin.params = list(fill = "treatment", color = "black", size = 0.2))





ggplot(gfprop_all_2015, aes(x=subplot, y=percentForb, color=treatment))+
  geom_boxplot()

gfprop_graph_2015<- gfprop_all_2015 %>% gather(func, percent, 9:13)%>%
  tbl_df()

ggplot(gfprop_graph_2015, aes(x=treatment, y=percent, fill=func))+
  facet_wrap(~subplot)+
  geom_boxplot()


ggplot(gfprop_all_2015, aes(x=AvCover, y=totcover, color=treatment))+
  #facet_grid(~treatment)+
  geom_point()+
  geom_smooth(method='lm', se=FALSE)

ggplot(gfprop_all_2015, aes(x=TaeCover, y=totcover, color=treatment))+
  #facet_grid(~treatment)+
  geom_point()+
  geom_smooth(method='lm', se=FALSE)

gfprop_graph2_2015<-gfprop_all_2015 %>% gather(func, percent, 11:13) %>% filter(subplot=="XC")
ggplot(gfprop_graph2_2015, aes(x=percent, y=totcover, group=func, color=func))+
  facet_grid(~treatment)+
  geom_point()+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method='lm', se=FALSE)

ggplot(gfprop_noF_2015, aes(x=treatment, y=AvCover))+
  geom_boxplot()

#does treatment affect % cover for Lolium?
ggplot(gfprop_noF_2015, aes(x=treatment, y=LolCover))+
  geom_boxplot()

g7<-lme(LolCover ~ treatment, random=~1|shelterBlock, gfprop_noF_2015, na.action=na.exclude)
summary(g7)
anova(g7)
r.squaredGLMM(g7) #16% of variation explained by fixed effects, 40% by whole model (spacial variation?)
qqnorm(residuals(g7))
qqline(residuals(g7))
shapiro.test(residuals(g7))
#not normally distributed, try log transformation
g8<-lme(log(LolCover+1) ~treatment, random=~1|shelterBlock, gfprop_noF_2015, na.action=na.exclude)
summary(g8)
anova(g8)
r.squaredGLMM(g8) #13% of variation explained by fixed effects, 60% by whole model (spatial variation?)
qqnorm(residuals(g8))
qqline(residuals(g8))
shapiro.test(residuals(g8))
#normal
gLS8<-lsmeans(g8, ~treatment)
contrast(gLS8, "pairwise") #differences between fall dry and spring dry

#does treatment affect medusahead
g5w<-lme(TaeCover ~ treatment, random=~1|shelterBlock, gfprop_noF_2015, na.action=na.exclude)
summary(g5w)
anova(g5w)
r.squaredGLMM(g5w) #15% of variation explained by fixed effects, 22% by whole model (interannual variation?)
qqnorm(residuals(g5w))
qqline(residuals(g5w))
shapiro.test(residuals(g5w))
#not normally distributed, try log transformation
g6w<-lme(log(TaeCover+1) ~treatment, random=~1|shelterBlock, gfprop_noF_2015, na.action=na.exclude)
summary(g6w)
anova(g6w)
r.squaredGLMM(g6w) #17% of variation explained by fixed effects, 33% by whole model (interannual variation?)
qqnorm(residuals(g6w))
qqline(residuals(g6w))
shapiro.test(residuals(g6w))
#still not normal, but slightly better
gLS6w<-lsmeans(g6w, ~treatment)
contrast(gLS6w, "pairwise") #differences in consistent/fall dry and spring dry

ggplot(gfprop_noF_2015, aes(x=treatment, y=TaeCover))+
  geom_boxplot()

#medusahead released from competition in consistent dry and fall dry?


gfprop_noF_season_2015 <- gfprop_graph_2015 %>% filter(treatment!="consistentDry") %>% filter(treatment!="controlRain") %>% filter(subplot!="F")
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

#does functional diversity/evenness/richness change with precipitation variability (drought vs. control) or seasonability (drought treatments)?
#how does functional diversity differ by block?
#does functional diversity stabilize the community?
#how does functional diversity relate to coefficient of variation for soil moisture?
