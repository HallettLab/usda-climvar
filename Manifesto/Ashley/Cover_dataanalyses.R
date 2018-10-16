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
  group_by(plot, subplot, func, treatment, shelterBlock, year) %>%
  summarize(cover=sum(cover)) %>%
  tbl_df()

gfproportion <- gf %>%
  group_by(plot, subplot, treatment, shelterBlock, year) %>%
  mutate(totcover = sum(cover)) %>%
  mutate(percentForb = (cover/totcover)*100) %>%
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

ggplot(gfproportion, aes(x=percentForb, y=totcover, color = treatment))+
  geom_point() + geom_smooth(method = "lm")

ggplot(gfproportion, aes(x=percentForb, y=totcover, color = treatment))+
  facet_grid(~subplot)+
  geom_point() + geom_smooth(method = "lm")
############################################

#calculate proportion grass
gfproportionG <- gf %>%
  group_by(plot, subplot, treatment, shelterBlock, year) %>%
  mutate(totcover = sum(cover)) %>%
  mutate(percentGrass = (cover/totcover)*100) %>%
  filter(func == "grass") %>%
  dplyr::select(-func)

#calc prop forb
gfproportionF <- gf %>%
  group_by(plot, subplot, treatment, shelterBlock, year) %>%
  mutate(totcover = sum(cover)) %>%
  mutate(percentForb = (cover/totcover)*100) %>%
  filter(func == "forb") %>%
  dplyr::select(-func) %>%
  dplyr::select(-cover)

gf2<-cover_noC %>%
  group_by(plot, year, subplot, species_name, func, treatment, shelterBlock) %>%
  mutate(cover=sum(cover)) %>%
  tbl_df()
  

#proportion avena
vegAv2 <- gf2 %>%
  group_by(plot, subplot, treatment, shelterBlock, year) %>%
  mutate(totcover = sum(cover)) %>%
  mutate(AvCover = cover/totcover * 100) %>%
  filter(species_name == "Avena barbata") %>%
  dplyr::select(-func,-species_name,-cover, -genus, -species, -status, -func2, -X)

#proportion lolium
vegLol <- gf2 %>%
  group_by(plot, subplot, treatment, shelterBlock, year)%>%
  mutate(totcover=sum(cover))%>%
  mutate(LolCover=cover/totcover*100) %>%
  filter(species_name == "Lolium multiflorum") %>%
  dplyr::select(-func,-species_name,-cover, -genus, -species, -status, -func2, -totcover, -X)

#proportion lolium
vegTae <- gf2 %>%
  group_by(plot, subplot, treatment, shelterBlock, year)%>%
  mutate(totcover=sum(cover))%>%
  mutate(TaeCover=cover/totcover*100) %>%
  filter(genus == "Taeniatherum") %>%
  dplyr::select(-func,-species_name,-cover, -genus, -species, -status, -func2, -totcover, -X)

# merge prop grass with forb prop and prop avena
gfprop_all<- merge(gfproportionG, gfproportionF, all.x = T) %>%
  tbl_df() 

gfprop_all<-merge(gfprop_all, vegAv2, all.x = T) %>%
  tbl_df() 

gfprop_all <-merge(gfprop_all,vegLol, all.x=T)%>%
  tbl_df() 

gfprop_all <-merge(gfprop_all,vegTae, all.x=T)%>% 
  tbl_df()

SE <- function(x){(sd(x)/sqrt(length(x)))}
MeanCov<-gfprop_all %>% group_by(treatment)%>%summarise(meanFCover=mean(percentForb))
CovFse<-aggregate(percentForb ~ treatment, data= gfprop_all, FUN = SE)
colnames(CovFse)[colnames(CovFse)=="percentForb"] <- "Fcov.SE"
#from ANNP_dataanalyses.R
MeanCov2<- merge(MeanFunc, MeanCov)
MeanCov2<-merge(MeanCov2,CovFse)

ggplot(MeanCov2, aes(x=sm_cv, y=meanFCover, color=treatment))+
  geom_errorbar(aes(ymin=meanFCover-Fcov.SE, ymax=meanFCover+Fcov.SE), width=1, position=pd)+
  geom_point(position=pd)+
  geom_smooth(method='lm')+
  labs(x="Soil Moisture CV", y="Percent Forb Cover")

gf3<-cover_noC %>%
  group_by(plot, year, subplot, treatment, shelterBlock) %>%
  mutate(totcover=sum(cover)) %>%
  tbl_df()


gf4 <- gf3 %>%
  group_by(year, subplot, treatment, shelterBlock, genus, species_name, func2) %>%
  summarize(percentCov = cover/totcover)%>%
  filter(subplot=="XC")%>%
  #filter(percentCov>0.10)%>%
  tbl_df()

gf5 <- gf3 %>%
  group_by(year, subplot, treatment, shelterBlock, genus, species_name, func2) %>%
  summarize(percentCov = cover/totcover)%>%
  #filter(percentCov>0.10)%>%
  tbl_df()

gf3_graph <- gf4 %>%
  group_by(treatment,genus, func2) %>%
  summarise(cover=mean(percentCov), secover=sd(percentCov)/sqrt(length(percentCov)))

gf5_graph <- gf5 %>%
  group_by(treatment, subplot, genus, func2)%>%
  summarise(cover=mean(percentCov), secover=sd(percentCov)/sqrt(length(percentCov)))

#create a stacked bar plot
ggplot(gf2_graph, aes(fill=species_name, y=cover, x=treatment)) +
geom_bar( stat="identity")

#now for control only using genera with > 1% cover
gf3_graph <- gf3_graph %>% arrange(func2, genus) %>% filter(cover>0.01)
gf3_graph$func2 <- factor(gf3_graph$func2, c("forb","Nfixer","grass"))
gf3_graph$genus<- factor(gf3_graph$genus, c("Centaurea", "Convolvulus", "Erodium","Hypochaeris","Trifolium","Vicia","Avena","Bromus","Cynodon","Hordeum","Lolium", "Taeniatherum"))
ggplot(gf3_graph, aes(fill=genus, colour=func2,  y=cover, x=treatment:func2)) +
  theme_bw()+
  scale_fill_manual(values = c("orange", "orangered", "firebrick","indianred4", "palegreen", "green4", "lightblue", "skyblue2", "skyblue4", "dodgerblue3", "royalblue3","navy"))+
  geom_bar( stat="identity", position='stack')+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#lets look at cover for the composition treatments
gf5_graph <- gf5_graph %>% arrange(func2, genus) %>% filter(cover>0.01) #%>%filter(genus!="NA")
levels(gf5_graph$genus)
gf5_graph$func2 <- factor(gf5_graph$func2, c("forb","Nfixer","grass"))
gf5_graph$genus<- factor(gf5_graph$genus, c("Anagalis", "Centaurea","Cerastium", "Convolvulus", "Erodium","Hypochaeris", "Rumex","Sherardia","Silene", "Trifolium","Vicia","Avena","Bromus","Cynodon","Hordeum","Lolium", "Taeniatherum", "Vulpia"))
ggplot(gf5_graph, aes(fill=genus, colour=func2,  y=cover, x=treatment:func2)) +
  theme_bw()+
  facet_wrap(~subplot)+
  scale_fill_manual(values = c("khaki1","yellow","goldenrod","orange", "darkorange2", "orangered", "firebrick","indianred4","saddlebrown", "palegreen", "green4", "lightblue", "skyblue2", "skyblue4", "dodgerblue3", "royalblue3","navy", "black"))+
  geom_bar( stat="identity", position='stack')+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


library(ggpubr)
#does proportion of Avena drive community structure?
ggplot(gfprop_all, aes(x=AvCover, y=percentForb, color=treatment))+
  geom_point()+
  geom_smooth(method='gam', formula= y~poly(x,2) )+
  labs(x="Percent Avena", y="Percent Forb Cover")

ggplot(gfprop_all, aes(x=AvCover, y=percentGrass, color=treatment))+
  geom_point()+
  geom_smooth(method='gam', formula= y~poly(x,2) )+
  labs(x="Percent Avena", y="Percent Grass Cover")

ggplot(gfprop_all, aes(x=AvCover, y=LolCover, color=treatment, shape=subplot))+
  geom_point()+
  #geom_smooth(method='lm' )+
  #facet_wrap(~treatment)+
  labs(x="Percent Avena", y="Percent Lolium")

gfprop_noF<-gfprop_all %>% filter(subplot!="F")
gfprop_noF_season <- gfprop_all %>% filter(treatment=="fallDry","springDry")

ggscatterhist(gfprop_noF, x = "AvCover", y = "LolCover",
  color = "treatment", size = 3, alpha = 0.6,
  margin.params = list(fill = "treatment", color = "black", size = 0.2))

ggscatterhist(gfprop_noF, x = "AvCover", y = "TaeCover",
              color = "treatment", size = 3, alpha = 0.6,
              margin.params = list(fill = "treatment", color = "black", size = 0.2))

#does grass cover increase with wet fall (springdry)? driven by avena?
ggscatterhist(gfprop_noF_wetyr, x = "AvCover", y = "percentGrass",
              color = "treatment", size = 3, alpha = 0.6,
              margin.params = list(fill = "treatment", color = "black", size = 0.2))

ggplot(gfprop_noF, aes(x=year, y=percentGrass, color=treatment))+
  geom_boxplot()+
  labs(x="year",y="Grass Cover")

g1<-lme(percentGrass ~treatment, random=~1|year/shelterBlock, gfprop_all, na.action=na.exclude)
summary(g1)
anova(g1)
r.squaredGLMM(g1) #1% of variation explained by fixed effects, 11% by whole model (interannual variation?)
qqnorm(residuals(g1))
qqline(residuals(g1))
shapiro.test(residuals(g1))
#not normally distributed, try sqrt transformation
g2<-lme(sqrt(percentGrass) ~treatment, random=~1|year/shelterBlock, gfprop_all, na.action=na.exclude)
summary(g2)
anova(g2)
r.squaredGLMM(g2) #1% of variation explained by fixed effects, 11% by whole model (interannual variation?)
qqnorm(residuals(g2))
qqline(residuals(g2))
shapiro.test(residuals(g2))
#log and sqrt are worse, stick with 1st model
LS1<-lsmeans(g1, ~treatment)
contrast(LS1, "pairwise")

g3<-lme(percentForb ~ year, random=~1|treatment/shelterBlock, gfprop_all, na.action=na.exclude)
summary(g3)
anova(g3)
r.squaredGLMM(g3) #7% of variation explained by fixed effects, 11% by whole model (interannual variation?)
qqnorm(residuals(g3))
qqline(residuals(g3))
shapiro.test(residuals(g3))
#not normally distributed, try sqrt transformation
g4<-lme(log(percentForb+1) ~year, random=~1|treatment/shelterBlock, gfprop_all, na.action=na.exclude)
summary(g4)
anova(g4)
r.squaredGLMM(g4) #7% of variation explained by fixed effects, 19% by whole model (interannual variation?)
qqnorm(residuals(g4))
qqline(residuals(g4))
shapiro.test(residuals(g4))
#normal
gLS4<-lsmeans(g4, ~year)
contrast(gLS4, "pairwise")

ggplot(gfprop_all, aes(x=year, y=percentForb, color=treatment))+
  geom_boxplot()

gfprop_graph<- gfprop_all %>% gather(func, percent, 9:12)%>%
tbl_df()

ggplot(gfprop_graph, aes(x=treatment, y=percent, fill=func))+
  facet_grid(~year)+
  geom_boxplot()


ggplot(gfprop_all, aes(x=percentGrass, y=AvCover, color=treatment))+
  #facet_grid(~treatment)+
  geom_point()+
  geom_smooth(method='lm', se=FALSE)

ggplot(gfprop_all, aes(x=percentGrass, y=LolCover, color=treatment))+
  #facet_grid(~treatment)+
  geom_point()+
  geom_smooth(method='lm', se=FALSE)

gfprop_graph2<-gfprop_all %>% gather(func, percent, 11:13) %>% filter(subplot=="XC")
ggplot(gfprop_graph2, aes(x=percentGrass, y=percent, color=func))+
  facet_grid(~treatment)+
  geom_point()+
  xlim(40,100)+
  ylim(0,100)+
  geom_smooth(method='lm', se=FALSE)

gfprop_noF_wetyr <- gfprop_noF %>% filter(year!="2015")

#first see if treatment effects % cover for Avena
g5<-lme(AvCover ~ treatment, random=~1|year/shelterBlock, gfprop_noF, na.action=na.exclude)
summary(g5)
anova(g5)
r.squaredGLMM(g5) #7% of variation explained by fixed effects, 36% by whole model (interannual variation?)
qqnorm(residuals(g5))
qqline(residuals(g5))
shapiro.test(residuals(g5))
#not normally distributed, try sqrt transformation
g6<-lme(log(AvCover+1) ~treatment, random=~1|year/shelterBlock, gfprop_noF, na.action=na.exclude)
summary(g6)
anova(g6)
r.squaredGLMM(g6) #9% of variation explained by fixed effects, 42% by whole model (interannual variation?)
qqnorm(residuals(g6))
qqline(residuals(g6))
shapiro.test(residuals(g6))
#normal
gLS6<-lsmeans(g6, ~treatment)
contrast(gLS6, "pairwise")

ggplot(gfprop_noF, aes(x=treatment, y=AvCover))+
  geom_boxplot()

#does treatment affect % cover for Lolium?
ggplot(gfprop_noF, aes(x=treatment, y=LolCover))+
  geom_boxplot()

g7<-lme(LolCover ~ treatment, random=~1|year/shelterBlock, gfprop_noF, na.action=na.exclude)
summary(g7)
anova(g7)
r.squaredGLMM(g7) #11% of variation explained by fixed effects, 51% by whole model (interannual variation?)
qqnorm(residuals(g7))
qqline(residuals(g7))
shapiro.test(residuals(g7))
#not normally distributed, try sqrt transformation
g8<-lme(log(LolCover+1) ~treatment, random=~1|year/shelterBlock, gfprop_noF, na.action=na.exclude)
summary(g8)
anova(g8)
r.squaredGLMM(g8) #7% of variation explained by fixed effects, 62% by whole model (interannual variation?)
qqnorm(residuals(g8))
qqline(residuals(g8))
shapiro.test(residuals(g8))
#normal
gLS8<-lsmeans(g8, ~treatment)
contrast(gLS8, "pairwise")

#are these effects stronger in the wet years?
#first see if treatment effects % cover for Avena
g5w<-lme(AvCover ~ treatment, random=~1|year/shelterBlock, gfprop_noF_wetyr, na.action=na.exclude)
summary(g5w)
anova(g5w)
r.squaredGLMM(g5w) #9% of variation explained by fixed effects, 42% by whole model (interannual variation?)
qqnorm(residuals(g5w))
qqline(residuals(g5w))
shapiro.test(residuals(g5w))
#not normally distributed, try sqrt transformation
g6w<-lme(log(AvCover+1) ~treatment, random=~1|year/shelterBlock, gfprop_noF_wetyr, na.action=na.exclude)
summary(g6w)
anova(g6w)
r.squaredGLMM(g6w) #14% of variation explained by fixed effects, 42% by whole model (interannual variation?)
qqnorm(residuals(g6w))
qqline(residuals(g6w))
shapiro.test(residuals(g6w))
#normal
gLS6w<-lsmeans(g6w, ~treatment)
contrast(gLS6w, "pairwise")

ggplot(gfprop_noF_wetyr, aes(x=treatment, y=AvCover))+
  geom_boxplot()

#does treatment affect % cover for Lolium?
ggplot(gfprop_noF_wetyr, aes(x=treatment, y=LolCover))+
  geom_boxplot()

g7w<-lme(LolCover ~ treatment, random=~1|year/shelterBlock, gfprop_noF_wetyr, na.action=na.exclude)
summary(g7w)
anova(g7w)
r.squaredGLMM(g7w) #10% of variation explained by fixed effects, 57% by whole model (interannual variation?)
qqnorm(residuals(g7w))
qqline(residuals(g7w))
shapiro.test(residuals(g7w))
#not normally distributed, try sqrt transformation
g8w<-lme(log(LolCover+1) ~treatment, random=~1|year/shelterBlock, gfprop_noF_wetyr, na.action=na.exclude)
summary(g8w)
anova(g8w)
r.squaredGLMM(g8w) #7% of variation explained by fixed effects, 69% by whole model (interannual variation?)
qqnorm(residuals(g8w))
qqline(residuals(g8w))
shapiro.test(residuals(g8w))
#normal
gLS8w<-lsmeans(g8w, ~treatment)
contrast(gLS8w, "pairwise")

gfprop_noF_season <- gfprop_graph %>% filter(treatment!="consistentDry") %>% filter(treatment!="controlRain") %>% filter(subplot!="F")
gfprop_noF_season_graph <- gfprop_noF_season %>%
  group_by(func, treatment) %>%
  summarize(meancover=mean(percent), secover=sd(percent)/sqrt(length(percent)))
gfprop_noF_season_graph<-gfprop_noF_season_graph%>%filter(func!="percentForb")

ggplot(gfprop_noF_season_graph, aes(x=treatment, y=meancover, color=func, group=func))+
  geom_line()+
  geom_point()+
  theme_bw()+
  geom_errorbar(aes(ymin=meancover-secover, ymax=meancover+secover), width=.2,position=position_dodge(0.05))+
  labs(x="Treatment", y="Mean Cover (%)")




#is cover specific to site (XC and control rain)?




#does functional diversity/evenness/richness change with precipitation variability (drought vs. control) or seasonability (drought treatments)?
#how does functional diversity relate to coefficient of variation for soil moisture?
