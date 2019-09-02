##The purpose of this script is to test treatment effects on ANPP and cover for 2015
#To compare with Lina's belowground results

library(tidyverse)
library(nlme)
library(ggplot2)
library(multcomp)
library(lsmeans)
library(MuMIn)
library(gridExtra)

setwd("~/Dropbox/ClimVar/DATA/")
May_ANPP<-read.csv("Plant_composition_data/ANPP/ANPP_CleanedData/ClimVar_ANPP-peak.csv")
#remove compost subplots for later
May_ANPP_noC <- May_ANPP %>% filter(subplot!="C")
#create subset with 2015 data only, remove compost subplot
May_ANPP_2015<-filter(May_ANPP, year=='2015', subplot !="C")

#####This first part is using control (XC) data to see effects of drought treatments on the background community
##1. Does variability of rainfall affect forage production (H1)?
#Using a mixed model to examine the effects of rainfall timing (fixed effect) with shelterbloc nested within year as random effect
#create subset with no species manipulations (control community) only
May_2015_XC<-filter(May_ANPP_2015, subplot=='XC')

m1a<-lme(weight_g_m ~treatment, random=~1|shelterBlock, May_2015_XC, na.action=na.exclude)
summary(m1a)
anova(m1a)
r.squaredGLMM(m1a) #6% of variation explained by fixed effects, 35% by whole model (spatial variation?)
qqnorm(residuals(m1a))
qqline(residuals(m1a))
shapiro.test(residuals(m1a))
#normally distributed, continue
LS1<-lsmeans(m1a, ~treatment)
contrast(LS1, "pairwise")
#no differences in ANPP for XC subplots across rainfall treatments
#let's see it
ggplot(d=May_2015_XC, aes(x=treatment, y=weight_g_m)) +
  theme_linedraw()+
  #facet_wrap(~shelterBlock)+
  labs(x="rainfall treatment", y="ANPP g/m2")+
  geom_boxplot(aes(y=weight_g_m), shape=16)
#lots of variation across treatments, althought control rain tends to have higher biomass than drought trts

##2. Does seasonality of rainfall affect forage production (H1)?
##Expect peak ANPP to be highest when rainfall occurs during peak season (spring dry) or consistently (controlRain)
##Expect peak ANPP to be lowest when rainfall occurs late in season (fall dry)
#try again with control rain removed to compare only treatments with same total rainfall
May_2015_drought<-filter(May_2015_XC, treatment!='controlRain')

m2a<-lme(weight_g_m ~treatment, random=~1|shelterBlock, May_2015_drought, na.action=na.exclude)
summary(m2a)
anova(m2a)
r.squaredGLMM(m2a) #only 2% of variation explained by fixed effects, 28% explained by whole model (spatial variation?)
qqnorm(residuals(m2a))
qqline(residuals(m2a))
shapiro.test(residuals(m2a))
#normally distributed, continue
LS2<-lsmeans(m2a, ~treatment)
contrast(LS2, "pairwise")
#no differences in total ANPP among drought treatments
#same plot as above but removes control rain
ggplot(d=May_2015_drought, aes(x=treatment, y=weight_g_m)) +
  theme_linedraw()+
  labs(x="rainfall treatment", y="ANPP g/m2")+
  geom_boxplot(aes(y=weight_g_m), shape=16)
#no differences in biomass, spring drought seems to have greater variation


###################
##In this second part,XC is removed to look at composition plots only
###################
#3. compensation (H2)
#H2 will be confirmed if mixed plots have greater ANPP compared to grass only or forb only
May_ANPP_2015_noXC<-filter(May_ANPP_2015, subplot !="XC")
m3a<-lme(weight_g_m ~treatment*subplot, random=~1|shelterBlock, May_ANPP_2015_noXC, na.action=na.exclude)
summary(m3a)
anova(m3a) #subplot significant
r.squaredGLMM(m3a)#32% of variation explained by fixed effects, 32% by whole model
qqnorm(residuals(m3a))
qqline(residuals(m3a))
shapiro.test(residuals(m3a)) #normal
LS3a<-lsmeans(m3a, ~subplot*treatment)
contrast(LS3a, "pairwise") #mixed plots have greater biomass than forb, but not grass

May_ANPP_2015_noXC$treatment <- factor(May_ANPP_2015_noXC$treatment, levels = c("controlRain","springDry", "fallDry","consistentDry"))
ggplot(d=May_ANPP_2015_noXC, aes(x=subplot, y=weight_g_m, fill=treatment)) +
  #facet_wrap(~treatment)+
  theme_classic()+
  labs(fill = "Drought\nTreatment", x="Community Treatment", y="ANPP g/m2")+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
  scale_fill_manual(values = c("#56B4E9","#999999","#E69F00","red")) +
  annotate("text", x=c("B", "F", "G"), y= c(925,925, 925), label = c("a", "b", "ab"))+
  geom_boxplot(aes(y=weight_g_m), shape=16)

ggplot(d=May_ANPP_2015_noXC, aes(x=subplot, y=weight_g_m, fill = subplot)) +
  theme_linedraw()+
  labs(x="composition treatment", y="ANPP g/m2")+
  annotate("text", x= c("B", "F","G"), y = c(925,925, 925), label = c("a", "b", "ab"), color = "#999999") +
  geom_boxplot(aes(y=weight_g_m), shape=16)

###################
##check Forb plots by themselves
###################
#3. Are there treatment effects on biomass?
May_ANPP_2015_G<-filter(May_ANPP_2015_noXC, subplot =="G")
m3b<-lme(weight_g_m ~treatment, random=~1|shelterBlock, May_ANPP_2015_G, na.action=na.exclude)
summary(m3b)
anova(m3b) #subplot significant
r.squaredGLMM(m3b)#32% of variation explained by fixed effects, 32% by whole model
qqnorm(residuals(m3b))
qqline(residuals(m3b))
shapiro.test(residuals(m3b)) #normal
LS3b<-lsmeans(m3b, ~treatment)
contrast(LS3b, "pairwise") #mixed plots have greater biomass than forb, but not grass

ggplot(d=May_ANPP_2015_noXC, aes(x=treatment, y=weight_g_m, fill=treatment)) +
  facet_wrap(~subplot)+
  theme_linedraw()+
  labs(x="drought treatment", y="ANPP g/m2")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  geom_boxplot(aes(y=weight_g_m), shape=16)

ggplot(d=May_ANPP_2015_noXC, aes(x=subplot, y=weight_g_m, fill = treatment)) +
  theme_linedraw()+
  labs(x="composition treatment", y="ANPP g/m2")+
  geom_boxplot(aes(y=weight_g_m), shape=16)

#make plots to match Lina's BNPP plots
#remove XC to match Lina's analyses
#note to see XC, use May_ANPP_2015
se <- function(x) sqrt(var(x)/length(x)) #create a function for SE
summary_ANPP_shelter <- May_ANPP_2015_noXC %>%
  filter(treatment == "consistentDry" | treatment == "controlRain") %>%
  group_by(shelter, subplot) %>% #group by shelter and functional groups
  summarise(mean = mean(weight_g_m), #summarise by mean and SE
            SE = se(weight_g_m))

summary_ANPP_fall <- May_ANPP_2015_noXC %>%
  filter(treatment == "controlRain" | treatment == "fallDry") %>%
  group_by(treatment, subplot) %>% #group by fall rain treatment and functional groups
  summarise(mean = mean(weight_g_m), #summarise by mean and SE
            SE = se(weight_g_m))

summary_ANPP_spring <- May_ANPP_2015_noXC %>%
  filter(treatment == "controlRain" | treatment == "springDry") %>%
  group_by(treatment, subplot) %>% #gropu by spring rain treatment and functional groups
  summarise(mean = mean(weight_g_m), #summarise by mean and SE
            SE = se(weight_g_m))


#interaction plot of shelter treatment and functional groups
ANPP_shelter <- ggplot(summary_ANPP_shelter, aes(x = as.factor(shelter), y = mean, group = subplot, color = subplot)) +
  geom_line(size = 1.0) + #add lines
  geom_point(size = 2.0) + #add points
  geom_errorbar(aes(ymin = mean-SE, ymax = mean+SE), width = 0.1, size = 1.0) + #add error bars
  theme_classic() + #simple, classic theme 
  labs(x = "Shelter", y = expression(paste("ANPP (g/m"^2,")"))) + #label the axis
  scale_color_manual(values = c("#999999","#E69F00","#56B4E9", "Red")) + #specify color 
  theme(legend.position = "none") + #remove the legend 
  scale_y_continuous(limits = c(100, 800)) +
  annotate("text", x= 1.5, y = 535, label = "Mixed", color = "#999999", angle = -40) +
  annotate("text", x= 1.35, y = 450, label = "Grass", color = "#56B4E9", angle = 40) +
  #annotate("text", x= 1.5, y = 485, label = "Control", color = "Red", angle = -10) +
  annotate("text", x= 1.5, y = 358, label = "Forb", color = "#E69F00", angle = 18)
ANPP_shelter

#interaction plot of fall rain treatment and functional groups
ANPP_fall <- 
  ggplot(summary_ANPP_fall, aes(x = treatment, y = mean, group = subplot, color = subplot)) +
  geom_line(size = 1.0) + #add lines
  geom_point(size = 2.0) + #add points
  geom_errorbar(aes(ymin = mean-SE, ymax = mean+SE), width = 0.1, size = 1.0) + #add error bars
  theme_classic() + #simple, classic theme 
  labs(x = "Fall Rain") + #label the axis
  scale_color_manual(values = c("#999999","#E69F00","#56B4E9", "Red")) + #specify color
  theme(legend.position = "none") + #remove the legend
  scale_y_continuous(labels = NULL, name = NULL, limits = c(100, 800)) + #remove y-axis label
  #scale_x_discrete(limits=c(1,0)) + #change order of discrete x scale 
  annotate("text", x= 1.5, y = 560, label = "Mixed", color = "#999999", angle = -10) +
  annotate("text", x= 1.5, y = 445, label = "Grass", color = "#56B4E9", angle = 28) +
  #annotate("text", x= 1.5, y = 515, label = "Control", color = "Red", angle = -31) +
  annotate("text", x= 1.5, y = 325, label = "Forb", color = "#E69F00", angle = -5)

ANPP_fall

#interaction plot of spring rain treatment and functional groups
ANPP_spring <- ggplot(summary_ANPP_spring, aes(x = treatment, y = mean, group = subplot, color = subplot)) +
  geom_line(size = 1.0) + #add lines
  geom_point(size = 2.0) + #add points
  geom_errorbar(aes(ymin = mean-SE, ymax = mean+SE), width = 0.1, size = 1.0) + #add error bars
  theme_classic() + #simple, classic theme 
  labs(x = "Spring Rain", y = expression(paste("Aggregated BNPP (g/m"^2,") in soil depth 0-30 cm"))) + #label the axis
  scale_color_manual(values = c("#999999","#E69F00","#56B4E9", "Red"), labels = c("Mixed", "Forb dominant", "Grass dominant")) + #legend colors and labels 
  theme(legend.position = "none") + #remove the legend
  scale_y_continuous(labels = NULL, name = NULL, limits = c(100, 800)) + #remove y-axis label
  annotate("text", x= 1.5, y = 650, label = "Mixed", color = "#999999", angle = 45) +
  annotate("text", x= 1.5, y = 425, label = "Grass", color = "#56B4E9", angle = 3) +
  #annotate("text", x= 1.5, y = 485, label = "Control", color = "Red", angle = -15) +
  annotate("text", x= 1.5, y = 365, label = "Forb", color = "#E69F00", angle = 30)
ANPP_spring

#compile interaction plots 
grid.arrange(ANPP_shelter, ANPP_fall, ANPP_spring, ncol = 3, widths = c(1.5,1.2,1.2))


############################
############################
##redo analyses with block A removed
#3. compensation (H2)
#H2 will be confirmed if mixed plots have greater ANPP compared to grass only or forb only
May_ANPP_2015_noA<-filter(May_ANPP_2015, shelterBlock !="A", subplot!="XC")
m4a<-lme(weight_g_m ~treatment*subplot, random=~1|shelterBlock, May_ANPP_2015_noA, na.action=na.exclude)
summary(m4a)
anova(m4a) #subplot significant
r.squaredGLMM(m4a)#39% of variation explained by fixed effects, 42% by whole model
qqnorm(residuals(m4a))
qqline(residuals(m4a))
shapiro.test(residuals(m4a)) #normal
LS4a<-lsmeans(m4a, ~subplot)
contrast(LS4a, "pairwise") #mixed plots have greater biomass forb plots and marginally greater than grass

ggplot(d=May_ANPP_2015_noA, aes(x=subplot, y=weight_g_m, fill=treatment)) +
  theme_linedraw()+
  labs(x="composition treatment", y="ANPP g/m2")+
  geom_boxplot(aes(y=weight_g_m), shape=16)

ggplot(d=May_ANPP_2015_noA, aes(x=subplot, y=weight_g_m, fill=subplot)) +
  theme_linedraw()+
  labs(x="composition treatment", y="ANPP g/m2")+
  geom_boxplot(aes(y=weight_g_m), shape=16)

##1. Does variability of rainfall affect forage production (H1)?
#Using a mixed model to examine the effects of rainfall timing (fixed effect) with shelterbloc nested within year as random effect

#create subset with no species manipulations (control community) only
May_2015_XC_noA<-filter(May_ANPP_2015, shelterBlock!="A", subplot=='XC')

m1A<-lme(weight_g_m ~treatment, random=~1|shelterBlock, May_2015_XC_noA, na.action=na.exclude)
summary(m1A)
anova(m1A)
r.squaredGLMM(m1A) #17% of variation explained by fixed effects, 33% by whole model (spatial variation?)
qqnorm(residuals(m1A))
qqline(residuals(m1A))
shapiro.test(residuals(m1A))
#normally distributed, continue
LS1<-lsmeans(m1a, ~treatment)
contrast(LS1, "pairwise")
#no differences in ANPP for XC subplots across rainfall treatments
#let's see it
ggplot(d=May_2015_XC_noA, aes(x=treatment, y=weight_g_m)) +
  theme_linedraw()+
  #facet_wrap(~shelterBlock)+
  labs(x="rainfall treatment", y="ANPP g/m2")+
  geom_boxplot(aes(y=weight_g_m), shape=16)

#make plots to match Lina's BNPP plots
summary_ANPP_shelter_noA <- May_ANPP_2015_noA %>%
  filter(treatment == "consistentDry" | treatment == "controlRain") %>%
  group_by(shelter, subplot) %>% #group by shelter and functional groups
  summarise(mean = mean(weight_g_m), #summarise by mean and SE
            SE = se(weight_g_m))

summary_ANPP_fall_noA <- May_ANPP_2015_noA %>%
  filter(treatment == "controlRain" | treatment == "fallDry") %>%
  group_by(treatment, subplot) %>% #group by fall rain treatment and functional groups
  summarise(mean = mean(weight_g_m), #summarise by mean and SE
            SE = se(weight_g_m))

summary_ANPP_spring_noA <- May_ANPP_2015_noA %>%
  filter(treatment == "controlRain" | treatment == "springDry") %>%
  group_by(treatment, subplot) %>% #gropu by spring rain treatment and functional groups
  summarise(mean = mean(weight_g_m), #summarise by mean and SE
            SE = se(weight_g_m))


#interaction plot of shelter treatment and functional groups
ANPP_shelter_noA <- ggplot(summary_ANPP_shelter_noA, aes(x = as.factor(shelter), y = mean, group = subplot, color = subplot)) +
  geom_line(size = 1.0) + #add lines
  geom_point(size = 2.0) + #add points
  geom_errorbar(aes(ymin = mean-SE, ymax = mean+SE), width = 0.1, size = 1.0) + #add error bars
  theme_classic() + #simple, classic theme 
  labs(x = "Shelter", y = expression(paste("ANPP (g/m"^2,")"))) + #label the axis
  scale_color_manual(values = c("#999999","#E69F00","#56B4E9", "Red")) + #specify color 
  theme(legend.position = "none") + #remove the legend 
  scale_y_continuous(limits = c(100, 1000)) +
  annotate("text", x= 1.5, y = 600, label = "Mixed", color = "#999999", angle = -40) +
  annotate("text", x= 1.35, y = 440, label = "Grass", color = "#56B4E9", angle = 60) +
  #annotate("text", x= 1.5, y = 518, label = "Control", color = "Red", angle = -40) +
  annotate("text", x= 1.5, y = 320, label = "Forb", color = "#E69F00", angle = 18)
ANPP_shelter_noA

#interaction plot of fall rain treatment and functional groups
ANPP_fall_noA <- 
  ggplot(summary_ANPP_fall_noA, aes(x = treatment, y = mean, group = subplot, color = subplot)) +
  geom_line(size = 1.0) + #add lines
  geom_point(size = 2.0) + #add points
  geom_errorbar(aes(ymin = mean-SE, ymax = mean+SE), width = 0.1, size = 1.0) + #add error bars
  theme_classic() + #simple, classic theme 
  labs(x = "Fall Rain") + #label the axis
  scale_color_manual(values = c("#999999","#E69F00","#56B4E9", "Red")) + #specify color
  theme(legend.position = "none") + #remove the legend
  scale_y_continuous(labels = NULL, name = NULL, limits = c(100, 1000)) + #remove y-axis label
  #scale_x_discrete(limits=c(1,0)) + #change order of discrete x scale 
  annotate("text", x= 1.5, y = 600, label = "Mixed", color = "#999999", angle = -40) +
  annotate("text", x= 1.35, y = 435, label = "Grass", color = "#56B4E9", angle = 55) +
  #annotate("text", x= 1.5, y = 500, label = "Control", color = "Red", angle = -45) +
  annotate("text", x= 1.5, y = 325, label = "Forb", color = "#E69F00", angle = 15)

ANPP_fall_noA

#interaction plot of spring rain treatment and functional groups
ANPP_spring_noA <- ggplot(summary_ANPP_spring_noA, aes(x = treatment, y = mean, group = subplot, color = subplot)) +
  geom_line(size = 1.0) + #add lines
  geom_point(size = 2.0) + #add points
  geom_errorbar(aes(ymin = mean-SE, ymax = mean+SE), width = 0.1, size = 1.0) + #add error bars
  theme_classic() + #simple, classic theme 
  labs(x = "Spring Rain", y = expression(paste("Aggregated BNPP (g/m"^2,") in soil depth 0-30 cm"))) + #label the axis
  scale_color_manual(values = c("#999999","#E69F00","#56B4E9", "Red"), labels = c("Mixed", "Forb dominant", "Grass dominant")) + #legend colors and labels 
  theme(legend.position = "none") + #remove the legend
  scale_y_continuous(labels = NULL, name = NULL, limits = c(100, 1000)) + #remove y-axis label
  annotate("text", x= 1.5, y = 695, label = "Mixed", color = "#999999", angle = 45) +
  annotate("text", x= 1.5, y = 405, label = "Grass", color = "#56B4E9", angle = 30) +
  #annotate("text", x= 1.5, y = 522, label = "Control", color = "Red", angle = -40) +
  annotate("text", x= 1.5, y = 340, label = "Forb", color = "#E69F00", angle = 30)
ANPP_spring_noA

#compile interaction plots 
grid.arrange(ANPP_shelter_noA, ANPP_fall_noA, ANPP_spring_noA, ncol = 3, widths = c(1.5,1.2,1.2))


##################
#################
#make plots for whole experiment (all years) to compare
se <- function(x) sqrt(var(x)/length(x)) #create a function for SE

summary_ANPP_shelter_all <- May_ANPP_noC %>%
  filter(treatment == "consistentDry" | treatment == "controlRain") %>%
  group_by(shelter, subplot) %>% #group by shelter and functional groups
  summarise(mean = mean(weight_g_m), #summarise by mean and SE
            SE = se(weight_g_m))

summary_ANPP_fall_all <- May_ANPP_noC %>%
  filter(treatment == "controlRain" | treatment == "fallDry") %>%
  group_by(treatment, subplot) %>% #group by fall rain treatment and functional groups
  summarise(mean = mean(weight_g_m), #summarise by mean and SE
            SE = se(weight_g_m))

summary_ANPP_spring_all <- May_ANPP_noC %>%
  filter(treatment == "controlRain" | treatment == "springDry") %>%
  group_by(treatment, subplot) %>% #gropu by spring rain treatment and functional groups
  summarise(mean = mean(weight_g_m), #summarise by mean and SE
            SE = se(weight_g_m))


#interaction plot of shelter treatment and functional groups
ANPP_shelter_all <- ggplot(summary_ANPP_shelter_all, aes(x = as.factor(shelter), y = mean, group = subplot, color = subplot)) +
  geom_line(size = 1.0) + #add lines
  geom_point(size = 2.0) + #add points
  geom_errorbar(aes(ymin = mean-SE, ymax = mean+SE), width = 0.1, size = 1.0) + #add error bars
  theme_classic() + #simple, classic theme 
  labs(x = "Shelter", y = expression(paste("ANPP (g/m"^2,")"))) + #label the axis
  scale_color_manual(values = c("#999999","#E69F00","#56B4E9", "Red")) + #specify color 
  theme(legend.position = "none") + #remove the legend 
  scale_y_continuous(limits = c(100, 800)) +
  annotate("text", x= 1.5, y = 535, label = "Mixed", color = "#999999", angle = -40) +
  annotate("text", x= 1.35, y = 450, label = "Grass", color = "#56B4E9", angle = 40) +
  annotate("text", x= 1.5, y = 485, label = "Control", color = "Red", angle = -10) +
  annotate("text", x= 1.5, y = 358, label = "Forb", color = "#E69F00", angle = 18)
ANPP_shelter_all

#interaction plot of fall rain treatment and functional groups
ANPP_fall_all <- 
  ggplot(summary_ANPP_fall_all, aes(x = treatment, y = mean, group = subplot, color = subplot)) +
  geom_line(size = 1.0) + #add lines
  geom_point(size = 2.0) + #add points
  geom_errorbar(aes(ymin = mean-SE, ymax = mean+SE), width = 0.1, size = 1.0) + #add error bars
  theme_classic() + #simple, classic theme 
  labs(x = "Fall Rain") + #label the axis
  scale_color_manual(values = c("#999999","#E69F00","#56B4E9", "Red")) + #specify color
  theme(legend.position = "none") + #remove the legend
  scale_y_continuous(labels = NULL, name = NULL, limits = c(100, 800)) + #remove y-axis label
  #scale_x_discrete(limits=c(1,0)) + #change order of discrete x scale 
  annotate("text", x= 1.5, y = 560, label = "Mixed", color = "#999999", angle = -10) +
  annotate("text", x= 1.5, y = 445, label = "Grass", color = "#56B4E9", angle = 28) +
  annotate("text", x= 1.5, y = 515, label = "Control", color = "Red", angle = -31) +
  annotate("text", x= 1.5, y = 325, label = "Forb", color = "#E69F00", angle = -5)

ANPP_fall_all

#interaction plot of spring rain treatment and functional groups
ANPP_spring_all <- ggplot(summary_ANPP_spring_all, aes(x = treatment, y = mean, group = subplot, color = subplot)) +
  geom_line(size = 1.0) + #add lines
  geom_point(size = 2.0) + #add points
  geom_errorbar(aes(ymin = mean-SE, ymax = mean+SE), width = 0.1, size = 1.0) + #add error bars
  theme_classic() + #simple, classic theme 
  labs(x = "Spring Rain", y = expression(paste("Aggregated BNPP (g/m"^2,") in soil depth 0-30 cm"))) + #label the axis
  scale_color_manual(values = c("#999999","#E69F00","#56B4E9", "Red"), labels = c("Mixed", "Forb dominant", "Grass dominant")) + #legend colors and labels 
  theme(legend.position = "none") + #remove the legend
  scale_y_continuous(labels = NULL, name = NULL, limits = c(100, 800)) + #remove y-axis label
  annotate("text", x= 1.5, y = 650, label = "Mixed", color = "#999999", angle = 45) +
  annotate("text", x= 1.5, y = 425, label = "Grass", color = "#56B4E9", angle = 3) +
  annotate("text", x= 1.5, y = 485, label = "Control", color = "Red", angle = -15) +
  annotate("text", x= 1.5, y = 365, label = "Forb", color = "#E69F00", angle = 30)
ANPP_spring_all

#compile interaction plots 
grid.arrange(ANPP_shelter_all, ANPP_fall_all, ANPP_spring_all, ncol = 3, widths = c(1.5,1.2,1.2))

