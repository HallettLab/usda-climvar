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
#create subset with 2015 data only, remove compost subplot
May_ANPP_2015<-filter(May_ANPP, year=='2015', subplot !="C")

##1. Does variability of rainfall affect forage production (H1)?
#Using a mixed model to examine the effects of rainfall timing (fixed effect) with shelterbloc nested within year as random effect

#create subset with no species manipulations (control community) only
May_2015_XC<-filter(May_ANPP_2015, subplot=='XC')

m1a<-lme(weight_g_m ~treatment, random=~1|shelterBlock, May_2015_XC, na.action=na.exclude)
summary(m1)
anova(m1)
r.squaredGLMM(m1) #6% of variation explained by fixed effects, 35% by whole model (interannual variation?)
qqnorm(residuals(m1))
qqline(residuals(m1))
shapiro.test(residuals(m1))
#normally distributed, continue
LS1<-lsmeans(m1, ~treatment)
contrast(LS1, "pairwise")
#no differences in ANPP for XC subplots across rainfall treatments
#let's see it
ggplot(d=May_2015_XC, aes(x=treatment, y=weight_g_m)) +
  theme_linedraw()+
  labs(x="rainfall treatment", y="ANPP g/m2")+
  geom_boxplot(aes(y=weight_g_m), shape=16)

##2. Does seasonality of rainfall affect forage production (H1)?
##Expect peak ANPP to be highest when rainfall occurs during peak season or consistently
##Expect peak ANPP to be lowest when rainfall occurs late in season
#try again with control rain removed to compare only treatments with same total rainfall
May_2015_drought<-filter(May_2015_XC, treatment!='controlRain')

m2a<-lme(weight_g_m ~treatment, random=~1|shelterBlock, May_2015_drought, na.action=na.exclude)
summary(m2a)
anova(m2a)
r.squaredGLMM(m2a) #only 2% of variation explained by fixed effects, 28% explained by whole model (interannual variation?)
qqnorm(residuals(m2a))
qqline(residuals(m2a))
shapiro.test(residuals(m2a))
#normally distributed, continue
LS2<-lsmeans(m2a, ~treatment)
contrast(LS2, "pairwise")
#no differences in total ANPP among drought treatments
#same plot as above but removes control rain
ggplot(d=May_XC_drought, aes(x=treatment, y=weight_g_m)) +
  theme_linedraw()+
  labs(x="rainfall treatment", y="ANPP g/m2")+
  geom_boxplot(aes(y=weight_g_m), shape=16)

#3. compensation (H2)
#H2 will be confirmed if mixed plots have greater ANPP compared to grass only or forb only
m3a<-lme(weight_g_m ~treatment*subplot, random=~1|shelterBlock, May_ANPP_2015, na.action=na.exclude)
summary(m2a)
anova(m2a)
r.squaredGLMM(m2a)#29% of variation explained by fixed effects, 35% by whole model
qqnorm(residuals(m2a))
qqline(residuals(m2a))
shapiro.test(residuals(m2a))
LS2a<-lsmeans(m2a, ~subplot*treatment)
contrast(LS2a, "pairwise")

ggplot(d=May_ANPP_2015, aes(x=treatment, y=weight_g_m, fill=subplot)) +
  theme_linedraw()+
  labs(x="composition treatment", y="ANPP g/m2")+
  geom_boxplot(aes(y=weight_g_m), shape=16)

#make plots to match Lina's BNPP plots
#note I kept control (XC) here for comparision, but can remove it
se <- function(x) sqrt(var(x)/length(x)) #create a function for SE
summary_ANPP_shelter <- May_ANPP_2015 %>%
  filter(treatment == "consistentDry" | treatment == "controlRain") %>%
  group_by(shelter, subplot) %>% #group by shelter and functional groups
  summarise(mean = mean(weight_g_m), #summarise by mean and SE
            SE = se(weight_g_m))

summary_ANPP_fall <- May_ANPP_2015 %>%
  filter(treatment == "controlRain" | treatment == "fallDry") %>%
  group_by(treatment, subplot) %>% #group by fall rain treatment and functional groups
  summarise(mean = mean(weight_g_m), #summarise by mean and SE
            SE = se(weight_g_m))

summary_ANPP_spring <- May_ANPP_2015 %>%
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
  annotate("text", x= 1.5, y = 485, label = "Control", color = "Red", angle = -10) +
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
  annotate("text", x= 1.5, y = 515, label = "Control", color = "Red", angle = -31) +
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
  annotate("text", x= 1.5, y = 485, label = "Control", color = "Red", angle = -15) +
  annotate("text", x= 1.5, y = 365, label = "Forb", color = "#E69F00", angle = 30)
ANPP_spring

#compile interaction plots 
grid.arrange(ANPP_shelter, ANPP_fall, ANPP_spring, ncol = 3, widths = c(1.5,1.2,1.2))

