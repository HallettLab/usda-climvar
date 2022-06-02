####DATA SET UP######

#load data
setwd("~/Dropbox/ClimVar/DATA/Plant_composition_data/BNPP/BNPP_CleanedData")
setwd("./Dropbox/ClimVar/DATA/Plant_composition_data/BNPP/BNPP_CleanedData")
BNPP0 <- read.csv("BNPP_MayHarvest_2015.csv", header = TRUE)

setwd("~/Dropbox/ClimVar/DATA/Plant_composition_data/ANPP/ANPP_CleanedData")
setwd("./Dropbox/ClimVar/DATA/Plant_composition_data/ANPP/ANPP_CleanedData")
ANPP0 <- read.csv("ClimVar_ANPP-peak.csv", header = TRUE)

#load packages
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(ggrepel)
library(nlme) 
library(multcomp)
library(MuMIn)

#inital QCQC
ANPP <- ANPP0 %>%
  unique() %>% #Remove potential duplicates
  filter(!is.na(weight_g_m)) #Check for any missing values

BNPP <- BNPP0 %>%
  unique() %>% #Remove potential duplicates
  filter(!is.na(bmass_g_m2)) #Check for any missing values

#Check for any spelling errors
str(ANPP)
levels(ANPP$subplot) #factor levels of subplot
levels(ANPP$treatment) #factor levels of treatment
levels(ANPP$shelterBlock) #factor levels of shelter block

str(BNPP)
levels(BNPP$subplot) #factor levels of subplot
levels(BNPP$depth) #factor levels of soil depth
levels(BNPP$treatment)  #factor levels of treatment
levels(BNPP$shelterBlock) #factor levels of shelter block

#Replace entry "20-Oct" to "10-20" in depth column in the BNPP dataset
BNPP$depth <- as.character(BNPP$depth) #set depth as character
BNPP <- BNPP %>%
  mutate(depth = replace(depth, depth == "20-Oct", "10-20")) %>%
  mutate(treatment_code = fct_recode(treatment, consDry = "consistentDry", control = "controlRain")) #make a column with shorter treatment names

#Calculate the aggregated BNPP in three levels of soil depth
BNPP1 <- BNPP %>%
  group_by(plot, subplot, treatment, treatment_code, shelterBlock, shelter, fall, spring) %>% #group data
  summarise(agg_BNPP = sum(bmass_g_m2))#sum BNPP 

#Join the aggregated BNPP dataset with ANPP
joined <- ANPP %>%
  filter(year == 2015) %>% #Filter the ANPP dataset by year 2015
  filter(subplot == "F" | subplot == "B" | subplot == "G") %>% #Filter the ANPP dataset by subplot B, F, or G
  inner_join(BNPP1, by = c( "plot", "subplot", "treatment", "shelterBlock")) %>% #Join datasets
  mutate(root_shoot = agg_BNPP/weight_g_m)  #Calculate the root:shoot ratio

####DATA VISUALIZATION####

#subset the data by depth
Top10 <- BNPP %>%
  filter(depth == "0-10")
Mid10 <- BNPP %>%
  filter(depth == "10-20")
Bottom10 <- BNPP %>%
  filter(depth == "20-30")

#name variables
Top10_trt <- Top10$treatment
Top10_biomass <- Top10$bmass_g_m2
Mid10_trt <- Mid10$treatment
Mid10_biomass <- Mid10$bmass_g_m2
Bottom10_trt <- Bottom10$treatment
Bottom10_biomass <- Bottom10$bmass_g_m2

#plot biomass vs. treatment 
par(mfrow = c(1,3))
boxplot(Top10_biomass ~ Top10_trt, main = "Depth 0-10")
boxplot(Mid10_biomass ~ Mid10_trt, main = "Depth 10-20")
boxplot(Bottom10_biomass ~ Bottom10_trt, main = "Depth 20-30")

#plot biomass vs. treatment (group by depth)
ggplot(BNPP, aes(x=treatment, y=bmass_g_m2, fill=depth, color=depth)) +
  geom_boxplot() +
  theme_classic()

#plot biomass vs. treatment (group by functional group)
Top10_plot <- ggplot(Top10, aes(x=treatment, y=bmass_g_m2, fill=subplot, color=subplot)) +
  geom_boxplot() +
  theme_classic()

Mid10_plot <- ggplot(Mid10, aes(x=treatment, y=bmass_g_m2, fill=subplot, color=subplot)) +
  geom_boxplot() +
  theme_classic()

Bottom10_plot <- ggplot(Bottom10, aes(x=treatment, y=bmass_g_m2, fill=subplot, color=subplot)) +
  geom_boxplot() +
  theme_classic()

ggarrange(Top10_plot, Mid10_plot, Bottom10_plot, 
          labels = c( " Depth 0-10" , " Depth 10-20", "Depth 20-30"),
          ncol = 1,
          nrow = 3)

#Plot a boxplot of BNPP by treatment, grouped by functional group and faceted by soil depth
ggplot(BNPP, aes(x=treatment_code, y=bmass_g_m2, fill = subplot)) +
  geom_boxplot() +
  labs(fill = "Functional\nGroups", x = "Treatment", y = expression(paste("BNPP (g/m"^2,")"))) + #label the axis and legend
  scale_fill_discrete(labels=c("Mixed", "Forb dominant", "Grass dominant")) + #rename the legend labels
  facet_wrap(~depth) + #facet by depth
  theme_bw() + #black and white background
  theme(legend.position=c(.9,.8)) #position the legend

#plot aggregated biomass vs. treatment by functional group
BNPP1$treatment <- factor(BNPP1$treatment, levels = c("controlRain","springDry", "fallDry","consistentDry"))
total_biomass <- ggplot(BNPP1, aes(x=subplot, y=agg_BNPP, fill = treatment)) +
  geom_boxplot() +
  theme_classic() + #simple classic theme
  labs(fill = "Drought\nTreatment", x = "Community Treatment", y = expression(paste("Aggregated BNPP (g/m"^2,") in soil depth 0-30 cm"))) + #label the legend title and axis
  scale_fill_manual(values = c("#56B4E9","#999999","#E69F00","red"))+
  annotate("text", x=c("B", "F", "G"), y= c(800,800, 800), label = c("a", "b", "b"))
total_biomass

#spatial variance?
block_biomass <- ggplot(BNPP1, aes(x=treatment, y=agg_BNPP, fill = shelterBlock)) +
  geom_boxplot()

block_biomass1 <- ggplot(BNPP1, aes(x=shelterBlock, y=agg_BNPP))+
  geom_boxplot()

block_biomass2 <- ggplot(BNPP1, aes(x=shelterBlock, y=agg_BNPP, fill = treatment))+
  geom_boxplot()

#plot root:shoot ratio vs. treatment by functional group 
root_shoot_plot <- ggplot(joined, aes(x=treatment, y=root_shoot, fill = subplot)) +
  geom_boxplot() +
  theme_classic() + #simple, classic theme 
  labs(fill = "Functional\nGroups", x = "Treatment", y = expression(paste("Root : shoot ratio"))) + #label the legend title and axis
  scale_fill_manual(values = c("#999999","#E69F00","#56B4E9"), labels = c("Mixed", "Forb dominant", "Grass dominant")) #legend labels
root_shoot_plot

#spatial variance?
block_root_shoot <- ggplot(joined, aes(x=treatment, y=root_shoot, fill = shelterBlock)) +
  geom_boxplot()

####ANALYSIS####
#Calculate mean and standard error of aggregated BNPP by treatment
se <- function(x) sqrt(var(x)/length(x)) #create a function for SE
summary_BNPP_shelter <- BNPP1 %>%
  filter(treatment == "consistentDry" | treatment == "controlRain") %>%
  group_by(shelter, subplot) %>% #group by shelter and functional groups
  summarise(mean = mean(agg_BNPP), #summarise by mean and SE
            SE = se(agg_BNPP))

summary_BNPP_fall <- BNPP1 %>%
  filter(treatment == "controlRain" | treatment == "fallDry") %>%
  group_by(fall, subplot) %>% #group by fall rain treatment and functional groups
  summarise(mean = mean(agg_BNPP), #summarise by mean and SE
            SE = se(agg_BNPP))

summary_BNPP_spring <- BNPP1 %>%
  filter(treatment == "controlRain" | treatment == "springDry") %>%
  group_by(spring, subplot) %>% #gropu by spring rain treatment and functional groups
  summarise(mean = mean(agg_BNPP), #summarise by mean and SE
            SE = se(agg_BNPP))

#Calculate mean and standard error of root:shoot by treatment
summary_RS_shelter <- joined %>%
  filter(treatment == "consistentDry" | treatment == "controlRain") %>%
  group_by(shelter.y, subplot) %>% #group by shelter and functional groups
  summarise(mean = mean(root_shoot), #summarise by mean and SE
            SE = se(root_shoot))

summary_RS_fall <- joined %>%
  filter(treatment == "controlRain" | treatment == "fallDry") %>%
  group_by(fall, subplot) %>% #group by fall rain treatment and functional groups
  summarise(mean = mean(root_shoot), #summarise by mean and SE
            SE = se(root_shoot))

summary_RS_spring <- joined %>%
  filter(treatment == "controlRain" | treatment == "springDry") %>%
  group_by(spring, subplot) %>% #gropu by spring rain treatment and functional groups
  summarise(mean = mean(root_shoot), #summarise by mean and SE
            SE = se(root_shoot))

#interaction plot of shelter treatment and functional groups
aggBNPP_shelter <- ggplot(summary_BNPP_shelter, aes(x = as.factor(shelter), y = mean, group = subplot, color = subplot)) +
  geom_line(size = 1.0) + #add lines
  geom_point(size = 2.0) + #add points
  geom_errorbar(aes(ymin = mean-SE, ymax = mean+SE), width = 0.1, size = 1.0) + #add error bars
  theme_classic() + #simple, classic theme 
  labs(x = "Shelter", y = expression(paste("Aggregated BNPP (g/m"^2,") in soil depth 0-30 cm"))) + #label the axis
  scale_color_manual(values = c("#999999","#E69F00","#56B4E9")) + #specify color 
  theme(legend.position = "none") + #remove the legend 
  scale_y_continuous(limits = c(100, 610)) +
  annotate("text", x= 1.5, y = 380, label = "Mixed", color = "#999999", angle = -60) +
  annotate("text", x= 1.5, y = 240, label = "Grass", color = "#56B4E9", angle = -30) +
  annotate("text", x= 1.5, y = 180, label = "Forb", color = "#E69F00", angle = -30)

#interaction plot of fall rain treatment and functional groups
summary_BNPP_fall$fall <- factor(summary_BNPP_fall$fall, levels = c(1,0))
aggBNPP_fall <- 
  ggplot(summary_BNPP_fall, aes(x = factor(fall), y = mean, group = subplot, color = subplot)) +
  geom_line(size = 1.0) + #add lines
  geom_point(size = 2.0) + #add points
  geom_errorbar(aes(ymin = mean-SE, ymax = mean+SE), width = 0.1, size = 1.0) + #add error bars
  theme_classic() + #simple, classic theme 
  labs(x = "Fall Rain") + #label the axis
  scale_color_manual(values = c("#999999","#E69F00","#56B4E9")) + #specify color
  theme(legend.position = "none") + #remove the legend
  scale_y_continuous(labels = NULL, name = NULL, limits = c(100, 610)) + #remove y-axis label
  #scale_x_discrete(limits=c(1,0)) + #change order of discrete x scale 
  annotate("text", x= 1.5, y = 385, label = "Mixed", color = "#999999", angle = -53) +
  annotate("text", x= 1.5, y = 275, label = "Grass", color = "#56B4E9", angle = 20) +
  annotate("text", x= 1.5, y = 190, label = "Forb", color = "#E69F00", angle = -20)

#interaction plot of spring rain treatment and functional groups
summary_BNPP_spring$spring <- factor(summary_BNPP_spring$spring, levels = c(1,0))
aggBNPP_spring <- ggplot(summary_BNPP_spring, aes(x = factor(spring), y = mean, group = subplot, color = subplot)) +
  geom_line(size = 1.0) + #add lines
  geom_point(size = 2.0) + #add points
  geom_errorbar(aes(ymin = mean-SE, ymax = mean+SE), width = 0.1, size = 1.0) + #add error bars
  theme_classic() + #simple, classic theme 
  labs(x = "Spring Rain", y = expression(paste("Aggregated BNPP (g/m"^2,") in soil depth 0-30 cm"))) + #label the axis
  scale_color_manual(values = c("#999999","#E69F00","#56B4E9"), labels = c("Mixed", "Forb dominant", "Grass dominant")) + #legend colors and labels 
  theme(legend.position = "none") + #remove the legend
  scale_y_continuous(labels = NULL, name = NULL, limits = c(100, 610)) + #remove y-axis label
  annotate("text", x= 1.5, y = 445, label = "Mixed", color = "#999999", angle = 9) +
  annotate("text", x= 1.5, y = 235, label = "Grass", color = "#56B4E9", angle = -43) +
  annotate("text", x= 1.5, y = 215, label = "Forb", color = "#E69F00", angle = 20)

#compile interaction plots (agg BNPP)
grid.arrange(aggBNPP_shelter, aggBNPP_fall, aggBNPP_spring, ncol = 3, widths = c(1.5,1.2,1.2))

#interaction plot of shelter treatment and functional groups
RS_shelter <- ggplot(summary_RS_shelter, aes(x = as.factor(shelter.y), y = mean, group = subplot, color = subplot)) +
  geom_line(size = 1.0) + #add lines
  geom_point(size = 2.0) + #add points
  geom_errorbar(aes(ymin = mean-SE, ymax = mean+SE), width = 0.1, size = 1.0) + #add error bars
  theme_classic() + #simple, classic theme 
  labs(x = "Shelter", y = expression(paste("Root: shoot ratio"))) + #label the axis
  scale_color_manual(values = c("#999999","#E69F00","#56B4E9")) + #specify color 
  theme(legend.position = "none") + #remove the legend 
  scale_y_continuous(limits = c(0.3, 1.5)) +
  annotate("text", x= 1.5, y = 0.71, label = "Mixed", color = "#999999", angle = -35) +
  annotate("text", x= 1.2, y = 0.70, label = "Grass", color = "#56B4E9", angle = -48) +
  annotate("text", x= 1.7, y = 0.56, label = "Forb", color = "#E69F00", angle = -48)

#interaction plot of fall rain treatment and functional groups
summary_RS_fall$fall <- factor(summary_RS_fall$fall, levels = c(1,0))
RS_fall <- ggplot(summary_RS_fall, aes(x = factor(fall), y = mean, group = subplot, color = subplot)) +
  geom_line(size = 1.0) + #add lines
  geom_point(size = 2.0) + #add points
  geom_errorbar(aes(ymin = mean-SE, ymax = mean+SE), width = 0.1, size = 1.0) + #add error bars
  theme_classic() + #simple, classic theme 
  labs(x = "Fall Rain") + #label the axis
  scale_color_manual(values = c("#999999","#E69F00","#56B4E9")) + #specify color
  theme(legend.position = "none") + #remove the legend
  scale_y_continuous(labels = NULL, name = NULL, limits = c(0.3, 1.5)) + #remove y-axis label
  annotate("text", x= 1.5, y = 0.75, label = "Mixed", color = "#999999", angle = -30) +
  annotate("text", x= 1.2, y = 0.72, label = "Grass", color = "#56B4E9", angle = -30) +
  annotate("text", x= 1.7, y = 0.66, label = "Forb", color = "#E69F00", angle = -30)


#interaction plot of spring rain treatment and functional groups
summary_RS_spring$spring <- factor(summary_RS_spring$spring, levels = c(1,0))
RS_spring <- ggplot(summary_RS_spring, aes(x = factor(spring), y = mean, group = subplot, color = subplot)) +
  geom_line(size = 1.0) + #add lines
  geom_point(size = 2.0) + #add points
  geom_errorbar(aes(ymin = mean-SE, ymax = mean+SE), width = 0.1, size = 1.0) + #add error bars
  theme_classic() + #simple, classic theme 
  labs(x = "Spring Rain", y = expression(paste("Aggregated BNPP (g/m"^2,") in soil depth 0-30 cm"))) + #label the axis
  scale_color_manual(values = c("#999999","#E69F00","#56B4E9"), labels = c("Mixed", "Forb dominant", "Grass dominant")) + #legend colors and labels 
  theme(legend.position = "none") + #remove the legend
  scale_y_continuous(labels = NULL, name = NULL, limits = c(0.3, 1.5)) + #remove y-axis label
  annotate("text", x= 1.5, y = 0.85, label = "Mixed", color = "#999999", angle = 25) +
  annotate("text", x= 1.5, y = 0.63, label = "Grass", color = "#56B4E9", angle = -50) +
  annotate("text", x= 1.5, y = 0.76, label = "Forb", color = "#E69F00", angle = 32)

#compile interaction plots (Root:shoot ratio)
grid.arrange(RS_shelter, RS_fall, RS_spring, ncol = 3, widths = c(1.5,1.2,1.2))

#compare CV for B, F, G
cv <- function(mean, sd){
  (sd/mean)*100
}
cv_all <- joined %>%
  group_by(subplot) %>%
  summarise(cv_BNPP = cv(mean = mean(agg_BNPP),sd = sd(agg_BNPP)), 
                cv_root_shoot = cv(mean= mean(root_shoot), sd = sd(root_shoot)))

#Remove Block A and reproduce the interaction plots
BNPP2 <- BNPP1 %>%
  filter(shelterBlock != "A")
  
summary_BNPP_shelter_minusA <- BNPP2 %>%
  filter(treatment == "consistentDry" | treatment == "controlRain") %>%
  group_by(shelter, subplot) %>% #group by shelter and functional groups
  summarise(mean = mean(agg_BNPP), #summarise by mean and SE
            SE = se(agg_BNPP))

summary_BNPP_fall_minusA <- BNPP2 %>%
  filter(treatment == "controlRain" | treatment == "fallDry") %>%
  group_by(fall, subplot) %>% #group by fall rain treatment and functional groups
  summarise(mean = mean(agg_BNPP), #summarise by mean and SE
            SE = se(agg_BNPP))

summary_BNPP_spring_minusA <- BNPP2 %>%
  filter(treatment == "controlRain" | treatment == "springDry") %>%
  group_by(spring, subplot) %>% #gropu by spring rain treatment and functional groups
  summarise(mean = mean(agg_BNPP), #summarise by mean and SE
            SE = se(agg_BNPP))


##Calculate mean and standard error of root:shoot by treatment
joined1 <- joined %>%
  filter(shelterBlock != "A")

summary_RS_shelter_minusA <- joined1 %>%
  filter(treatment == "consistentDry" | treatment == "controlRain") %>%
  group_by(shelter.y, subplot) %>% #group by shelter and functional groups
  summarise(mean = mean(root_shoot), #summarise by mean and SE
            SE = se(root_shoot))

summary_RS_fall_minusA <- joined1 %>%
  filter(treatment == "controlRain" | treatment == "fallDry") %>%
  group_by(fall, subplot) %>% #group by fall rain treatment and functional groups
  summarise(mean = mean(root_shoot), #summarise by mean and SE
            SE = se(root_shoot))

summary_RS_spring_minusA <- joined1 %>%
  filter(treatment == "controlRain" | treatment == "springDry") %>%
  group_by(spring, subplot) %>% #gropu by spring rain treatment and functional groups
  summarise(mean = mean(root_shoot), #summarise by mean and SE
            SE = se(root_shoot))

#interaction plot of shelter treatment and functional groups minus Block A
aggBNPP_shelter_minusA <- ggplot(summary_BNPP_shelter_minusA, aes(x = as.factor(shelter), y = mean, group = subplot, color = subplot)) +
  geom_line(size = 1.0) + #add lines
  geom_point(size = 2.0) + #add points
  geom_errorbar(aes(ymin = mean-SE, ymax = mean+SE), width = 0.1, size = 1.0) + #add error bars
  theme_classic() + #simple, classic theme 
  labs(x = "Shelter", y = expression(paste("Aggregated BNPP (g/m"^2,") in soil depth 0-30 cm"))) + #label the axis
  scale_color_manual(values = c("#999999","#E69F00","#56B4E9")) + #specify color 
  theme(legend.position = "none") + #remove the legend 
  scale_y_continuous(limits = c(100, 610)) +
  annotate("text", x= 1.5, y = 358, label = "Mixed", color = "#999999", angle = -43) +
  annotate("text", x= 1.5, y = 255, label = "Grass", color = "#56B4E9", angle = -39) +
  annotate("text", x= 1.5, y = 205, label = "Forb", color = "#E69F00", angle = -38)

#interaction plot of fall rain treatment and functional groups minus Block A
summary_BNPP_fall_minusA$fall <- factor(summary_BNPP_fall_minusA$fall, levels = c(1,0))
aggBNPP_fall_minusA <- 
  ggplot(summary_BNPP_fall_minusA, aes(x = factor(fall), y = mean, group = subplot, color = subplot)) +
  geom_line(size = 1.0) + #add lines
  geom_point(size = 2.0) + #add points
  geom_errorbar(aes(ymin = mean-SE, ymax = mean+SE), width = 0.1, size = 1.0) + #add error bars
  theme_classic() + #simple, classic theme 
  labs(x = "Fall Rain") + #label the axis
  scale_color_manual(values = c("#999999","#E69F00","#56B4E9")) + #specify color
  theme(legend.position = "none") + #remove the legend
  scale_y_continuous(labels = NULL, name = NULL, limits = c(100, 610)) + #remove y-axis label
  #scale_x_discrete(limits=c(1,0)) + #change order of discrete x scale 
  annotate("text", x= 1.5, y = 355, label = "Mixed", color = "#999999", angle = -51) +
  annotate("text", x= 1.5, y = 310, label = "Grass", color = "#56B4E9", angle = 20) +
  annotate("text", x= 1.5, y = 225, label = "Forb", color = "#E69F00", angle = -20)

#interaction plot of spring rain treatment and functional groups minus Block A
summary_BNPP_spring_minusA$spring <- factor(summary_BNPP_spring_minusA$spring, levels = c(1,0))
aggBNPP_spring_minusA <- ggplot(summary_BNPP_spring_minusA, aes(x = factor(spring), y = mean, group = subplot, color = subplot)) +
  geom_line(size = 1.0) + #add lines
  geom_point(size = 2.0) + #add points
  geom_errorbar(aes(ymin = mean-SE, ymax = mean+SE), width = 0.1, size = 1.0) + #add error bars
  theme_classic() + #simple, classic theme 
  labs(x = "Spring Rain", y = expression(paste("Aggregated BNPP (g/m"^2,") in soil depth 0-30 cm"))) + #label the axis
  scale_color_manual(values = c("#999999","#E69F00","#56B4E9"), labels = c("Mixed", "Forb dominant", "Grass dominant")) + #legend colors and labels 
  theme(legend.position = "none") + #remove the legend
  scale_y_continuous(labels = NULL, name = NULL, limits = c(100, 610)) + #remove y-axis label
  annotate("text", x= 1.5, y = 382, label = "Mixed", color = "#999999", angle = -12) +
  annotate("text", x= 1.5, y = 247, label = "Grass", color = "#56B4E9", angle = -47) +
  annotate("text", x= 1.3, y = 234, label = "Forb", color = "#E69F00", angle = -10)

#compile interaction plots (agg BNPP) minus Block A
grid.arrange(aggBNPP_shelter_minusA, aggBNPP_fall_minusA, aggBNPP_spring_minusA, ncol = 3, widths = c(1.5,1.2,1.2))

#interaction plot of shelter treatment and functional groups minus Block A
RS_shelter_minusA <- ggplot(summary_RS_shelter_minusA, aes(x = as.factor(shelter.y), y = mean, group = subplot, color = subplot)) +
  geom_line(size = 1.0) + #add lines
  geom_point(size = 2.0) + #add points
  geom_errorbar(aes(ymin = mean-SE, ymax = mean+SE), width = 0.1, size = 1.0) + #add error bars
  theme_classic() + #simple, classic theme 
  labs(x = "Shelter", y = expression(paste("Root: shoot ratio"))) + #label the axis
  scale_color_manual(values = c("#999999","#E69F00","#56B4E9")) + #specify color 
  theme(legend.position = "none") + #remove the legend 
  scale_y_continuous(limits = c(0.3, 1.8)) +
  annotate("text", x= 1.5, y = 0.62, label = "Mixed", color = "#999999", angle = 15) +
  annotate("text", x= 1.2, y = 0.73, label = "Grass", color = "#56B4E9", angle = -52) +
  annotate("text", x= 1.5, y = 0.75, label = "Forb", color = "#E69F00", angle = -50)

#interaction plot of fall rain treatment and functional groups minus Block A
summary_RS_fall_minusA$fall <- factor(summary_RS_fall_minusA$fall, levels = c(1,0))
RS_fall_minusA <- ggplot(summary_RS_fall_minusA, aes(x = factor(fall), y = mean, group = subplot, color = subplot)) +
  geom_line(size = 1.0) + #add lines
  geom_point(size = 2.0) + #add points
  geom_errorbar(aes(ymin = mean-SE, ymax = mean+SE), width = 0.1, size = 1.0) + #add error bars
  theme_classic() + #simple, classic theme 
  labs(x = "Fall Rain") + #label the axis
  scale_color_manual(values = c("#999999","#E69F00","#56B4E9")) + #specify color
  theme(legend.position = "none") + #remove the legend
  scale_y_continuous(labels = NULL, name = NULL, limits = c(0.3, 1.8)) + #remove y-axis label
  annotate("text", x= 1.5, y = 0.67, label = "Mixed", color = "#999999", angle = 30) +
  annotate("text", x= 1.2, y = 0.80, label = "Grass", color = "#56B4E9", angle = -40) +
  annotate("text", x= 1.7, y = 0.80, label = "Forb", color = "#E69F00", angle = -30)


#interaction plot of spring rain treatment and functional groups minus Block A
summary_RS_spring_minusA$spring <- factor(summary_RS_spring_minusA$spring, levels = c(1,0))
RS_spring_minusA <- ggplot(summary_RS_spring_minusA, aes(x = factor(spring), y = mean, group = subplot, color = subplot)) +
  geom_line(size = 1.0) + #add lines
  geom_point(size = 2.0) + #add points
  geom_errorbar(aes(ymin = mean-SE, ymax = mean+SE), width = 0.1, size = 1.0) + #add error bars
  theme_classic() + #simple, classic theme 
  labs(x = "Spring Rain", y = expression(paste("Aggregated BNPP (g/m"^2,") in soil depth 0-30 cm"))) + #label the axis
  scale_color_manual(values = c("#999999","#E69F00","#56B4E9"), labels = c("Mixed", "Forb dominant", "Grass dominant")) + #legend colors and labels 
  theme(legend.position = "none") + #remove the legend
  scale_y_continuous(labels = NULL, name = NULL, limits = c(0.3, 1.8)) + #remove y-axis label
  annotate("text", x= 1.2, y = 0.65, label = "Mixed", color = "#999999", angle = 35) +
  annotate("text", x= 1.5, y = 0.73, label = "Grass", color = "#56B4E9", angle = -55) +
  annotate("text", x= 1.5, y = 0.98, label = "Forb", color = "#E69F00", angle = 27)

#compile interaction plots (Root:shoot ratio) minus Block A
grid.arrange(RS_shelter_minusA, RS_fall_minusA, RS_spring_minusA, ncol = 3, widths = c(1.5,1.2,1.2))


###Randomized block ANOVA biomass ~ functional group * treatment (block as random)
#check assumptions for randomized block ANOVA

#Normality of the response variable at each level of the factor
BNPP.agg <- with(BNPP, aggregate(data.frame(bmass_g_m2), 
                                    by = list(A = treatment, B = subplot), mean))
boxplot(bmass_g_m2 ~ A, BNPP.agg) 
#Homogeneity of variance
with(BNPP.agg, plot(tapply(bmass_g_m2, A, var), 
                      tapply(bmass_g_m2, A, mean)))

#determine whether the design is balanced
!is.list(replications(bmass_g_m2 ~ treatment + subplot + depth, BNPP)) #returns FALSE
replications(bmass_g_m2 ~ treatment + subplot + depth, BNPP) #unbalanced 

#fit mixed model using lme
Top10_fit <- aov(data = Top10, Top10_biomass ~ Top10_trt * subplot + Error(shelterBlock))
summary(Top10_fit)
TukeyHSD(aov(data = Top10, Top10_biomass ~ Top10_trt * subplot))

Mid10_fit <- aov(data = Mid10, Mid10_biomass ~ Mid10_trt * subplot + Error(shelterBlock))
summary(Mid10_fit)
TukeyHSD(aov(data = Mid10, Mid10_biomass ~ Mid10_trt * subplot))

Bottom10_fit <- aov(data = Bottom10, Bottom10_biomass ~ Bottom10_trt * subplot + Error(shelterBlock))
summary(Bottom10_fit)
TukeyHSD(aov(data = Bottom10, Bottom10_biomass ~ Bottom10_trt * subplot))

total.lme <- lme(data = joined, agg_BNPP ~ treatment * subplot, random = ~1| shelterBlock)
anova(total.lme)
summary(total.lme)
r.squaredGLMM(total.lme) #30% of variation explained by fixed effects, 30% of variation explained by whole model
qqnorm(residuals(total.lme))
qqline(residuals(total.lme))
shapiro.test(residuals(total.lme)) #normally distributed

total <- aov(data = joined, agg_BNPP ~ treatment * subplot + Error(shelterBlock))
summary(total)
TukeyHSD(aov(data= joined, agg_BNPP ~ treatment * subplot ))
###Randomized block ANOVA root:shoot ~ functional group * treatment (block as random)
root_shoot.lme <- lme(data = joined, root_shoot ~ treatment * subplot, random = ~1| shelterBlock)
anova(root_shoot.lme)
summary(root_shoot.lme)
r.squaredGLMM(root_shoot.lme) #9% of variation explained by fixed effects, 10% of variation explained by whole model
qqnorm(residuals(root_shoot.lme))
qqline(residuals(root_shoot.lme))
shapiro.test(residuals(root_shoot.lme)) #not normally distributed

###Interaction of treatment and functional groups
#interaction aggBNPP shelter
joined$shelter.y <- as.factor(joined$shelter.y)
ft1 <- aov(data= joined, agg_BNPP ~ subplot * shelter.y + Error(shelterBlock))
summary(ft1) 
TukeyHSD(aov(data= joined, agg_BNPP ~ subplot * shelter.y))

#interaction aggBNPP fall rain
ft2 <- aov(data= joined, agg_BNPP ~ subplot * fall + Error(shelterBlock))
summary(ft2) 
TukeyHSD(aov(data= joined, agg_BNPP ~ subplot * fall))

#interaction aggBNPP spring rain
ft3 <- aov(data= joined, agg_BNPP ~ subplot * spring + Error(shelterBlock))
summary(ft3) 
TukeyHSD(aov(data= joined, agg_BNPP ~ subplot * spring))

#interaction root_shoot shelter
ft4 <- aov(data= joined, root_shoot ~ subplot * shelter.y + Error(shelterBlock))
summary(ft4) 
TukeyHSD(aov(data= joined,root_shoot ~ subplot * shelter.y))

#interaction root_shoot fall rain
ft5 <- aov(data= joined, root_shoot  ~ subplot * fall + Error(shelterBlock))
summary(ft5) 
TukeyHSD(aov(data= joined, root_shoot  ~ subplot * fall))

#interaction root_shoot spring rain
ft6 <- aov(data= joined, root_shoot  ~ subplot * spring + Error(shelterBlock))
summary(ft6)
TukeyHSD(aov(data= joined, root_shoot ~ subplot * spring))

#subset data by treatment
Amount_rain <- joined %>%
  filter(treatment == "consistentDry" | treatment == "controlRain")
Fall_rain <- joined %>%
  filter(treatment == "controlRain" | treatment == "fallDry")
Spring_rain <- joined %>%
  filter(treatment == "controlRain" | treatment == "springDry")

#aggBNPP controlRain vs. consistentlyDry 
ft7 <- aov(data = Amount_rain, agg_BNPP ~ subplot * shelter.y + Error(shelterBlock))
summary(ft7)
TukeyHSD(aov(data = Amount_rain, agg_BNPP ~ subplot * shelter.y))

#root_shoot controlRain vs. consistentDry
ft8 <- aov(data = Amount_rain, root_shoot ~ subplot * shelter.y + Error(shelterBlock))
summary(ft8)

#aggBNPP controlRain vs. fallDry
ft9 <- aov(data = Fall_rain, agg_BNPP ~ subplot * fall + Error(shelterBlock))
summary(ft9)

#root_shoot controlRain vs. fallDry
ft10 <- aov(data = Fall_rain, root_shoot ~ subplot * fall + Error(shelterBlock))
summary(ft10)

#aggBNPP controlRain vs. springDry
ft11 <- aov(data = Spring_rain, agg_BNPP ~ subplot * spring + Error(shelterBlock))
summary(ft11)
TukeyHSD(aov(data = Spring_rain, agg_BNPP ~ subplot * spring))

#root_shoot controlRain vs. springDry
ft12 <- aov(data = Spring_rain, root_shoot ~ subplot * spring + Error(shelterBlock))
summary(ft12)

#Is the decline of grass BNPP in spring drought significant?
Spring_rain_grass <- Spring_rain %>%
  filter(subplot == "G")
ft12 <- lme(data = Spring_rain_grass, agg_BNPP ~ spring, random=~1|shelterBlock)
summary(ft12) #Not significant p = 0.3


