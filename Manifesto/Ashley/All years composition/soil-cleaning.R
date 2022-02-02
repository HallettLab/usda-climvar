
#the purpose of this script is to clean up soil extracts data 
library(tidyverse)
library(stringr)
setwd("~/Dropbox/ClimVar/DATA/")

## Import soil extracts
soil <- read.csv("~/Dropbox/ClimVar/DATA/Soil Extracts/EXTRACTS_out_HEADERS2.csv", na.strings="99999")
soil$Treatment<-as.factor(soil$Treatment)

#extract notes & remove notes from df
soil_notes <- soil[,14:15] 
soil<-soil %>% dplyr::select(-(12:19)) %>% 
  dplyr::rename(plot=Plot) %>% dplyr::rename(subplot = Treatment)%>%
  mutate(subplot = recode(subplot, "1"= "XC", "2"="C", "3"="L", "4"="G", "5"="F", "6"="B"))

## Import import shelter key
shelterkey <- read.csv("~/Dropbox/ClimVar/DATA/Plant_composition_data/Cover/Cover_EnteredData/Shelter_key.csv")

#merge data with shelter key
soil <- merge(soil, shelterkey) %>%
  tbl_df() 

#write.csv(soil, "~/Dropbox/ClimVar/DATA/Soil Extracts/Soil_CleanedData/ClimVar_soil_all.csv")

soil[,'plot'] <- as.factor(as.character(soil[,'plot']))
soil[,'subplot'] <- as.factor(as.character(soil[,'subplot']))
soil[,'shelter'] <- as.factor(as.character(soil[,'shelter']))

#explore data
ggplot(soil, aes(x=subplot, y=SOC))+
  geom_point(aes(color=treatment))+
  facet_wrap(~Year*DOY)

soil_comp<-subset(soil, subplot %in% c("B","G","F"))
soil_comp<-soil_comp %>% group_by(Year, DOY, subplot, Depth, treatment)%>%
  mutate(TN=NH4+NO3.NO2)

#rename trt groups
soil_comp$treatment2[soil_comp$treatment=="controlRain"] <- "Ambient"
soil_comp$treatment2[soil_comp$treatment=="springDry"] <- "Late"
soil_comp$treatment2[soil_comp$treatment=="fallDry"] <- "Early"
soil_comp$treatment2[soil_comp$treatment=="consistentDry"] <- "Consistent"

fs1b<-ggplot(subset(soil_comp, Depth==1), aes(x=treatment2, y=SOC, fill=subplot))+
  #theme_linedraw()+
  geom_boxplot()+
  labs(x="  ", y="SOC (ug/g)")+
  #scale_fill_manual(values = c("royalblue2","sienna","lightsteelblue1", "peachpuff"), guide = guide_legend(title = "Treatment")) +
  scale_fill_manual(values = c("white","black","gray"), guide = guide_legend(title = "Composition"),
                    labels=c("Grass", "Forb", "Mixed")) +
  facet_wrap(~Year)
fs1b

ggplot(subset(soil_comp, Depth==1), aes(x=treatment2, y=SIC, fill=subplot))+
  #theme_linedraw()+
  geom_boxplot()+
  #scale_fill_manual(values = c("royalblue2","sienna","lightsteelblue1", "peachpuff"), guide = guide_legend(title = "Treatment")) +
  scale_fill_manual(values = c("white","black","gray"), guide = guide_legend(title = "Composition"),
                    labels=c("Grass-only", "Forb-only", "Mixed")) +
  facet_wrap(~Year)

ggplot(subset(soil_comp, Depth==1), aes(x=treatment2, y=MBC, fill=subplot))+
  #theme_linedraw()+
  geom_boxplot()+
  #scale_fill_manual(values = c("royalblue2","sienna","lightsteelblue1", "peachpuff"), guide = guide_legend(title = "Treatment")) +
  #scale_fill_manual(values = c("white","black","gray"), guide = guide_legend(title = "Composition"),
                    #labels=c("Grass-only", "Forb-only", "Mixed")) +
  facet_wrap(~Year)

ggplot(subset(soil_comp, Depth==1), aes(x=treatment2, y=MBIC, fill=subplot))+
  theme_linedraw()+
  #scale_fill_manual(values = c("royalblue2","sienna","lightsteelblue1", "peachpuff"), guide = guide_legend(title = "Treatment")) +
  scale_fill_manual(values = c("white","black","gray"), guide = guide_legend(title = "Composition"),
                    labels=c("Grass", "Forb", "Mixed")) +
  geom_boxplot()+
  facet_wrap(~Year)

fs1a<-ggplot(subset(soil_comp, Depth==1), aes(x=treatment2, y=TN, fill=subplot))+
  #theme_linedraw()+
  #scale_fill_manual(values = c("royalblue2","sienna","lightsteelblue1", "peachpuff"), guide = guide_legend(title = "Treatment")) +
  scale_fill_manual(values = c("white","black","gray"), guide = guide_legend(title = "Composition"),
                    labels=c("Grass", "Forb", "Mixed")) +
  labs(x=" ", y="Total N (ug/g)")+
  geom_boxplot()+
  facet_wrap(~Year)
fs1a

ggarrange(fs1a,fs1b,ncol=1, common.legend = TRUE, legend = "bottom",
          labels=c("a)","b)"))

ggplot(subset(soil_comp, Depth==1), aes(x=treatment2, y=NO3.NO2,fill=subplot))+
  #theme_linedraw()+
  #scale_fill_manual(values = c("royalblue2","sienna","lightsteelblue1", "peachpuff"), guide = guide_legend(title = "Treatment")) +
  scale_fill_manual(values = c("white","black","gray"), guide = guide_legend(title = "Composition"),
                    labels=c("Grass-only", "Forb-only", "Mixed")) +
   geom_boxplot()+
  facet_wrap(~Year)

ggplot(subset(soil_comp, Depth==1), aes(x=treatment2, y=TN,fill=subplot))+
  #theme_linedraw()+
  #scale_fill_manual(values = c("royalblue2","sienna","lightsteelblue1", "peachpuff"), guide = guide_legend(title = "Treatment")) +
  scale_fill_manual(values = c("white","black","gray"), guide = guide_legend(title = "Composition"),
                    labels=c("Grass-only", "Forb-only", "Mixed")) +
  geom_boxplot()+
  facet_wrap(~Year)


m.soc<-lme(SOC ~subplot*treatment, random=~1|shelterBlock, data=subset(soil_comp, Year==2015&Depth==1), na.action=na.exclude)
summary(m.soc)
anova(m.soc) #treatment is significant
r.squaredGLMM(m.soc) #35% of variation explained by fixed effects, 43% by whole model 
qqnorm(residuals(m.soc))
qqline(residuals(m.soc))
shapiro.test(residuals(m.soc))
LS.g<-lsmeans(m.soc, ~treatment)
contrast(LS.g, "pairwise")

m.mbc<-lme(MBC ~subplot*treatment, random=~1|shelterBlock, data=subset(soil_comp, Year==2015&Depth==1), na.action=na.exclude)
summary(m.mbc)
anova(m.mbc) #treatment is significant
r.squaredGLMM(m.mbc) #35% of variation explained by fixed effects, 43% by whole model 
qqnorm(residuals(m.mbc))
qqline(residuals(m.mbc))
shapiro.test(residuals(m.mbc))
LS.mbc<-lsmeans(m.mbc, ~treatment)
contrast(LS.mbc, "pairwise")

m.nh4<-lme(NH4 ~treatment*subplot, random=~1|shelterBlock, subset(soil_comp, Year==2015&Depth==1), na.action=na.exclude)
summary(m.nh4)
anova(m.nh4) #treatment is significant
r.squaredGLMM(m.nh4) #25% of variation explained by fixed effects, 38% by whole model (interannual variation?)
qqnorm(residuals(m.nh4))
qqline(residuals(m.nh4))
shapiro.test(residuals(m.nh4))
LS.nh4<-lsmeans(m.nh4, ~subplot, by="treatment")
contrast(LS.nh4, "pairwise")

m.no2<-lme(log(NO3.NO2+1) ~subplot*treatment, random=~1|shelterBlock, subset(soil_comp, Year==2015&Depth==1), na.action=na.exclude)
summary(m.no2)
anova(m.no2) #treatment is significant
r.squaredGLMM(m.no2) #36% of variation explained by fixed effects, 65% by whole model (interannual variation?)
qqnorm(residuals(m.no2))
qqline(residuals(m.no2))
shapiro.test(residuals(m.no2))
LS.no2<-lsmeans(m.no2, ~subplot, by="treatment")
contrast(LS.no2, "pairwise")

soil_comp2<-subset(soil_comp, Year==2015&Depth==1&DOY==168)%>%select(-Depth,-Year)%>%
  filter(subplot=="B")
part_soil<-merge(FG_B,soil_comp2)

f4c<-ggplot(subset(soil_comp, Depth==1), aes(x=treatment2, y=TN, fill=subplot))+
  #theme_linedraw()+
  #scale_fill_manual(values = c("royalblue2","sienna","lightsteelblue1", "peachpuff"), guide = guide_legend(title = "Treatment")) +
  scale_fill_manual(values = c("white","black","gray"), guide = guide_legend(title = "Composition"),
                    labels=c("Grass-only", "Forb-only", "Mixed")) +
  labs(x="Drought Treatment", y="NH4 (ug/g)")+
  geom_boxplot()#+
  #facet_wrap(~Year)
f4c

soil_comp$subplot <- ordered(soil_comp$subplot, levels = c("G","F","B"))
f4c<-ggplot(subset(soil_comp, Depth==1 & Year=="2015"), aes(x=subplot, y=TN, fill=subplot))+
  labs(x="", y=expression(Total~N~(ug/g)))+
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  scale_fill_manual(values = c("white","black", "gray")) +
  theme(legend.position="none")+
  annotate("text", x= c("B", "F","G"), y = c(10, 10, 10), label = c("a", "a", "b"), color = "black", size=5) +
  scale_x_discrete(labels=c("F" = "Forb", "G" = "Grass", "B" = "Mixed"), limits=levels(soil_comp$subplot))+
  geom_boxplot()#+
#facet_wrap(~Year)
f4c

f4d<-ggplot(subset(soil_comp, Depth==1 & subplot=="B"), aes(x=as.factor(Year), y=TN))+
  labs(x="", y=expression(Total~N~(ug/g)))+
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  annotate("text", x= c("2015", "2016"), y = c(9, 9), label = c("a", "b"), color = "black", size=5) +
  #theme_linedraw()+
  #scale_fill_manual(values = c("royalblue2","sienna","lightsteelblue1", "peachpuff"), guide = guide_legend(title = "Treatment")) +
  #scale_fill_manual(values = c("white","black","gray"), guide = guide_legend(title = "Composition"),
                    #labels=c("Grass-only", "Forb-only", "Mixed")) +
  #labs(x="Drought Treatment", y="NH4 (ug/g)")+
  geom_boxplot(fill="gray")#+
#facet_wrap(~Year)
f4d

f4a<-ggplot(subset(soil_comp, Depth==1), aes(x=treatment2, y=TN, fill=subplot))+
  ylim(0,15)+
  #theme_linedraw()+
  #scale_fill_manual(values = c("royalblue2","sienna","lightsteelblue1", "peachpuff"), guide = guide_legend(title = "Treatment")) +
  scale_fill_manual(values = c("white","black","gray"), guide = guide_legend(title = "Composition"),
                    labels=c("Grass-only", "Forb-only", "Mixed")) +
  labs(x="", y="Total N (ug/g)")+
  geom_boxplot()
  #facet_wrap(~Year)
f4a

f4a.v2<-ggplot(subset(soil_comp, Depth==1), aes(x=treatment2, y=TN, fill=subplot))+
  ylim(0,10)+
  #theme_linedraw()+
  #scale_fill_manual(values = c("royalblue2","sienna","lightsteelblue1", "peachpuff"), guide = guide_legend(title = "Treatment")) +
  scale_fill_manual(values = c("white","black","gray"), guide = guide_legend(title = "Composition"),
                    labels=c("Grass-only", "Forb-only", "Mixed")) +
  labs(x="", y="Total N (ug/g)")+
  geom_boxplot()
#facet_wrap(~Year)
f4a.v2

ggdraw(f4a) +
  draw_plot(f4c, .15, .65, .3, .3) +
  draw_plot(f4d, 0.5,0.65, 0.3,0.3)+
  draw_plot_label(
    c("A", "B", "C"),
    c(0, 0.125, 0.475),
    c(1, 0.95, 0.95),
    size = 14
  )

#plot big panel plus small panels on right
right_col2<-plot_grid(f4c, f4d, labels = c('E', 'F'), label_size = 14, ncol=1)
right_col2

F1D<-plot_grid(f4a.v2 + theme(legend.position="none"), right_col2, labels=c("D",""), rel_widths=c(2,1.5), label_size=14, ncol=2)
F1D

f5c<-ggplot(subset(soil_comp, Depth==1 & Year=="2015"), aes(x=subplot, y=SOC, fill=subplot))+
  labs(x="", y=expression(SOC~(ug/g)))+
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  scale_fill_manual(values = c("white","black", "gray")) +
  theme(legend.position="none")+
  annotate("text", x= c( "F","G","B"), y = c(300, 300, 300), label = c("a", "a", "a"), color = "black", size=5) +
  scale_x_discrete(labels=c("F" = "Forb", "G" = "Grass", "B" = "Mixed"), limits=levels(soil_comp$subplot))+
  geom_boxplot()#+
#facet_wrap(~Year)
f5c

f5d<-ggplot(subset(soil_comp, Depth==1 & subplot=="B"), aes(x=as.factor(Year), y=SOC))+
  labs(x="", y=expression(SOC~(ug/g)))+
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  annotate("text", x= c("2015", "2016"), y = c(300, 300), label = c("a", "b"), color = "black", size=5) +
  #theme_linedraw()+
  #scale_fill_manual(values = c("royalblue2","sienna","lightsteelblue1", "peachpuff"), guide = guide_legend(title = "Treatment")) +
  #scale_fill_manual(values = c("white","black","gray"), guide = guide_legend(title = "Composition"),
  #labels=c("Grass-only", "Forb-only", "Mixed")) +
  #labs(x="Drought Treatment", y="NH4 (ug/g)")+
  geom_boxplot(fill="gray")#+
#facet_wrap(~Year)
f5d

f5a<-ggplot(subset(soil_comp, Depth==1), aes(x=treatment2, y=SOC, fill=subplot))+
  ylim(0,400)+
  #theme_linedraw()+
  #scale_fill_manual(values = c("royalblue2","sienna","lightsteelblue1", "peachpuff"), guide = guide_legend(title = "Treatment")) +
  scale_fill_manual(values = c("white","black","gray"), guide = guide_legend(title = "Composition"),
                    labels=c("Grass-only", "Forb-only", "Mixed")) +
  labs(x="", y="SOC (ug/g)")+
  geom_boxplot()
#facet_wrap(~Year)
f5a

f5a.v2<-ggplot(subset(soil_comp, Depth==1), aes(x=treatment2, y=SOC, fill=subplot))+
  ylim(0,300)+
  #theme_linedraw()+
  #scale_fill_manual(values = c("royalblue2","sienna","lightsteelblue1", "peachpuff"), guide = guide_legend(title = "Treatment")) +
  scale_fill_manual(values = c("white","black","gray"), guide = guide_legend(title = "Composition"),
                    labels=c("Grass-only", "Forb-only", "Mixed")) +
  labs(x="", y="SOC (ug/g)")+
  geom_boxplot()
#facet_wrap(~Year)
f5a.v2

ggdraw(f5a) +
  draw_plot(f5c, .15, .65, .3, .3) +
  draw_plot(f5d, 0.5,0.65, 0.3,0.3)+
  draw_plot_label(
    c("D", "E", "F"),
    c(0, 0.125, 0.475),
    c(1, 0.95, 0.95),
    size = 14
  )

#plot big panel plus small panels on right
right_col3<-plot_grid(f5c, f5d, labels = c('H', 'I'), label_size = 14, ncol=1)
right_col3

F1G<-plot_grid(f5a.v2 + theme(legend.position="none"), right_col3, labels=c("G",""), rel_widths=c(2,1.5), label_size=14, ncol=2)
F1G

ggplot(subset(part_soil, year==2016), aes(x=TN, y=dY))+
  geom_point(aes(shape=treatment2),size=4)+
  geom_smooth(method="lm", se=F, aes(group=1), size=2, color="black")+
  scale_shape_manual(values = c(15,16,17,18), guide = guide_legend(title = "Treatment")) +
  stat_cor()+ #for p-values and Rho
  stat_regline_equation(label.y=500)+# regresion line equation will be shown
  labs(x="Total nitrogen 2015", y="Deviation from expected yield (ΔY) 2016")
  #facet_wrap(~treatment2)

soil_comp.sum<-soil_comp %>% group_by(Year, subplot, DOY, Depth, treatment2)%>% 
  summarize(TN.m=mean(TN), NH4.m=mean(NH4), NO3.m=mean(NO3.NO2), SOC.m=mean(SOC),
            TN.se=calcSE(TN),NH4.se=calcSE(NH4), NO3.se=calcSE(NO3.NO2), SOC.se=calcSE(SOC))

nfig<-ggplot(subset(soil_comp.sum, Depth==1), aes(x=treatment2, color=subplot))+
  geom_point(aes(y=TN.m),size=4, position=position_dodge(0.5))+
  geom_errorbar(aes(ymin=TN.m-TN.se,ymax=TN.m+TN.se),
                position = position_dodge(width = 0.5), width=0)+
  #theme_linedraw()+
  #scale_fill_manual(values = c("royalblue2","sienna","lightsteelblue1", "peachpuff"), guide = guide_legend(title = "Treatment")) +
  scale_color_manual(values = c("blue","black","gray"), guide = guide_legend(title = "Composition"),
                    labels=c("Grass-only", "Forb-only", "Mixed")) +
  #labs(x="Drought Treatment", y="NH4 (ug/g)")+
  #geom_boxplot()+
  facet_wrap(~Year)
nfig<-nfig + geom_point(aes(y=NH4.m),shape=17, position=position_dodge(0.5))  
nfig<-nfig +  geom_point(aes(y=NO3.m),shape=15, position=position_dodge(0.5))
nfig
