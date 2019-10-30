library(tidyverse)
library(nlme)
library(ggplot2)
library(multcomp)
library(lsmeans)
library(MuMIn)
library(RColorBrewer)
library(gridExtra)
setwd("~/Dropbox/ClimVar/DATA/")
FG_ANPP<-read.csv("Plant_composition_data/ANPP/ANPP_CleanedData/ClimVar_ANPP-func-groups.csv")
levels(FG_ANPP$plot)
levels(FG_ANPP$year)
levels(FG_ANPP$func) #note FG names are not consistent
levels(FG_ANPP$subplot)
levels(FG_ANPP$treatment)
levels(FG_ANPP$shelter)
levels(FG_ANPP$shelterBlock)
head(FG_ANPP)

#change plots, years, shelter to factors
FG_ANPP[,'plot'] <- as.factor(as.character(FG_ANPP[,'plot']))
FG_ANPP[,'year'] <- as.factor(as.character(FG_ANPP[,'year']))
FG_ANPP[,'shelter'] <- as.factor(as.character(FG_ANPP[,'shelter']))

#rename func groups
FG_ANPP$func[FG_ANPP$func=="N"] <- "N-fixer"
FG_ANPP$func[FG_ANPP$func=="G"] <- "Grass"
FG_ANPP$func[FG_ANPP$func=="F"] <- "Forb"

##1. Does variability of rainfall affect forage production (H1)?
#Using a mixed model to examine the effects of rainfall timing (fixed effect) with shelterblock nested within year as random effect

#create subset with no species manipulations (control community) only
FG_XC<-filter(FG_ANPP, harvest=="Second", subplot=="XC") %>% dplyr::select(-harvest) %>% arrange(treatment,shelterBlock)

m1<-lme(log(weight_g_m+1) ~treatment*func, random=~1|year/shelterBlock, FG_XC, na.action=na.exclude)
summary(m1)
anova(m1)
r.squaredGLMM(m1) #73% of variation explained by fixed effects, 77% by whole model
qqnorm(residuals(m1))
qqline(residuals(m1))
shapiro.test(residuals(m1))
#normally distributed, continue
LS1<-lsmeans(m1, ~treatment*func)
contrast(LS1, "pairwise")
#control rain ANPP is significantly greater than all the drought except spring dry, no surprise
#control rain is most similar to spring dry
#let's see it
ggplot(d=FG_XC, aes(x=treatment, y=weight_g_m, fill=func)) +
  theme_linedraw()+
  facet_wrap(~subplot*year, ncol=3)+
  labs(x="rainfall treatment", y="ANPP g/m2")+
  #annotate("text", x= c("consistentDry", "controlRain","fallDry","springDry"), y = c(900, 975, 900,900), label = c("a", "b", "ab", "ab"), color = "#999999") +
  geom_boxplot(aes(y=weight_g_m), shape=16)

FG_XC$treatment2 <- as.character(FG_XC$treatment)
#Then turn it back into a factor with the levels in the correct order
FG_XC$treatment2 <- factor(FG_XC$treatment2, levels = c("controlRain", "springDry","fallDry","consistentDry"))
ggplot(d=FG_XC, aes(x=treatment2, y=weight_g_m, fill=func)) +
  theme_linedraw()+
  labs(x="rainfall treatment", y="ANPP g/m2")+
  #annotate("text", x= c("controlRain","consistentDry", "fallDry","springDry"), y = c( 975, 900, 900,900), label = c("a", "b", "ab", "ab"), color = "#999999") +
  geom_boxplot(aes(y=weight_g_m), shape=16)

#let's do this again but test functional groups individually
##GRASS##
#create subset with no species manipulations (control community) only
FG_XC_grass<-filter(FG_ANPP, subplot=='XC', harvest=="Second", func=="Grass") %>% dplyr::select(-harvest) %>% arrange(treatment,shelterBlock)

m2<-lme(weight_g_m ~treatment*year, random=~1|shelterBlock, FG_XC_grass, na.action=na.exclude)
summary(m2)
anova(m2) #treatment is significant
r.squaredGLMM(m2) #10% of variation explained by fixed effects, 59% by whole model (interannual variation?)
qqnorm(residuals(m2))
qqline(residuals(m2))
shapiro.test(residuals(m2))
#normally distributed, continue
LS2<-lsmeans(m2, ~treatment*year)
contrast(LS2, "pairwise")

FG_XC_grass$treatment2 <- as.character(FG_XC_grass$treatment)
#Then turn it back into a factor with the levels in the correct order
FG_XC_grass$treatment2 <- factor(FG_XC_grass$treatment2, levels = c("controlRain", "fallDry","springDry","consistentDry"))

fg_G<-ggplot(FG_XC_grass, aes(x = treatment2, y = weight_g_m)) + 
  geom_boxplot(aes(fill = treatment2)) + 
  scale_fill_manual(values = c("royalblue2", "peachpuff","lightsteelblue1", "sienna"), guide = guide_legend(title = "Treatment")) +
  geom_jitter(position=position_jitter(0.2), aes(color=year)) +
  #geom_jitter(position=position_jitter(0.2), aes(color=May_XC$shelterBlock, shape=as.factor(year))) +
  scale_color_manual(values = c("gray80","gray50", "black"), guide = guide_legend(title = "Year")) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")+
  theme(axis.text.x=element_text(angle=90, hjust=1))+
  #facet_wrap(~year)+
  annotate("text", x= c("controlRain","consistentDry", "fallDry","springDry"), y = c( 990, 900, 900,900), label = c("a", "b", "b", "ab"), color = "black") +
  xlab("Rainfall Treatment") +
  ylab("Grass ANPP g/m2")
fg_G

ggplot(d=FG_XC_grass, aes(x=year, y=weight_g_m)) +
  geom_boxplot(aes(y=weight_g_m, fill=year), shape=16)+
  scale_fill_manual(values = c("gray99","gray80", "gray50"), guide = guide_legend(title = "Year")) +
  geom_jitter(position=position_jitter(0.2), aes(color=treatment2)) +
  scale_color_manual(values = c("royalblue2", "peachpuff","lightsteelblue1", "sienna"), guide = guide_legend(title = "Treatment")) +
  labs(x="", y="Grass ANPP g/m2")+
  #annotate("text", x= c("2015", "2016","2017"), y = c(900, 975, 975), label = c("a", "b", "b"), color = "black") +
  theme_linedraw()

##FORB##
#create subset with no species manipulations (control community) only
FG_XC_forb<-filter(FG_ANPP, subplot=='XC', harvest=="Second", func=="Forb") %>% dplyr::select(-harvest)

m3<-lme(log(weight_g_m+1) ~treatment, random=~1|year/shelterBlock, FG_XC_forb, na.action=na.exclude)
summary(m3)
anova(m3)#no treatment effect
r.squaredGLMM(m3) #<1% of variation explained by fixed effects, 35% by whole model (interannual variation?)
qqnorm(residuals(m3))
qqline(residuals(m3))
shapiro.test(residuals(m3))
#normally distributed, continue
LS3<-lsmeans(m3, ~treatment)
contrast(LS3, "pairwise") #no differences

FG_XC_forb$treatment2 <- as.character(FG_XC_forb$treatment)
#Then turn it back into a factor with the levels in the correct order
FG_XC_forb$treatment2 <- factor(FG_XC_forb$treatment2, levels = c("controlRain", "fallDry", "springDry", "consistentDry"))

fg_F<-ggplot(FG_XC_forb, aes(x = treatment2, y = weight_g_m)) + 
  geom_boxplot(aes(fill = treatment2)) + 
  scale_fill_manual(values = c("royalblue2", "peachpuff","lightsteelblue1", "sienna"), guide = guide_legend(title = "Treatment")) +
  geom_jitter(position=position_jitter(0.2), aes(color=year)) +
  #geom_jitter(position=position_jitter(0.2), aes(color=May_XC$shelterBlock, shape=as.factor(year))) +
  scale_color_manual(values = c("gray80","gray50", "black"), guide = guide_legend(title = "Year")) +
  theme_bw(base_size = 14) +
  theme(axis.text.x=element_text(angle=90, hjust=1))+
  theme(legend.position = "none")+
  #facet_wrap(~year)+
  annotate("text", x= c("controlRain","consistentDry", "fallDry","springDry"), y = c( 175, 175, 175, 175), label = c("a", "a", "a", "a"), color = "black") +
  xlab("Rainfall Treatment") +
  ylab("Forb ANPP g/m2")
fg_F

ggplot(d=FG_XC_forb, aes(x=year, y=weight_g_m)) +
  geom_boxplot(aes(y=weight_g_m, fill=year), shape=16)+
  scale_fill_manual(values = c("gray99","gray80", "gray50"), guide = guide_legend(title = "Year")) +
  geom_jitter(position=position_jitter(0.2), aes(color=treatment2)) +
  scale_color_manual(values = c("royalblue2", "peachpuff","lightsteelblue1", "sienna"), guide = guide_legend(title = "Treatment")) +
  labs(x="", y="Forb ANPP g/m2")+
  #annotate("text", x= c("2015", "2016","2017"), y = c(900, 975, 975), label = c("a", "b", "b"), color = "black") +
  theme_linedraw()

##N-FIXER##
#create subset with no species manipulations (control community) only
FG_XC_nfix<-filter(FG_ANPP, subplot=='XC', harvest=="Second", func=="N-fixer") %>% dplyr::select(-harvest)

m4<-lme(log(weight_g_m+1) ~treatment, random=~1|year/shelterBlock, FG_XC_nfix, na.action=na.exclude)
summary(m4)
anova(m4) #treatment is significant
r.squaredGLMM(m4) #30% of variation explained by fixed effects, 39% by whole model (interannual variation?)
qqnorm(residuals(m4))
qqline(residuals(m4))
shapiro.test(residuals(m4))
#normally distributed, continue
LS4<-lsmeans(m4, ~treatment)
contrast(LS4, "pairwise") #consistent drought is different from fall and spring, but not control rain

FG_XC_nfix$treatment2 <- as.character(FG_XC_nfix$treatment)
#Then turn it back into a factor with the levels in the correct order
FG_XC_nfix$treatment2 <- factor(FG_XC_nfix$treatment2, levels = c("controlRain", "fallDry", "springDry","consistentDry"))

fg_N<-ggplot(FG_XC_nfix, aes(x = treatment2, y = weight_g_m)) + 
  geom_boxplot(aes(fill = treatment2)) + 
  scale_fill_manual(values = c("royalblue2", "peachpuff", "lightsteelblue1", "sienna"), guide = guide_legend(title = "Treatment")) +
  geom_jitter(position=position_jitter(0.2), aes(color=year)) +
  #geom_jitter(position=position_jitter(0.2), aes(color=May_XC$shelterBlock, shape=as.factor(year))) +
  scale_color_manual(values = c("gray80","gray50", "black"), guide = guide_legend(title = "Year")) +
  theme_bw(base_size = 14) +
  theme(axis.text.x=element_text(angle=90, hjust=1))+
  #theme(legend.position = "none")+
  #facet_wrap(~year)+
  annotate("text", x= c("controlRain","consistentDry", "fallDry","springDry"), y = c( 175, 175, 175, 175), label = c("ab", "b", "a", "a"), color = "black") +
  xlab("Rainfall Treatment") +
  ylab("N-Fixer ANPP g/m2")
fg_N

ggplot(d=FG_XC_nfix, aes(x=year, y=weight_g_m)) +
  geom_boxplot(aes(y=weight_g_m, fill=year), shape=16)+
  scale_fill_manual(values = c("gray99","gray80", "gray50"), guide = guide_legend(title = "Year")) +
  geom_jitter(position=position_jitter(0.2), aes(color=treatment2)) +
  scale_color_manual(values = c("royalblue2", "peachpuff", "lightsteelblue1","sienna"), guide = guide_legend(title = "Treatment")) +
  labs(x="", y="Forb ANPP g/m2")+
  #annotate("text", x= c("2015", "2016","2017"), y = c(900, 975, 975), label = c("a", "b", "b"), color = "black") +
  theme_linedraw()

#compile FG plots 
grid.arrange(fg_G, fg_F, fg_N, ncol = 3, widths = c(1.1,1.1,1.5))

#create subset with no species manipulations (control community) only for both harvests
FG_XC<-filter(FG_ANPP) 

fgabund <- FG_ANPP %>% filter(subplot=="XC")%>%
  group_by (func, harvest) %>%
  summarize(tot=mean(weight_g_m), setot=sd(weight_g_m)/sqrt(length(weight_g_m))) %>%
  tbl_df() 

ggplot(data=fgabund, aes(x=harvest, y=tot, fill=harvest))+
  geom_bar(stat="identity")+
  geom_errorbar(aes(ymax = tot+setot, ymin = tot-setot), width=.25) + 
  facet_wrap(~func)

m1<-lme(log(weight_g_m+1) ~treatment*func*harvest, random=~1|year/shelterBlock, FG_XC, na.action=na.exclude)
summary(m1)
anova(m1)
r.squaredGLMM(m1) #73% of variation explained by fixed effects, 77% by whole model
qqnorm(residuals(m1))
qqline(residuals(m1))
shapiro.test(residuals(m1))
#normally distributed, continue
LS1<-lsmeans(m1, ~treatment)
contrast(LS1, "pairwise")
#control rain ANPP is significantly greater than all the drought except spring dry, no surprise
#control rain is most similar to spring dry
#let's see it
ggplot(d=FG_XC, aes(x=treatment, y=weight_g_m, fill=func, color=harvest)) +
  theme_linedraw()+
  facet_wrap(~subplot*year, ncol=3)+
  labs(x="rainfall treatment", y="ANPP g/m2")+
  #annotate("text", x= c("consistentDry", "controlRain","fallDry","springDry"), y = c(900, 975, 900,900), label = c("a", "b", "ab", "ab"), color = "#999999") +
  geom_boxplot(aes(y=weight_g_m), shape=16)

FG_XC$treatment2 <- as.character(FG_XC$treatment)
#Then turn it back into a factor with the levels in the correct order
FG_XC$treatment2 <- factor(FG_XC$treatment2, levels = c("controlRain","fallDry", "springDry","consistentDry"))
ggplot(d=FG_XC, aes(x=treatment2, y=weight_g_m, fill=func)) +
  theme_linedraw()+
  labs(x="rainfall treatment", y="ANPP g/m2")+
  #annotate("text", x= c("controlRain","consistentDry", "fallDry","springDry"), y = c( 975, 900, 900,900), label = c("a", "b", "ab", "ab"), color = "#999999") +
  geom_boxplot(aes(y=weight_g_m), shape=16)

#let's do this again but test functional groups individually
##GRASS##
#create subset with no species manipulations (control community) only, 2nd harvest is when grass is most abundant
FG_XC_grass<-filter(FG_ANPP, subplot=='XC', harvest=="First", func=="Grass") %>% dplyr::select(-harvest)

m2<-lme(weight_g_m ~treatment, random=~1|year/shelterBlock, FG_XC_grass, na.action=na.exclude)
summary(m2)
anova(m2) #treatment is significant
r.squaredGLMM(m2) #10% of variation explained by fixed effects, 59% by whole model (interannual variation?)
qqnorm(residuals(m2))
qqline(residuals(m2))
shapiro.test(residuals(m2))
#normally distributed, continue
LS2<-lsmeans(m2, ~treatment)
contrast(LS2, "pairwise")

FG_XC_grass$treatment2 <- as.character(FG_XC_grass$treatment)
#Then turn it back into a factor with the levels in the correct order
FG_XC_grass$treatment2 <- factor(FG_XC_grass$treatment2, levels = c("controlRain","springDry", "fallDry","consistentDry"))

fg_G<-ggplot(FG_XC_grass, aes(x = treatment2, y = weight_g_m)) + 
  geom_boxplot(aes(fill = treatment2)) + 
  scale_fill_manual(values = c("royalblue2","lightsteelblue1", "peachpuff", "sienna"), guide = guide_legend(title = "Treatment")) +
  geom_jitter(position=position_jitter(0.2), aes(color=year)) +
  #geom_jitter(position=position_jitter(0.2), aes(color=May_XC$shelterBlock, shape=as.factor(year))) +
  scale_color_manual(values = c("gray80","gray50", "black"), guide = guide_legend(title = "Year")) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")+
  theme(axis.text.x=element_text(angle=90, hjust=1))+
  #facet_wrap(~year)+
  annotate("text", x= c("controlRain","consistentDry", "fallDry","springDry"), y = c( 990, 920, 920,920), label = c("a", "b", "ab", "ab"), color = "black") +
  xlab("Rainfall Treatment") +
  ylab("Grass ANPP g/m2")

ggplot(d=FG_XC_grass, aes(x=year, y=weight_g_m)) +
  geom_boxplot(aes(y=weight_g_m, fill=year), shape=16)+
  scale_fill_manual(values = c("gray99","gray80", "gray50"), guide = guide_legend(title = "Year")) +
  geom_jitter(position=position_jitter(0.2), aes(color=treatment2)) +
  scale_color_manual(values = c("royalblue2","lightsteelblue1", "peachpuff", "sienna"), guide = guide_legend(title = "Treatment")) +
  labs(x="", y="Grass ANPP g/m2")+
  #annotate("text", x= c("2015", "2016","2017"), y = c(900, 975, 975), label = c("a", "b", "b"), color = "black") +
  theme_linedraw()

##FORB##
#create subset with no species manipulations (control community) only, first harvest b/c forbs most abundant then
FG_XC_forb<-filter(FG_ANPP, subplot=='XC', harvest=="First", func=="Forb") %>% dplyr::select(-harvest)

m3<-lme(log(weight_g_m+1) ~treatment, random=~1|year/shelterBlock, FG_XC_forb, na.action=na.exclude)
summary(m3)
anova(m3)#treatment significant
r.squaredGLMM(m3) #15% of variation explained by fixed effects, 50% by whole model (interannual variation?)
qqnorm(residuals(m3))
qqline(residuals(m3))
shapiro.test(residuals(m3))
#normally distributed, continue
LS3<-lsmeans(m3, ~treatment)
contrast(LS3, "pairwise") #no differences

FG_XC_forb$treatment2 <- as.character(FG_XC_forb$treatment)
#Then turn it back into a factor with the levels in the correct order
FG_XC_forb$treatment2 <- factor(FG_XC_forb$treatment2, levels = c("controlRain","springDry", "fallDry","consistentDry"))

fg_F<-ggplot(FG_XC_forb, aes(x = treatment2, y = weight_g_m)) + 
  geom_boxplot(aes(fill = treatment2)) + 
  scale_fill_manual(values = c("royalblue2","lightsteelblue1", "peachpuff", "sienna"), guide = guide_legend(title = "Treatment")) +
  geom_jitter(position=position_jitter(0.2), aes(color=year)) +
  #geom_jitter(position=position_jitter(0.2), aes(color=May_XC$shelterBlock, shape=as.factor(year))) +
  scale_color_manual(values = c("gray80","gray50", "black"), guide = guide_legend(title = "Year")) +
  theme_bw(base_size = 14) +
  theme(axis.text.x=element_text(angle=90, hjust=1))+
  theme(legend.position = "none")+
  #facet_wrap(~year)+
  annotate("text", x= c("controlRain","consistentDry", "fallDry","springDry"), y = c( 185, 185, 185, 185), label = c("a", "a", "b", "a"), color = "black") +
  xlab("Rainfall Treatment") +
  ylab("Forb ANPP g/m2")

ggplot(d=FG_XC_forb, aes(x=year, y=weight_g_m)) +
  geom_boxplot(aes(y=weight_g_m, fill=year), shape=16)+
  scale_fill_manual(values = c("gray99","gray80", "gray50"), guide = guide_legend(title = "Year")) +
  geom_jitter(position=position_jitter(0.2), aes(color=treatment2)) +
  scale_color_manual(values = c("royalblue2","lightsteelblue1", "peachpuff", "sienna"), guide = guide_legend(title = "Treatment")) +
  labs(x="", y="Forb ANPP g/m2")+
  #annotate("text", x= c("2015", "2016","2017"), y = c(900, 975, 975), label = c("a", "b", "b"), color = "black") +
  theme_linedraw()

##N-FIXER##
#create subset with no species manipulations (control community) only, first harvest bc n-fixers most abundant then
FG_XC_nfix<-filter(FG_ANPP, subplot=='XC', harvest=="First", func=="N-fixer") %>% dplyr::select(-harvest)

m4<-lme(log(weight_g_m+1) ~treatment, random=~1|year/shelterBlock, FG_XC_nfix, na.action=na.exclude)
summary(m4)
anova(m4) #treatment is not significant
r.squaredGLMM(m4) #30% of variation explained by fixed effects, 39% by whole model (interannual variation?)
qqnorm(residuals(m4))
qqline(residuals(m4))
shapiro.test(residuals(m4))
#normally distributed, continue
LS4<-lsmeans(m4, ~treatment)
contrast(LS4, "pairwise") #consistent drought is different from fall and spring, but not control rain

FG_XC_nfix$treatment2 <- as.character(FG_XC_nfix$treatment)
#Then turn it back into a factor with the levels in the correct order
FG_XC_nfix$treatment2 <- factor(FG_XC_nfix$treatment2, levels = c("controlRain","springDry", "fallDry","consistentDry"))

fg_N<-ggplot(FG_XC_nfix, aes(x = treatment2, y = weight_g_m)) + 
  geom_boxplot(aes(fill = treatment2)) + 
  scale_fill_manual(values = c("royalblue2","lightsteelblue1", "peachpuff", "sienna"), guide = guide_legend(title = "Treatment")) +
  geom_jitter(position=position_jitter(0.2), aes(color=year)) +
  #geom_jitter(position=position_jitter(0.2), aes(color=May_XC$shelterBlock, shape=as.factor(year))) +
  scale_color_manual(values = c("gray80","gray50", "black"), guide = guide_legend(title = "Year")) +
  theme_bw(base_size = 14) +
  theme(axis.text.x=element_text(angle=90, hjust=1))+
  #theme(legend.position = "none")+
  #facet_wrap(~year)+
  annotate("text", x= c("controlRain","consistentDry", "fallDry","springDry"), y = c( 185, 185, 185, 185), label = c("a", "a", "a", "a"), color = "black") +
  xlab("Rainfall Treatment") +
  ylab("N-Fixer ANPP g/m2")

ggplot(d=FG_XC_nfix, aes(x=year, y=weight_g_m)) +
  geom_boxplot(aes(y=weight_g_m, fill=year), shape=16)+
  scale_fill_manual(values = c("gray99","gray80", "gray50"), guide = guide_legend(title = "Year")) +
  geom_jitter(position=position_jitter(0.2), aes(color=treatment2)) +
  scale_color_manual(values = c("royalblue2","lightsteelblue1", "peachpuff", "sienna"), guide = guide_legend(title = "Treatment")) +
  labs(x="", y="Forb ANPP g/m2")+
  #annotate("text", x= c("2015", "2016","2017"), y = c(900, 975, 975), label = c("a", "b", "b"), color = "black") +
  theme_linedraw()

#compile FG plots 
grid.arrange(fg_G, fg_F, fg_N, ncol = 3, widths = c(1.1,1.1,1.7))
grid.arrange(fg_G, fg_F, ncol = 2)

#total ANPP at each harvest for XC
fgabund2 <- FG_ANPP %>% filter(subplot=="XC")%>%
  group_by (plot, year, harvest, treatment, shelterBlock) %>%
  summarize(sum=sum(weight_g_m)) %>%
  tbl_df() 

m_tot<-lme(sum ~treatment, random=~1|year/shelterBlock, subset(fgabund2, harvest=="First"), na.action=na.exclude)
summary(m_tot)
anova(m_tot)#treatment significant
r.squaredGLMM(m_tot) #12% of variation explained by fixed effects, 50% by whole model (interannual variation?)
qqnorm(residuals(m_tot))
qqline(residuals(m_tot))
shapiro.test(residuals(m_tot))
#normally distributed, continue
LS3<-lsmeans(m_tot, ~treatment)
contrast(LS3, "pairwise") #no differences

fgabund2$treatment2 <- as.character(fgabund2$treatment)
#Then turn it back into a factor with the levels in the correct order
fgabund2$treatment2 <- factor(fgabund2$treatment2, levels = c("controlRain","springDry", "fallDry","consistentDry"))
ggplot(subset(fgabund2, harvest=="First"), aes(x = treatment2, y = sum)) + 
  geom_boxplot(aes(fill = treatment2)) + 
  scale_fill_manual(values = c("royalblue2","lightsteelblue1", "peachpuff", "sienna"), guide = guide_legend(title = "Treatment")) +
  geom_jitter(position=position_jitter(0.2), aes(color=year)) +
  #geom_jitter(position=position_jitter(0.2), aes(color=May_XC$shelterBlock, shape=as.factor(year))) +
  scale_color_manual(values = c("gray80","gray50", "black"), guide = guide_legend(title = "Year")) +
  theme_bw(base_size = 14) +
  #facet_wrap(~year)+
  annotate("text", x= c("controlRain","consistentDry", "fallDry","springDry"), y = c( 1050, 1050, 1050,1050), label = c("a", "b", "ab", "ab"), color = "black") +
  xlab("Rainfall Treatment") +
  ylab("ANPP g/m2")+
  ylim(200,1100)

m_tot2<-lme(sum ~treatment, random=~1|year/shelterBlock, subset(fgabund2, harvest=="Second"), na.action=na.exclude)
summary(m_tot2)
anova(m_tot2)#treatment significant
r.squaredGLMM(m_tot2) #12% of variation explained by fixed effects, 50% by whole model (interannual variation?)
qqnorm(residuals(m_tot2))
qqline(residuals(m_tot2))
shapiro.test(residuals(m_tot2))
#normally distributed, continue
LS3<-lsmeans(m_tot2, ~treatment)
contrast(LS3, "pairwise") #no differences

ggplot(subset(fgabund2, harvest=="Second"), aes(x = treatment2, y = sum)) + 
  geom_boxplot(aes(fill = treatment2)) + 
  scale_fill_manual(values = c("royalblue2","lightsteelblue1", "peachpuff", "sienna"), guide = guide_legend(title = "Treatment")) +
  geom_jitter(position=position_jitter(0.2), aes(color=year)) +
  #geom_jitter(position=position_jitter(0.2), aes(color=May_XC$shelterBlock, shape=as.factor(year))) +
  scale_color_manual(values = c("gray80","gray50", "black"), guide = guide_legend(title = "Year")) +
  theme_bw(base_size = 14) +
  #facet_wrap(~year)+
  annotate("text", x= c("controlRain","consistentDry", "fallDry","springDry"), y = c( 1050, 1050, 1050,1050), label = c("a", "b", "b", "ab"), color = "black") +
  xlab("Rainfall Treatment") +
  ylab("ANPP g/m2")+
  ylim(200,1100)
