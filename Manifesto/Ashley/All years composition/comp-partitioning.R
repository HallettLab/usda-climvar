#NOTE: mistakes in partitioning calculations here, corrected in comp_part2.R which is the code for the manuscript
##the purpose of this code is additive partitioning of biodiversity effects in the ClimVar mixed composition plots
#code below to be deleted due to errors in calculatations, but leaving up for now so steps of what was done are clear 
##Author: Ashley Shaw

library(tidyverse)
library(nlme)
library(ggplot2)
library(multcomp)
library(lsmeans)
library(MuMIn)
library(RColorBrewer)
library(ggpubr)
theme_set(theme_bw())
theme_update( panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
              strip.background = element_blank(),
              text = element_text(size = 14),
              strip.text= element_text(size = 12), 
              axis.text = element_text(size = 12))

setwd("~/Dropbox/ClimVar/DATA/")
May_ANPP<-read.csv("Plant_composition_data/ANPP/ANPP_CleanedData/ClimVar_ANPP-peak.csv")

head(May_ANPP)

calcSE<-function(x){
  x <- x[!is.na(x)]
  sd(x)/sqrt(length(x))
}

#change plots, years, shelter to factors
May_ANPP[,'plot'] <- as.factor(as.character(May_ANPP[,'plot']))
May_ANPP[,'year'] <- as.factor(as.character(May_ANPP[,'year']))
May_ANPP[,'shelter'] <- as.factor(as.character(May_ANPP[,'shelter']))

compensation (H2)
#H2 will be confirmed if mixed plots have greater ANPP compared to grass only or forb only
#first remove compost treatment 
May_ANPP_noC<-filter(May_ANPP, subplot!='C', subplot!='XC', year!=2017)
May_ANPP_noC <- May_ANPP_noC %>% mutate(treatment=ordered(treatment, levels = c(controlRain="controlRain", springDry="springDry", fallDry="fallDry", consistentDry="consistentDry")))

m3<-lme(weight_g_m ~treatment*subplot, random=~1|year/shelterBlock, May_ANPP_noC, na.action=na.exclude)
summary(m3)
anova(m3)
r.squaredGLMM(m3)
qqnorm(residuals(m3))
qqline(residuals(m3))
shapiro.test(residuals(m3))
LS3<-lsmeans(m3, ~treatment*subplot, by="treatment")
contrast(LS3, "pairwise")
#not normally distributed, try log transform

m4<-lme(log(weight_g_m+1) ~treatment*subplot*year, random=~1|shelterBlock, May_ANPP_noC, na.action=na.exclude)
summary(m4)
anova(m4)
r.squaredGLMM(m4)#7% of variation explained by fixed effects, 47% explained by entire model (lots of interannual variation?)
qqnorm(residuals(m4))
qqline(residuals(m4))
shapiro.test(residuals(m4))
#barely normal
LS4<-lsmeans(m4, ~subplot, by="year")
contrast(LS4, "pairwise")
#ANOVA: overall sign effect of subplot on ANPP
#lsmeans: B (mixed plots) have greater ANPP compared to F (forb-only), but not G (grass-only)
#no sig difference in ANPP between B (mixed plots) and XC (no manipulation)
#plot the effect of species composition on ANPP:
ggplot(d=May_ANPP_noC, aes(x=treatment, y=weight_g_m, fill=subplot)) +
  facet_wrap(~year)+
  #theme_linedraw()+
  labs(x="composition treatment", y="ANPP g/m2")+
  geom_boxplot(aes(y=weight_g_m), shape=16)

ggplot(d=May_ANPP_noC, aes(x=subplot, y=weight_g_m)) +
  #theme_linedraw()+
  labs(x="", y="ANPP g/m2")+
  annotate("text", x= c("B", "F","G"), y = c(1200, 1250, 1200), label = c("a", "b", "ab"), color = "black") +
  geom_boxplot(aes(y=weight_g_m), shape=16)


m3.b<-lme(weight_g_m ~treatment*year, random=~1|shelterBlock, subset(May_ANPP_noC, subplot=="B"), na.action=na.exclude)
summary(m3.b)
anova(m3.b)
r.squaredGLMM(m3.b)
qqnorm(residuals(m3.b))
qqline(residuals(m3.b))
shapiro.test(residuals(m3.b))
LS3<-lsmeans(m3.b, ~treatment*year, by="year")
contrast(LS3, "pairwise")

#rename trt groups
May_ANPP_noC$treatment2[May_ANPP_noC$treatment=="controlRain"] <- "Ambient"
May_ANPP_noC$treatment2[May_ANPP_noC$treatment=="springDry"] <- "Late"
May_ANPP_noC$treatment2[May_ANPP_noC$treatment=="fallDry"] <- "Early"
May_ANPP_noC$treatment2[May_ANPP_noC$treatment=="consistentDry"] <- "Consistent"

ggplot(d=May_ANPP_noC, aes(x=treatment2, y=weight_g_m, fill=treatment2)) +
  #theme_linedraw()+
  facet_wrap(~year)+
  labs(x="", y="ANPP g/m2")+
  scale_fill_manual(values = c("royalblue2","sienna","lightsteelblue1", "peachpuff"), guide = guide_legend(title = "Treatment")) +
  #annotate("text", x= c("B", "F","G"), y = c(1200, 1250, 1200), label = c("a", "b", "ab"), color = "black") +
  geom_boxplot(aes(y=weight_g_m), shape=16)

ggplot(d=subset(May_ANPP_noC, subplot=='B'), aes(x=treatment2, y=weight_g_m, fill=treatment2)) +
  facet_wrap(~year)+
  ggtitle("Observed Yield (Yo)")+
  #theme_linedraw()+
  theme(legend.position="none")+
  scale_fill_manual(values = c("royalblue2","sienna","lightsteelblue1", "peachpuff"), guide = guide_legend(title = "Treatment")) +
  labs(x="Drought Treatment", y="ANPP g/m2")+
  geom_boxplot(aes(y=weight_g_m), shape=16)

ggplot(d=subset(May_ANPP_noC), aes(x=treatment2, y=weight_g_m, fill=subplot)) +
  facet_wrap(~year)+
  ggtitle("a) Observed Yield (Yo)")+
  #theme_linedraw()+
  #theme(legend.position="none")+
  scale_fill_manual(values = c("gray","black","white"), guide = guide_legend(title = "Composition"),
                    labels=c("Mixed", "Forb-only", "Grass-only")) +
  labs(x="Drought Treatment", y="ANPP g/m2")+
  geom_boxplot(aes(y=weight_g_m), shape=16)



FG_ANPP<-read.csv("Plant_composition_data/ANPP/ANPP_CleanedData/ClimVar_ANPP-func-groups.csv")
FG_ANPP<-FG_ANPP %>% dplyr::select(-totweight,-propweight)%>%
  group_by(year, subplot, plot, treatment, date) %>%
  mutate(totweight = sum(weight_g_m), propweight = weight_g_m/totweight) %>%
  tbl_df() %>% ungroup()%>%
  mutate(func = ifelse(func == "F", "Forb", func),
         func = ifelse(func == "N", "N-fixer", func),
         func = ifelse(func == "G", "Grass", func)) %>%
  mutate(treatment=ordered(treatment, levels = c(consistentDry="consistentDry", springDry="springDry", fallDry="fallDry", controlRain="controlRain"))) 


##FORB##
#create subset with no species manipulations (control community) only, first harvest b/c forbs most abundant then
FG_forb<-filter(FG_ANPP, subplot=='F', harvest=="First", func=="Forb", year!="2017") %>% dplyr::select(-harvest, -func, -date) %>% 
  mutate(treatment=ordered(treatment, levels = c(controlRain="controlRain", springDry="springDry", fallDry="fallDry", consistentDry="consistentDry")))

FG_B.f<-filter(FG_ANPP, subplot=='B', harvest=="First", func=="Forb", year!="2017") %>% dplyr::select(-harvest, -func, -date)
FG_B.f$M<-FG_forb$weight_g_m
FG_B.f<-FG_B.f %>% group_by(treatment, shelterBlock, year)%>%mutate(RYof=propweight*M)%>%mutate(RYef=0.20*M)%>%mutate(dYf=RYof-RYef)%>%
  mutate(dRY.f=propweight-0.20)

##GRASS##
FG_grass<-filter(FG_ANPP, subplot=='G', harvest=="Second", func=="Grass", year!="2017") %>% dplyr::select(-harvest, -func, -date)%>% 
  mutate(treatment=ordered(treatment, levels = c(controlRain="controlRain", springDry="springDry", fallDry="fallDry", consistentDry="consistentDry")))
FG_B.g<-filter(FG_ANPP, subplot=='B', harvest=="Second", func=="Grass", year!="2017") %>% dplyr::select(-harvest, -func, -date)
FG_B.g$M<-FG_grass$weight_g_m
FG_B.g<-FG_B.g %>% group_by(treatment, shelterBlock, year)%>%mutate(RYog=propweight*M)%>%mutate(RYeg=0.80*M)%>%mutate(dYg=RYog-RYeg)%>%
  mutate(dRY.g=(propweight-0.80))

#calculate complementarity and selection effects
FG_B.g$dYf<-FG_B.f$dYf
FG_B.g$dRY.f<-FG_B.f$dRY.f
FG_B.g <- FG_B.g %>% ungroup() %>% mutate(complement=2*(((dRY.f+dRY.g))*((M+FG_B.f$M))))%>%mutate(x=dRY.f+dRY.g)%>% mutate(y=(M+FG_B.f$M))
FG_B.g<-FG_B.g%>%group_by(treatment,year)%>%mutate(select=2*cov(x,y))


FG_B.g<-FG_B.g %>% group_by(treatment, shelterBlock, year)%>%mutate(dY=dYg+dYf)%>% ungroup()%>%
  mutate(treatment=ordered(treatment, levels = c(controlRain="controlRain", springDry="springDry", fallDry="fallDry", consistentDry="consistentDry")))


#rename trt groups
FG_B.g$treatment2[FG_B.g$treatment=="controlRain"] <- "Ambient"
FG_B.g$treatment2[FG_B.g$treatment=="springDry"] <- "Late"
FG_B.g$treatment2[FG_B.g$treatment=="fallDry"] <- "Early"
FG_B.g$treatment2[FG_B.g$treatment=="consistentDry"] <- "Consistent"


#plot biomass of forb, grass, mixed
FG_B.g$Fweight<-FG_forb$weight_g_m
FG_B.g$Gweight<-FG_grass$weight_g_m
FG_all<-FG_B.g%>%gather("Fweight","Gweight","totweight",key=comp, value=biomass)

pcomp<-ggplot(d=FG_all, aes(x=comp, y=biomass, fill=comp)) +
  #theme_linedraw()+
  labs(x="", y=expression(ANPP~(g/m^2)))+
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  scale_fill_manual(values = c("black","white", "gray")) +
  theme(legend.position="none")+
  annotate("text", x= c("Fweight", "Gweight","totweight"), y = c(1000, 1000, 1000), label = c("a", "b", "c"), color = "black") +
  scale_x_discrete(labels=c("Fweight" = "Forb", "Gweight" = "Grass", "totweight" = "Mixed"))+
  geom_boxplot(aes(y=biomass), shape=16)
pcomp

pyr<-ggplot(d=subset(FG_all, comp=="totweight"), aes(x=as.factor(year), y=biomass)) +
  #theme_linedraw()+
  labs(x="", y=expression(ANPP~(g/m^2)))+
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  annotate("text", x= c("2015", "2016"), y = c(1000, 1000), label = c("a", "b"), color = "black") +
  #scale_x_discrete(labels=c("Fweight" = "Forb-only", "Gweight" = "Grass-only", "totweight" = "Mixed"))+
  geom_boxplot(aes(y=biomass), fill="gray")
pyr

bm<-lme(biomass ~treatment*comp*year, random=~1|shelterBlock, FG_all, na.action=na.exclude)
summary(bm)
anova(bm)
r.squaredGLMM(bm)
qqnorm(residuals(bm))
qqline(residuals(bm))
shapiro.test(residuals(bm))
LS3<-lsmeans(bm, ~treatment, by="comp")
contrast(LS3, "pairwise")
#not normally distributed, try log transform



#summarize for graphs by treatment within year
FG_B.sum_yr<- FG_B.g %>% group_by(treatment2, year)%>%
  summarize(sel.m=mean(select), sel.se=calcSE(select), 
         comp.m=mean(complement), comp.se=calcSE(complement),
         net.m=mean(dY), net.se=calcSE(dY))%>%
  gather(net.m, sel.m, comp.m,net.se,comp.se,sel.se, key="type", value="value")%>%
  separate(type, c("type", "x"))%>%
  spread(x,value)

#summarize for graphs by treatment
FG_B.sum<- FG_B.g %>% group_by(treatment2)%>%
  summarize(sel.m=mean(select), sel.se=calcSE(select), 
            comp.m=mean(complement), comp.se=calcSE(complement),
            net.m=mean(dY), net.se=calcSE(dY))%>%
  gather(net.m, sel.m, comp.m,net.se,comp.se,sel.se, key="type", value="value")%>%
  separate(type, c("type", "x"))%>%
  spread(x,value)

f1a<-ggplot(d=subset(FG_all), aes(x=treatment2, y=biomass, fill=comp)) +
  #facet_wrap(~year)+
  ggtitle(expression(Observed~Yield~(Y[o])))+
  #theme_linedraw()+
  #theme(legend.position="none")+
  scale_fill_manual(values = c("black","white", "gray"), guide = guide_legend(title = "Composition"),
      labels=c("Forb-only", "Grass-only", "Mixed")) +
  labs(x="", y=expression(ANPP~(g/m^2)))+
  geom_boxplot(aes(y=biomass), shape=16)+
  ylim(0,1500)
f1a

library(cowplot)

ggdraw(f1a) +
  draw_plot(pcomp, .15, .6, .3, .3) +
  draw_plot(pyr, 0.5,0.6, 0.3,0.3)+
  draw_plot_label(
    c("A", "B", "C"),
    c(0, 0.15, 0.5),
    c(1, 0.92, 0.92),
    size = 14
  )



f1a2<-ggplot(d=subset(FG_all), aes(x=treatment2, y=biomass, fill=comp)) +
  facet_wrap(~year)+
  #ggtitle(expression(Observed~Yield~(Y[o])))+
  #theme_linedraw()+
  #theme(legend.position="none")+
  scale_fill_manual(values = c("black","white", "gray"), guide = guide_legend(title = "Composition"),
                    labels=c("Forb-only", "Grass-only", "Mixed")) +
  labs(x="Drought Treatment", y=expression(ANPP~(g/m^2)))+
  geom_boxplot(aes(y=biomass), shape=16)
f1a2

#plot dY, deviation from expected yield
f1b<-ggplot(d=FG_B.g, aes(x=treatment2, y=dY, fill=treatment2)) +
  facet_wrap(~year)+
  ggtitle("b) Net Effect ")+
  scale_fill_manual(values = c("royalblue2","sienna","lightsteelblue1", "peachpuff"), guide = guide_legend(title = "Treatment")) +
  #theme_linedraw()+
  theme(legend.position="none")+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  labs(x="Drought Treatment", y="Net Effect (ΔY)")+
  geom_boxplot(aes(y=dY), shape=16)
f1b

ggarrange(f1a, f1b, ncol = 1, nrow = 2)

#test for treatment effects on selection
m.sel<-lme(select ~treatment, random=~1|shelterBlock, FG_B.g, na.action=na.exclude)
summary(m.sel)
anova(m.sel) #no sign. effects
r.squaredGLMM(m.sel) #67% of variation explained by fixed effects, 67% by whole model (interannual variation?)
qqnorm(residuals(m.sel))
qqline(residuals(m.sel))
shapiro.test(residuals(m.sel))
LS.sel<-lsmeans(m.sel, ~treatment)
contrast(LS.sel, "pairwise")

f2a<-ggplot(d=FG_B.g, aes(x=treatment2, y=select, fill=treatment2)) +
  #facet_wrap(~year)+
  ggtitle("a) Selection effect")+
  scale_fill_manual(values = c("royalblue2","sienna","lightsteelblue1", "peachpuff"), guide = guide_legend(title = "Treatment")) +
  #theme_classic()+
  theme(legend.position="none")+
  labs(x="Drought treatment", y="Selection Effect")+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  geom_boxplot(aes(y=select), shape=16)
f2a

#test for treatment effects on complementarity
m.comp<-lme(complement ~treatment*year, random=~1|shelterBlock, FG_B.g, na.action=na.exclude)
summary(m.comp)
anova(m.comp) #treatment is significant
r.squaredGLMM(m.comp) #36% of variation explained by fixed effects, 65% by whole model (interannual variation?)
qqnorm(residuals(m.comp))
qqline(residuals(m.comp))
shapiro.test(residuals(m.comp))
LS.all<-lsmeans(m.comp, ~treatment, by="year")
contrast(LS.all, "pairwise")

f2b<-ggplot(d=FG_B.g, aes(x=treatment2, y=complement, fill=treatment2)) +
  #facet_wrap(~year)+
  ggtitle("b) Complementarity Effect")+
  scale_fill_manual(values = c("royalblue2","sienna","lightsteelblue1", "peachpuff"), guide = guide_legend(title = "Treatment")) +
  #theme_classic()+
  #annotate("text", x= c("Ambient", "Consistent Drought","Early Drought", "Late Drought"), y = c(50,50,50,50), label = c("a", "b", "a", "a"), color = "black") +
  theme(legend.position="none")+
  labs(x="Drought Treatment", y="Complementarity Effect")+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  geom_boxplot(aes(y=complement), shape=16)
f2b

f2v2<-ggplot(FG_B.sum, aes(x=treatment2, y=m))+
  geom_bar(aes(fill=type), stat="identity", position = "dodge")+
  geom_errorbar(aes(ymin=m-se,ymax=m+se, group=type), position="dodge")
f2v2

f2v3<-ggplot(FG_B.sum_yr, aes(x=treatment2, y=m))+
  facet_wrap(~year)+
  geom_bar(aes(fill=type), stat="identity", position = "dodge")+
  geom_errorbar(aes(ymin=m-se,ymax=m+se, group=type), position="dodge")
f2v3

ggarrange(f2a,f2b,ncol=2)

ggplot(d=FG_B.g, aes(x=treatment2, y=complement, fill=treatment2)) +
  facet_wrap(~year)+
  #ggtitle("b) Complementarity Effect")+
  scale_fill_manual(values = c("royalblue2","sienna","lightsteelblue1", "peachpuff"), guide = guide_legend(title = "Treatment")) +
  #theme_classic()+
  #annotate("text", x= c("Ambient", "Consistent Drought","Early Drought", "Late Drought"), y = c(50,50,50,50), label = c("a", "b", "a", "a"), color = "black") +
  theme(legend.position="none")+
  labs(x="Drought Treatment", y="Complementarity Effect")+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  geom_boxplot(aes(y=complement), shape=16)

ggplot(d=FG_B.g, aes(x=treatment2, y=select, fill=treatment2)) +
  facet_wrap(~year)+
  #ggtitle("b) Complementarity Effect")+
  scale_fill_manual(values = c("royalblue2","sienna","lightsteelblue1", "peachpuff"), guide = guide_legend(title = "Treatment")) +
  #theme_classic()+
  #annotate("text", x= c("Ambient", "Consistent Drought","Early Drought", "Late Drought"), y = c(50,50,50,50), label = c("a", "b", "a", "a"), color = "black") +
  theme(legend.position="none")+
  labs(x="Drought Treatment", y="Selection Effect")+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  geom_boxplot(aes(y=select), shape=16)

#do forbs do better or worse than expected?
m.f<-lme(dYf ~treatment*year, random=~1|shelterBlock, FG_B.g, na.action=na.exclude)
summary(m.f)
anova(m.f) #treatment is not significant
r.squaredGLMM(m.f) #36% of variation explained by fixed effects, 65% by whole model (interannual variation?)
qqnorm(residuals(m.f))
qqline(residuals(m.f))
shapiro.test(residuals(m.f))
LS.f<-lsmeans(m.f, ~year)
contrast(LS.f, "pairwise")

f4b<-ggplot(d=FG_B.g, aes(x=treatment2, y=dYf, fill=treatment2)) +
  facet_wrap(~year,ncol=1)+
  ylim(-80,200)+
  ggtitle("Forbs")+
  scale_fill_manual(values = c("royalblue2","sienna","lightsteelblue1", "peachpuff"), guide = guide_legend(title = "Treatment")) +
  #theme_linedraw()+
  theme(legend.position="none")+
  labs(x="Rainfall Treatment", y="Deviation from Expected Yield (ΔYf)")+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  geom_boxplot(aes(y=dYf), shape=16)
f4b

#do grasses do better or worse than expected?
m.g<-lme(dYg ~treatment*year, random=~1|shelterBlock, FG_B.g, na.action=na.exclude)
summary(m.g)
anova(m.g) #treatment is not significant
r.squaredGLMM(m.g) #36% of variation explained by fixed effects, 65% by whole model (interannual variation?)
qqnorm(residuals(m.g))
qqline(residuals(m.g))
shapiro.test(residuals(m.g))
LS.g<-lsmeans(m.g, ~year)
LS.g2<-emmeans(m.g, "treatment")
plot(LS.g2, comparisons=T)
pairs(LS.g2)
contrast(LS.g, "pairwise")

f4a<-ggplot(d=FG_B.g, aes(x=treatment2, y=dYg, fill=treatment2)) +
  facet_wrap(~year,ncol=1)+
  ylim(-80,200)+
  ggtitle("Grasses")+
  scale_fill_manual(values = c("royalblue2","sienna","lightsteelblue1", "peachpuff"), guide = guide_legend(title = "Treatment")) +
  #theme_linedraw()+
  theme(legend.position="none")+
  labs(x="Rainfall Treatment", y="Deviation from Expected Yield (ΔYg)")+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  geom_boxplot(aes(y=dYg), shape=16)
f4a

ggarrange(f4a,f4b,labels=c("a)","b)"))

#rename trt groups
gf_graphic2$treatment2[gf_graphic2$treatment=="controlRain"] <- "Ambient"
gf_graphic2$treatment2[gf_graphic2$treatment=="springDry"] <- "Late Drought"
gf_graphic2$treatment2[gf_graphic2$treatment=="fallDry"] <- "Early Drought"
gf_graphic2$treatment2[gf_graphic2$treatment=="consistentDry"] <- "Consistent Drought"

ggplot(data = subset(gf_graphic2, subplot %in% c("B","G","F")), aes(x=treatment, y=meancover, color=func)) + 
  geom_point() + 
  facet_wrap(~subplot*year, ncol=3) +
  theme_bw() + 
  geom_errorbar(aes(ymax = meancover+secover, ymin = meancover-secover), width=.25) + 
  labs(x="Treatment", y="Percent cover", fill="Functional
       group") +
  theme(text = element_text(size=20))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

library(FD)
setwd("~/Dropbox/ClimVar/DATA/Plant_composition_data")
cover<-read.csv("Cover/Cover_CleanedData/ClimVar_species-cover.csv")
cover<- cover %>% dplyr::select(-X, -species, -genus, -status, -func, -func2) %>% spread(species_name, cover) %>% filter(subplot!='C')%>%filter(subplot!="XC")
trait.dat1<-read.csv("Traits/Traits_GHscreening_Brad/BradTraits_Cleaned.csv") 
trait.dat<-trait.dat1%>%dplyr::select(-Taxon, -Origin, -GF, -Trt, -Ht, -LDMC, -SLA)
abovetr15<-read.csv("Traits/Above_Traits_Cleaned_2015_2.csv")
all_trait<-left_join(abovetr15,trait.dat, by="ID")
levels(abovetr15$Taxon)
names(cover)<-str_replace_all(names(cover), c("\\ " = "_"))
names(cover)<-str_replace_all(names(cover), c("\\-" = "_"))
cover_2<- cover %>% dplyr::select(-Anagalis_arvensis, -Bromus_sp., -Bromus_sterilis, -Convolvulus_arvensis, 
                                      -Gastridium_phleoides, -Hordeum_sp., -Juncus_bufonius, -Kickxia_spuria,
                                      -Linum_bienne, -Lythrum_hyssopifolia, -Medicago_arabica, -Medicago_polymorpha, -Rumex_pulcher,
                                      -Sherardia_arvensis, -Sonchus_oleraceus, -Zeltnera_muehlenbergii)
cover_fd<- cover_2 %>% dplyr::select(-plot, -subplot, -treatment, -shelterBlock, -shelter, -year)
abovetr_fd<-abovetr15 %>% filter(Ht > 0) %>% dplyr::select(-Taxon, -source, -GF, -Origin) 
cover_fd2 <- cover_fd %>% rename(ACHMIL=Achillea_millefolium, AVEBAR=Avena_barbata, AVEFAT=Avena_fatua, BRADIS=Brachypodium_distachyon, BRIMIN=Briza_minor, BRODIA=Bromus_diandrus, BROHOR=Bromus_hordeaceus, BROMAD=Bromus_madritensis_madritensis, CARPYC=Carduus_pycnocephalus, CENSOL=Centaurea_solstitialis, CERGLO=Cerastium_glomeratum, CLAPUR=Clarkia_amoena, CYNDAC=Cynodon_dactylon, CYNECH=Cynosaurus_echinatus,
                                     EROBOT=Erodium_botrys, EROCIC=Erodium_cicutarium, EROMOS=Erodium_moschatum, FILGAL=Fillago_gallica, GALPAR=Galium_parisiense, GERMOL=Geranium_sp., HORMAR=Hordeum_marinum, HORMUR=Hordeum_murinum, HYPGLA=Hypochaeris_glabra, HYPRAD=Hypochaeris_radicata, LACSER=Lactuca_serriola, LOLMUL=Lolium_multiflorum, LUPBIC=Lupinus_bicolor, SENVUL=Senecio_vulgaris,
                                     SILGAL=Silene_gallica, TAECAP=Taeniatherum_caput_medusae, TORARV=Torilis_arvensis, TRIDUB=Trifolium_dubium, TRIGLO=Trifolium_glomeratum, TRIHIR=Trifolium_hirtum, TRISP=Trifolium_sp.,  TRISUB=Trifolium_subterraneum, TRIWIL=Trifolium_wildenovii, TRIHYA=Triteleia_hyacintha, VICSAT=Vicia_sativa, VULBRO=Vulpia_bromoides, VULMYU=Vulpia_myuros)

#remove species that do not occur in any community
colSums(cover_fd2)
cover_fd3 <- cover_fd2 %>% dplyr::select(-LUPBIC, -SENVUL, -TRIWIL, -Unknown)
cover_fd3<-data.matrix(cover_fd3)
#remove species from traits so matches cover data
abovetr15_fd2 <- abovetr15_fd[-c(27,28,37), ]
abovetr15_fd2 <- abovetr15_fd2 %>% dplyr::select(-ID)



aboveFD_2<-dbFD (abovetr15_fd2, cover_fd3, w.abun = T, stand.x = F,
                 calc.FRic = TRUE, m = "max", stand.FRic = FALSE,
                 scale.RaoQ = FALSE, calc.FGR = FALSE, clust.type = "ward",
                 km.inf.gr = 2, km.sup.gr = nrow(x) - 1, km.iter = 100, calc.CWM = TRUE,
                 calc.FDiv = TRUE, print.pco = FALSE, messages = TRUE)

aboveFD_2<-as.data.frame(aboveFD_2)
aboveFD.z <- data.frame(scale(aboveFD_2)) #convert to z-scores

May_ANPP_noC<- arrange(May_ANPP_noC, as.numeric(as.character(plot)), treatment, subplot)
aboveFD_2b<-bind_cols(May_ANPP_noC, aboveFD_2)
rao_anpp<-ggplot(aboveFD_2b, aes(x=RaoQ, y=weight_g_m, group=treatment2, color=treatment2))+
  #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
  #annotate("text", x= c(1500, 6000, 9000), y = c(750, 450, 315), label = c("R2=0.31", "R2=0.04", "R2=0.09"), color = c("brown2", "dodgerblue3","forestgreen")) +
  geom_point()+
  scale_color_manual(values = c("royalblue2","sienna","lightsteelblue1", "peachpuff"), guide = guide_legend(title = "Treatment")) +
  #xlim(40,100)+
  #ylim(0,100)
  theme_bw()
  #geom_smooth(method='lm', se=FALSE, color="black", aes(group=4))
rao_anpp

rao_anppB<-ggplot(subset(aboveFD_2b, subplot=="B"), aes(x=RaoQ, y=weight_g_m, group=treatment, color=treatment))+
  #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
  #annotate("text", x= c(1500, 6000, 9000), y = c(750, 450, 315), label = c("R2=0.31", "R2=0.04", "R2=0.09"), color = c("brown2", "dodgerblue3","forestgreen")) +
  geom_point()+
  #xlim(40,100)+
  #ylim(0,100)
  theme_bw()
#geom_smooth(method='lm', se=FALSE, color="black", aes(group=4))
rao_anppB

#make bray-curtis dissimilarity matrix
spp.bcd <- vegdist(cover_fd3)
spp.mds<-metaMDS(cover_fd3, trace = FALSE, autotransform=T, trymax=100, k=6) #runs several with different starting configurations
#trace= TRUE will give output for step by step what its doing
#default is 2 dimensions, can put k=4 for 4 dimensions
spp.mds #solution converged after 20 tries, stress = 9.04
summary(spp.mds)

stressplot(spp.mds, spp.bcd) #stressplot to show fit, fit is decent
ordiplot(spp.mds)
spscores1<-scores(spp.mds,display="sites",choices=1)
spscores2<-scores(spp.mds,display="sites",choices=2)
tplots<-dat_may[,3]
tplot_levels<-levels(tplots)
spscoresall<-data.frame(tplots,spscores1,spscores2)

bio.plot <- ordiplot(spp.mds,choices=c(1,2), type = "none")   #Set up the plot
colsT1 <- rep(c("indianred4",  "dodgerblue1", "darkgoldenrod", "purple"), each = 9) #color based on trt
Lcols <- rep(c("indianred4",  "dodgerblue1", "darkgoldenrod", "purple"))
shapes <- rep(c(15, 5, 17), each=3) #shapes on comp treatment
Lshapes <- rep(c(15,5,17))
points(spscoresall$NMDS1,spscoresall$NMDS2,col=colsT1,pch=shapes) 
text(spp.mds, display = "species", cex=0.5, col="grey30") #label species
legend("bottomright",legend=levels(as.factor(cover_2$treatment)), col=Lcols, pch=15, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)
# add legend for treatment
legend("topright",legend=levels(as.factor(cover_2$subplot)), col="black", pch=Lshapes, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)
#help(ordiplot)

bio.plot2 <- ordiplot(spp.mds,choices=c(1,2), type = "none")   #Set up the plot
colsT2 <- rep(c("grey70",  "black"), each = 18) #color based on site
Lcols2 <- rep(c("grey70",  "black", "darkgoldenrod"))
shapes <- rep(c(15, 17), each=18) #shapes on site
Lshapes <- rep(c(15,17))
points(spscoresall$NMDS1,spscoresall$NMDS2,col=colsT2,pch=shapes) 
text(spp.mds, display = "species", cex=0.5, col="grey30") #label species
legend("bottomright",legend=levels(as.factor(dat$site)), col=Lcols2, pch=15, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)
# add legend for treatment
legend("topright",legend=levels(dat$ppt_trt), col="black", pch=Lshapes, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)
#help(ordiplot)


## PCA for trait groups by plant species; visualized by orgin and functional group ##
## Followed by CMW of PCA scores, regressed against ANPP and BNPP
## NOTE: Uses all_trait, abovetr15_fd2, cover15_fd2 from traitxANPP.R and BNPP1 from Lina's BNPP_exploratory_analysis.R


# Select out correlated traits(MD, actual_area, Total) and those I don't have as much faith in (RMF, RGR)
all_trait2 <- all_trait %>%
  dplyr::select( -Seed.mass.grams, -C.N.Ratio)
#all_trait2 <- all_trait2[-c(42:77), ]

# Remove legumes from analysis if desired
#trait.dat2 <- subset(trait.dat2, GF != "L")
#trait.dat2$GF <- as.character(trait.dat2$GF)
#trait.dat2$GF <- as.factor(trait.dat2$GF)

# First a PCA for ALL TRAITS (ABOVE AND BELOW)
# matrix for PCA
traits <- as.data.frame(all_trait2[,c(6:ncol(all_trait2))]) 
traits <- traits %>% dplyr::select(-Total,-MD, -actual_area, -RGR)
row.names(traits) <- all_trait2$ID

# run PCA
myrda <- rda(na.omit(traits), scale = TRUE)

# extract values
siteout <- as.data.frame(scores(myrda, choices=c(1,2), display=c("sites")))
siteout$ID<-rownames(siteout)
siteout$name <- siteout$ID

enviroout<-as.data.frame(scores(myrda, choices=c(1,2), display=c("species")))
enviroout$type<-"traits"
enviroout$name<-rownames(enviroout)

# merge PC axes with trait data
tog <- left_join(all_trait2, siteout) %>% 
  mutate(func = paste(Origin, GF, sep = "_"))

# Remove legumes from legend key (if running PCA without legumes)
#tog <- subset(tog, GF != "L")

#pdf("TraitPCA.pdf", width = 9, height = 7.5)
#pdf("TraitPCA_noLegumes.pdf", width = 9, height = 7.5)

ggplot(tog, aes(x=PC1, y=PC2))+ 
  geom_hline(aes(yintercept=0), color="grey") + 
  geom_vline(aes(xintercept=0), color="grey") +
  geom_text(aes(label = name, color = func), size = 5) +
  # scale_color_manual(values = c("grey20", "grey70")) +
  geom_segment(data = enviroout,
               aes(x = 0, xend =  PC1,
                   y = 0, yend =  PC2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "black") + #grid is required for arrow to work.
  geom_text(data = enviroout,
            aes(x=  PC1*1.2, y =  PC2*1.2, #we add 10% to the text to push it slightly out from arrows
                label = name), #otherwise you could use hjust and vjust
            size = 6,
            hjust = 0.5, 
            color="black") + 
  theme_bw() +theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                    text=element_text(size = 20))+ 
  xlab(paste("Axis 1 (",sprintf("%.1f",myrda$CA$eig["PC1"]/myrda$tot.chi*100,3),"%)",sep="")) +
  ylab(paste("Axis 2 (",sprintf("%.1f",myrda$CA$eig["PC2"]/myrda$tot.chi*100,3),"%)",sep="")) 

#dev.off()

#keep species that we have trait data for
cover_fd4 <- cover_fd2 %>% dplyr::select(ACHMIL, AVEFAT, BRADIS, BRODIA, BROHOR, BROMAD, CENSOL, CYNDAC, CYNECH, EROBOT, HORMUR, LACSER, LOLMUL, TAECAP,TRIHIR,VULMYU)
cover_fd4<-data.matrix(cover_fd4)

siteout<- siteout %>% dplyr::select(-ID, -name) 
siteout <- siteout[-c(14), ] #drop LUPBIC from PCA scores bc it has 0 % cover in any plot 

siteout_fd<-dbFD (siteout, cover_fd4, w.abun = T, stand.x = F,
                  calc.FRic = TRUE, m = "max", stand.FRic = FALSE,
                  scale.RaoQ = FALSE, calc.FGR = FALSE, clust.type = "ward",
                  km.inf.gr = 2, km.sup.gr = nrow(x) - 1, km.iter = 100, calc.CWM = TRUE,
                  calc.FDiv = TRUE, print.pco = FALSE, messages = TRUE)
siteout_fd<-as.data.frame(siteout_fd)

May_ANPP_noC<- arrange(May_ANPP_noC, as.numeric(as.character(plot)), treatment, subplot)
May_ANPP_noC2<-bind_cols(May_ANPP_noC, siteout_fd)

ggplot(data=May_ANPP_noC2, aes(x=CWM.PC1, y=weight_g_m, group=treatment, color=treatment))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
  geom_point()+
  theme_bw()+
  #scale_color_manual(values=c("tomato", "green3", "dodgerblue"), guide = guide_legend(title = "Treatment"), #change legend title
                     #labels=c("Mixed", "Forb", "Grass"))+ #change labels in the legend)+
  xlab("Community Weighted Means of PC1 Scores")+
  ylab("ANPP g/m2")+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method="lm", formula= y ~ x, se=FALSE,  aes(group=1))

ggplot(data=May_ANPP_noC2, aes(x=CWM.PC2, y=weight_g_m, group=treatment, color=treatment))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
  geom_point()+
  theme_bw()+
  xlab("Community Weighted Means of PC2 Scores")+
  ylab("ANPP g/m2")+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method="lm", formula= y ~ x, se=FALSE, color="black", aes(group=1))


#############################################
##run PCA again but for aboveground traits only
traits2 <- as.data.frame(all_trait2[,c(6:8)])
row.names(traits2) <- all_trait2$ID

myrda2 <- rda(na.omit(traits2), scale = TRUE)

# extract values
siteout2 <- as.data.frame(scores(myrda2, choices=c(1,2), display=c("sites")))
siteout2$ID<-rownames(siteout2)
siteout2$name <- siteout2$ID

enviroout2<-as.data.frame(scores(myrda2, choices=c(1,2), display=c("species")))
enviroout2$type<-"traits"
enviroout2$name<-rownames(enviroout2)

# merge PC axes with trait data
tog2 <- left_join(all_trait2, siteout2) %>%
  mutate(func = paste(Origin, GF, sep = "_"))

a.pca<-ggplot(tog2, aes(x=PC1, y=PC2))+ 
  #ggtitle("a) PCA on aboveground traits")+
  geom_hline(aes(yintercept=0), color="grey") + 
  geom_vline(aes(xintercept=0), color="grey") +
  geom_text(aes(label = name, color = GF), size = 5) +
  scale_color_manual(values = c("mediumpurple", "seagreen3","sienna1"), labels=c("Forb","Grass", "Legume"), guide = guide_legend(title = "Functional Group")) +
  geom_segment(data = enviroout2,
               aes(x = 0, xend =  PC1,
                   y = 0, yend =  PC2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "black") + #grid is required for arrow to work.
  geom_text(data = enviroout2,
            aes(x=  PC1*1.2, y =  PC2*1.2, #we add 10% to the text to push it slightly out from arrows
                label = name), #otherwise you could use hjust and vjust
            size = 6,
            hjust = 0.5, 
            color="black") + 
  theme_bw() +theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                    text=element_text(size = 20))+ 
  xlab(paste("Axis 1 (",sprintf("%.1f",myrda2$CA$eig["PC1"]/myrda2$tot.chi*100,3),"%)",sep="")) +
  ylab(paste("Axis 2 (",sprintf("%.1f",myrda2$CA$eig["PC2"]/myrda2$tot.chi*100,3),"%)",sep="")) 
a.pca

aboveFD_2b$CWM.PC1<-May_ANPP_noC2$CWM.PC1
aboveFD_2b$CWM.PC2<-May_ANPP_noC2$CWM.PC2

a1<-ggplot(data=aboveFD_2b, aes(x=CWM.PC1, y=RaoQ, group=treatment, color=treatment))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=treatment)) +
  #ggtitle("a)")+
  geom_point()+
  facet_wrap(~year, ncol=1)+
  theme_bw()+
  #theme(legend.position="none")+
  #scale_color_manual(values=c("tomato", "green3", "dodgerblue"), guide = guide_legend(title = "Treatment"), #change legend title
  #labels=c("Mixed", "Forb", "Grass"))+ #change labels in the legend)+
  xlab("Community Weighted Means of Aboveground PC1 Scores")+
  ylab("RaoQ")+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method="lm", formula= y ~ x, se=FALSE, aes(group=treatment))
a1

#keep species that we have trait data for
cover_fd5 <- cover_fd2 %>% dplyr::select(ACHMIL, AVEBAR, AVEFAT, BRADIS, BRIMIN, BRODIA, BROHOR, BROMAD, CARPYC, CENSOL, CERGLO, CLAPUR, CYNDAC, CYNECH, EROBOT, EROCIC, EROMOS, FILGAL, GALPAR,GERMOL,HORMAR, HORMUR,HYPGLA, HYPRAD, LACSER, LOLMUL, SILGAL, TAECAP,TORARV, TRIDUB,TRIGLO, TRIHIR,TRISP, TRISUB,TRIHYA,VICSAT,VULBRO, VULMYU)
cover_fd5<-data.matrix(cover15_fd5)

#remove excess rows from PCA scores
#siteout2 <- siteout2[-c(15), ]
siteout2<- siteout2 %>% dplyr::select(-ID, -name)
siteout2 <- siteout2[-c(27,28,37), ] #drop species from PCA scores w/ 0 % cover in any plot 

siteout2_fd<-dbFD (siteout2, cover_fd5, w.abun = T, stand.x = F,
                   calc.FRic = TRUE, m = "max", stand.FRic = FALSE,
                   scale.RaoQ = FALSE, calc.FGR = FALSE, clust.type = "ward",
                   km.inf.gr = 2, km.sup.gr = nrow(x) - 1, km.iter = 100, calc.CWM = TRUE,
                   calc.FDiv = TRUE, print.pco = FALSE, messages = TRUE)
siteout2_fd<-as.data.frame(siteout2_fd)


May_ANPP_noC<- arrange(May_ANPP_noC, as.numeric(as.character(plot)), treatment, subplot)
May_ANPP_noC3<-bind_cols(May_ANPP_noC, siteout2_fd)

a.1<-ggplot(data=May_ANPP_noC3, aes(x=CWM.PC1, y=weight_g_m, group=treatment, color=treatment))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
  #ggtitle("a)")+
  geom_point()+
  theme_bw()+
  theme(legend.position="none")+
  #scale_color_manual(values=c("tomato", "green3", "dodgerblue"), guide = guide_legend(title = "Treatment"), #change legend title
                     #labels=c("Mixed", "Forb", "Grass"))+ #change labels in the legend)+
  xlab("Community Weighted Means of Aboveground PC1 Scores")+
  ylab("ANPP g/m2")+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method="lm", formula= y ~ x, se=FALSE, color="black", aes(group=1))
a.1

a.2<-ggplot(data=May_ANPP_noC3, aes(x=CWM.PC2, y=weight_g_m, group=treatment, color=treatment))+
  #ggtitle("b)")+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
  geom_point()+
  theme_bw()+
  theme(legend.position="none")+
  #scale_color_manual(values=c("tomato", "green3", "dodgerblue"), guide = guide_legend(title = "Treatment"), #change legend title
                     #labels=c("Mixed", "Forb", "Grass"))+ #change labels in the legend)+
  xlab("Community Weighted Means of Aboveground PC2 Scores")+
  ylab("")+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method="lm", formula= y ~ x, se=FALSE, color="black", aes(group=1))
a.2

aboveFD_2b$CWM.PC1<-May_ANPP_noC3$CWM.PC1
aboveFD_2b$CWM.PC2<-May_ANPP_noC3$CWM.PC2

b.1<-ggplot(data=subset(aboveFD_2b, subplot=="B"), aes(x=CWM.PC1, y=RaoQ, group=treatment, color=treatment))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=treatment)) +
  #ggtitle("a)")+
  geom_point()+
  facet_wrap(~year, ncol=1)+
  theme_bw()+
  #theme(legend.position="none")+
  #scale_color_manual(values=c("tomato", "green3", "dodgerblue"), guide = guide_legend(title = "Treatment"), #change legend title
  #labels=c("Mixed", "Forb", "Grass"))+ #change labels in the legend)+
  xlab("Community Weighted Means of Aboveground PC1 Scores")+
  ylab("RaoQ")+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method="lm", formula= y ~ x, se=FALSE, aes(group=treatment))
b.1

FD_plot <- aboveFD_2b %>% 
  group_by (treatment2,year, subplot) %>%
  summarize(Rao=mean(RaoQ), seR=sd(RaoQ)/sqrt(length(RaoQ)), 
            CWM.PC1=mean(CWM.PC1),
            ANPP=mean(weight_g_m)) %>%
  tbl_df() 

b.2<-ggplot(subset(FD_plot, year!="2017"), aes(x=CWM.PC1, y=Rao, group=treatment2, color=treatment2))+
  #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=treatment)) +
  #ggtitle("a)")+
  #geom_point()+
  facet_wrap(~year, ncol=1)+
  theme_classic()+
  geom_density_2d()+
  #theme(legend.position="none")+
  scale_color_manual(values = c("royalblue2","sienna","lightsteelblue1", "peachpuff"), guide = guide_legend(title = "Treatment")) +
  #labels=c("Mixed", "Forb", "Grass"))+ #change labels in the legend)+
  xlab("Community Weighted Means of Aboveground PC1 Scores")+
  ylab("RaoQ")
  #xlim(40,100)+
  #ylim(0,100)
  #geom_smooth(method="lm", formula= y ~ x, se=FALSE, aes(group=treatment))
b.2

b.3<-ggplot(data=aboveFD_2b, aes(x=RaoQ, y=weight_g_m, group=treatment2, color=treatment2))+
  #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=treatment)) +
  #ggtitle("a)")+
  #geom_point()+
  facet_wrap(~year, ncol=1)+
  theme_bw()+
  geom_density_2d()+
  #theme(legend.position="none")+
  scale_color_manual(values = c("royalblue2","sienna","lightsteelblue1", "peachpuff"), guide = guide_legend(title = "Treatment")) +  #labels=c("Mixed", "Forb", "Grass"))+ #change labels in the legend)+
  xlab("RaoQ")+
  ylab("ANPP")
#xlim(40,100)+
#ylim(0,100)
#geom_smooth(method="lm", formula= y ~ x, se=FALSE, aes(group=treatment))
b.3

FD_plot2<-FD_plot%>%filter(subplot=="B")

b.4<-ggplot(FD_plot2, aes(x=treatment2, y=Rao, group=treatment2, color=treatment2))+
  #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=treatment)) +
  #ggtitle("a)")+
  geom_point()+
  geom_errorbar(aes(ymin=Rao-seR, ymax=Rao+seR), width=.25)+
  facet_wrap(~year, ncol=4)+
  theme_classic()+
  #geom_density_2d()+
  #theme(legend.position="none")+
  scale_color_manual(values = c("royalblue2","sienna","lightsteelblue1", "peachpuff"), guide = guide_legend(title = "Treatment"))+
  scale_x_discrete(labels=c("Ambient", "Consistent\nDrought", "Early\nDrought", "Late\nDrought"))+
  xlab("Rainfall Treatment")+
  ylab("RaoQ")+
  theme(legend.position="none")
#xlim(40,100)+
#ylim(0,100)
#geom_smooth(method="lm", formula= y ~ x, se=FALSE, aes(group=treatment))
b.4



