#Biodiversity effects across rainfall treatments in CA annual grasslands
#the purpose of this code is additive partitioning of biodiversity effects in the ClimVar mixed composition plots
#this code accompanies the partitioning manuscript and reproducible, publication figures are in the final section 
#Date of last update: 1/February/2021
#Author: Ashley Shaw

library(tidyverse)
library(nlme)
library(ggplot2)
library(multcomp)
library(lsmeans)
library(MuMIn)
library(RColorBrewer)
library(ggpubr)
library(cowplot)

#Setup----
theme_set(theme_bw())
theme_update( panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
              strip.background = element_blank(),
              text = element_text(size = 14),
              strip.text= element_text(size = 12), 
              axis.text = element_text(size = 12))

#function for standard error calculations, for use later
calcSE<-function(x){
  x <- x[!is.na(x)]
  sd(x)/sqrt(length(x))
}

setwd("~/Dropbox/ClimVar/DATA/")
#Read in data for total biomass and biomass by fuctional groups
May_ANPP<-read.csv("Plant_composition_data/ANPP/ANPP_CleanedData/ClimVar_ANPP-peak.csv")
FG_ANPP<-read.csv("Plant_composition_data/ANPP/ANPP_CleanedData/ClimVar_ANPP-func-groups.csv")

#check calculations of total biomass and proportional biomass for each group
FG_ANPP<-FG_ANPP %>% dplyr::select(-totweight,-propweight)%>%
  group_by(year, subplot, plot, treatment, date) %>%
  mutate(totweight = sum(weight_g_m), propweight = weight_g_m/totweight) %>%
  tbl_df() %>% ungroup()%>%
  mutate(func = ifelse(func == "F", "Forb", func), # rename functional groups, check that names are consistent
         func = ifelse(func == "N", "N-fixer", func),
         func = ifelse(func == "G", "Grass", func)) %>%
  mutate(treatment=ordered(treatment, levels = c(consistentDry="consistentDry", springDry="springDry", fallDry="fallDry", controlRain="controlRain"))) 


#FORB contributions to mixed plot production ----
#create subset with forb monocultures only, first harvest b/c forbs most abundant then
#Here we include N-fixers 
FG_forb<-filter(FG_ANPP, subplot=='F', harvest=="First", func=="Forb"|func=="N-fixer", year!="2017") %>% dplyr::select(-harvest, -date) %>% 
  mutate(treatment=ordered(treatment, levels = c(controlRain="controlRain", springDry="springDry", fallDry="fallDry", consistentDry="consistentDry")))

#summarise forb productivity, for use later (see outliers, below)
FG_forb.sum<- FG_forb %>% 
  #filter(weight_g_m>1)%>%
  group_by(treatment, year,func)%>%
  summarise(mean(weight_g_m))

#Now filter mixed plots only, 2015 and 2016, focusing on forb production in the mixed plots
#Here we include N-fixers
FG_B.f2<-filter(FG_ANPP, subplot=='B', harvest=="First", func=="Forb"|func=="N-fixer", year!="2017") %>% dplyr::select(-harvest,  -date)
FG_B.f2$M<-FG_forb$weight_g_m #add a column for 'M' or total forb production in monoculture

#Add N-fixers and forb biomass together, since forb monocultures had both
FG_B.f2<-FG_B.f2 %>% group_by(year, treatment, plot,subplot, shelterBlock)%>%summarize(weight_g_m=sum(weight_g_m), 
                                                                                  totweight=unique(totweight), 
                                                                                  propweight=sum(propweight), M=sum(M))
FG_B.f2$M[FG_B.f2$M == 0.8] <- 177 #replace outlier low monoculture value with group average

#components of additive partitioning for forbs
FG_B.f2<-FG_B.f2 %>% group_by(treatment, shelterBlock, year)%>%
  mutate(RYo=weight_g_m/M)%>% #calculate the observed relative yield of forbs in the mixture, based on the monoculture (RYo)
  mutate(RYe=0.30)%>% #set the expected relative yield of forbs in the mixture to 30%, approximately 1/3
  mutate(dYi=weight_g_m-(RYe*M))%>% #calculate the deviation of forbs from their expected yield (dYi)
  mutate(dRY=RYo-0.30) #calculate the deviation in the relative yield observed (RYo) from the relative yield expected (RYe)


#GRASS contributions to mixed plot production ----
#create subset with grass monocultures only, second harvest b/c grasses most abundant then
FG_grass<-filter(FG_ANPP, subplot=='G', harvest=="Second", func=="Grass", year!="2017") %>% dplyr::select(-harvest,  -date)%>% 
  mutate(treatment=ordered(treatment, levels = c(controlRain="controlRain", springDry="springDry", fallDry="fallDry", consistentDry="consistentDry")))

#Now filter mixed plots only, 2015 and 2016, focusing on grass production in the mixed plots
FG_B.g2<-filter(FG_ANPP, subplot=='B', harvest=="Second", func=="Grass", year!="2017") %>% dplyr::select(-harvest, -date)

#components of additive partitioning for grasses
FG_B.g2$M<-FG_grass$weight_g_m #add a column for 'M' or total grass production in monoculture
FG_B.g2<-FG_B.g2 %>% group_by(treatment, shelterBlock, year)%>%
  mutate(RYo=weight_g_m/M)%>%#calculate the observed relative yield of grass in the mixture, based on the monoculture (RYo)
  mutate(RYe=0.70)%>% #set the expected relative yield of grass in the mixture to 70%, approximately 2/3
  mutate(dYi=weight_g_m-(RYe*M))%>% #calculate the deviation of grasses from their expected yield (dYi)
  mutate(dRY=(RYo-0.70)) #calculate the deviation in the relative yield observed (RYo) from the relative yield expected (RYe)

#Biodiversity effects----
FG_B<-bind_rows(FG_B.f2,FG_B.g2) #combine grass and forb data frames

#calculate complementarity and selection effects
FG_B<-FG_B%>% group_by(year,treatment, shelterBlock)%>% 
  mutate(complement=2*(mean(dRY))*mean(M))%>% #complementarity
  mutate(select=2*cov(dRY,M))%>% #selection
  mutate(dY=sum(dYi), RYT=sum(RYo)) #net effect (dY) and relative yield total (RYT)
FG_B$func[is.na(FG_B$func)] <- "Forb" #name forbs again as a functional group, was deleted when combining

#Exploratory visualizations ----
#rename trt groups for plots
FG_B$treatment2[FG_B$treatment=="controlRain"] <- "Ambient"
FG_B$treatment2[FG_B$treatment=="springDry"] <- "Late"
FG_B$treatment2[FG_B$treatment=="fallDry"] <- "Early"
FG_B$treatment2[FG_B$treatment=="consistentDry"] <- "Consistent"

#to create a plot of biomass of forb, grass, mixed, we include forb biomass in monoculture, grass biomass in monoculture
FG_B.g2$Fweight<-FG_B.f2$M
FG_B.g2$Gweight<-FG_grass$weight_g_m
FG_all<-FG_B.g2%>%gather("Fweight","Gweight","totweight",key=comp, value=biomass)

#rename trt groups for plots
FG_all$treatment2[FG_all$treatment=="controlRain"] <- "Ambient"
FG_all$treatment2[FG_all$treatment=="springDry"] <- "Late"
FG_all$treatment2[FG_all$treatment=="fallDry"] <- "Early"
FG_all$treatment2[FG_all$treatment=="consistentDry"] <- "Consistent"

#summarize for graphs by treatment within year
FG_B.sum_yr2<- FG_B %>% group_by(treatment2, year)%>%
  summarize(sel.m=mean(select), sel.se=calcSE(select), 
            comp.m=mean(complement), comp.se=calcSE(complement),
            net.m=mean(dY), net.se=calcSE(dY))%>%
  gather(net.m, sel.m, comp.m,net.se,comp.se,sel.se, key="type", value="value")%>%
  separate(type, c("type", "x"))%>%
  spread(x,value)

#summarize for graphs by treatment
FG_B.sum2<- FG_B %>% group_by(treatment2)%>%
  summarize(sel.m=mean(select), sel.se=calcSE(select), 
            comp.m=mean(complement), comp.se=calcSE(complement),
            net.m=mean(dY), net.se=calcSE(dY))%>%
  gather(net.m, sel.m, comp.m,net.se,comp.se,sel.se, key="type", value="value")%>%
  separate(type, c("type", "x"))%>%
  spread(x,value)

#summarize for RYT for graphs by treatment
FG_B.sum.RYT<- FG_B %>% group_by(treatment2,year)%>% filter (func=="Grass")%>%
  summarize(RYT.m=mean(RYT), ryt.se=calcSE(RYT))
 

#plot dY, deviation from expected yield
p.dY<-ggplot(d=FG_B, aes(x=treatment2, y=dY, fill=treatment2)) +
  facet_wrap(~year)+
  ggtitle("b) Net Effect ")+
  scale_fill_manual(values = c("royalblue2","sienna","lightsteelblue1", "peachpuff"), guide = guide_legend(title = "Treatment")) +
  #theme_linedraw()+
  theme(legend.position="none")+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  labs(x="Drought Treatment", y="Net Effect (ΔY)")+
  geom_boxplot(aes(y=dY), shape=16)
p.dY

#treament effects on selection
ggplot(d=FG_B, aes(x=treatment2, y=select, fill=treatment2)) +
  #facet_wrap(~year)+
  ggtitle("a) Selection effect")+
  scale_fill_manual(values = c("royalblue2","sienna","lightsteelblue1", "peachpuff"), guide = guide_legend(title = "Treatment")) +
  #theme_classic()+
  theme(legend.position="none")+
  labs(x="Drought treatment", y="Selection Effect")+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  geom_boxplot(aes(y=select), shape=16)

#selection effects by year
ggplot(d=FG_B, aes(x=treatment2, y=select, fill=treatment2)) +
  facet_wrap(~year)+
  #ggtitle("b) Complementarity Effect")+
  scale_fill_manual(values = c("royalblue2","sienna","lightsteelblue1", "peachpuff"), guide = guide_legend(title = "Treatment")) +
  #theme_classic()+
  #annotate("text", x= c("Ambient", "Consistent Drought","Early Drought", "Late Drought"), y = c(50,50,50,50), label = c("a", "b", "a", "a"), color = "black") +
  theme(legend.position="none")+
  labs(x="Drought Treatment", y="Selection Effect")+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  geom_boxplot(aes(y=select), shape=16)

#treatment effects on complementarity
ggplot(d=FG_B, aes(x=treatment2, y=complement, fill=treatment2)) +
  #facet_wrap(~year)+
  ggtitle("b) Complementarity Effect")+
  scale_fill_manual(values = c("royalblue2","sienna","lightsteelblue1", "peachpuff"), guide = guide_legend(title = "Treatment")) +
  #theme_classic()+
  #annotate("text", x= c("Ambient", "Consistent Drought","Early Drought", "Late Drought"), y = c(50,50,50,50), label = c("a", "b", "a", "a"), color = "black") +
  theme(legend.position="none")+
  labs(x="Drought Treatment", y="Complementarity Effect")+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  geom_boxplot(aes(y=complement), shape=16)

#complementarity by year
ggplot(d=FG_B, aes(x=treatment2, y=complement, fill=treatment2)) +
  facet_wrap(~year)+
  #ggtitle("b) Complementarity Effect")+
  scale_fill_manual(values = c("royalblue2","sienna","lightsteelblue1", "peachpuff"), guide = guide_legend(title = "Treatment")) +
  #theme_classic()+
  #annotate("text", x= c("Ambient", "Consistent Drought","Early Drought", "Late Drought"), y = c(50,50,50,50), label = c("a", "b", "a", "a"), color = "black") +
  theme(legend.position="none")+
  labs(x="Drought Treatment", y="Complementarity Effect")+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  geom_boxplot(aes(y=complement), shape=16)

#deviation from expected yields of each functional group by treatment
FG_B.sum3<- FG_B %>% group_by(func,year,treatment2)%>%
  summarize(dYi.m=mean(dYi), dYi.se=calcSE(dYi))

#Deviation in expected yield by functional group (dYi)
ggplot(d=FG_B, aes(x=treatment2, y=dYi, fill=treatment2)) +
  facet_wrap(~func*year, scales="free")+
  #ylim(-80,200)+
  scale_fill_manual(values = c("royalblue2","sienna","lightsteelblue1", "peachpuff"), guide = guide_legend(title = "Treatment")) +
  #theme_linedraw()+
  theme(legend.position="none")+
  labs(x="Rainfall Treatment", y="Deviation from Expected Yield (ΔYi)")+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  geom_boxplot(aes(y=dYi), shape=16)

#Relative yield totals for the mixed plots by rainfall treatment and year
ggplot(d=FG_B.sum.RYT, aes(x=treatment2, y=RYT.m, fill=as.factor(year))) +
  #facet_wrap(~func,ncol=2)+
  #ylim(-80,200)+
  #ggtitle()+
  scale_fill_manual(values = c("white","lightgray"), guide = guide_legend(title = "")) +
  #theme_linedraw()+
  theme(legend.position=c(0.05,0.9), legend.key=element_blank())+
  guides(fill=guide_legend(title=""))+
  labs(x="Rainfall Treatment", y=expression(Relative~yield~total))+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  geom_bar(aes(y=RYT.m, x=treatment2), color="black",stat="identity", position = position_dodge(width = 0.9))+
  geom_errorbar(aes(ymin=RYT.m-ryt.se,ymax=RYT.m+ryt.se),position = position_dodge(width = 0.9), width=0)


#Test treatment effects ----
#Here we use mixed models where treatment and year are fixed effects and block is a random effect
#For total biomass, we also compare composition plots, with composition as a fixed effect

#1. total biomass
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

#2. net effects
m.dy<-lme(dY ~treatment2*year, random=~1|shelterBlock, FG_B, na.action=na.exclude)
summary(m.dy)
anova(m.dy) #effect of treatment, treatment*year
r.squaredGLMM(m.dy) #32% of variation explained by fixed effects, 33% by whole model (interannual variation?)
qqnorm(residuals(m.dy))
qqline(residuals(m.dy))
shapiro.test(residuals(m.dy)) #normal
LS.dy<-lsmeans(m.dy, ~treatment2, by="year")
contrast(LS.dy, "pairwise")

#3. selection
m.sel<-lme(select ~treatment*year, random=~1|shelterBlock, FG_B, na.action=na.exclude)
summary(m.sel)
anova(m.sel) #no sign. effects
r.squaredGLMM(m.sel) #13% of variation explained by fixed effects, 15% by whole model (interannual variation?)
qqnorm(residuals(m.sel))
qqline(residuals(m.sel))
shapiro.test(residuals(m.sel)) #note, not normally distributed
LS.sel<-lsmeans(m.sel, ~treatment)
contrast(LS.sel, "pairwise")

#4. complementarity
m.comp<-lme(complement ~treatment*year, random=~1|shelterBlock, FG_B, na.action=na.exclude)
summary(m.comp)
anova(m.comp) #treatment and treatment*year are significant
r.squaredGLMM(m.comp) #27% of variation explained by fixed effects, 30% by whole model (interannual variation?)
qqnorm(residuals(m.comp))
qqline(residuals(m.comp))
shapiro.test(residuals(m.comp)) #normal
LS.all<-lsmeans(m.comp, ~treatment, by="year")
contrast(LS.all, "pairwise")

#5. do forbs do better or worse than expected?
m.f<-lme(dYi ~treatment*year, random=~1|shelterBlock, FG_B.f2, na.action=na.exclude)
summary(m.f)
anova(m.f) #treatment is not significant
r.squaredGLMM(m.f) #8% of variation explained by fixed effects, 14% by whole model (interannual variation?)
qqnorm(residuals(m.f))
qqline(residuals(m.f))
shapiro.test(residuals(m.f)) #normal
LS.f<-lsmeans(m.f, ~treatment, by="year")
contrast(LS.f, "pairwise")

#6. do grasses do better or worse than expected?
m.g<-lme(dYi ~treatment*year, random=~1|shelterBlock, subset(FG_B, func=="Grass"), na.action=na.exclude)
summary(m.g)
anova(m.g) #treatment is not significant
r.squaredGLMM(m.g) #27% of variation explained by fixed effects, 27% by whole model (interannual variation?)
qqnorm(residuals(m.g))
qqline(residuals(m.g))
shapiro.test(residuals(m.g)) #normal
LS.g<-lsmeans(m.g, ~treatment, by="year") 
contrast(LS.g, "pairwise")


#Manuscript figures----

# * FIGURE 1 ----
#real biomass by composition treatment
#total biomass of the mixed, forb, and grass plots by treatment 
FG_all$comp <- ordered(FG_all$comp, levels = c("Gweight","Fweight","totweight"))
f1a<-ggplot(d=subset(FG_all), aes(x=treatment2, y=biomass, fill=comp)) +
  #facet_wrap(~year)+
  ggtitle("Observed Yield (Yo)")+
  #theme_linedraw()+
  #theme(legend.position="none")+
  scale_fill_manual(values = c("white","black", "gray"), guide = guide_legend(title = "Composition"),
                    labels=c( "Grass-only", "Forb-only","Mixed")) +
  labs(x="Drought Treatment", y="ANPP g/m2")+
  geom_boxplot(aes(y=biomass))+
  ylim(0,1500)
f1a

f1a.v2<-ggplot(d=subset(FG_all), aes(x=treatment2, y=biomass, fill=comp)) +
  #facet_wrap(~year)+
  #ggtitle("Observed Yield (Yo)")+
  #theme_linedraw()+
  #theme(legend.position="none")+
  scale_fill_manual(values = c("white","black", "gray"), guide = guide_legend(title = "Composition"),
                    labels=c("Grass","Forb",  "Mixed")) +
  labs(x="", y=expression(ANPP~(g/m^2)))+
  geom_boxplot(aes(y=biomass))+
  ylim(0,1000)
f1a.v2

f1a.v3<-ggplot(d=subset(FG_all, comp=="totweight"), aes(x=treatment2, y=biomass, fill="gray")) +
  #facet_wrap(~year)+
  #ggtitle("Observed Yield (Yo)")+
  #theme_linedraw()+
  theme(legend.position="none")+
  scale_fill_manual(values = c("gray")) +
  #scale_fill_manual(values = c("white","black", "gray"), guide = guide_legend(title = "Composition"),
  #                  labels=c("Grass","Forb",  "Mixed")) +
  labs(x="", y=expression(ANPP~(g/m^2)))+
  geom_boxplot(aes(y=biomass))+
  ylim(0,1000)
f1a.v3

FG_all$comp2 <- ordered(FG_all$comp, levels = c("Gweight","Fweight","totweight"))
pcomp<-ggplot(d=FG_all, aes(x=comp2, y=biomass, fill=comp2)) +
  #theme_linedraw()+
  labs(x="", y=expression(ANPP~(g/m^2)))+
  annotate("text", x= c( "Gweight","Fweight","totweight"), y = c(1000, 1000, 1000), label = c("a", "b", "c"), color = "black") +
  scale_fill_manual(values = c("white","black", "gray"), guide = guide_legend(title = "Composition"),
                    labels=c("Grass","Forb",  "Mixed")) +  theme(legend.position="none")+
  scale_x_discrete(labels=c("Gweight" = "Grass", "Fweight" = "Forb",  "totweight" = "Mixed"), limits=levels(FG_all$comp2))+
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

#cowplot package to create a plot with inset figures
ggdraw(f1a) +
  draw_plot(pcomp, .15, .6, .3, .3) +
  draw_plot(pyr, 0.5,0.6, 0.3,0.3)+
  draw_plot_label(
    c("A", "B", "C"),
    c(0, 0.15, 0.5),
    c(1, 0.92, 0.92),
    size = 14
  )

#plot big panel plus small panels on right
right_col1<-plot_grid(pcomp, pyr, labels = c('B', 'C'), label_size = 14, ncol=1)
right_col1

#plot big panel plus small panels on right
right_col2<-plot_grid(f1a.v3, pyr, labels = c('B', 'C'), label_size = 14, ncol=1)
right_col2

F1A<-plot_grid(f1a.v2 + theme(legend.position="none"), right_col1, labels=c("A",""), rel_widths=c(2,1.5), label_size=14, ncol=2)
F1A

F1A.2<-plot_grid(f1a.v3 + theme(legend.position="none"), right_col1, labels=c("A",""), rel_widths=c(2,1.5), label_size=14, ncol=2)
F1A.2

F1A.3<-plot_grid(pcomp + theme(legend.position="none"), right_col2, labels=c("A",""), rel_widths=c(2,1.5), label_size=14, ncol=2)
F1A.3

#run soil-cleaning.R
legend_b <- get_legend(
  f1a.v2 + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)

Figure1<- plot_grid(F1A, F1D, F1G,legend_b, ncol=1, rel_heights = c(1, 1,1,.1))
Figure1

Figure1.2<- plot_grid(F1A.2, F1D.2, F1G.2,legend_b, ncol=1, rel_heights = c(1, 1,1,.1))
Figure1.2

Figure1.3<- plot_grid(F1A.3, F1D.3, F1G.3,legend_b, ncol=1, rel_heights = c(1, 1,1,.1))
Figure1.3
#save as 850 x 1300

#alternatively, ggpubr for panels
#ggarrange(pcomp, pyr, f1a, ncol = 2, nrow = 2, labels=c("a)","b)", "c)"), heights=c(1,2))

# * FIGURE 2 ----
#Figure 2a - biodiversity effects 
f2a<-ggplot(FG_B.sum2, aes(x=treatment2, y=m))+
  geom_abline(intercept = 0, slope = 0)+
  geom_bar(aes(fill=type), stat="identity", position = "dodge", color="black")+
  geom_errorbar(aes(ymin=m-se,ymax=m+se, group=type), position="dodge", color="black")+
  labs(x="", y="Biodiversity Effect on Productivity")+
  theme(legend.position=c(0.35,0.85))+
  theme(legend.background = element_blank())+
  scale_fill_manual(values = c("white","gray75","gray25"), 
                    guide = guide_legend(title = "Effect"),
                    labels= c("Complementarity", "Net", "Selection")) 
f2a

#Figure 2b - biodiversity effects by year 
f2b<-ggplot(FG_B.sum_yr2, aes(x=treatment2, y=m))+
  geom_abline(intercept = 0, slope = 0)+
  facet_wrap(~year)+
  geom_bar(aes(fill=type), stat="identity", position = "dodge", color="black")+
  geom_errorbar(aes(ymin=m-se,ymax=m+se, group=type), position="dodge", color="black")+
  labs(x="Drought Treatment", y="Biodiversity Effect on Productivity")+
  theme(legend.position="none")+
  scale_fill_manual(values = c("white","gray75","gray25"), 
                    guide = guide_legend(title = "Effect"),
                    labels= c("Complementarity", "Net", "Selection")) 
f2b

#compile figure 2
ggarrange(f2a, f2b, ncol = 1, nrow = 2, labels=c("a)","b)"), heights=c(2,1))


# * FIGURE 3 ---- 
#deviation in expected yield by functional group
f3<-ggplot(d=FG_B.sum3, aes(x=treatment2, y=dYi.m, shape=as.factor(year))) +
  facet_wrap(~func,ncol=2)+
  #ylim(-80,200)+
  #ggtitle()+
  #scale_color_manual(values = c("royalblue2","sienna","lightsteelblue1", "peachpuff"), guide = guide_legend(title = "Treatment")) +
  #theme_linedraw()+
  theme(legend.position=c(0.05,0.9), legend.background=element_blank())+
  guides(shape=guide_legend(title=""))+
  labs(x="Rainfall Treatment", y=expression(Deviation~from~Expected~Yield~(ΔY[i])))+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  geom_point(aes(y=dYi.m), position = position_dodge(width = 0.50), size=3)+
  geom_errorbar(aes(ymin=dYi.m-dYi.se,ymax=dYi.m+dYi.se),position = position_dodge(width = 0.50), width=0)
f3

# * FIGURE S1 ---- 
#biomass by rainfall and composition treatment by year
fs1<-ggplot(d=subset(FG_all), aes(x=treatment2, y=biomass, fill=comp)) +
  facet_wrap(~year)+
  #ggtitle(expression(Observed~Yield~(Y[o])))+
  #theme_linedraw()+
  #theme(legend.position="none")+
  scale_fill_manual(values = c("black","white", "gray"), guide = guide_legend(title = "Composition"),
                    labels=c("Forb-only", "Grass-only", "Mixed")) +
  labs(x="Drought Treatment", y=expression(ANPP~(g/m^2)))+
  geom_boxplot(aes(y=biomass), shape=16)
fs1

# * FIGURE S2 ----
#dRY vs M to help understand selection effects
ggplot(FG_B, aes(y=dRY, x=M))+
  geom_point(aes(color=treatment2,shape=as.factor(func)))+
  geom_smooth(method="lm",(aes(color=treatment2)),se=F)+
  facet_wrap(~year)

##note that figure 4, figure 5, and figure S4 are currently in soil_cleaning.R
#figure S3 is in rainfall.R
##TO-DO: move code for remaining figures here











