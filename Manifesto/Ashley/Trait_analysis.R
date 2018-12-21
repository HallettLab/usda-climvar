library(vegan)

###check trait differences
#does functional diversity change with precipitation variability (drought vs. control) or seasonability (drought treatments)?
#uses objects created in ANPP_dataanalyses
traits<-read.csv("Traits/Traits_ProcessedData-GH/ClimVar_trait-diversity-GH.csv")
May_ANPP<-merge(May_ANPP,traits) 
May_ANPP <- May_ANPP %>%filter( subplot!="C")
May_all_XC <- May_ANPP %>% filter (subplot=="XC")

ggplot(May_all_XC, aes(y=RaoQ, x=forbCover, color = treatment))+
  facet_wrap(~as.factor(as.character(year)))+
  geom_point()+
  geom_smooth(method="lm",se=F)
#scale_colour_gradient(low="magenta", high="plum4")

gf_graph_grass <- May_all_XC %>%
  group_by(treatment, year) %>%
  summarize(Grass=mean(grassCover), GrassSE=sd(grassCover)/sqrt(length(grassCover)))

gf_graph_forb <- May_all_XC %>%
  group_by(treatment, year) %>%
  summarize(Forb=mean(forbCover), ForbSE=sd(forbCover)/sqrt(length(forbCover)))

gf_graph_all <- merge(gf_graph_grass, gf_graph_forb)
gf_graph_all1 <- gf_graph_all %>% dplyr::select(-ForbSE,-GrassSE) %>% gather(func, cover, -treatment, -year)
gf_graph_all2 <- gf_graph_all %>% dplyr::select(-Forb,-Grass) %>% gather(func2, SE, -treatment, -year)
gf_graph_all <- merge(gf_graph_all1, gf_graph_all2)

ggplot(gf_graph_all, aes(x=as.factor(year), y=cover, group = func, color=func)) + 
  geom_errorbar(aes(ymax = cover+SE, ymin = cover-SE), width=.25) + 
  geom_point() + 
  theme_bw() + 
  facet_wrap(~treatment)+
  labs(x="Year", y="Cover %") +
  theme(text = element_text(size=20))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#does functional diversity differ by subplot or treatment?
Fd<-lme(FDis ~ treatment*year, random=~1|shelterBlock, May_all_XC, na.action=na.exclude)
summary(Fd)
anova(Fd)#no effects of treatment on fdis
r.squaredGLMM(Fd) #7% of variation explained by fixed effects, 25% by whole model (spatial variation?)
qqnorm(residuals(Fd))
qqline(residuals(Fd))
shapiro.test(residuals(Fd))
#close to normal
FdLS<-lsmeans(Fd, ~as.factor(year)*treatment)
contrast(FdLS, "pairwise") #no dif

ggplot(May_all_XC, aes(x=as.factor(year), y=FDis, fill=treatment))+
  #facet_wrap(~treatment)+
  geom_boxplot()

Rao<-lme(RaoQ ~ treatment*as.factor(year), random=~1|shelterBlock, May_all_XC, na.action=na.exclude)
summary(Rao)
anova(Rao)
r.squaredGLMM(Rao) #12% of variation explained by fixed effects, 26% by whole model (spatial variation?)
qqnorm(residuals(Rao))
qqline(residuals(Rao))
shapiro.test(residuals(Rao))
#normal
rLS<-lsmeans(Rao, ~treatment*year)
contrast(rLS, "pairwise") #no dif

ggplot(May_all_XC, aes(x=treatment, y=RaoQ, fill=as.factor(year)))+
  geom_boxplot()

#how does functional diversity differ by treatment for 2015?
May_2015_XC2<-May_all_XC %>% filter(year=="2015")
Rao2<-lme(RaoQ ~ treatment, random=~1|shelterBlock, May_2015_XC2, na.action=na.exclude)
summary(Rao2)
anova(Rao2)
r.squaredGLMM(Rao2) #22% of variation explained by fixed effects, 39% by whole model (spatial variation?)
qqnorm(residuals(Rao2))
qqline(residuals(Rao2))
shapiro.test(residuals(Rao2))
#normal
r2LS<-lsmeans(Rao, ~treatment)
contrast(r2LS, "pairwise") 

Fd2<-lme(FDis ~ treatment, random=~1|shelterBlock, May_2015_XC2, na.action=na.exclude)
summary(Fd2)
anova(Fd2)
r.squaredGLMM(Fd2) #18% of variation explained by fixed effects, 33% by whole model (spatial variation?)
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

#which traits are prevalent in each community?
trait_graph_all1 <- traits %>%
  gather(trait, CWM, 14:24)%>%
  group_by(subplot, year, treatment, trait) %>%
  filter(subplot=="XC")%>%
  summarize(CWM.M=mean(CWM), CWM.SE=sd(CWM)/sqrt(length(CWM)))

#let's check traits one by one

#coarse root diameter
mDiamC<-lme(CWM.DiamC ~ treatment*as.factor(year), random=~1|shelterBlock, May_all_XC, na.action=na.exclude)
summary(mDiamC)
anova(mDiamC)#no treatment effect
r.squaredGLMM(mDiamC) #14% of variation explained by fixed effects, 19% by whole model
qqnorm(residuals(mDiamC))
qqline(residuals(mDiamC))
shapiro.test(residuals(mDiamC))
#not normal
mDiamC2<-lme(log(CWM.DiamC+1) ~ treatment*as.factor(year), random=~1|shelterBlock, May_all_XC, na.action=na.exclude)
summary(mDiamC2)
anova(mDiamC2)#no treatment effect
r.squaredGLMM(mDiamC2) #14% of variation explained by fixed effects, 19% by whole model
qqnorm(residuals(mDiamC2))
qqline(residuals(mDiamC2))
shapiro.test(residuals(mDiamC2)) # better 
dcLS<-lsmeans(mDiamC2, ~treatment*year)
contrast(dcLS, "pairwise") #no dif

DC_graph<-ggplot(subset(trait_graph_all1, trait=="CWM.DiamC"), aes(x=year, y=CWM.M, color=treatment))+
  ggtitle("Coarse Root Diameter")+
  #geom_bar(stat="identity", position="dodge")+
  theme_bw()+
  facet_wrap(~treatment, ncol=4)+
  labs(x="", y="mm") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1))  + #remove the legend
  geom_errorbar(aes(ymax=CWM.M+CWM.SE, ymin=CWM.M-CWM.SE), color="black")+
  geom_point(cex=4)
DC_graph

#Height
mHt<-lme(CWM.Ht ~ treatment*year, random=~1|shelterBlock, May_all_XC, na.action=na.exclude)
summary(mHt)
anova(mHt)#no effects
r.squaredGLMM(mHt) #5% of variation explained by fixed effects, 34% by whole model (spatial variation?)
qqnorm(residuals(mHt))
qqline(residuals(mHt))
shapiro.test(residuals(mHt))
#normal
htLS<-lsmeans(mHt, ~treatment)
contrast(htLS, "pairwise") #no effects

Ht_graph<-ggplot(subset(trait_graph_all1, trait=="CWM.Ht"), aes(x=year, y=CWM.M, color=treatment))+
  ggtitle("Height")+
  #geom_bar(stat="identity", position="dodge")+
  theme_bw()+
  labs(x="", y="cm") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1))  + #remove the legend
  facet_wrap(~treatment, ncol=4)+
  geom_errorbar(aes(ymax=CWM.M+CWM.SE, ymin=CWM.M-CWM.SE), color="black")+
  geom_point(cex=4)
Ht_graph

#Leaf Dry Matter Content
mLDMC<-lme(CWM.LDMC ~ treatment, random=~1|year/shelterBlock, May_all_XC, na.action=na.exclude)
summary(mLDMC)
anova(mLDMC)#treatment is significant!!!!!!!!!!!!!!!
r.squaredGLMM(mLDMC) #15% of variation explained by fixed effects, 27% by whole model (spatial variation?)
qqnorm(residuals(mLDMC))
qqline(residuals(mLDMC))
shapiro.test(residuals(mLDMC))
#normal
ldmcLS<-lsmeans(mLDMC, ~treatment)
contrast(ldmcLS, "pairwise") #no differences

LDMC_graph<-ggplot(subset(trait_graph_all1, trait=="CWM.LDMC"), aes(x=year, y=CWM.M, color=treatment))+
  #geom_bar(stat="identity", position="dodge")+
  ggtitle("Leaf Dry Matter Content")+
  theme_bw()+
  labs(x="Treatment", y="proportion (dry:fresh") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1))  + #remove the legend
  facet_wrap(~treatment, ncol=4)+
  geom_errorbar(aes(ymax=CWM.M+CWM.SE, ymin=CWM.M-CWM.SE), color="black")+
  geom_point(cex=4)
LDMC_graph

#Relative growth rate
mRGR<-lme(CWM.RGR ~ treatment*year, random=~1|shelterBlock, May_all_XC, na.action=na.exclude)
summary(mRGR)
anova(mRGR)#trt significant
r.squaredGLMM(mRGR) #14% of variation explained by fixed effects, 50% by whole model (spatial variation?)
qqnorm(residuals(mRGR))
qqline(residuals(mRGR))
shapiro.test(residuals(mRGR))
#normal
rgrLS<-lsmeans(mRGR, ~treatment)
contrast(rgrLS, "pairwise") #no differences

RGR_graph<-ggplot(subset(trait_graph_all1, trait=="CWM.RGR"), aes(x=year, y=CWM.M, color=treatment))+
  #geom_bar(stat="identity", position="dodge")+
  ggtitle("Relative Growth Rate")+
  theme_bw()+
  labs(x="", y="1/t") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) + #remove the legend
  facet_wrap(~treatment, ncol=4)+
  geom_errorbar(aes(ymax=CWM.M+CWM.SE, ymin=CWM.M-CWM.SE), color="black")+
  geom_point(cex=4)
RGR_graph

#Root mass fraction
mRMF<-lme(CWM.RMF ~ treatment*year, random=~1|shelterBlock, May_all_XC, na.action=na.exclude)
summary(mRMF)
anova(mRMF)#nothing significant
r.squaredGLMM(mRMF) #12% of variation explained by fixed effects, 28% by whole model (spatial variation?)
qqnorm(residuals(mRMF))
qqline(residuals(mRMF))
shapiro.test(residuals(mRMF))
#normal
rmfLS<-lsmeans(mRMF, ~treatment)
contrast(rmfLS, "pairwise") #no differences

RMF_graph<-ggplot(subset(trait_graph_all1, trait=="CWM.RMF"), aes(x=year, y=CWM.M, color=treatment))+
  #geom_bar(stat="identity", position="dodge")+
  ggtitle("Root Mass Fraction")+
  theme_bw()+
  labs(x="", y="Proportion") +
  theme(legend.position = "none") + #remove the legend
  facet_wrap(~treatment, ncol=4)+
  geom_errorbar(aes(ymax=CWM.M+CWM.SE, ymin=CWM.M-CWM.SE), color="black")+
  geom_point(cex=4)
RMF_graph

#Specific Leaf Area
mSLA<-lme(CWM.SLA ~ treatment*as.factor(year), random=~1|shelterBlock, May_all_XC, na.action=na.exclude)
summary(mSLA)
anova(mSLA)#nothing significant
r.squaredGLMM(mSLA) #11% of variation explained by fixed effects, 43% by whole model (spatial variation?)
qqnorm(residuals(mSLA))
qqline(residuals(mSLA))
shapiro.test(residuals(mSLA))
#v. close to normal
slaLS<-lsmeans(mSLA, ~treatment)
contrast(slaLS, "pairwise") 

SLA_graph<-ggplot(subset(trait_graph_all1, trait=="CWM.SLA"), aes(x=year, y=CWM.M, color=treatment))+
  #geom_bar(stat="identity", position="dodge")+
  ggtitle("Specific Leaf Area")+
  theme_bw()+
  labs(x="", y="cm^2/g") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1))  + #remove the legend
  facet_wrap(~treatment, ncol=4)+
  geom_errorbar(aes(ymax=CWM.M+CWM.SE, ymin=CWM.M-CWM.SE), color="black")+
  geom_point(cex=4)
SLA_graph

#Specific root length, coarse roots
mSRLC<-lme(CWM.SRLC ~ treatment*year, random=~1|shelterBlock, May_all_XC, na.action=na.exclude)
summary(mSRLC)
anova(mSRLC)#none significant
r.squaredGLMM(mSRLC) #10% of variation explained by fixed effects, 14% by whole model (spatial variation?)
qqnorm(residuals(mSRLC))
qqline(residuals(mSRLC))
shapiro.test(residuals(mSRLC))
#close to normal, log is worse
srlcLS<-lsmeans(mSRLC, ~treatment)
contrast(srlcLS, "pairwise") #both and grass > forb
#grass in control rain > forb in control rain


SRLC_graph<-ggplot(subset(trait_graph_all1, trait=="CWM.SRLC"), aes(x=year, y=CWM.M, color=treatment))+
  #geom_bar(stat="identity", position="dodge")+
  ggtitle("SRL: COARSE")+
  theme_bw()+
  facet_wrap(~treatment, ncol=4)+
  labs(x="", y="cm/g") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1))  + #remove the legend
  geom_errorbar(aes(ymax=CWM.M+CWM.SE, ymin=CWM.M-CWM.SE), color="black")+
  geom_point(cex=4)
SRLC_graph

#Specific root length, fine roots
mSRLF<-lme(CWM.SRLF ~ treatment*as.factor(year), random=~1|shelterBlock, May_all_XC, na.action=na.exclude)
summary(mSRLF)
anova(mSRLF)#nothing significant
r.squaredGLMM(mSRLF) #11% of variation explained by fixed effects, 38% by whole model (spatial variation?)
qqnorm(residuals(mSRLF))
qqline(residuals(mSRLF))
shapiro.test(residuals(mSRLF))
#normal
srlfLS<-lsmeans(mSRLF, ~treatment)
contrast(srlfLS, "pairwise")

SRLF_graph<-ggplot(subset(trait_graph_all1, trait=="CWM.SRLF"), aes(x=year, y=CWM.M, color=treatment))+
  #geom_bar(stat="identity", position="dodge")+
  ggtitle("SRL: FINE")+
  theme_bw()+
  facet_wrap(~treatment, ncol=4)+
  labs(x="", y="cm/g") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1))  + #remove the legend
  geom_errorbar(aes(ymax=CWM.M+CWM.SE, ymin=CWM.M-CWM.SE), color="black")+
  geom_point(cex=4)
SRLF_graph

#density of roots
mDens<-lme(CWM.Dens ~ treatment*year, random=~1|shelterBlock, May_all_XC, na.action=na.exclude)
summary(mDens)
anova(mDens)#nothing significant
r.squaredGLMM(mDens) #10% of variation explained by fixed effects, 24% by whole model (spatial variation?)
qqnorm(residuals(mDens))
qqline(residuals(mDens))
shapiro.test(residuals(mDens))
#not normal
mDens2<-lme(log(CWM.Dens+1) ~ treatment*year, random=~1|shelterBlock, May_all_XC, na.action=na.exclude)
summary(mDens2)
anova(mDens2)#nothing significant
r.squaredGLMM(mDens2) #10% of variation explained by fixed effects, 24% by whole model (spatial variation?)
qqnorm(residuals(mDens2))
qqline(residuals(mDens2))
shapiro.test(residuals(mDens2))
densLS2<-lsmeans(mDens2, ~treatment)
contrast(densLS2, "pairwise")

Dens_graph<-ggplot(subset(trait_graph_all1, trait=="CWM.Dens"), aes(x=year, y=CWM.M, color=treatment))+
  #geom_bar(stat="identity", position="dodge")+
  theme_bw()+
  ggtitle("Root Density")+
  labs(x="", y="Density (g/cm^3)") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1))  + #remove the legend
  facet_wrap(~treatment, ncol=4)+
  geom_errorbar(aes(ymax=CWM.M+CWM.SE, ymin=CWM.M-CWM.SE), color="black")+
  geom_point(cex=4)
Dens_graph

#total biomass
#I'm confused why the total biomass is so low - just ~0.15g?
mtot<-lme(CWM.Total ~ treatment*year, random=~1|shelterBlock, May_all_XC, na.action=na.exclude)
summary(mtot)
anova(mtot)#nothing significant
r.squaredGLMM(mtot) #6% of variation explained by fixed effects, 36% by whole model (spatial variation?)
qqnorm(residuals(mtot))
qqline(residuals(mtot))
shapiro.test(residuals(mtot))
#normal
densLS2<-lsmeans(mDens2, ~treatment)
contrast(densLS2, "pairwise")

Total_graph<-ggplot(subset(trait_graph_all1, trait=="CWM.Total"), aes(x=year, y=CWM.M, color=treatment))+
  #geom_bar(stat="identity", position="dodge")+
  theme_bw()+
  labs(x="", y="Total Biomass (g)") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Total Biomass")+
  facet_wrap(~treatment, ncol=4)+
  geom_errorbar(aes(ymax=CWM.M+CWM.SE, ymin=CWM.M-CWM.SE), color="black")+
  geom_point(cex=4)
Total_graph

#Proportion fine roots
mpf<-lme(CWM.PropF ~ treatment, random=~1|year/shelterBlock, May_all_XC, na.action=na.exclude)
summary(mpf)
anova(mpf)#treatment sign.
r.squaredGLMM(mpf) #14% of variation explained by fixed effects, 45% by whole model (spatial variation?)
qqnorm(residuals(mpf))
qqline(residuals(mpf))
shapiro.test(residuals(mpf))
#normal
pfLS2<-lsmeans(mpf, ~treatment)
contrast(pfLS2, "pairwise")

PF_graph<-ggplot(subset(trait_graph_all1, trait=="CWM.PropF"), aes(x=year, y=CWM.M, color=treatment))+
  ggtitle("Proportion fine roots")+
  #geom_bar(stat="identity", position="dodge")+
  theme_bw()+
  facet_wrap(~treatment, ncol=4)+
  labs(x="", y="Proportion") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) + #remove the legend
  geom_errorbar(aes(ymax=CWM.M+CWM.SE, ymin=CWM.M-CWM.SE), color="black")+
  geom_point(cex=4)
PF_graph

#compile root trait plots 
grid.arrange(PF_graph, RMF_graph, SRLF_graph, SRLC_graph, DC_graph, Dens_graph,  ncol = 2, widths = c(4,4))

#compile aboveground traits
grid.arrange(SLA_graph, LDMC_graph, Ht_graph, ncol=1)



## PCA for traits by treatment and year; visualized by orgin and functional group ##
## Uses "Traits_2015"

# Select out correlated traits(MD, actual_area, Total) and those I don't have as much faith in (RMF, RGR)
trait.dat2 <- traits_2015 %>% dplyr::select(- CWM.Total, -CWM.RMF, - CWM.RGR) %>%
  mutate(ID = paste(subplot, treatment, sep = "_"))

# matrix for PCA
tr <- as.matrix(trait.dat2[,c(14:21)])
row.names(tr) <- trait.dat2$ID

# run PCA
myrda <- rda(tr, scale = TRUE)

# extract values
siteout <- as.data.frame(scores(myrda, choices=c(1,2), display=c("sites")))
siteout$ID<-rownames(siteout)
siteout$name <- siteout$ID

enviroout<-as.data.frame(scores(myrda, choices=c(1,2), display=c("species")))
enviroout$type<-"traits"
enviroout$name<-rownames(enviroout)

# merge PC axes with trait data
tog <- left_join(trait.dat2, siteout) 

pdf("CompTraitPCA_2015.pdf", width = 9, height = 7.5)

ggplot(tog, aes(x=PC1, y=PC2))+ 
  geom_hline(aes(yintercept=0), color="grey") + 
  geom_vline(aes(xintercept=0), color="grey") +
  geom_text(aes(label = name, color = subplot), size = 5) +
  # scale_color_manual(values = c("grey20", "grey70")) +
  geom_segment(data = enviroout,
               aes(x = 0, xend =  PC1,
                   y = 0, yend =  PC2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "black") + #grid is required for arrow to work.
  geom_text(data = enviroout,
            aes(x=  PC1*1.1, y =  PC2*1.1, #we add 10% to the text to push it slightly out from arrows
                label = name), #otherwise you could use hjust and vjust. I prefer this option
            size = 3,
            hjust = 0.5, 
            color="black") + 
  theme_bw() +theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                    text=element_text(size = 15))+ 
  xlab(paste("Axis 1 (",sprintf("%.1f",myrda$CA$eig["PC1"]/myrda$tot.chi*100,3),"%)",sep="")) +
  ylab(paste("Axis 2 (",sprintf("%.1f",myrda$CA$eig["PC2"]/myrda$tot.chi*100,3),"%)",sep="")) 

dev.off()
