library(tidyverse)
library(nlme)
library(dplyr)
library(multcomp)
library(lsmeans)
library(MuMIn)
library(vegan)
library(RColorBrewer)
library(gridExtra)
library(FD)

##Run Aboveground_2015_2.R and Lina's BNPP_exploratory_analysis.R scripts first for ANPP and BNPP objects

setwd("~/Dropbox/ClimVar/DATA/Plant_composition_data")
cover15<-read.csv("Cover/Cover_CleanedData/ClimVar_2015_species-cover_wide.csv")
trait.dat1<-read.csv("Traits/Traits_GHscreening_Brad/BradTraits_Cleaned.csv") 
trait.dat<-trait.dat1%>%dplyr::select(-Taxon, -Origin, -GF, -Trt, -Ht, -LDMC, -SLA)
abovetr15<-read.csv("Traits/Above_Traits_Cleaned_2015_2.csv")
all_trait<-left_join(abovetr15,trait.dat, by="ID")
levels(abovetr15$Taxon)
names(cover15)<-str_replace_all(names(cover15), c("\\." = "_"))
cover15_2<- cover15 %>% dplyr::select(-Anagalis_arvensis, -Bromus_sp_, -Bromus_sterilis, -Convolvulus_arvensis, 
                          -Gastridium_phleoides, -Hordeum_sp_, -Juncus_bufonius, -Kickxia_spuria,
                          -Linum_bienne, -Lythrum_hyssopifolia, -Medicago_arabica, -Medicago_polymorpha, -Rumex_pulcher,
                          -Sherardia_arvensis, -Sonchus_oleraceus, -Zeltnera_muehlenbergii)
cover15_fd<- cover15_2 %>% dplyr::select(-plot, -subplot, -treatment, -shelterBlock, -shelter, -year, -X)
abovetr15_fd<-abovetr15 %>% filter(Ht > 0) %>% dplyr::select(-Taxon, -source, -GF, -Origin) 
cover15_fd2 <- cover15_fd %>% rename(ACHMIL=Achillea_millefolium, AVEBAR=Avena_barbata, AVEFAT=Avena_fatua, BRADIS=Brachypodium_distachyon, BRIMIN=Briza_minor, BRODIA=Bromus_diandrus, BROHOR=Bromus_hordeaceus, BROMAD=Bromus_madritensis_madritensis, CARPYC=Carduus_pycnocephalus, CENSOL=Centaurea_solstitialis, CERGLO=Cerastium_glomeratum, CLAPUR=Clarkia_amoena, CYNDAC=Cynodon_dactylon, CYNECH=Cynosaurus_echinatus,
                                     EROBOT=Erodium_botrys, EROCIC=Erodium_cicutarium, EROMOS=Erodium_moschatum, FILGAL=Fillago_gallica, GALPAR=Galium_parisiense, GERMOL=Geranium_sp_, HORMAR=Hordeum_marinum, HORMUR=Hordeum_murinum, HYPGLA=Hypochaeris_glabra, HYPRAD=Hypochaeris_radicata, LACSER=Lactuca_serriola, LOLMUL=Lolium_multiflorum, LUPBIC=Lupinus_bicolor, SENVUL=Senecio_vulgaris,
                                     SILGAL=Silene_gallica, TAECAP=Taeniatherum_caput_medusae, TORARV=Torilis_arvensis, TRIDUB=Trifolium_dubium, TRIGLO=Trifolium_glomeratum, TRIHIR=Trifolium_hirtum, TRISP=Trifolium_sp_,  TRISUB=Trifolium_subterraneum, TRIWIL=Trifolium_wildenovii, TRIHYA=Triteleia_hyacintha, VICSAT=Vicia_sativa, VULBRO=Vulpia_bromoides, VULMYU=Vulpia_myuros)

#remove species that do not occur in any community
cover15_fd3 <- cover15_fd2 %>% dplyr::select(-EROCIC, -GERMOL, -LUPBIC, -SENVUL, -TRIWIL, -TRIHYA)
cover15_fd3<-data.matrix(cover15_fd3)
#remove species from traits so matches cover data
abovetr15_fd2 <- abovetr15_fd[-c(16,20,27,28,37,38), ]
row.names(abovetr15_fd2)<-abovetr15_fd2$ID 
abovetr15_fd2 <- abovetr15_fd2 %>% dplyr::select(-ID, -Seed.mass.grams, -C.N.Ratio)



aboveFD_2<-dbFD (abovetr15_fd2, cover15_fd3, w.abun = T, stand.x = F,
  calc.FRic = TRUE, m = "max", stand.FRic = FALSE,
  scale.RaoQ = FALSE, calc.FGR = FALSE, clust.type = "ward",
  km.inf.gr = 2, km.sup.gr = nrow(x) - 1, km.iter = 100, calc.CWM = TRUE,
  calc.FDiv = TRUE, print.pco = FALSE, messages = TRUE)

aboveFD_2<-as.data.frame(aboveFD_2)

aboveFD_3<-dbFD (abovetr15_fd2, cover15_fd3, w.abun = T, stand.x = F,
                 calc.FRic = TRUE, m = "max", stand.FRic = T,
                 scale.RaoQ = FALSE, calc.FGR = F, clust.type = "ward",
                 km.inf.gr = 2, km.sup.gr = nrow(x) - 1, km.iter = 100, calc.CWM = TRUE,
                 calc.FDiv = TRUE, print.pco = FALSE, messages = TRUE)
aboveFD_3<-as.data.frame(aboveFD_3)

#RaoQ of height
abovetr15_fd_ht<-abovetr15_fd2 %>% dplyr::select(Ht)
aboveFD_2_ht<-dbFD (abovetr15_fd_ht, cover15_fd3, w.abun = T, stand.x = F,
                    calc.FRic = TRUE, m = "max", stand.FRic = FALSE,
                    scale.RaoQ = FALSE, calc.FGR = FALSE, clust.type = "ward",
                    km.inf.gr = 2, km.sup.gr = nrow(x) - 1, km.iter = 100, calc.CWM = TRUE,
                    calc.FDiv = TRUE, print.pco = FALSE, messages = TRUE)
aboveFD_2_ht<-as.data.frame(aboveFD_2_ht)
aboveFD_2$Rht<-aboveFD_2_ht$RaoQ

#RaoQ of SLA
abovetr15_fd_sla<-abovetr15_fd2 %>% dplyr::select(SLA)
aboveFD_2_sla<-dbFD (abovetr15_fd_sla, cover15_fd3, w.abun = T, stand.x = F,
                     calc.FRic = TRUE, m = "max", stand.FRic = FALSE,
                     scale.RaoQ = FALSE, calc.FGR = FALSE, clust.type = "ward",
                     km.inf.gr = 2, km.sup.gr = nrow(x) - 1, km.iter = 100, calc.CWM = TRUE,
                     calc.FDiv = TRUE, print.pco = FALSE, messages = TRUE)
aboveFD_2_sla<-as.data.frame(aboveFD_2_sla)
aboveFD_2$Rsla<-aboveFD_2_sla$RaoQ

#RaoQ of LDMC
abovetr15_fd_ldmc<-abovetr15_fd2 %>% dplyr::select(LDMC)
aboveFD_2_ldmc<-dbFD (abovetr15_fd_ldmc, cover15_fd3, w.abun = T, stand.x = F,
                      calc.FRic = TRUE, m = "max", stand.FRic = FALSE,
                      scale.RaoQ = FALSE, calc.FGR = FALSE, clust.type = "ward",
                      km.inf.gr = 2, km.sup.gr = nrow(x) - 1, km.iter = 100, calc.CWM = TRUE,
                      calc.FDiv = TRUE, print.pco = FALSE, messages = TRUE)
aboveFD_2_ldmc<-as.data.frame(aboveFD_2_ldmc)
aboveFD_2$Rldmc<-aboveFD_2_ldmc$RaoQ

###standardize the variables (scale function converts to z-scores)
aboveFD.z <- data.frame(scale(aboveFD_2))

traits_2015_2<- arrange(traits_2015, treatment, shelterBlock)
aboveFD_2$ANPPgm<-traits_2015_2$ANPPgm
aboveFD_2$shelterBlock<-traits_2015_2$shelterBlock
aboveFD_2$subplot<-traits_2015_2$subplot
aboveFD_2$treatment<-traits_2015_2$treatment

aboveFD.z$ANPPgm<-traits_2015_2$ANPPgm
aboveFD.z$shelterBlock<-traits_2015_2$shelterBlock
aboveFD.z$subplot<-traits_2015_2$subplot
aboveFD.z$treatment<-traits_2015_2$treatment

ggplot(data=aboveFD_2, aes(x=log(CWM.Ht)))+
geom_density(aes(y = ..count..)) +
  geom_vline(aes(xintercept = mean(log(CWM.Ht)), group=subplot), 
             linetype = "dashed", size = 0.6)

ggplot(data=aboveFD_2, aes(x=log(nbsp), group=subplot, color=subplot))+
  geom_density() 
  #geom_vline(aes(xintercept = grp.mean, group=subplot), 
             #linetype = "dashed", size = 0.6)

#calculate shannon-weiner
aboveFD_2$H <- diversity(cover15_fd3)

ggplot(data=May_2015, aes(x=Avena, y=ANPPgm, group=subplot, color=subplot))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  geom_point()+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method="lm", formula= y ~ x, se=FALSE)

ggplot(aboveFD_2, aes(x=nbsp, y=ANPPgm, group=subplot, color=subplot))+
  #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  geom_point()+
  annotate("text", x= c(7.5, 6.25, 14), y = c(800,525, 260), label = c("R2=0.30", "R2=0.07", "R2=0.36"), color = c("brown2", "dodgerblue3","forestgreen")) +
  #xlim(40,100)+
  #ylim(0,100)
  theme_bw()+
  labs(x="Species Richness", y="ANPP g/m2")+
  geom_smooth(method='lm', se=FALSE)

ggplot(aboveFD_2, aes(x=H, y=ANPPgm, group=subplot, color=subplot, shape=treatment))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  geom_point()+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method='lm', se=FALSE)

May_2015<-arrange(May_2015, treatment, shelterBlock, subplot)
aboveFD_2$Avena<-May_2015$Avena

ggplot(aboveFD_2, aes(x=Avena, y=ANPPgm, group=subplot, color=subplot))+
  #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  geom_point()+
  annotate("text", x= c(50, 60, 30), y = c(775, 525, 260), label = c("R2=0.46", "R2=0.09", "R2=0.01"), color = c("brown2", "dodgerblue3","forestgreen")) +
  #xlim(40,100)+
  #ylim(0,100)
  theme_bw()+
  labs(x="Avena % Cover", y="ANPP g/m2")+
  geom_smooth(method='lm', se=FALSE)

m_H<-lme(ANPPgm ~ nbsp*subplot, random=~1|shelterBlock, aboveFD_2, na.action=na.exclude)
summary(m_H)
anova(m_H)
r.squaredGLMM(m_H) #18% of variation explained by fixed effects, 30% by whole model (spatial variation?)
qqnorm(residuals(m_H))
qqline(residuals(m_H))
shapiro.test(residuals(m_H))
#normal
residuals.H<-m_H$residuals
residuals.H<-as.data.frame(residuals.H)
aboveFD_2$residuals.Avena<-residuals.H$fixed

ggplot(aboveFD_2, aes(x=FDis, y=residuals.Avena, group=subplot, color=subplot))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  geom_point()+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method='lm', se=FALSE)


m_ht<-lme(ANPPgm ~ CWM.Ht, random=~1|shelterBlock, aboveFD_2, na.action=na.exclude)
summary(m_ht)
anova(m_ht)
r.squaredGLMM(m_ht) #18% of variation explained by fixed effects, 30% by whole model (spatial variation?)
qqnorm(residuals(m_ht))
qqline(residuals(m_ht))
shapiro.test(residuals(m_ht))
#normal

residuals.ht<-m_ht$residuals
residuals.ht<-as.data.frame(residuals.ht)
aboveFD_2$residuals.ht<-residuals.ht$fixed

ht<-ggplot(aboveFD_2, aes(x=CWM.Ht, y=ANPPgm, group=subplot, color=subplot))+
  #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  #annotate("text", x= c(50, 40, 15), y = c(750, 490, 400), label = c("R2=0.38", "R2=0.04", "R2=0.02"), color = c("brown2", "dodgerblue3","forestgreen")) +
  geom_point()+
  #ggtitle("a)")+
  theme_bw()+
  ylab("")+
  xlab("CWM Height")+
  scale_color_manual(values=c("tomato", "green3", "dodgerblue"), guide = guide_legend(title = "Treatment"), #change legend title
                     labels=c("Mixed", "Forb", "Grass"))+ #change labels in the legend)+
  theme(legend.position = "none")+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method='lm', se=FALSE, color="black", aes(group=1))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) 
ht

ggplot(aboveFD.z, aes(x=CWM.Ht, y=ANPPgm))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  geom_point()+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method='lm', se=FALSE)

m_sla<-lme(ANPPgm ~ CWM.SLA, random=~1|shelterBlock, aboveFD.z, na.action=na.exclude)
summary(m_sla)
anova(m_sla) #significant
r.squaredGLMM(m_sla) #24% of variation explained by fixed effects, 42% by whole model (spatial variation?)
qqnorm(residuals(m_sla))
qqline(residuals(m_sla))
shapiro.test(residuals(m_sla))
#normal

residuals<-m_sla$residuals
residuals<-as.data.frame(residuals)
aboveFD_2$residuals.sla<-residuals$fixed

sla<- ggplot(aboveFD.z, aes(x=CWM.SLA, y=ANPPgm, group=subplot, color=subplot))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
  #annotate("text", x= c(225, 225, 275), y = c(750, 515, 250), label = c("R2=0.24", "R2=0.04", "R2=0.04"), color = c("brown2", "dodgerblue3","forestgreen")) +
  geom_point()+
  #ggtitle("b)")+
  theme_bw()+
  ylab("")+
  xlab("CWM Specific Leaf Area")+
  scale_color_manual(values=c("tomato", "green3", "dodgerblue"), guide = guide_legend(title = "Treatment"), #change legend title
                     labels=c("Mixed", "Forb", "Grass"))+ #change labels in the legend)+
  theme(legend.position = "none")+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method='lm', se=FALSE, color="black", aes(group=1))
sla

ggplot(aboveFD_2, aes(x=CWM.SLA, y=residuals.ht, group=subplot, color=subplot))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  geom_point()+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method='lm', se=FALSE)

ggplot(aboveFD_2, aes(x=CWM.SLA, y=ANPPgm))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  geom_point()+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method='lm', se=FALSE)

m_ldmc<-lme(ANPPgm ~ CWM.LDMC, random=~1|shelterBlock, aboveFD.z, na.action=na.exclude)
summary(m_ldmc)
anova(m_ldmc) #not significant
r.squaredGLMM(m_ldmc) #30% of variation explained by fixed effects, 38% by whole model (spatial variation?)
qqnorm(residuals(m_ldmc))
qqline(residuals(m_ldmc))
shapiro.test(residuals(m_ldmc))
#not normal

ldmc<-ggplot(aboveFD.z, aes(x=CWM.LDMC, y=ANPPgm, group=subplot, color=subplot))+
  #annotate("text", x= c(0.3, 0.3, 0.2), y = c(625, 250, 375), label = c("R2=0.19", "R2=0.05", "R2=0.03"), color = c("brown2", "dodgerblue3","forestgreen")) +
  geom_point()+
  #ggtitle("c)")+
  theme_bw()+
  ylab("")+
  xlab("CWM Leaf Dry Matter Content")+
  scale_color_manual(values=c("tomato", "green3", "dodgerblue"), guide = guide_legend(title = "Treatment"), #change legend title
                     labels=c("Mixed", "Forb", "Grass"))+ #change labels in the legend)+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method='lm', se=FALSE, color="black", aes(group=1))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) 
ldmc

ggplot(aboveFD_2, aes(x=CWM.LDMC, y=ANPPgm))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  geom_point()+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method='lm', se=FALSE)

#compile FG plots 
grid.arrange(ht, sla, ldmc, ncol = 3, widths = c(1.1,1.1,1.35))
#save as 1200wX600h for PNG

figureS4 <- ggarrange(ht, sla, ldmc,
                    ncol =3, nrow =1, common.legend = TRUE, legend = "bottom",
                    align = "v",labels = c("a)", "b)", "c)"))
figureS4<-annotate_figure(figureS4, 
                left = text_grob("ANPP g/m2", rot = 90))
figureS4

m_rq<-lme(ANPPgm ~ RaoQ, random=~1|shelterBlock, aboveFD_2, na.action=na.exclude)
summary(m_rq)
anova(m_rq)
r.squaredGLMM(m_rq) #30% of variation explained by fixed effects, 38% by whole model (spatial variation?)
qqnorm(residuals(m_rq))
qqline(residuals(m_rq))
shapiro.test(residuals(m_rq))
#close to normal

rao_anpp<-ggplot(aboveFD.z, aes(x=RaoQ, y=ANPPgm, group=subplot, color=subplot))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
  #annotate("text", x= c(1500, 6000, 9000), y = c(750, 450, 315), label = c("R2=0.31", "R2=0.04", "R2=0.09"), color = c("brown2", "dodgerblue3","forestgreen")) +
  geom_point()+
  #xlim(40,100)+
  #ylim(0,100)
  theme_bw()+
  geom_smooth(method='lm', se=FALSE, color="black", aes(group=1))
rao_anpp
  
ggplot(aboveFD.z, aes(x=RaoQ, y=ANPPgm))+
  stat_smooth_func(geom="text",method="lm", formula= y ~ x + I(x^2), hjust=0,parse=TRUE) +
  geom_point()+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method="lm", formula= y ~ x + I(x^2), se=FALSE)

#BNPP1 from Lina's BNPP_exploratory_analysis.R code:
aboveFD.z<-merge(aboveFD.z, BNPP1)

rao_bnpp<-ggplot(aboveFD.z, aes(x=RaoQ, y=agg_BNPP, group=subplot, color=subplot))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
  #annotate("text", x= c(1500, 6000, 9000), y = c(750, 450, 315), label = c("R2=0.31", "R2=0.04", "R2=0.09"), color = c("brown2", "dodgerblue3","forestgreen")) +
  geom_point()+
  #xlim(40,100)+
  #ylim(0,100)
  theme_bw()+
  geom_smooth(method='lm', se=FALSE, color="black", aes(group=1))
rao_bnpp

m_rht<-lme(ANPPgm ~ Rht, random=~1|shelterBlock, aboveFD_2, na.action=na.exclude)
summary(m_rht)
anova(m_rht)
r.squaredGLMM(m_rht) #30% of variation explained by fixed effects, 38% by whole model (spatial variation?)
qqnorm(residuals(m_rht))
qqline(residuals(m_rht))
shapiro.test(residuals(m_rht))
#close to normal

rht<-ggplot(aboveFD.z, aes(x=Rht, y=ANPPgm, group=subplot, color=subplot))+
  #annotate("text", x= c(100,400,200), y = c(625, 375, 250), label = c("R2=0.05", "R2=0.01", "R2=0.04"), color = c("brown2", "dodgerblue3","forestgreen")) +
  geom_point()+
  theme_bw()+
  ylab("")+
  xlab("Rao's Q Height")+
  scale_color_manual(values=c("tomato", "green3", "dodgerblue"), guide = guide_legend(title = "Treatment"), #change legend title
                     labels=c("Mixed", "Forb", "Grass"))+ #change labels in the legend)+
  theme(legend.position = "none")+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method='lm', se=FALSE, color="black", aes(group=1))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1))
rht

ggplot(aboveFD_2, aes(x=Rht, y=ANPPgm))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  geom_point()+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method="lm", formula= y ~ x, se=FALSE)

m_rsla<-lme(ANPPgm ~ Rsla, random=~1|shelterBlock, aboveFD_2, na.action=na.exclude)
summary(m_rsla)
anova(m_rsla)
r.squaredGLMM(m_rsla) #30% of variation explained by fixed effects, 38% by whole model (spatial variation?)
qqnorm(residuals(m_rsla))
qqline(residuals(m_rsla))
shapiro.test(residuals(m_rsla))
#close to normal

rsla<-ggplot(aboveFD.z, aes(x=Rsla, y=ANPPgm, group=subplot, color=subplot))+
  #annotate("text", x= c(2500, 7500, 6250), y = c(750, 500, 225), label = c("R2=0.31", "R2=0.04", "R2=0.09"), color = c("brown2", "dodgerblue3","forestgreen")) +
  geom_point()+
  theme_bw()+
  ylab("")+
  xlab("Rao's Q Specific Leaf Area")+
  scale_color_manual(values=c("tomato", "green3", "dodgerblue"), guide = guide_legend(title = "Treatment"), #change legend title
                     labels=c("Mixed", "Forb", "Grass"))+ #change labels in the legend)+
  theme(legend.position = "none")+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method='lm', se=FALSE, color="black", aes(group=1))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) 
rsla

ggplot(aboveFD_2, aes(x=Rsla, y=ANPPgm))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  geom_point()+
  theme_bw()+
  theme(legend.position = "none")+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method="lm", formula= y ~ x, se=FALSE)

m_rldmc<-lme(ANPPgm ~ Rldmc, random=~1|shelterBlock, aboveFD_2, na.action=na.exclude)
summary(m_rldmc)
anova(m_rldmc)
r.squaredGLMM(m_rldmc) #30% of variation explained by fixed effects, 38% by whole model (spatial variation?)
qqnorm(residuals(m_rldmc))
qqline(residuals(m_rldmc))
shapiro.test(residuals(m_rldmc))
#close to normal

rldmc<-ggplot(aboveFD.z, aes(x=Rldmc, y=ANPPgm, group=subplot, color=subplot))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
  #annotate("text", x= c(0.0025, 0.0012, 0.008), y = c(750, 600, 420), label = c("R2=0.51", "R2=0.05", "R2=0.03"), color = c("brown2", "dodgerblue3","forestgreen")) +
  geom_point()+
  theme_bw()+
  ylab("")+
  xlab("Rao's Q Leaf Dry Matter Content")+
  scale_color_manual(values=c("tomato", "green3", "dodgerblue"), guide = guide_legend(title = "Treatment"), #change legend title
                     labels=c("Mixed", "Forb", "Grass"))+ #change labels in the legend)+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method='lm', se=FALSE, color="black", aes(group=1))
rldmc

ggplot(aboveFD.z, aes(x=Rldmc, y=ANPPgm))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  geom_point()+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method="lm", formula= y ~ x, se=FALSE)

grid.arrange(rht,rsla,rldmc, ncol = 3, widths = c(1.1,1.1,1.35))

figureS2 <- ggarrange(rht, rsla, rldmc,
                      ncol =3, nrow =1, common.legend = TRUE, legend = "bottom",
                      align = "v",labels = c("a)", "b)", "c)"))
figureS2<-annotate_figure(figureS2, 
                          left = text_grob("ANPP g/m2", rot = 90))
figureS2

library(MASS)
fit <- lm(ANPPgm~CWM.Ht+CWM.SLA+CWM.LDMC+Rht+Rsla+Rldmc,data=aboveFD_2)
step <- stepAIC(fit, direction="backward")
step$anova # display results

m_fit<-lme(ANPPgm ~ CWM.SLA + CWM.LDMC + Rsla, random=~1|shelterBlock, aboveFD_2, na.action=na.exclude)
summary(m_fit)
anova(m_fit)
r.squaredGLMM(m_fit) #38% of variation explained by fixed effects, 44% by whole model (spatial variation?)
qqnorm(residuals(m_fit))
qqline(residuals(m_fit))
shapiro.test(residuals(m_fit))
#normal

library(leaps)
leaps<-regsubsets(ANPPgm~CWM.Ht+CWM.SLA+CWM.LDMC+Rht+Rsla+Rldmc,aboveFD_2,nbest=5)
# view results 
summary(leaps)
# plot a table of models showing variables in each model.
# models are ordered by the selection statistic.
plot(leaps,scale="r2")
# plot statistic by subset size 
library(car)
subsets(leaps, statistic="adjr")

library(relaimpo)
calc.relimp(fit,type=c("lmg","last","first","pratt"),
            rela=TRUE)

# Bootstrap Measures of Relative Importance (1000 samples) 
boot <- boot.relimp(fit, b = 1000, type = c("lmg", 
                                            "last", "first", "pratt"), rank = TRUE, 
                    diff = TRUE, rela = TRUE)
booteval.relimp(boot) # print result
plot(booteval.relimp(boot,sort=TRUE)) # plot result

fit <- lm(ANPPgm~CWM.Ht+CWM.SLA+CWM.LDMC,data=subset(aboveFD_2))
step <- stepAIC(fit, direction="backward")
summary(step)
step$anova # display results

fit_b <- lm(ANPPgm~CWM.Ht+CWM.SLA+CWM.LDMC,data=subset(aboveFD_2, (subplot=="B")))
step_b <- stepAIC(fit_b, direction="backward")
summary(step_b)
step_b$anova # display results

aboveFD_2_b<-subset(aboveFD_2, subplot=="B")
m_fit<-lm(ANPPgm ~ CWM.Ht, aboveFD_2_b, na.action=na.exclude)
summary(m_fit)
anova(m_fit)
r.squaredGLMM(m_fit) #38% of variation explained by fixed effects, 44% by whole model (spatial variation?)
qqnorm(residuals(m_fit))
qqline(residuals(m_fit))
shapiro.test(residuals(m_fit))
#normal

fit_f <- lm(ANPPgm~CWM.Ht+CWM.SLA+CWM.LDMC,data=subset(aboveFD_2, (subplot=="F")))
step_f <- stepAIC(fit_f, direction="backward")
summary(step_f)
step_f$anova # display results

fit_g <- lm(ANPPgm~CWM.Ht+CWM.SLA+CWM.LDMC,data=subset(aboveFD_2, (subplot=="G")))
step_g <- stepAIC(fit_g, direction="backward")
summary(step_g)
step_g$anova # display results

fit <- lm(ANPPgm~Rsla+Rldmc+Rht,data=subset(aboveFD_2))
step <- stepAIC(fit, direction="backward")
summary(step)
step$anova # display results

fit_b <- lm(ANPPgm~Rsla+Rldmc+Rht,data=subset(aboveFD_2, (subplot=="B")))
step_b <- stepAIC(fit_b, direction="backward")
summary(step_b)
step_b$anova # display results

fit_f <- lm(ANPPgm~Rsla+Rldmc+Rht,data=subset(aboveFD_2, (subplot=="F")))
step_f <- stepAIC(fit_f, direction="backward")
summary(step_f)
step_f$anova # display results

fit_g <- lm(ANPPgm~Rsla+Rldmc+Rht,data=subset(aboveFD_2, (subplot=="G")))
step_g <- stepAIC(fit_g, direction="backward")
summary(step_g)
step_g$anova # display results

library(leaps)
leaps<-regsubsets(ANPPgm~RaoQ+CWM.Ht+CWM.SLA+CWM.LDMC+Rht+Rsla+Rldmc,aboveFD_2_b,nbest=5)
# view results 
summary(leaps)
# plot a table of models showing variables in each model.
# models are ordered by the selection statistic.
plot(leaps,scale="r2")
# plot statistic by subset size 
library(car)
subsets(leaps, statistic="adjr")

##try analyses with z-scores
ggplot(data=aboveFD.z, aes(x=log(ANPPgm), group=subplot, color=subplot))+
  geom_density(aes(y = ..count..)) +
  geom_vline(aes(xintercept = mean(log(ANPPgm)), group=subplot), 
             linetype = "dashed", size = 0.6)

ggplot(data=aboveFD.z, aes(x=nbsp, group=subplot, color=subplot))+
  geom_density() +
  geom_vline(aes(xintercept = mean(nbsp), group=subplot), 
             linetype = "dashed", size = 0.6)

#calculate shannon-weiner & simpson
aboveFD.z$H <- diversity(cover15_fd3, index="shannon")
aboveFD.z$simpson <- diversity(cover15_fd3, index="simpson")

ggplot(aboveFD.z, aes(x=nbsp, y=ANPPgm, group=subplot, color=subplot))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  geom_point()+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method='lm', se=FALSE)

ggplot(aboveFD.z, aes(x=H, y=ANPPgm, group=subplot, color=subplot))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  geom_point()+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method='lm', se=FALSE)

ggplot(aboveFD.z, aes(x=simpson, y=ANPPgm, group=subplot, color=subplot))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  geom_point()+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method='lm', se=FALSE)

m_ht<-lme(ANPPgm ~ CWM.Ht, random=~1|shelterBlock, aboveFD.z, na.action=na.exclude)
summary(m_ht)
anova(m_ht)
r.squaredGLMM(m_ht) #28% of variation explained by fixed effects, 40% by whole model (spatial variation?)
qqnorm(residuals(m_ht))
qqline(residuals(m_ht))
shapiro.test(residuals(m_ht))
#normal

ggplot(aboveFD.z, aes(x=CWM.Ht, y=ANPPgm, group=subplot, color=subplot))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  geom_point()+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method='lm', se=FALSE)

ggplot(aboveFD.z, aes(x=CWM.Ht, y=ANPPgm))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  geom_point()+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method='lm', se=FALSE)

m_sla<-lme(ANPPgm ~ CWM.SLA, random=~1|shelterBlock, aboveFD.z, na.action=na.exclude)
summary(m_sla)
anova(m_sla) #significant
r.squaredGLMM(m_sla) #24% of variation explained by fixed effects, 37% by whole model (spatial variation?)
qqnorm(residuals(m_sla))
qqline(residuals(m_sla))
shapiro.test(residuals(m_sla))
#normal

ggplot(aboveFD.z, aes(x=CWM.SLA, y=ANPPgm, group=subplot, color=subplot))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  geom_point()+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method='lm', se=FALSE)

ggplot(aboveFD.z, aes(x=CWM.SLA, y=ANPPgm))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  geom_point()+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method='lm', se=FALSE)

m_ldmc<-lme(ANPPgm ~ CWM.LDMC, random=~1|shelterBlock, aboveFD.z, na.action=na.exclude)
summary(m_ldmc)
anova(m_ldmc) #not significant
r.squaredGLMM(m_ldmc) #26% of variation explained by fixed effects, 31% by whole model (spatial variation?)
qqnorm(residuals(m_ldmc))
qqline(residuals(m_ldmc))
shapiro.test(residuals(m_ldmc))
#not normal

ggplot(aboveFD.z, aes(x=CWM.LDMC, y=ANPPgm, group=subplot, color=subplot))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  geom_point()+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method='lm', se=FALSE)

ggplot(aboveFD.z, aes(x=CWM.LDMC, y=ANPPgm))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  geom_point()+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method='lm', se=FALSE)

m_rq<-lme(ANPPgm ~ RaoQ, random=~1|shelterBlock, aboveFD.z, na.action=na.exclude)
summary(m_rq)
anova(m_rq)
r.squaredGLMM(m_rq) #15% of variation explained by fixed effects, 15% by whole model (spatial variation?)
qqnorm(residuals(m_rq))
qqline(residuals(m_rq))
shapiro.test(residuals(m_rq))
#close to normal

ggplot(aboveFD.z, aes(x=RaoQ, y=ANPPgm, group=subplot, color=subplot))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  geom_point()+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method='lm', se=FALSE)

ggplot(aboveFD.z, aes(x=RaoQ, y=ANPPgm))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  geom_point()+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method="lm", se=FALSE)

abovetr15_fd_ht<-abovetr15_fd2 %>% dplyr::select(Ht)
aboveFD_2_ht<-dbFD (abovetr15_fd_ht, cover15_fd3, w.abun = T, stand.x = F,
                    calc.FRic = TRUE, m = "max", stand.FRic = FALSE,
                    scale.RaoQ = FALSE, calc.FGR = FALSE, clust.type = "ward",
                    km.inf.gr = 2, km.sup.gr = nrow(x) - 1, km.iter = 100, calc.CWM = TRUE,
                    calc.FDiv = TRUE, print.pco = FALSE, messages = TRUE)
aboveFD_z_ht<-as.data.frame(scale(aboveFD_2_ht))
aboveFD.z$Rht<-aboveFD_z_ht$RaoQ

m_rht<-lme(ANPPgm ~ Rht, random=~1|shelterBlock, aboveFD.z, na.action=na.exclude)
summary(m_rht)
anova(m_rht)
r.squaredGLMM(m_rht) #1% of variation explained by fixed effects, 1% by whole model (spatial variation?)
qqnorm(residuals(m_rht))
qqline(residuals(m_rht))
shapiro.test(residuals(m_rht))
#not normal

ggplot(aboveFD.z, aes(x=Rht, y=ANPPgm, group=subplot, color=subplot))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  geom_point()+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method='lm', se=FALSE)

ggplot(aboveFD.z, aes(x=Rht, y=ANPPgm))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  geom_point()+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method="lm", formula= y ~ x, se=FALSE)

aboveFD_z_sla<-as.data.frame(scale(aboveFD_2_sla))
aboveFD.z$Rsla<-aboveFD_z_sla$RaoQ

m_rsla<-lme(ANPPgm ~ Rsla, random=~1|shelterBlock, aboveFD.z, na.action=na.exclude)
summary(m_rsla)
anova(m_rsla)
r.squaredGLMM(m_rsla) #15% of variation explained by fixed effects, 15% by whole model (spatial variation?)
qqnorm(residuals(m_rsla))
qqline(residuals(m_rsla))
shapiro.test(residuals(m_rsla))
#close to normal

ggplot(aboveFD.z, aes(x=Rsla, y=ANPPgm, group=subplot, color=subplot))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  geom_point()+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method='lm', se=FALSE)

ggplot(aboveFD.z, aes(x=Rsla, y=ANPPgm))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  geom_point()+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method="lm", formula= y ~ x, se=FALSE)

aboveFD_z_ldmc<-as.data.frame(scale(aboveFD_2_ldmc))
aboveFD.z$Rldmc<-aboveFD_z_ldmc$RaoQ

m_rldmc<-lme(ANPPgm ~ Rldmc, random=~1|shelterBlock, aboveFD.z, na.action=na.exclude)
summary(m_rldmc)
anova(m_rldmc)
r.squaredGLMM(m_rldmc) #14% of variation explained by fixed effects, 22% by whole model (spatial variation?)
qqnorm(residuals(m_rldmc))
qqline(residuals(m_rldmc))
shapiro.test(residuals(m_rldmc))
#close to normal

ggplot(aboveFD.z, aes(x=Rldmc, y=ANPPgm, group=subplot, color=subplot))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  geom_point()+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method='lm', se=FALSE)

ggplot(aboveFD.z, aes(x=Rldmc, y=ANPPgm))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  geom_point()+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method="lm", formula= y ~ x, se=FALSE)

library(MASS)
#CWM
fit <- lm(ANPPgm~CWM.Ht+CWM.SLA+CWM.LDMC,data=aboveFD.z)
step <- stepAIC(fit, direction="both")
summary(step)
step$anova # display results

m_fit<-lme(ANPPgm ~ RaoQ + CWM.SLA + CWM.LDMC, random=~1|shelterBlock, aboveFD.z, na.action=na.exclude)
summary(m_fit)
anova(m_fit)
r.squaredGLMM(m_fit) #38% of variation explained by fixed effects, 45% by whole model (spatial variation?)
qqnorm(residuals(m_fit))
qqline(residuals(m_fit))
shapiro.test(residuals(m_fit))
#normal

library(leaps)
leaps<-regsubsets(ANPPgm~RaoQ+CWM.Ht+CWM.SLA+CWM.LDMC+Rht+Rsla+Rldmc,aboveFD.z,nbest=5)
# view results 
summary(leaps)
# plot a table of models showing variables in each model.
# models are ordered by the selection statistic.
plot(leaps,scale="r2")
# plot statistic by subset size 
library(car)
subsets(leaps, statistic="adjr")

#full model, both only
fit_b <- lm(ANPPgm~CWM.Ht+CWM.SLA+CWM.LDMC,data=subset(aboveFD.z, (subplot=="B")))
step_b <- stepAIC(fit_b, direction="backward")
summary(step_b)
step_b$anova # display results

aboveFD_z_b<-subset(aboveFD.z, subplot=="B")
m_fit<-lm(ANPPgm ~ CWM.Ht, aboveFD_2_b, na.action=na.exclude)
summary(m_fit)
anova(m_fit)
r.squaredGLMM(m_fit) #38% of variation explained by fixed effects, 44% by whole model (spatial variation?)
qqnorm(residuals(m_fit))
qqline(residuals(m_fit))
shapiro.test(residuals(m_fit))
#normal

fit <- lm(ANPPgm~RaoQ+Rsla+Rldmc+Rht,data=aboveFD.z)
step <- stepAIC(fit, direction="backward")
summary(step)
step$anova # display results

fit_b <- lm(ANPPgm~RaoQ+Rsla+Rldmc+Rht,data=subset(aboveFD.z, (subplot=="B")))
step_b <- stepAIC(fit_b, direction="backward")
summary(step_b)
step_b$anova # display results

#full model, grass only
fit_g <- lm(ANPPgm~RaoQ+Rsla+Rldmc+Rht,data=subset(aboveFD.z, (subplot=="G")))
step_g <- stepAIC(fit_g, direction="backward")
summary(step_g)
step_g$anova # display results

#full model, forb only
fit_f <- lm(ANPPgm~RaoQ+Rsla+Rldmc+Rht,data=subset(aboveFD.z, (subplot=="F")))
step_f <- stepAIC(fit_f, direction="backward")
summary(step_f)
step_f$anova # display results

library(leaps)
leaps<-regsubsets(ANPPgm~RaoQ+CWM.Ht+CWM.SLA+CWM.LDMC+Rht+Rsla+Rldmc,aboveFD_z_b,nbest=5)
# view results 
summary(leaps)
# plot a table of models showing variables in each model.
# models are ordered by the selection statistic.
plot(leaps,scale="r2")
# plot statistic by subset size 
library(car)
subsets(leaps, statistic="adjr")

#full model, grass only
fit_g <- lm(ANPPgm~CWM.Ht+CWM.SLA+CWM.LDMC,data=subset(aboveFD.z, (subplot=="G")))
step_g <- stepAIC(fit_g, direction="backward")
summary(step_g)
step_g$anova # display results

#full model, forb only
fit_f <- lm(ANPPgm~CWM.Ht+CWM.SLA+CWM.LDMC,data=subset(aboveFD.z, (subplot=="F")))
step_f <- stepAIC(fit_f, direction="backward")
summary(step_f)
step_f$anova # display results


##now include environmental variables in the model
aboveFD.z<-cbind(aboveFD.z, env)

#full model, all subplots
fit <- lm(ANPPgm~Percent_litter+avg_sm+SOC+MBC+NH4+NO3.NO2,data=aboveFD.z)
step <- stepAIC(fit, direction="both")
summary(step)
step$anova # display results

#full model, both only
fit_b <- lm(ANPPgm~CWM.Ht+CWM.SLA+CWM.LDMC+RaoQ+Rsla+Rldmc+Rht+Litter_depth_cm+avg_sm+SOC+MBC+NH4+NO3.NO2,data=subset(aboveFD.z, (subplot=="B")))
step_b <- stepAIC(fit_b, direction="backward")
summary(step_b)
step_b$anova # display results

#full model, grass only
fit_g <- lm(ANPPgm~CWM.Ht+CWM.SLA+CWM.LDMC+RaoQ+Rsla+Rldmc+Rht+Litter_depth_cm+avg_sm+SOC+MBC+NH4+NO3.NO2,data=subset(aboveFD.z, (subplot=="G")))
step_g <- stepAIC(fit_g, direction="both")
summary(step_g)
step_g$anova # display results

#full model, forb only
fit_f <- lm(ANPPgm~CWM.Ht+CWM.SLA+CWM.LDMC+RaoQ+Rsla+Rldmc+Rht+Litter_depth_cm+avg_sm+SOC+MBC+NH4+NO3.NO2,data=subset(aboveFD.z, (subplot=="F")))
step_f <- stepAIC(fit_f, direction="backward")
summary(step_f)
step_f$anova # display results

ggplot(data=May_2015, aes(x=Avena, y=ANPPgm, group=subplot, color=subplot))+
stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  geom_point()+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method="lm", formula= y ~ x, se=FALSE)

