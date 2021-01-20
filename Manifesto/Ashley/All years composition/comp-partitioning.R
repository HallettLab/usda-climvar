library(tidyverse)
library(nlme)
library(ggplot2)
library(multcomp)
library(lsmeans)
library(MuMIn)
library(RColorBrewer)
setwd("~/Dropbox/ClimVar/DATA/")
May_ANPP<-read.csv("Plant_composition_data/ANPP/ANPP_CleanedData/ClimVar_ANPP-peak.csv")

head(May_ANPP)

#change plots, years, shelter to factors
May_ANPP[,'plot'] <- as.factor(as.character(May_ANPP[,'plot']))
May_ANPP[,'year'] <- as.factor(as.character(May_ANPP[,'year']))
May_ANPP[,'shelter'] <- as.factor(as.character(May_ANPP[,'shelter']))

compensation (H2)
#H2 will be confirmed if mixed plots have greater ANPP compared to grass only or forb only
#first remove compost treatment 
May_ANPP_noC<-filter(May_ANPP, subplot!='C', subplot!='XC')
May_ANPP_noC <- May_ANPP_noC %>% mutate(treatment=ordered(treatment, levels = c(controlRain="controlRain", springDry="springDry", fallDry="fallDry", consistentDry="consistentDry")))

m3<-lme(weight_g_m ~treatment*subplot, random=~1|year/shelterBlock, May_ANPP_noC, na.action=na.exclude)
summary(m3)
anova(m3)
r.squaredGLMM(m3)
qqnorm(residuals(m3))
qqline(residuals(m3))
shapiro.test(residuals(m3))
LS3<-lsmeans(m3, ~subplot)
contrast(LS3, "pairwise")
#not normally distributed, try log transform

m4<-lme(log(weight_g_m+1) ~treatment*subplot, random=~1|year/shelterBlock, May_ANPP_noC, na.action=na.exclude)
summary(m4)
anova(m4)
r.squaredGLMM(m4)#7% of variation explained by fixed effects, 47% explained by entire model (lots of interannual variation?)
qqnorm(residuals(m4))
qqline(residuals(m4))
shapiro.test(residuals(m4))
#barely normal
LS4<-lsmeans(m4, ~subplot)
contrast(LS4, "pairwise")
#ANOVA: overall sign effect of subplot on ANPP
#lsmeans: B (mixed plots) have greater ANPP compared to F (forb-only), but not G (grass-only)
#no sig difference in ANPP between B (mixed plots) and XC (no manipulation)
#plot the effect of species composition on ANPP:
ggplot(d=May_ANPP_noC, aes(x=treatment, y=weight_g_m, fill=subplot)) +
  facet_wrap(~year)+
  theme_linedraw()+
  labs(x="composition treatment", y="ANPP g/m2")+
  geom_boxplot(aes(y=weight_g_m), shape=16)

ggplot(d=May_ANPP_noC, aes(x=subplot, y=weight_g_m)) +
  theme_linedraw()+
  labs(x="", y="ANPP g/m2")+
  annotate("text", x= c("B", "F","G"), y = c(1200, 1250, 1200), label = c("a", "b", "ab"), color = "black") +
  geom_boxplot(aes(y=weight_g_m), shape=16)


#rename trt groups
May_ANPP_noC$treatment2[May_ANPP_noC$treatment=="controlRain"] <- "Ambient"
May_ANPP_noC$treatment2[May_ANPP_noC$treatment=="springDry"] <- "Late Drought"
May_ANPP_noC$treatment2[May_ANPP_noC$treatment=="fallDry"] <- "Early Drought"
May_ANPP_noC$treatment2[May_ANPP_noC$treatment=="consistentDry"] <- "Consistent Drought"

ggplot(d=subset(May_ANPP_noC, subplot=='B'), aes(x=treatment2, y=weight_g_m, fill=treatment2)) +
  #facet_wrap(~year)+
  ggtitle("Observed Yield (Yo)")+
  theme_linedraw()+
  theme(legend.position="none")+
  scale_fill_manual(values = c("royalblue2","sienna","lightsteelblue1", "peachpuff"), guide = guide_legend(title = "Treatment")) +
  labs(x="Rainfall Treatment", y="ANPP g/m2")+
  geom_boxplot(aes(y=weight_g_m), shape=16)

FG_ANPP<-read.csv("Plant_composition_data/ANPP/ANPP_CleanedData/ClimVar_ANPP-func-groups.csv")
FG_ANPP<-FG_ANPP %>% dplyr::select(-totweight,-propweight)%>%
  group_by(year, subplot, plot, treatment, date) %>%
  mutate(totweight = sum(weight_g_m), propweight = weight_g_m/totweight) %>%
  tbl_df() %>% ungroup()%>%
  mutate(treatment=ordered(treatment, levels = c(consistentDry="consistentDry", springDry="springDry", fallDry="fallDry", controlRain="controlRain"))) 


##FORB##
#create subset with no species manipulations (control community) only, first harvest b/c forbs most abundant then
FG_forb<-filter(FG_ANPP, subplot=='F', harvest=="First", func=="Forb") %>% dplyr::select(-harvest, -func, -date) %>% 
  mutate(treatment=ordered(treatment, levels = c(controlRain="controlRain", springDry="springDry", fallDry="fallDry", consistentDry="consistentDry")))

FG_B.f<-filter(FG_ANPP, subplot=='B', harvest=="First", func=="Forb") %>% dplyr::select(-harvest, -func, -date)
FG_B.f$M<-FG_forb$weight_g_m
FG_B.f<-FG_B.f %>% group_by(treatment, shelterBlock, year)%>%mutate(RYof=propweight*M)%>%mutate(RYef=0.25*M)%>%mutate(dYf=RYof-RYef)%>%
  mutate(dRY.f=propweight-0.25)

FG_grass<-filter(FG_ANPP, subplot=='G', harvest=="Second", func=="Grass") %>% dplyr::select(-harvest, -func, -date)%>% 
  mutate(treatment=ordered(treatment, levels = c(controlRain="controlRain", springDry="springDry", fallDry="fallDry", consistentDry="consistentDry")))
FG_B.g<-filter(FG_ANPP, subplot=='B', harvest=="Second", func=="Grass") %>% dplyr::select(-harvest, -func, -date)
FG_B.g$M<-FG_grass$weight_g_m
FG_B.g<-FG_B.g %>% group_by(treatment, shelterBlock, year)%>%mutate(RYog=propweight*M)%>%mutate(RYeg=0.75*M)%>%mutate(dYg=RYog-RYeg)%>%
  mutate(dRY.g=(propweight-0.75))
FG_B.g$dYf<-FG_B.f$dYf
FG_B.g$dRY.f<-FG_B.f$dRY.f
FG_B.g <- FG_B.g %>% mutate (select=2*(((dRY.f+dRY.g))*((M+FG_B.f$M))))%>%mutate(x=dRY.f+dRY.g)%>% mutate(y=(M+FG_B.f$M))
FG_B.g<-FG_B.g%>%group_by(treatment, shelterBlock)%>%mutate(complement=2*cov(x,y))


FG_B.g<-FG_B.g %>% group_by(treatment, shelterBlock, year)%>%mutate(dY=dYg+dYf)%>% ungroup()%>%
  mutate(treatment=ordered(treatment, levels = c(controlRain="controlRain", springDry="springDry", fallDry="fallDry", consistentDry="consistentDry")))


#rename trt groups
FG_B.g$treatment2[FG_B.g$treatment=="controlRain"] <- "Ambient"
FG_B.g$treatment2[FG_B.g$treatment=="springDry"] <- "Late Drought"
FG_B.g$treatment2[FG_B.g$treatment=="fallDry"] <- "Early Drought"
FG_B.g$treatment2[FG_B.g$treatment=="consistentDry"] <- "Consistent Drought"

ggplot(d=FG_B.g, aes(x=treatment2, y=dY, fill=treatment2)) +
  #facet_wrap(~year)+
  ggtitle("Net Effect ")+
  scale_fill_manual(values = c("royalblue2","sienna","lightsteelblue1", "peachpuff"), guide = guide_legend(title = "Treatment")) +
  theme_linedraw()+
  theme(legend.position="none")+
  labs(x="Rainfall Treatment", y="Net Effect (ΔY)")+
  geom_boxplot(aes(y=dY), shape=16)

ggplot(d=FG_B.g, aes(x=treatment2, y=select, fill=treatment2)) +
  #facet_wrap(~year)+
  ggtitle("Selection Effect")+
  scale_fill_manual(values = c("royalblue2","sienna","lightsteelblue1", "peachpuff"), guide = guide_legend(title = "Treatment")) +
  theme_linedraw()+
  theme(legend.position="none")+
  labs(x="Rainfall treatment", y="Selection Effect (NΔRYM)")+
  geom_boxplot(aes(y=select), shape=16)

ggplot(d=FG_B.g, aes(x=treatment2, y=complement, fill=treatment2)) +
  #facet_wrap(~year)+
  ggtitle("Complementarity Effect")+
  scale_fill_manual(values = c("royalblue2","sienna","lightsteelblue1", "peachpuff"), guide = guide_legend(title = "Treatment")) +
  theme_linedraw()+
  theme(legend.position="none")+
  labs(x="Rainfall Treatment", y="Complementarity Effect (N*cov(ΔRY, M)")+
  geom_boxplot(aes(y=complement), shape=16)

#rename trt groups
gf_graphic2$treatment2[gf_graphic2$treatment=="controlRain"] <- "Ambient"
gf_graphic2$treatment2[gf_graphic2$treatment=="springDry"] <- "Late Drought"
gf_graphic2$treatment2[gf_graphic2$treatment=="fallDry"] <- "Early Drought"
gf_graphic2$treatment2[gf_graphic2$treatment=="consistentDry"] <- "Consistent Drought"

ggplot(data = subset(gf_graphic2, subplot %in% c("B","G","F")), aes(x=treatment2, y=meancover, color=func)) + 
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
rao_anpp<-ggplot(aboveFD_2b, aes(x=RaoQ, y=weight_g_m, group=treatment, color=treatment))+
  #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
  #annotate("text", x= c(1500, 6000, 9000), y = c(750, 450, 315), label = c("R2=0.31", "R2=0.04", "R2=0.09"), color = c("brown2", "dodgerblue3","forestgreen")) +
  geom_point()+
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
Lcols <- rep(c("indianred4",  "dodgerblue1", "darkgoldenrod"))
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