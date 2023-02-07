veg <- read.csv("./Dropbox/ClimVar/DATA/Plant_composition_data/Cover/Cover_CleanedData/ClimVar_2015_species-cover_wide.csv", header = TRUE)
traits <- read.csv("./Dropbox/ClimVar/DATA/Plant_composition_data/Traits/Traits_GHscreening_Lina_Ashley/Combined_GHtraits.csv", header = TRUE)
BNPP <- read.csv("./Dropbox/ClimVar/DATA/Plant_composition_data/BNPP/BNPP_CleanedData/BNPP_MayHarvest_2015.csv", header = TRUE)
ANPP <- read.csv("./Dropbox/ClimVar/DATA/Plant_composition_data/ANPP/ANPP_CleanedData/ClimVar_ANPP-peak.csv", header = TRUE)

library(tidyverse)
library(FD)
library(vegan)
library(ggplot2)
library(ggpubr)
library(MASS)

###Set up data for FD package###
#Species names must be in same order in both files
veg_keys <- veg[,1:7]
vegdat <- veg[,8:64]
vegdat1 <- vegdat[, -c(10,15:16,20,23:25,28,31,32,34,36:39,41:44,46,53,57)] #remove all columns that we do not have trait data for and any species that has no cover in 2015
#group some species
vegdat1$BRODIA <- vegdat1$Bromus.diandrus + vegdat1$Bromus.sterilis
vegdat1$HORMUR <- vegdat1$Hordeum.marinum + vegdat1$Hordeum.murinum
vegdat1$TRISP <- vegdat1$Trifolium.sp. + vegdat1$Trifolium.wildenovii
vegdat2 <- vegdat1[ , -c(7,10,19,20,30,32)] #remove grouped columns 
colnames(vegdat2) <- c("ACHMIL","ANAARV", "AVEBAR", "AVEFAT", "BRADIS", "BRIMIN", "BROHOR", "BROMAD","CARPYC", "CENSOL","CERGLO", 
                       "CYNDAC", "CYNECH", "EROBOT", "EROMOS", "FILGAL", "HYPGLA", "HYPRAD", 
                       "LACSER", "LOLMUL","RUMPUL", "TAECAP", "TRIDUB", "TRIGLO", "TRIHIR", "TRISUB", "VICSAT", "VULBRO", "VULMYU",
                       "BRODIA", "HORMUR", "TRISP") #rename columns
vegdat3 <- vegdat2[c("ACHMIL","ANAARV", "AVEBAR", "AVEFAT", "BRADIS","BRIMIN", "BRODIA","BROHOR", "BROMAD","CARPYC", "CENSOL","CERGLO", 
                     "CYNDAC", "CYNECH", "EROBOT", "EROMOS",  "FILGAL", "HORMUR","HYPGLA", "HYPRAD", 
                     "LACSER", "LOLMUL","RUMPUL", "TAECAP","TRIDUB", "TRIGLO", "TRIHIR", "TRISP", "TRISUB", "VICSAT", "VULBRO", "VULMYU"
                      )] #reorder columns
dim(vegdat3) #check the dimension of veg data
comp <- as.matrix(vegdat3[1:nrow(vegdat3), 1:ncol(vegdat3)])

traits$ID <- as.character(traits$ID) #set ID column as character
traits$ID <- ifelse(traits$ID == "TRIREP", "TRISP", traits$ID) #use Trifolium repens traits for Trifolium sp. 
all_traits <- traits %>%
  filter(ID %in% c("ACHMIL","ANAARV", "AVEBAR", "AVEFAT", "BRADIS","BRIMIN", "BRODIA","BROHOR", "BROMAD","CARPYC", "CENSOL","CERGLO", 
                   "CYNDAC", "CYNECH", "EROBOT", "EROMOS",  "FILGAL", "HORMUR","HYPGLA", "HYPRAD", 
                   "LACSER", "LOLMUL","RUMPUL", "TAECAP","TRIDUB", "TRIGLO", "TRIHIR", "TRISP", "TRISUB", "VICSAT", "VULBRO", "VULMYU")) %>%
  dplyr::select(ID, SLA, LDMC, Ht, Dens, DiamC, SRLC, SRLF, PropF)
dim(all_traits) #check the dimension of trait data
tr <- as.matrix(all_traits[1:nrow(all_traits), 2:ncol(all_traits)])
row.names(tr) <- all_traits$ID

above_traits <- traits %>%
  filter(ID %in% c("ACHMIL","ANAARV", "AVEBAR", "AVEFAT", "BRADIS","BRIMIN", "BRODIA","BROHOR", "BROMAD","CARPYC", "CENSOL","CERGLO", 
                   "CYNDAC", "CYNECH", "EROBOT", "EROMOS",  "FILGAL", "HORMUR","HYPGLA", "HYPRAD", 
                   "LACSER", "LOLMUL","RUMPUL", "TAECAP","TRIDUB", "TRIGLO", "TRIHIR", "TRISP", "TRISUB", "VICSAT", "VULBRO", "VULMYU")) %>%
  dplyr::select(ID, SLA, LDMC, Ht)
dim(above_traits) #check the dimension of trait data
above_tr <- as.matrix(above_traits[1:nrow(above_traits), 2:ncol(above_traits)])
row.names(above_tr) <- above_traits$ID

below_traits <- traits %>%
  filter(ID %in% c("ACHMIL","ANAARV", "AVEBAR", "AVEFAT", "BRADIS","BRIMIN", "BRODIA","BROHOR", "BROMAD","CARPYC", "CENSOL","CERGLO", 
                   "CYNDAC", "CYNECH", "EROBOT", "EROMOS",  "FILGAL", "HORMUR","HYPGLA", "HYPRAD", 
                   "LACSER", "LOLMUL","RUMPUL", "TAECAP","TRIDUB", "TRIGLO", "TRIHIR", "TRISP", "TRISUB", "VICSAT", "VULBRO", "VULMYU")) %>%
  dplyr::select(ID, Dens, DiamC, SRLC, SRLF, PropF)
dim(below_traits) #check the dimension of trait data
below_tr <- as.matrix(below_traits[1:nrow(below_traits), 2:ncol(below_traits)])
row.names(below_tr) <- below_traits$ID

#Calculate ind RaoQ and ind FEve
SLA_results <- dbFD(tr[,1], comp, corr = "cailliez")
SLA_results <- as.data.frame(SLA_results) %>%
  tbl_df()

LDMC_results <- dbFD(tr[,2], comp, corr = "cailliez", calc.CWM = TRUE)
LDMC_results <- as.data.frame(LDMC_results) %>%
  tbl_df()

Ht_results <- dbFD(tr[,3], comp, corr = "cailliez")
Ht_results <- as.data.frame(Ht_results) %>%
  tbl_df()

Dens_results <- dbFD(tr[,4], comp, corr="cailliez")
Dens_results <- as.data.frame(Dens_results) %>%
  tbl_df()

DiamC_results <- dbFD(tr[,5], comp, corr="cailliez")
DiamC_results <- as.data.frame(DiamC_results) %>%
  tbl_df()

SRLC_results <- dbFD(tr[,6], comp, corr="cailliez")
SRLC_results <- as.data.frame(SRLC_results) %>%
  tbl_df()

SRLF_results <- dbFD(tr[,7], comp, corr="cailliez")
SRLF_results <- as.data.frame(SRLF_results) %>%
  tbl_df()

PropF_results <- dbFD(tr[,8], comp, corr="cailliez")
PropF_results <- as.data.frame(PropF_results)%>%
  tbl_df()

overall_rao <- dbFD(tr, comp, corr="cailliez")
overall_rao <- as.data.frame(overall_rao) %>%
  tbl_df()

#Above- and belowground trait RaoQ and FEve
above_fd <- dbFD(above_tr, comp, corr="cailliez")
above_fd <- as.data.frame(above_fd) %>%
  tbl_df()

below_fd <- dbFD(below_tr, comp, corr="cailliez")
below_fd <- as.data.frame(below_fd) %>%
  tbl_df()

##Join ind RaoQ 
rao_table <- as.data.frame(cbind(veg_keys[,2:5],SLA_results[,7], LDMC_results[,7],
                                 Ht_results[,7], Dens_results[,7], DiamC_results[,7],
                                 SRLC_results[,7], SRLF_results[,7], PropF_results[,7], overall_rao[,8]))
colnames(rao_table) <- c("plot", "subplot", "year", "treatment", "SLARaoQ", "LDMCRaoQ", "HtRaoQ", "DensRaoQ", 
                          "DiamCRaoQ", "SRLCRaoQ", "SRLFRaoQ", "PropFRaoQ", "RaoQ")

##Join ind FEve
FEve_table <- as.data.frame(cbind(veg_keys[,2:5],SLA_results[,5], LDMC_results[,5],
                                 Ht_results[,5], Dens_results[,5], DiamC_results[,5],
                                 SRLC_results[,5], SRLF_results[,5], PropF_results[,5], overall_rao[,5],
                                 overall_rao[,9:16]))
colnames(FEve_table) <- c("plot", "subplot", "year", "treatment", "SLAFEve", "LDMCFEve", "HtFEve", "DensFEve", 
                         "DiamCFEve", "SRLCFEve", "SRLFFEve", "PropFFEve", "FEve", "CWM.SLA", "CWM.LDMC", "CWM.Ht",
                         "CWM.Dens", "CWM.DiamC", "CWM.SRLC", "CWM.SRLF", "CWM.PropF")

#Join Above- and belowground trait RaoQ and FEve
FD_above_below <- as.data.frame(cbind(veg_keys[,2:5],above_fd[,3], above_fd[,5:8], below_fd[,3], below_fd[,5:8]))
colnames(FD_above_below) <- c("plot", "subplot", "year", "treatment",  
                              "aboveFRic", "aboveFEve", "aboveFDiv", "aboveFDis", "aboveRaoQ",
                              "belowFRic", "belowFEve", "belowFDiv", "belowFDis", "belowRaoQ")

#Filter 2015 ANPP
ANPP1 <- ANPP %>%
  filter(year == 2015) %>%
  filter(subplot %in% c("B", "F", "G")) %>%
  dplyr::select(plot, subplot, treatment, weight_g_m)

#Replace entry "20-Oct" to "10-20" in depth column in the BNPP dataset
BNPP$depth <- as.character(BNPP$depth) #set depth as character
BNPP <- BNPP %>%
  mutate(depth = replace(depth, depth == "20-Oct", "10-20")) %>%
  mutate(treatment_code = fct_recode(treatment, consDry = "consistentDry", control = "controlRain")) #make a column with shorter treatment names

#Calculate the aggregated BNPP in three levels of soil depth
BNPP1 <- BNPP %>%
  group_by(plot, subplot, treatment) %>% #group data
  summarise(agg_BNPP = sum(bmass_g_m2)) %>% #sum BNPP
  dplyr::select(plot, subplot, treatment, agg_BNPP)

#standardize ind RaoQ 
stand_rao_num <- decostand(rao_table[,5:13], "standardize") 
stand_rao <- as.data.frame(cbind(rao_table[,1:4], stand_rao_num))

#isolate CENSOL and TAECAP
CENSOL <- vegdat3$CENSOL
TAECAP <- vegdat3$TAECAP

#Join RaoQ and biomass
joined_rao <- stand_rao %>%
  left_join(ANPP1, by = c("plot", "subplot", "treatment")) %>%
  left_join(BNPP1, by = c("plot", "subplot", "treatment")) %>%
  mutate(total = weight_g_m + agg_BNPP) %>%
  cbind(CENSOL) %>%
  cbind(TAECAP) 

#Join FEve and biomass
joined_FEve <- FEve_table %>%
  left_join(ANPP1, by = c("plot", "subplot", "treatment")) %>%
  left_join(BNPP1, by = c("plot", "subplot", "treatment")) %>%
  mutate(total = weight_g_m + agg_BNPP) %>%
  cbind(CENSOL) %>%
  cbind(TAECAP) 

#Join community FD metrics and biomass
FD_table <- as.data.frame(cbind(veg_keys[,2:5], overall_rao))
joined_FD <- FD_table %>%
  left_join(ANPP1, by = c("plot", "subplot", "treatment")) %>%
  left_join(BNPP1, by = c("plot", "subplot", "treatment")) %>%
  mutate(total = weight_g_m + agg_BNPP)%>%
  cbind(CENSOL) %>%
  cbind(TAECAP) 

#Join community FD metrics and biomass
joined_FD_above_below <- FD_above_below %>%
  left_join(ANPP1, by = c("plot", "subplot", "treatment")) %>%
  left_join(BNPP1, by = c("plot", "subplot", "treatment")) %>%
  mutate(total = weight_g_m + agg_BNPP)%>%
  cbind(CENSOL) %>%
  cbind(TAECAP)
levels(joined_FD_above_below$subplot) <- list(Mixed = "B", Forb = "F", Grass = "G")

#Correlation FD indices and CWM traits
pairs(~FRic + FEve + FDis + RaoQ , data = FD_table)
pairs(~FEve + CWM.SLA + CWM.LDMC + CWM.Ht + CWM.Dens + CWM.DiamC + CWM.SRLC + CWM.SRLF + CWM.PropF, data = FD_table)
pairs(~RaoQ + CWM.SLA + CWM.LDMC + CWM.Ht + CWM.Dens + CWM.DiamC + CWM.SRLC + CWM.SRLF + CWM.PropF, data = FD_table)

#Heatmap of productivity in FEve and RaoQ space
f_heat_below <- ggplot(joined_FD_above_below, aes(x = belowRaoQ, y = belowFEve))+
  geom_point(aes(size = agg_BNPP, col = agg_BNPP))+
  theme_classic()+
  stat_density_2d()+
  xlab("Belowground RaoQ")+
  ylab("Belowground FEve")+
  guides( size = "none")+
  labs(col = "BNPP g/m2")
  
f_heat_above <- ggplot(joined_FD_above_below, aes(x = aboveRaoQ, y = aboveFEve))+
  geom_point(aes(col = weight_g_m, size = weight_g_m))+
  theme_classic()+
  stat_density_2d()+
  xlab("Aboveground RaoQ")+
  ylab("Aboveground FEve")+
  guides( size = "none")+
  labs(col = "ANPP g/m2")

f_heat_map <- ggarrange(f_heat_above, f_heat_below)

#Fig 2: ANPP and BNPP ~ FEve and RaoQ
p_FEve_above <- ggplot(joined_FD_above_below, aes(x = aboveFEve, y = weight_g_m)) +
  geom_point(aes(color = treatment, shape = subplot)) +
  theme_classic() +
  #ylim(50,900)+
  labs(y = bquote('ANPP'~(g/m^2)), x = "Aboveground FEve", color = "Rainfall Treatment", shape = "Composition Treatment") +
  #geom_smooth(aes(color = subplot, linetype = subplot), method = lm, size = 1, se = FALSE, fullrange = FALSE) +
  stat_cor(aes(group=1,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), show.legend = FALSE)+
  scale_color_manual(name = "Rainfall Treatment", labels = c("Control", "Spring Dry", "Fall Dry","Consistent Dry" ), values= c("#0070b8", "#b2c7e4", "#fccaaf", "#c85b23"))+
  #scale_linetype_manual( values = c("solid", "dashed", "dashed"), guide = "none")+
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE, color = "black", linetype = "solid") 
p_FEve_below <- ggplot(joined_FD_above_below, aes(x = belowFEve, y = agg_BNPP)) +
  geom_point(aes(color = treatment, shape = subplot)) +
  theme_classic() +
  #ylim(50,900)+
  labs(y = bquote('BNPP'~(g/m^2)), x = "Belowground FEve", color = "Rainfall Treatment", shape = "Composition Treatment") +
  #geom_smooth(aes(color = subplot, linetype = subplot), method = lm, size = 1, se = FALSE, fullrange = FALSE) +
  stat_cor(aes(group=1,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), show.legend = FALSE)+
  scale_color_manual(name = "Rainfall Treatment", labels = c("Control", "Spring Dry", "Fall Dry","Consistent Dry" ), values= c("#0070b8", "#b2c7e4", "#fccaaf", "#c85b23"))+
  #scale_linetype_manual( values = c("solid", "dashed", "dashed"), guide = "none")+
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE, color = "black", linetype = "solid") 
p_RaoQ_above <- ggplot(joined_FD_above_below, aes(x = aboveRaoQ, y = weight_g_m)) +
  geom_point(aes(color = treatment, shape = subplot)) +
  theme_classic() +
  #ylim(50,900)+
  labs(y = bquote('ANPP'~(g/m^2)), x = "Aboveground Rao's Q", color = "Rainfall Treatment", shape = "Composition Treatment") +
  #geom_smooth(aes(color = subplot, linetype = subplot), method = lm, size = 1, se = FALSE, fullrange = FALSE) +
  stat_cor(aes(group=1,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), show.legend = FALSE)+
  scale_color_manual(name = "Rainfall Treatment", labels = c("Control", "Spring Dry", "Fall Dry","Consistent Dry" ), values= c("#0070b8", "#b2c7e4", "#fccaaf", "#c85b23"))+
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE, color = "black", linetype = "solid")
p_RaoQ_below <- ggplot(joined_FD_above_below, aes(x = belowRaoQ, y = agg_BNPP)) +
  geom_point(aes(color = treatment, shape = subplot)) +
  theme_classic() +
  #ylim(50,900)+
  labs(y = bquote('BNPP'~(g/m^2)), x = "Belowground Rao's Q", color = "Rainfall Treatment", shape = "Composition Treatment") +
  #geom_smooth(aes(color = subplot, linetype = subplot), method = lm, size = 1, se = FALSE, fullrange = FALSE) +
  stat_cor(aes(group=1,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), show.legend = FALSE)+
  scale_color_manual(name = "Rainfall Treatment", labels = c("Control", "Spring Dry", "Fall Dry","Consistent Dry" ), values= c("#0070b8", "#b2c7e4", "#fccaaf", "#c85b23"))+
  #scale_linetype_manual( values = c("solid", "dashed", "dashed"), guide = "none")+
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE, color = "black", linetype = "solid")

Figure2 <- ggarrange(p_FEve_above, p_FEve_below, p_RaoQ_above, p_RaoQ_below, ncol = 2, nrow = 2, common.legend = TRUE, legend = "right", labels = c("a)", "b)", "c)", "d)"))

#Treatment effects on FD
FD_seeding_tr <- joined_FD_above_below %>% 
  group_by(subplot) %>% 
  summarize(aboveRaoQ.m=mean(aboveRaoQ), belowRaoQ.m=mean(belowRaoQ), 
            se.aboveRaoQ=calcSE(aboveRaoQ),se.belowRaoQ=calcSE(belowRaoQ),
            aboveFEve.m=mean(aboveFEve), belowFEve.m=mean(belowFEve),
            se.aboveFEve=calcSE(aboveFEve), se.belowFEve=calcSE(belowFEve))

FD_rain_tr <- joined_FD_above_below %>% 
  group_by(treatment) %>% 
  summarize(aboveRaoQ.m=mean(aboveRaoQ), belowRaoQ.m=mean(belowRaoQ), 
            se.aboveRaoQ=calcSE(aboveRaoQ),se.belowRaoQ=calcSE(belowRaoQ),
            aboveFEve.m=mean(aboveFEve), belowFEve.m=mean(belowFEve),
            se.aboveFEve=calcSE(aboveFEve), se.belowFEve=calcSE(belowFEve))

fit_seeding_a <- aov(aboveFEve~subplot, joined_FD_above_below)
TukeyHSD(aov(fit_seeding_a))
fit_seeding_b <- aov(belowFEve~subplot, joined_FD_above_below)
TukeyHSD(aov(fit_seeding_b))
fit_seeding_c <- aov(aboveRaoQ~subplot, joined_FD_above_below)
TukeyHSD(aov(fit_seeding_c))
fit_seeding_d <- aov(belowRaoQ~subplot, joined_FD_above_below)
TukeyHSD(aov(fit_seeding_d))

fit_rain_a <- aov(aboveFEve~treatment, joined_FD_above_below)
TukeyHSD(aov(fit_rain_a))
fit_rain_b <- aov(belowFEve~treatment, joined_FD_above_below)
TukeyHSD(aov(fit_rain_b))
fit_rain_c <- aov(aboveRaoQ~treatment, joined_FD_above_below)
TukeyHSD(aov(fit_rain_c))
fit_rain_d <- aov(belowRaoQ~treatment, joined_FD_above_below)
TukeyHSD(aov(fit_rain_d))

summary(aov(aboveFEve~subplot*treatment, joined_FD_above_below))
summary(aov(belowFEve~subplot*treatment, joined_FD_above_below))
summary(aov(aboveRaoQ~subplot*treatment, joined_FD_above_below))
summary(aov(belowRaoQ~subplot*treatment, joined_FD_above_below))

f_fd_a<-ggplot(data=FD_seeding_tr, aes(x=subplot, y=aboveFEve.m, col = as.factor(subplot)))+
     geom_point(position=position_dodge(0.9),size=4)+
     labs(y=expression(atop("Aboveground FEve")), color = "Treatment")+
     theme(axis.title.x=element_blank(),
                     text = element_text(size=14),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(),
                     axis.line = element_line(colour = "black"),
                     panel.border = element_rect(colour = "black", fill = NA, size = 1.2))+
     scale_x_discrete(labels = c("Mixed","Forb", "Grass"))+
     scale_color_manual(labels = c("Mixed","Forb", "Grass"),  values = c("#fc8d62", "#66c2a5", "#8da0cb"))+
     geom_errorbar(aes(ymin = aboveFEve.m-se.aboveFEve, ymax = aboveFEve.m+se.aboveFEve), size = 0.5, width = 0,position=position_dodge(0.9))+
     annotate("text", x = 1, y = 0.6, label = "a", size = 6) +
     annotate("text", x = 2, y = 0.6, label = "a", size = 6) +
     annotate("text", x = 3, y = 0.6, label = "b", size = 6)
f_fd_b<-ggplot(data=FD_seeding_tr, aes(x=subplot, y=belowFEve.m, col = as.factor(subplot)))+
  geom_point(position=position_dodge(0.9),size=4)+
  labs(y=expression(atop("Belowground FEve")), color = "Treatment")+
  theme(axis.title.x=element_blank(),
                 text = element_text(size=14),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 axis.line = element_line(colour = "black"),
                 panel.border = element_rect(colour = "black", fill = NA, size = 1.2))+
  scale_x_discrete(labels = c("Mixed","Forb", "Grass"))+
  scale_color_manual(labels = c("Mixed","Forb", "Grass"),  values = c("#fc8d62", "#66c2a5", "#8da0cb"))+
  geom_errorbar(aes(ymin = belowFEve.m-se.belowFEve, ymax = belowFEve.m+se.belowFEve), size = 0.5, width = 0,position=position_dodge(0.9))+
  annotate("text", x = 1, y = 0.55, label = "ab", size = 6) +
  annotate("text", x = 2, y = 0.55, label = "a", size = 6) +
  annotate("text", x = 3, y = 0.55, label = "b", size = 6)
f_fd_c<-ggplot(data=FD_seeding_tr, aes(x=subplot, y=aboveRaoQ.m, col = as.factor(subplot)))+
  geom_point(position=position_dodge(0.9),size=4)+
  labs(y=expression(atop("Aboveground Rao's Q")), color = "Treatment")+
  theme(axis.title.x=element_blank(),
        text = element_text(size=14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2))+
  scale_x_discrete(labels = c("Mixed","Forb", "Grass"))+
  scale_color_manual(labels = c("Mixed","Forb", "Grass"),  values = c("#fc8d62", "#66c2a5", "#8da0cb"))+
  geom_errorbar(aes(ymin = aboveRaoQ.m-se.aboveRaoQ, ymax = aboveRaoQ.m+se.aboveRaoQ), size = 0.5, width = 0,position=position_dodge(0.9))+
  annotate("text", x = 1, y = 1.8, label = "b", size = 6) +
  annotate("text", x = 2, y = 1.8, label = "a", size = 6) +
  annotate("text", x = 3, y = 1.8, label = "b", size = 6)
f_fd_d<-ggplot(data=FD_seeding_tr, aes(x=subplot, y=belowRaoQ.m, col = as.factor(subplot)))+
  geom_point(position=position_dodge(0.9),size=4)+
  labs(y=expression(atop("Belowground Rao's Q")), color = "Treatment")+
  theme(axis.title.x=element_blank(),
        text = element_text(size=14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2))+
  scale_x_discrete(labels = c("Mixed","Forb", "Grass"))+
  scale_color_manual(labels = c("Mixed","Forb", "Grass"),  values = c("#fc8d62", "#66c2a5", "#8da0cb"))+
  geom_errorbar(aes(ymin = belowRaoQ.m-se.belowRaoQ, ymax = belowRaoQ.m+se.belowRaoQ), size = 0.5, width = 0,position=position_dodge(0.9))+
  annotate("text", x = 1, y = 3.1, label = "b", size = 6) +
  annotate("text", x = 2, y = 3.1, label = "a", size = 6) +
  annotate("text", x = 3, y = 3.1, label = "b", size = 6)
f_fd_e<-ggplot(data=FD_rain_tr, aes(x=treatment, y=aboveFEve.m, col = treatment))+
  geom_point(position=position_dodge(0.9), size=4)+
  labs(y=expression(atop("Aboveground FEve")))+
  scale_x_discrete(labels = c("Control\n ","Spring\n Dry\n", "Fall\n Dry\n", "Consistent\n Dry\n"))+
  scale_color_manual(name = "Treatment", labels = c("Control", "Spring Dry", "Fall Dry","Consistent Dry" ), values= c("#0070b8", "#b2c7e4", "#fccaaf", "#c85b23"))+
  geom_errorbar(aes(ymin = aboveFEve.m-se.aboveFEve, ymax = aboveFEve.m+se.aboveFEve), size = 0.5, width = 0,position=position_dodge(0.9))+
  theme(axis.title.x=element_blank(),
        text = element_text(size=14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2))+
  annotate("text", x = 1, y = 0.6, label = "a", size = 6) +
  annotate("text", x = 2, y = 0.6, label = "a", size = 6) +
  annotate("text", x = 3, y = 0.6, label = "a", size = 6) +
  annotate("text", x = 4, y = 0.6, label = "a", size = 6)

f_fd_e<-ggplot(data=FD_rain_tr, aes(x=treatment, y=aboveFEve.m, col = treatment))+
  geom_point(position=position_dodge(0.9), size=4)+
  labs(y=expression(atop("Aboveground FEve")))+
  scale_x_discrete(labels = c("Control\n ","Spring\n Dry\n", "Fall\n Dry\n", "Consistent\n Dry\n"))+
  scale_color_manual(name = "Treatment", labels = c("Control", "Spring Dry", "Fall Dry","Consistent Dry" ), values= c("#0070b8", "#b2c7e4", "#fccaaf", "#c85b23"))+
  geom_errorbar(aes(ymin = aboveFEve.m-se.aboveFEve, ymax = aboveFEve.m+se.aboveFEve), size = 0.5, width = 0,position=position_dodge(0.9))+
  theme(axis.title.x=element_blank(),
        text = element_text(size=14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2))+
  annotate("text", x = 1, y = 0.6, label = "a", size = 6) +
  annotate("text", x = 2, y = 0.6, label = "a", size = 6) +
  annotate("text", x = 3, y = 0.6, label = "a", size = 6) +
  annotate("text", x = 4, y = 0.6, label = "a", size = 6)
f_fd_f<-ggplot(data=FD_rain_tr, aes(x=treatment, y=belowFEve.m, col = treatment))+
  geom_point(position=position_dodge(0.9), size=4)+
  labs(y=expression(atop("Belowground FEve")))+
  scale_x_discrete(labels = c("Control\n ","Spring\n Dry\n", "Fall\n Dry\n", "Consistent\n Dry\n"))+
  scale_color_manual(name = "Treatment", labels = c("Control", "Spring Dry", "Fall Dry","Consistent Dry" ), values= c("#0070b8", "#b2c7e4", "#fccaaf", "#c85b23"))+
  geom_errorbar(aes(ymin = belowFEve.m-se.belowFEve, ymax = belowFEve.m+se.belowFEve), size = 0.5, width = 0,position=position_dodge(0.9))+
  theme(axis.title.x=element_blank(),
        text = element_text(size=14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2))+
  annotate("text", x = 1, y = 0.53, label = "a", size = 6) +
  annotate("text", x = 2, y = 0.53, label = "a", size = 6) +
  annotate("text", x = 3, y = 0.53, label = "a", size = 6) +
  annotate("text", x = 4, y = 0.53, label = "a", size = 6)
f_fd_g<-ggplot(data=FD_rain_tr, aes(x=treatment, y=aboveRaoQ.m, col = treatment))+
  geom_point(position=position_dodge(0.9), size=4)+
  labs(y=expression(atop("Aboveground Rao's Q")))+
  scale_x_discrete(labels = c("Control\n ","Spring\n Dry\n", "Fall\n Dry\n", "Consistent\n Dry\n"))+
  scale_color_manual(name = "Treatment", labels = c("Control", "Spring Dry", "Fall Dry","Consistent Dry" ), values= c("#0070b8", "#b2c7e4", "#fccaaf", "#c85b23"))+
  geom_errorbar(aes(ymin = aboveRaoQ.m-se.aboveRaoQ, ymax = aboveRaoQ.m+se.aboveRaoQ), size = 0.5, width = 0,position=position_dodge(0.9))+
  theme(axis.title.x=element_blank(),
        text = element_text(size=14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2))+
  annotate("text", x = 1, y = 1.4, label = "a", size = 6) +
  annotate("text", x = 2, y = 1.4, label = "a", size = 6) +
  annotate("text", x = 3, y = 1.4, label = "a", size = 6) +
  annotate("text", x = 4, y = 1.4, label = "a", size = 6)
f_fd_h<-ggplot(data=FD_rain_tr, aes(x=treatment, y=belowRaoQ.m, col = treatment))+
  geom_point(position=position_dodge(0.9), size=4)+
  labs(y=expression(atop("Belowground Rao's Q")))+
  scale_x_discrete(labels = c("Control\n ","Spring\n Dry\n", "Fall\n Dry\n", "Consistent\n Dry\n"))+
  scale_color_manual(name = "Treatment", labels = c("Control", "Spring Dry", "Fall Dry","Consistent Dry" ), values= c("#0070b8", "#b2c7e4", "#fccaaf", "#c85b23"))+
  geom_errorbar(aes(ymin = belowRaoQ.m-se.belowRaoQ, ymax = belowRaoQ.m+se.belowRaoQ), size = 0.5, width = 0,position=position_dodge(0.9))+
  theme(axis.title.x=element_blank(),
        text = element_text(size=14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2))+
  annotate("text", x = 1, y = 2.2, label = "a", size = 6) +
  annotate("text", x = 2, y = 2.2, label = "a", size = 6) +
  annotate("text", x = 3, y = 2.2, label = "a", size = 6) +
  annotate("text", x = 4, y = 2.2, label = "a", size = 6)

fig_FD_subplot <- ggarrange(f_fd_a, f_fd_b, f_fd_c, f_fd_d,
                       ncol =2, nrow =2, 
                       align = "v",labels = c("a)", "b)", "c)", "d)"), common.legend = TRUE, legend = "top")
fig_FD_tr <- ggarrange(f_fd_e, f_fd_f, f_fd_g, f_fd_h,
                     ncol =2, nrow =2, 
                     align = "v",labels = c("e)", "f)", "g)", "h)"), common.legend = TRUE, legend = "top")
ggarrange(fig_FD_subplot, fig_FD_tr, ncol = 1, nrow = 2, heights = c(1, 1.2))

# #Stepwise regression of BNPP ~ Ind RaoQ
# model8 <- lm(total ~ CENSOL + TAECAP + SLARaoQ + LDMCRaoQ + HtRaoQ + DensRaoQ + DiamCRaoQ + SRLCRaoQ + SRLFRaoQ + PropFRaoQ, joined_rao)
# step_all <- stepAIC(model8, direction = "backward", trace = FALSE)
# step_all$anova
# model9 <- lm(total ~ SLARaoQ , joined_rao)
# summary(model9)
# 
# model10 <- lm(weight_g_m ~ CENSOL + TAECAP +SLARaoQ + LDMCRaoQ + HtRaoQ , joined_rao)
# step_all <- stepAIC(model10, direction = "backward", trace = FALSE)
# step_all$anova
# model11 <- lm(weight_g_m ~ SLARaoQ + LDMCRaoQ , joined_rao)
# summary(model11)
# 
# model12 <- lm(agg_BNPP ~ CENSOL + TAECAP +DensRaoQ + DiamCRaoQ + SRLCRaoQ + SRLFRaoQ + PropFRaoQ, joined_rao)
# step_all <- stepAIC(model12, direction = "backward", trace = FALSE)
# step_all$anova
# model13 <- lm(agg_BNPP ~ DiamCRaoQ , joined_rao)
# summary(model13)

#Stepwise regression of BNPP ~ Ind FEve
model14 <- lm(total ~ SLAFEve + LDMCFEve + HtFEve + DensFEve + DiamCFEve + SRLCFEve + SRLFFEve + PropFFEve, joined_FEve)
summary(model14)
model15 <- lm(total ~ CWM.SLA + CWM.LDMC + CWM.Ht + CWM.Dens + CWM.DiamC + CWM.SRLC + CWM.SRLF + CWM.PropF, joined_FEve)
summary(model15)
model16 <- lm(total ~ CENSOL + TAECAP, joined_FEve)
summary(model16)
model17 <- lm(total ~ SLAFEve + LDMCFEve + HtFEve + DensFEve + DiamCFEve + SRLCFEve + SRLFFEve + PropFFEve + CWM.SLA + CWM.LDMC + CWM.Ht + CWM.Dens + CWM.DiamC + CWM.SRLC + CWM.SRLF + CWM.PropF + TAECAP + CENSOL, joined_FEve)
step_all <- stepAIC(model17, direction = "backward", trace = FALSE)
step_all$anova
model_total <- lm(total ~ SLAFEve + CWM.SLA + CWM.Ht + CWM.DiamC + CENSOL, joined_FEve)
model18 <- lm(total ~ CWM.SRLC + CWM.SRLF, joined_FEve)
step_all <- stepAIC(model18, direction = "backward", trace = FALSE)
step_all$anova
summary(model_total)

model14 <- lm(weight_g_m ~ SLAFEve + LDMCFEve + HtFEve + DensFEve + DiamCFEve + SRLCFEve + SRLFFEve + PropFFEve, joined_FEve)
summary(model14)
model15 <- lm(weight_g_m ~ CWM.SLA + CWM.LDMC + CWM.Ht + CWM.Dens + CWM.DiamC + CWM.SRLC + CWM.SRLF + CWM.PropF, joined_FEve)
summary(model15)
model16 <- lm(weight_g_m ~ CENSOL + TAECAP, joined_FEve)
summary(model16)
model17 <- lm(weight_g_m ~ DiamCFEve +  PropFFEve + CWM.SRLF, joined_FEve)
step_all <- stepAIC(model17, direction = "backward", trace = FALSE)
step_all$anova


model13 <- lm(agg_BNPP ~ DiamCRaoQ , joined_rao)
summary(model13)

#community RaoQ with ANPP and BNPP 
fitoverall <- lm(weight_g_m ~ RaoQ, joined_rao) 
summary(fitoverall) 
fitB <- lm(weight_g_m ~ RaoQ, joined_rao%>%filter(subplot == "B")) 
summary(fitB) 
fitF <- lm(weight_g_m ~ RaoQ, joined_rao%>%filter(subplot == "F"))
summary(fitF)
fitG <- lm(weight_g_m ~ RaoQ, joined_rao%>%filter(subplot == "G"))
summary(fitG)

fitoverall <- lm(agg_BNPP ~ RaoQ, joined_rao) 
summary(fitoverall)
fitB <- lm(agg_BNPP ~ RaoQ, joined_rao%>%filter(subplot == "B")) 
summary(fitB) 
fitF <- lm(agg_BNPP ~ RaoQ, joined_rao%>%filter(subplot == "F"))
summary(fitF)
fitG <- lm(agg_BNPP ~ RaoQ, joined_rao%>%filter(subplot == "G"))
summary(fitG)

#plot ANPP BNPP by RaoQ 
#By functional composition treatment Fig 2
p1 <- ggplot(joined_rao, aes(x = RaoQ, y = weight_g_m, color = subplot, linetype = subplot)) +
  geom_point() +
  theme_classic() +
  ylim(50,900)+
  labs(y = bquote('ANPP'~(g/m^2)), x = "Rao's Q Aboveground Traits", color = "Treatment") +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE) +
  #stat_cor(aes(group=subplot,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x.npc = 0.5)+
  scale_color_manual(name = "Treatment", labels = c("Mixed", "Forb", "Grass"), values= c("#fc8d62", "#66c2a5", "#8da0cb")) +
  scale_linetype_manual( values = c("solid", "dashed", "dashed"), guide = "none")

p2 <- ggplot(joined_rao, aes(x = RaoQ, y = agg_BNPP, color = subplot, linetype = subplot)) +
  geom_point() +
  theme_classic() +
  ylim(50, 900)+
  theme(legend.position = "none") +
  labs(y = bquote('BNPP'~(g/m^2)), x = "Rao's Q Belowground Traits", color = "Treatment") +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE) +
  #stat_cor(aes(group=subplot,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x.npc = 0.5)+
  scale_color_manual(name = "Treatment", labels = c("Mixed", "Forb", "Grass"), values= c("#fc8d62", "#66c2a5", "#8da0cb")) +
  scale_linetype_manual( values = c("dashed", "dashed", "dashed"), guide = "none")

joined_rao$treatment <- factor(joined_rao$treatment, levels =c("controlRain",  "springDry", "fallDry","consistentDry"))

p3 <- ggplot(joined_rao, aes(x = RaoQ, y = weight_g_m, color = treatment, linetype = treatment)) +
  geom_point() +
  theme_classic() +
  ylim(50,900)+
  labs(y = bquote('ANPP'~(g/m^2)), x = "Rao's Q Aboveground Traits", color = "Treatment") +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE) +
  #stat_cor(aes(group=treatment,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x.npc = 0.5)+
  scale_color_manual(name = "Treatment", labels = c("Control", "Spring Dry", "Fall Dry","Consistent Dry" ), values= c("#0070b8", "#b2c7e4", "#fccaaf", "#c85b23"))+
  scale_linetype_manual( values = c("dashed", "dashed", "dashed", "solid"), guide = "none")
p4 <- ggplot(joined_rao, aes(x = RaoQ, y = agg_BNPP, color = treatment, linetype = treatment)) +
  geom_point() +
  theme_classic() +
  ylim(50, 900)+
  theme(legend.position = "none") +
  labs(y = bquote('BNPP'~(g/m^2)), x = "Rao's Q Belowground Traits", color = "Treatment") +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE) +
  #stat_cor(aes(group=treatment,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x.npc = 0.5)+
  scale_color_manual(name = "Treatment", labels = c("Control", "Spring Dry", "Fall Dry","Consistent Dry" ), values= c("#0070b8", "#b2c7e4", "#fccaaf", "#c85b23"))+
  scale_linetype_manual( values = c("solid", "dashed", "dashed", "dashed"), guide = "none")

legend_comp <- get_legend(p1)
legend_rain <- get_legend(p3)
p1 <- p1 + theme(legend.position = "none")
p3 <- p3 + theme(legend.position = "none")

ggarrange(p1, p2, legend_comp, p3, p4, legend_rain, ncol = 3, nrow = 2, 
          align = "v",labels = c("a)", "b)", "", "c)", "d)", ""), 
          widths=c(1, 1, 0.3))

#community FEve with ANPP and BNPP 
fitoverall <- lm(weight_g_m ~ FEve, joined_FEve) 
summary(fitoverall) 
fitB <- lm(weight_g_m ~ FEve, joined_FEve%>%filter(subplot == "B")) 
summary(fitB) 
fitF <- lm(weight_g_m ~ FEve, joined_FEve%>%filter(subplot == "F"))
summary(fitF)
fitG <- lm(weight_g_m ~ FEve, joined_FEve%>%filter(subplot == "G"))
summary(fitG)

fitoverall <- lm(agg_BNPP ~ FEve, joined_FEve) 
summary(fitoverall)
fitB <- lm(agg_BNPP ~ FEve, joined_FEve%>%filter(subplot == "B")) 
summary(fitB) 
fitF <- lm(agg_BNPP ~ FEve, joined_FEve%>%filter(subplot == "F"))
summary(fitF)
fitG <- lm(agg_BNPP ~ FEve, joined_FEve%>%filter(subplot == "G"))
summary(fitG)

fitconstDry <- lm(weight_g_m ~ FEve, joined_FEve%>%filter(treatment == "consistentDry")) 
summary(fitconstDry) 
fitcontrol <- lm(weight_g_m ~ FEve, joined_FEve%>%filter(treatment == "controlRain"))
summary(fitcontrol)
fitfallDry <- lm(weight_g_m ~ FEve, joined_FEve%>%filter(treatment == "fallDry"))
summary(fitfallDry)
fitspringDry <- lm(weight_g_m ~ FEve, joined_FEve%>%filter(treatment == "springDry"))
summary(fitspringDry)

fitconstDry <- lm(agg_BNPP ~ FEve, joined_FEve%>%filter(treatment == "consistentDry")) 
summary(fitconstDry)
fitcontrol <- lm(agg_BNPP ~ FEve, joined_FEve%>%filter(treatment == "controlRain")) 
summary(fitcontrol) 
fitfallDry  <- lm(agg_BNPP ~ FEve, joined_FEve%>%filter(treatment == "fallDry"))
summary(fitfallDry)
fitspringDry <- lm(agg_BNPP ~ FEve, joined_FEve%>%filter(treatment == "springDry"))
summary(fitspringDry)

#plot ANPP BNPP by FEve 
#By functional composition treatment Fig 2 revision
p1_FEve <- ggplot(joined_FEve, aes(x = FEve, y = weight_g_m, color = subplot, linetype = subplot)) +
  geom_point() +
  theme_classic() +
  ylim(50,900)+
  labs(y = bquote('ANPP'~(g/m^2)), x = "Functional Evenness", color = "Treatment") +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE) +
  #stat_cor(aes(group=subplot,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x.npc = 0.5)+
  scale_color_manual(name = "Treatment", labels = c("Mixed", "Forb", "Grass"), values= c("#fc8d62", "#66c2a5", "#8da0cb")) +
  scale_linetype_manual( values = c("dashed", "dashed", "dashed"), guide = "none")

p2_FEve <- ggplot(joined_FEve, aes(x = FEve, y = agg_BNPP, color = subplot, linetype = subplot)) +
  geom_point() +
  theme_classic() +
  ylim(50, 900)+
  theme(legend.position = "none") +
  labs(y = bquote('BNPP'~(g/m^2)), x = "Functional Evenness", color = "Treatment") +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE) +
  #stat_cor(aes(group=subplot,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x.npc = 0.5)+
  scale_color_manual(name = "Treatment", labels = c("Mixed", "Forb", "Grass"), values= c("#fc8d62", "#66c2a5", "#8da0cb")) +
  scale_linetype_manual( values = c("dashed", "dashed", "dashed"), guide = "none")

joined_FEve$treatment <- factor(joined_FEve$treatment, levels =c("controlRain",  "springDry", "fallDry","consistentDry"))

p3_FEve <- ggplot(joined_FEve, aes(x = FEve, y = weight_g_m, color = treatment, linetype = treatment)) +
  geom_point() +
  theme_classic() +
  ylim(50,900)+
  labs(y = bquote('ANPP'~(g/m^2)), x = "Functional Evenness", color = "Treatment") +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE) +
  #stat_cor(aes(group=treatment,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x.npc = 0.5)+
  scale_color_manual(name = "Treatment", labels = c("Control", "Spring Dry", "Fall Dry","Consistent Dry" ), values= c("#0070b8", "#b2c7e4", "#fccaaf", "#c85b23"))+
  scale_linetype_manual( values = c("dashed", "dashed", "dashed", "dashed"), guide = "none")
p4_FEve <- ggplot(joined_FEve, aes(x = FEve, y = agg_BNPP, color = treatment, linetype = treatment)) +
  geom_point() +
  theme_classic() +
  ylim(50, 900)+
  theme(legend.position = "none") +
  labs(y = bquote('BNPP'~(g/m^2)), x = "Functional Evenness", color = "Treatment") +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE) +
  #stat_cor(aes(group=treatment,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x.npc = 0.5)+
  scale_color_manual(name = "Treatment", labels = c("Control", "Spring Dry", "Fall Dry","Consistent Dry" ), values= c("#0070b8", "#b2c7e4", "#fccaaf", "#c85b23"))+
  scale_linetype_manual( values = c("dashed", "dashed", "dashed", "dashed"), guide = "none")

legend_comp_FEve <- get_legend(p1_FEve)
legend_rain_FEve <- get_legend(p3_FEve)
p1_FEve <- p1_FEve + theme(legend.position = "none")
p3_FEve <- p3_FEve + theme(legend.position = "none")

ggarrange(p1_FEve, p2_FEve, legend_comp_FEve, p3_FEve, p4_FEve, legend_rain_FEve, ncol = 3, nrow = 2, 
          align = "v",labels = c("a)", "b)", "", "c)", "d)", ""), 
          widths=c(1, 1, 0.3))

#Fig2 revision Total biomass and FD metrics
fitoverall <- lm(total ~ FRic, joined_FD) 
summary(fitoverall) 
fitB <- lm(total ~ FRic, joined_FD%>%filter(subplot == "B")) 
summary(fitB) 
fitF <- lm(total ~ FRic, joined_FD%>%filter(subplot == "F"))
summary(fitF)
fitG <- lm(total ~ FRic, joined_FD%>%filter(subplot == "G"))
summary(fitG)
p_FRic <- ggplot(joined_FD, aes(x = FRic, y = total)) +
  geom_point(aes(color = subplot)) +
  theme_classic() +
  #ylim(50,900)+
  labs(y = bquote('Total Biomass'~(g/m^2)), x = "FRic", color = "Treatment") +
  geom_smooth(aes(color = subplot, linetype = subplot), method = lm, size = 1, se = FALSE, fullrange = FALSE) +
  #stat_cor(aes(group=subplot,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x.npc = 0.5)+
  scale_color_manual(name = "Treatment", labels = c("Mixed", "Forb", "Grass"), values= c("#fc8d62", "#66c2a5", "#8da0cb")) +
  scale_linetype_manual( values = c("solid", "dashed", "dashed"), guide = "none")+
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE, color = "black", linetype = "dashed") 

fitoverall <- lm(total ~ FEve, joined_FD) 
summary(fitoverall) 
fitB <- lm(total ~ FEve, joined_FD%>%filter(subplot == "B")) 
summary(fitB) 
fitF <- lm(total ~ FEve, joined_FD%>%filter(subplot == "F"))
summary(fitF)
fitG <- lm(total ~ FEve, joined_FD%>%filter(subplot == "G"))
summary(fitG)
p_FEve <- ggplot(joined_FD, aes(x = FEve, y = total)) +
  geom_point(aes(color = subplot)) +
  theme_classic() +
  #ylim(50,900)+
  theme(legend.position = "none") +
  labs(y = "", x = "FEve", color = "Treatment") +
  geom_smooth(aes(color = subplot, linetype = subplot), method = lm, size = 1, se = FALSE, fullrange = FALSE) +
  #stat_cor(aes(group=subplot,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x.npc = 0.5)+
  scale_color_manual(name = "Treatment", labels = c("Mixed", "Forb", "Grass"), values= c("#fc8d62", "#66c2a5", "#8da0cb")) +
  scale_linetype_manual( values = c("dashed", "dashed", "dashed"), guide = "none")+
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE, color = "black", linetype = "dashed") 

fitoverall <- lm(total ~ FDiv, joined_FD) 
summary(fitoverall) 
fitB <- lm(total ~ FDiv, joined_FD%>%filter(subplot == "B")) 
summary(fitB) 
fitF <- lm(total ~ FDiv, joined_FD%>%filter(subplot == "F"))
summary(fitF)
fitG <- lm(total ~ FDiv, joined_FD%>%filter(subplot == "G"))
summary(fitG)
p_FDiv <- ggplot(joined_FD, aes(x = FDiv, y = total)) +
  geom_point(aes(color = subplot)) +
  theme_classic() +
  #ylim(50,900)+
  theme(legend.position = "none") +
  labs(y = "", x = "FDiv", color = "Treatment") +
  geom_smooth(aes(color = subplot, linetype = subplot), method = lm, size = 1, se = FALSE, fullrange = FALSE) +
  #stat_cor(aes(group=subplot,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x.npc = 0.5)+
  scale_color_manual(name = "Treatment", labels = c("Mixed", "Forb", "Grass"), values= c("#fc8d62", "#66c2a5", "#8da0cb")) +
  scale_linetype_manual( values = c("dashed", "dashed", "dashed"), guide = "none")+
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE, color = "black", linetype = "dashed") 

fitoverall <- lm(total ~ RaoQ, joined_FD) 
summary(fitoverall) 
fitB <- lm(total ~ RaoQ, joined_FD%>%filter(subplot == "B")) 
summary(fitB) 
fitF <- lm(total ~ RaoQ, joined_FD%>%filter(subplot == "F"))
summary(fitF)
fitG <- lm(total ~ RaoQ, joined_FD%>%filter(subplot == "G"))
summary(fitG)
p_raoQ <- ggplot(joined_FD, aes(x = RaoQ, y = total)) +
  geom_point(aes(color = subplot)) +
  theme_classic() +
  #ylim(50,900)+
  theme(legend.position = "none") +
  labs(y = "", x = "Rao's Q", color = "Treatment") +
  geom_smooth(aes(color = subplot, linetype = subplot), method = lm, size = 1, se = FALSE, fullrange = FALSE) +
  #stat_cor(aes(group=subplot,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x.npc = 0.5)+
  scale_color_manual(name = "Treatment", labels = c("Mixed", "Forb", "Grass"), values= c("#fc8d62", "#66c2a5", "#8da0cb")) +
  scale_linetype_manual( values = c("solid", "dashed", "dashed"), guide = "none")+
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE, color = "black")

joined_FD$treatment <- factor(joined_FD$treatment, levels =c("controlRain",  "springDry", "fallDry","consistentDry"))

fitconstDry <- lm(total ~ FRic, joined_FD%>%filter(treatment == "consistentDry")) 
summary(fitconstDry) 
fitcontrol <- lm(total ~ FRic, joined_FD%>%filter(treatment == "controlRain"))
summary(fitcontrol)
fitfallDry <- lm(total ~ FRic, joined_FD%>%filter(treatment == "fallDry"))
summary(fitfallDry)
fitspringDry <- lm(total ~ FRic, joined_FD%>%filter(treatment == "springDry"))
summary(fitspringDry)
p_FRic_PPT <- ggplot(joined_FD, aes(x = FRic, y = total)) +
  geom_point(aes(color = treatment)) +
  theme_classic() +
  #ylim(50,900)+
  labs(y = bquote('Total Biomass'~(g/m^2)), x = "FRic", color = "Treatment") +
  geom_smooth(aes(color = treatment, linetype = treatment), method = lm, size = 1, se = FALSE, fullrange = FALSE) +
  #stat_cor(aes(group=subplot,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x.npc = 0.5)+
  scale_color_manual(name = "Treatment", labels = c("Control", "Spring Dry", "Fall Dry","Consistent Dry" ), values= c("#0070b8", "#b2c7e4", "#fccaaf", "#c85b23")) +
  scale_linetype_manual( values = c("dashed", "dashed", "dashed", "dashed"), guide = "none")+
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE, color = "black", linetype = "dashed") 

fitconstDry <- lm(total ~ FEve, joined_FD%>%filter(treatment == "consistentDry")) 
summary(fitconstDry) 
fitcontrol <- lm(total ~ FEve, joined_FD%>%filter(treatment == "controlRain"))
summary(fitcontrol)
fitfallDry <- lm(total ~ FEve, joined_FD%>%filter(treatment == "fallDry"))
summary(fitfallDry)
fitspringDry <- lm(total ~ FEve, joined_FD%>%filter(treatment == "springDry"))
summary(fitspringDry)
p_FEve_PPT <- ggplot(joined_FD, aes(x = FEve, y = total)) +
  geom_point(aes(color = treatment)) +
  theme_classic() +
  #ylim(50,900)+
  theme(legend.position = "none") +
  labs(y = "", x = "FEve", color = "Treatment") +
  geom_smooth(aes(color = treatment, linetype = treatment), method = lm, size = 1, se = FALSE, fullrange = FALSE) +
  #stat_cor(aes(group=subplot,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x.npc = 0.5)+
  scale_color_manual(name = "Treatment", labels = c("Control", "Spring Dry", "Fall Dry","Consistent Dry" ), values= c("#0070b8", "#b2c7e4", "#fccaaf", "#c85b23")) +
  scale_linetype_manual( values = c("dashed", "dashed", "dashed", "dashed"), guide = "none")+
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE, color = "black", linetype = "dashed") 

fitconstDry <- lm(total ~ FDiv, joined_FD%>%filter(treatment == "consistentDry")) 
summary(fitconstDry) 
fitcontrol <- lm(total ~ FDiv, joined_FD%>%filter(treatment == "controlRain"))
summary(fitcontrol)
fitfallDry <- lm(total ~ FDiv, joined_FD%>%filter(treatment == "fallDry"))
summary(fitfallDry)
fitspringDry <- lm(total ~ FDiv, joined_FD%>%filter(treatment == "springDry"))
summary(fitspringDry)
p_FDiv_PPT <- ggplot(joined_FD, aes(x = FDiv, y = total)) +
  geom_point(aes(color = treatment)) +
  theme_classic() +
  #ylim(50,900)+
  theme(legend.position = "none") +
  labs(y = "", x = "FDiv", color = "Treatment") +
  geom_smooth(aes(color = treatment, linetype = treatment), method = lm, size = 1, se = FALSE, fullrange = FALSE) +
  #stat_cor(aes(group=subplot,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x.npc = 0.5)+
  scale_color_manual(name = "Treatment", labels = c("Control", "Spring Dry", "Fall Dry","Consistent Dry" ), values= c("#0070b8", "#b2c7e4", "#fccaaf", "#c85b23")) +
  scale_linetype_manual( values = c("dashed", "dashed", "dashed", "dashed"), guide = "none")+
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE, color = "black", linetype = "dashed") 

fitconstDry <- lm(total ~ RaoQ, joined_FD%>%filter(treatment == "consistentDry")) 
summary(fitconstDry) 
fitcontrol <- lm(total ~ RaoQ, joined_FD%>%filter(treatment == "controlRain"))
summary(fitcontrol)
fitfallDry <- lm(total ~ RaoQ, joined_FD%>%filter(treatment == "fallDry"))
summary(fitfallDry)
fitspringDry <- lm(total ~ RaoQ, joined_FD%>%filter(treatment == "springDry"))
summary(fitspringDry)
p_RaoQ_PPT <- ggplot(joined_FD, aes(x = RaoQ, y = total)) +
  geom_point(aes(color = treatment)) +
  theme_classic() +
  #ylim(50,900)+
  theme(legend.position = "none") +
  labs(y = "", x = "RaoQ", color = "Treatment") +
  geom_smooth(aes(color = treatment, linetype = treatment), method = lm, size = 1, se = FALSE, fullrange = FALSE) +
  #stat_cor(aes(group=subplot,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x.npc = 0.5)+
  scale_color_manual(name = "Treatment", labels = c("Control", "Spring Dry", "Fall Dry","Consistent Dry" ), values= c("#0070b8", "#b2c7e4", "#fccaaf", "#c85b23")) +
  scale_linetype_manual( values = c("dashed", "solid", "dashed", "dashed"), guide = "none")+
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE, color = "black", linetype = "solid") 

legend_comp_FD <- get_legend(p_FRic)
legend_rain_FD <- get_legend(p_FRic_PPT)
p_FRic <- p_FRic + theme(legend.position = "none")
p_FRic_PPT <- p_FRic_PPT + theme(legend.position = "none")

panel1 <- ggarrange(p_FRic, p_FEve, p_FDiv, p_raoQ,  ncol = 4, nrow = 1, 
          align = "v",labels = c("a)", "b)", "c)", "d)"), common.legend = FALSE)

panel2 <- ggarrange(p_FRic_PPT, p_FEve_PPT, p_FDiv_PPT, p_RaoQ_PPT,  ncol = 4, nrow = 1, 
          align = "v",labels = c("e)", "f)", "g)", "h)"), common.legend = FALSE)

ggarrange(panel1, legend_comp_FD, panel2,legend_rain_FD, ncol = 2, nrow = 2, widths=c(1,0.1))

# CENSOL and TAECAP vs biomass
ggplot(joined_FD, aes(x = log(CENSOL), y = agg_BNPP)) +
  geom_point(aes(color = subplot)) +
  theme_classic() +
  ylim(50, 900)+
  #theme(legend.position = "none") +
  labs(y = bquote('BNPP'~(g/m^2)), x = "CENSOL", color = "Treatment") +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE, color = "black") +
  #stat_cor(aes(group=treatment,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x.npc = 0.5)+
  scale_color_manual(name = "Treatment", labels = c("Mixed", "Forb", "Grass"), values= c("#fc8d62", "#66c2a5", "#8da0cb")) 
  #scale_linetype_manual( values = c("dashed", "dashed", "dashed", "dashed"), guide = "none")
ggplot(joined_FD, aes(x = log(CENSOL), y = total)) +
  geom_point(aes(color = subplot)) +
  theme_classic() +
  ylim(50, 900)+
  #theme(legend.position = "none") +
  labs(y = bquote('Total Biomass'~(g/m^2)), x = "CENSOL", color = "Treatment") +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE, color = "black") +
  #stat_cor(aes(group=treatment,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x.npc = 0.5)+
  scale_color_manual(name = "Treatment", labels = c("Mixed", "Forb", "Grass"), values= c("#fc8d62", "#66c2a5", "#8da0cb")) 
  #scale_linetype_manual( values = c("dashed", "dashed", "dashed", "dashed"), guide = "none")
ggplot(joined_FD, aes(x = log(CENSOL), y = weight_g_m)) +
  geom_point(aes(color = subplot)) +
  theme_classic() +
  ylim(50, 900)+
  #theme(legend.position = "none") +
  labs(y = bquote('ANPP'~(g/m^2)), x = "CENSOL", color = "Treatment") +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE, color = "black") +
  #stat_cor(aes(group=treatment,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x.npc = 0.5)+
  scale_color_manual(name = "Treatment", labels = c("Mixed", "Forb", "Grass"), values= c("#fc8d62", "#66c2a5", "#8da0cb")) 
  #scale_linetype_manual( values = c("dashed", "dashed", "dashed", "dashed"), guide = "none")

ggplot(joined_FD, aes(x = log(TAECAP), y = agg_BNPP, color = treatment, linetype = treatment)) +
  geom_point() +
  theme_classic() +
  ylim(50, 900)+
  #theme(legend.position = "none") +
  labs(y = bquote('BNPP'~(g/m^2)), x = "TAECAP", color = "Treatment") +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE) +
  #stat_cor(aes(group=treatment,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x.npc = 0.5)+
  scale_color_manual(name = "Treatment", labels = c("Control", "Spring Dry", "Fall Dry","Consistent Dry" ), values= c("#0070b8", "#b2c7e4", "#fccaaf", "#c85b23"))+
  scale_linetype_manual( values = c("dashed", "dashed", "dashed", "dashed"), guide = "none")
ggplot(joined_FD, aes(x = log(TAECAP), y = total, color = treatment, linetype = treatment)) +
  geom_point() +
  theme_classic() +
  ylim(50, 900)+
  #theme(legend.position = "none") +
  labs(y = bquote('Total Biomass'~(g/m^2)), x = "TAECAP", color = "Treatment") +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE) +
  #stat_cor(aes(group=treatment,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x.npc = 0.5)+
  scale_color_manual(name = "Treatment", labels = c("Control", "Spring Dry", "Fall Dry","Consistent Dry" ), values= c("#0070b8", "#b2c7e4", "#fccaaf", "#c85b23"))+
  scale_linetype_manual( values = c("dashed", "dashed", "dashed", "dashed"), guide = "none")
ggplot(joined_FD, aes(x = log(TAECAP), y = weight_g_m, color = treatment, linetype = treatment)) +
  geom_point() +
  theme_classic() +
  ylim(50, 900)+
  #theme(legend.position = "none") +
  labs(y = bquote('ANPP'~(g/m^2)), x = "TAECAP", color = "Treatment") +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE) +
  #stat_cor(aes(group=treatment,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x.npc = 0.5)+
  scale_color_manual(name = "Treatment", labels = c("Control", "Spring Dry", "Fall Dry","Consistent Dry" ), values= c("#0070b8", "#b2c7e4", "#fccaaf", "#c85b23"))+
  scale_linetype_manual( values = c("dashed", "dashed", "dashed", "dashed"), guide = "none")

ggplot(joined_FD, aes(x = log(AVEBAR), y = agg_BNPP, color = treatment, linetype = treatment)) +
  geom_point() +
  theme_classic() +
  ylim(50, 900)+
  #theme(legend.position = "none") +
  labs(y = bquote('BNPP'~(g/m^2)), x = "AVEBAR", color = "Treatment") +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE) +
  #stat_cor(aes(group=treatment,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x.npc = 0.5)+
  scale_color_manual(name = "Treatment", labels = c("Control", "Spring Dry", "Fall Dry","Consistent Dry" ), values= c("#0070b8", "#b2c7e4", "#fccaaf", "#c85b23"))+
  scale_linetype_manual( values = c("dashed", "dashed", "dashed", "dashed"), guide = "none")
ggplot(joined_FD, aes(x = log(AVEBAR), y = total, color = treatment, linetype = treatment)) +
  geom_point() +
  theme_classic() +
  ylim(50, 900)+
  #theme(legend.position = "none") +
  labs(y = bquote('Total Biomass'~(g/m^2)), x = "AVEBAR", color = "Treatment") +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE) +
  #stat_cor(aes(group=treatment,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x.npc = 0.5)+
  scale_color_manual(name = "Treatment", labels = c("Control", "Spring Dry", "Fall Dry","Consistent Dry" ), values= c("#0070b8", "#b2c7e4", "#fccaaf", "#c85b23"))+
  scale_linetype_manual( values = c("dashed", "dashed", "dashed", "dashed"), guide = "none")
ggplot(joined_FD, aes(x = log(AVEBAR), y = weight_g_m, color = treatment, linetype = treatment)) +
  geom_point() +
  theme_classic() +
  ylim(50, 900)+
  #theme(legend.position = "none") +
  labs(y = bquote('ANPP'~(g/m^2)), x = "TAECAP", color = "Treatment") +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE) +
  #stat_cor(aes(group=treatment,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x.npc = 0.5)+
  scale_color_manual(name = "Treatment", labels = c("Control", "Spring Dry", "Fall Dry","Consistent Dry" ), values= c("#0070b8", "#b2c7e4", "#fccaaf", "#c85b23"))+
  scale_linetype_manual( values = c("dashed", "dashed", "dashed", "dashed"), guide = "none")


#ind trait RaoQ with ANPP and BNPP Table S6
#Regression BNPP ~ Rao's Q for each trait
FDSLA_all <- lm(weight_g_m ~ SLARaoQ, joined_rao)
summary(FDSLA_all) #significant
FDSLA_both<- lm(weight_g_m ~ SLARaoQ, joined_rao%>%filter(subplot == "B")) 
summary(FDSLA_both) #significant
FDSLA_forb <- lm(weight_g_m ~ SLARaoQ, joined_rao%>%filter(subplot == "F"))
summary(FDSLA_forb) #not significant
FDSLA_grass <- lm(weight_g_m ~ SLARaoQ, joined_rao%>%filter(subplot == "G"))
summary(FDSLA_grass ) #not significant

FDLDMC_all <- lm(weight_g_m ~ LDMCRaoQ, joined_rao)
summary(FDLDMC_all) #significant
FDLDMC_both<- lm(weight_g_m ~ LDMCRaoQ, joined_rao%>%filter(subplot == "B")) 
summary(FDLDMC_both) #significant
FDLDMC_forb <- lm(weight_g_m ~ LDMCRaoQ, joined_rao%>%filter(subplot == "F"))
summary(FDLDMC_forb) #not significant
FDLDMC_grass <- lm(weight_g_m ~ LDMCRaoQ, joined_rao%>%filter(subplot == "G"))
summary(FDLDMC_grass ) #not significant

FDHt_all <- lm(weight_g_m ~ HtRaoQ, joined_rao)
summary(FDHt_all) #not significant
FDHt_both<- lm(weight_g_m ~ HtRaoQ, joined_rao%>%filter(subplot == "B")) 
summary(FDHt_both) #not significant
FDHt_forb <- lm(weight_g_m ~ HtRaoQ, joined_rao%>%filter(subplot == "F"))
summary(FDHt_forb) #not significant
FDHt_grass <- lm(weight_g_m ~ HtRaoQ, joined_rao%>%filter(subplot == "G"))
summary(FDHt_grass ) #not significant

FDDens_all <- lm(agg_BNPP ~ DensRaoQ, joined_rao)
summary(FDDens_all) #not significant
FDDens_both<- lm(agg_BNPP ~ DensRaoQ, joined_rao%>%filter(subplot == "B")) 
summary(FDDens_both) #not significant
FDDens_forb <- lm(agg_BNPP ~ DensRaoQ, joined_rao%>%filter(subplot == "F"))
summary(FDDens_forb) #not significant
FDDens_grass <- lm(agg_BNPP ~ DensRaoQ, joined_rao%>%filter(subplot == "G"))
summary(FDDens_grass ) #not significant

FDSRLF_all <- lm(agg_BNPP ~ SRLFRaoQ, joined_rao)
summary(FDSRLF_all) #not significant
FDSRLF_both<- lm(agg_BNPP ~ SRLFRaoQ, joined_rao%>%filter(subplot == "B")) 
summary(FDSRLF_both) #not significant
FDSRLF_forb <- lm(agg_BNPP ~ SRLFRaoQ, joined_rao%>%filter(subplot == "F"))
summary(FDSRLF_forb) #not significant
FDSRLF_grass <- lm(agg_BNPP ~ SRLFRaoQ, joined_rao%>%filter(subplot == "G"))
summary(FDSRLF_grass) #not significant

FDSRLC_all <- lm(agg_BNPP ~ SRLCRaoQ, joined_rao)
summary(FDSRLC_all) #not significant
FDSRLC_both<- lm(agg_BNPP ~ SRLCRaoQ, joined_rao%>%filter(subplot == "B")) 
summary(FDSRLC_both) #not significant
FDSRLC_forb <- lm(agg_BNPP ~ SRLCRaoQ, joined_rao%>%filter(subplot == "F"))
summary(FDSRLC_forb) #not significant
FDSRLC_grass <- lm(agg_BNPP ~ SRLCRaoQ, joined_rao%>%filter(subplot == "G"))
summary(FDSRLC_grass ) #not significant

FDDiamC_all <- lm(agg_BNPP ~ DiamCRaoQ, joined_rao)
summary(FDDiamC_all) #not significant
FDDiamC_both<- lm(agg_BNPP ~ DiamCRaoQ, joined_rao%>%filter(subplot == "B")) 
summary(FDDiamC_both) #not significant
FDDiamC_forb <- lm(agg_BNPP ~ DiamCRaoQ, joined_rao%>%filter(subplot == "F"))
summary(FDDiamC_forb) #not significant
FDDiamC_grass <- lm(agg_BNPP ~ DiamCRaoQ, joined_rao%>%filter(subplot == "G"))
summary(FDDiamC_grass ) #not significant

FDPropF_all <- lm(agg_BNPP ~ PropFRaoQ, joined_rao)
summary(FDPropF_all) #not significant
FDPropF_both<- lm(agg_BNPP ~ PropFRaoQ, joined_rao%>%filter(subplot == "B")) 
summary(FDPropF_both) #not significant
FDPropF_forb <- lm(agg_BNPP ~ PropFRaoQ, joined_rao%>%filter(subplot == "F"))
summary(FDPropF_forb) #not significant
FDPropF_grass <- lm(agg_BNPP ~ PropFRaoQ, joined_rao%>%filter(subplot == "G"))
summary(FDPropF_grass ) #not significant

# Biomass ANPP and BNPP by overall Rao's Q
RaoANPP_all <- lm(weight_g_m ~ RaoQ, joined_rao)
summary(RaoANPP_all)
RaoANPP_both <- lm(weight_g_m ~ RaoQ, joined_rao%>%filter(subplot == "B"))
summary(RaoANPP_both)
RaoANPP_forb <- lm(weight_g_m ~ RaoQ, joined_rao%>%filter(subplot == "F"))
summary(RaoANPP_forb)
RaoANPP_grass <- lm(weight_g_m ~ RaoQ, joined_rao%>%filter(subplot == "G"))
summary(RaoANPP_grass)

RaoBNPP_all <- lm(agg_BNPP ~ RaoQ, joined_rao)
summary(RaoBNPP_all)
RaoBNPP_both <- lm(agg_BNPP ~ RaoQ, joined_rao%>%filter(subplot == "B"))
summary(RaoBNPP_both)
RaoBNPP_forb <- lm(agg_BNPP ~ RaoQ, joined_rao%>%filter(subplot == "F"))
summary(RaoBNPP_forb)
RaoBNPP_grass <- lm(agg_BNPP ~ RaoQ, joined_rao%>%filter(subplot == "G"))
summary(RaoBNPP_grass)

###Relationships bw RaoQ of three aboveground traits and ANPP
p1ar <- ggplot(data = joined_rao, aes(x = HtRaoQ, y = weight_g_m)) +
  geom_point(aes(color = as.factor(subplot))) +
  theme_classic() +
  xlab("Rao's Q Height") +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE, color = "black")+
  labs(y= NULL) +
  theme(legend.position="none") +
  scale_color_manual(name = "treatment", labels = c("Mixed", "Forb", "Grass"), values = c("#fc8d62", "#66c2a5", "#8da0cb") )+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))

p2ar <- ggplot(data = joined_rao, aes(x = SLARaoQ, y = weight_g_m)) +
  geom_point(aes(color = as.factor(subplot))) +
  theme_classic() +
  xlab("Rao's Q Specific Leaf Area") +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE, color = "black")+
  labs(y= NULL) +
  theme(legend.position= "none")+
  scale_color_manual(name = "treatment", labels = c("Mixed", "Forb", "Grass"), values = c("#fc8d62", "#66c2a5", "#8da0cb") )+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))

p3ar <- ggplot(data = joined_rao, aes(x = LDMCRaoQ, y = weight_g_m)) +
  geom_point(aes(color = as.factor(subplot))) +
  theme_classic() +
  xlab("Rao's Q Leaf Dry Matter Content") +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE, color = "black")+
  labs(y= NULL) +
  theme(legend.position = "none")+
  scale_color_manual(name = "treatment", labels = c("Mixed", "Forb", "Grass"), values = c("#fc8d62", "#66c2a5", "#8da0cb") )+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))


figurear <- ggarrange(p1ar, p2ar, p3ar, 
                     ncol =3, nrow =1, common.legend = TRUE, legend = "bottom",
                     align = "v",labels = c("a)", "b)", "c)"))
annotate_figure(figurear, 
                left = text_grob("ANPP (g/m2)", rot = 90))

###Relationships bw RaoQ of five root traits and BNPP
p1br <- ggplot(data = joined_rao, aes(x = DensRaoQ, y = agg_BNPP)) +
  geom_point(aes(color = as.factor(subplot))) +
  theme_classic() +
  xlab("Rao's Q Root Density") +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE, color = "black")+
  labs(y= NULL) +
  theme(legend.position="none") +
  scale_color_manual(name = "treatment", labels = c("Mixed", "Forb", "Grass"), values = c("#fc8d62", "#66c2a5", "#8da0cb") )+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))

p2br <- ggplot(data = joined_rao, aes(x = SRLFRaoQ, y = agg_BNPP)) +
  geom_point(aes(color = as.factor(subplot))) +
  theme_classic() +
  xlab("Rao's Q SRLF") +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE, color = "black")+
  labs(y= NULL) +
  theme(legend.position= "none")+
  scale_color_manual(name = "treatment", labels = c("Mixed", "Forb", "Grass"), values = c("#fc8d62", "#66c2a5", "#8da0cb") )+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))

p3br <- ggplot(data = joined_rao, aes(x = SRLCRaoQ, y = agg_BNPP)) +
  geom_point(aes(color = as.factor(subplot))) +
  theme_classic() +
  xlab("Rao's Q SRLC") +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE, color = "black")+
  labs(y= NULL) +
  theme(legend.position = "none")+
  scale_color_manual(name = "treatment", labels = c("Mixed", "Forb", "Grass"), values = c("#fc8d62", "#66c2a5", "#8da0cb") )+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))

p4br <- ggplot(data = joined_rao, aes(x = DiamCRaoQ, y = agg_BNPP)) +
  geom_point(aes(color = as.factor(subplot))) +
  theme_classic() +
  xlab("Rao's Q Diameter Coarse") +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE, color = "black")+
  labs(y= NULL) +
  theme(legend.position = "none")+
  scale_color_manual(name = "treatment", labels = c("Mixed", "Forb", "Grass"), values = c("#fc8d62", "#66c2a5", "#8da0cb") )+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))

p5br <- ggplot(data = joined_rao, aes(x = PropFRaoQ, y = agg_BNPP)) +
  geom_point(aes(color = as.factor(subplot))) +
  theme_classic() +
  xlab("Rao's Q Proportion of Fine") +
  labs(y= NULL) +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE, color = "black")+
  theme(legend.position = "none")+
  scale_color_manual(name = "treatment", labels = c("Mixed", "Forb", "Grass"), values = c("#fc8d62", "#66c2a5", "#8da0cb") )+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))

figurebr <- ggarrange(p1br, p2br, p3br, p4br, p5br,
                    ncol =3, nrow =2, common.legend = TRUE, legend = "bottom",
                    align = "v",labels = c("a)", "b)", "c)", "d)", "e)"))
annotate_figure(figurebr, 
                left = text_grob("BNPP (g/m2) depth 0-30 cm", rot = 90))
