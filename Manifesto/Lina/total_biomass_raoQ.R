veg <- read.csv("./Dropbox/ClimVar/DATA/Plant_composition_data/Cover/Cover_CleanedData/ClimVar_2015_species-cover_wide.csv", header = TRUE)
traits <- read.csv("./Dropbox/ClimVar/DATA/Plant_composition_data/Traits/Traits_GHscreening_Brad/BradTraits_Cleaned.csv", header = TRUE)
BNPP <- read.csv("./Dropbox/ClimVar/DATA/Plant_composition_data/BNPP/BNPP_CleanedData/BNPP_MayHarvest_2015.csv", header = TRUE)
ANPP <- read.csv("./Dropbox/ClimVar/DATA/Plant_composition_data/ANPP/ANPP_CleanedData/ClimVar_ANPP-peak.csv", header = TRUE)

library(tidyverse)
library(FD)

###Set up data for FD package###
#Species names must be in same order in both files
veg_keys <- veg[,1:7]
vegdat <- veg[,8:64]
vegdat1 <- vegdat[, -c(2,6,10,12,14:16,22:25,28:32,34,37:44,46,53:54,57)] #remove all columns that we do not have trait data on
vegdat1$AVESP <- vegdat1$Avena.barbata + vegdat1$Avena.fatua #group some species
vegdat1$BRODIA <- vegdat1$Bromus.diandrus + vegdat1$Bromus.sterilis
vegdat1$EROBOT <- vegdat1$Erodium.botrys + vegdat1$Erodium.cicutarium + vegdat1$Erodium.moschatum
vegdat1$HORMUR <- vegdat1$Hordeum.marinum + vegdat1$Hordeum.murinum
vegdat1$TRISP <- vegdat1$Trifolium.dubium + vegdat1$Trifolium.glomeratum + vegdat1$Trifolium.sp. + vegdat1$Trifolium.subterraneum + vegdat1$Trifolium.wildenovii
vegdat1$VULMYU <- vegdat1$Vulpia.bromoides + vegdat1$Vulpia.myuros
vegdat2 <- vegdat1[ , -c(2,3,5,8,12:16,19,21,22,24:28)] #remove grouped columns and Lupinus bicolor because it does not occur in any community
colnames(vegdat2) <- c("ACHMIL", "BRADIS", "BROHOR", "BROMAD", "CENSOL", "CYNDAC", "CYNECH",
                       "LACSER", "LOLMUL", "TAECAP", "TRIHIR", "AVESP", "BRODIA", 
                       "EROBOT", "HORMUR", "TRISP", "VULMYU") #rename columns
vegdat3 <- vegdat2[c("ACHMIL", "AVESP", "BRADIS", "BRODIA", "BROHOR", "BROMAD", 
                     "CENSOL", "CYNDAC", "CYNECH", "EROBOT", "HORMUR", "LACSER",
                     "LOLMUL", "TAECAP", "TRIHIR", "TRISP", "VULMYU")] #reorder columns
dim(vegdat3) #check the dimension of veg data
comp <- as.matrix(vegdat3[1:nrow(vegdat3), 1:ncol(vegdat3)])

traits$ID <- as.character(traits$ID) #set ID column as character
traits$ID <- ifelse(traits$ID == "AVEFAT", "AVESP", traits$ID) #use Avena fatua traits for Avena sp.
traits$ID <- ifelse(traits$ID == "TRIREP", "TRISP", traits$ID) #use Trifolium repens traits for Trifolium sp. 
all_traits <- traits %>%
  filter(ID %in% c("ACHMIL", "AVESP", "BRADIS", "BRODIA", "BROHOR", "BROMAD", 
                   "CENSOL", "CYNDAC", "CYNECH", "EROBOT", "HORMUR", "LACSER",
                   "LOLMUL", "TAECAP", "TRIHIR","TRISP", "VULMYU")) %>%
  select(ID, SLA, LDMC, Ht, Dens, DiamC, SRLC, SRLF, PropF)
dim(all_traits) #check the dimension of trait data
tr <- as.matrix(all_traits[1:nrow(all_traits), 2:ncol(all_traits)])
row.names(tr) <- all_traits$ID

#Calculate ind RaoQ
SLA_results <- dbFD(tr[,1], comp, corr = "cailliez")
SLA_results <- as.data.frame(SLA_results) %>%
  tbl_df()

LDMC_results <- dbFD(tr[,2], comp, corr = "cailliez")
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

##Join ind RaoQ with biomass
rao_table <- as.data.frame(cbind(veg_keys[,2:3],SLA_results[,7], LDMC_results[,7],
                                 Ht_results[,7], Dens_results[,7], DiamC_results[,7],
                                 SRLC_results[,7], SRLF_results[,7], PropF_results[,7], overall_rao[,8]))
colnames(rao_table) <- c("plot", "subplot", "SLARaoQ", "LDMCRaoQ", "HtRaoQ", "DensRaoQ", 
                          "DiamCRaoQ", "SRLCRaoQ", "SRLFRaoQ", "PropFRaoQ", "RaoQ")
#Filter 2015 ANPP
ANPP1 <- ANPP %>%
  filter(year == 2015) %>%
  filter(subplot %in% c("B", "F", "G")) %>%
  select(plot, subplot, weight_g_m)

#Replace entry "20-Oct" to "10-20" in depth column in the BNPP dataset
BNPP$depth <- as.character(BNPP$depth) #set depth as character
BNPP <- BNPP %>%
  mutate(depth = replace(depth, depth == "20-Oct", "10-20")) %>%
  mutate(treatment_code = fct_recode(treatment, consDry = "consistentDry", control = "controlRain")) #make a column with shorter treatment names

#Calculate the aggregated BNPP in three levels of soil depth
BNPP1 <- BNPP %>%
  group_by(plot, subplot) %>% #group data
  summarise(agg_BNPP = sum(bmass_g_m2)) %>% #sum BNPP
  select(plot, subplot, agg_BNPP)

joined_rao <- rao_table %>%
  left_join(ANPP1, by = c("plot", "subplot")) %>%
  left_join(BNPP1, by = c("plot", "subplot")) %>%
  mutate(total = weight_g_m + agg_BNPP)

#standardize ind RaoQ and Biomass
stand_joined_rao <- decostand(joined_rao[,3:14], "standardize") 
stand_rao <- as.data.frame(cbind(joined_rao[,1:2], stand_joined_rao))

#Stepwise regression of BNPP ~ Ind RaoQ
rao_both <- stand_rao %>%
  filter(subplot == "B")
rao_forb <- stand_rao %>%
  filter(subplot == "F")
rao_grass <- stand_rao %>%
  filter(subplot == "G")

library(MASS)
model0 <- lm(total ~ SLARaoQ + LDMCRaoQ + HtRaoQ + DensRaoQ + DiamCRaoQ + SRLCRaoQ + SRLFRaoQ + PropFRaoQ, rao_both)
step_both <- stepAIC(model0, direction = "backward", trace = FALSE)
step_both$anova
model1 <- lm(total ~ SLARaoQ + LDMCRaoQ + HtRaoQ + DensRaoQ + DiamCRaoQ + SRLCRaoQ + SRLFRaoQ + PropFRaoQ, rao_both)
summary(model1)

model2 <- lm(total ~ SLARaoQ + LDMCRaoQ + HtRaoQ + DensRaoQ + DiamCRaoQ + SRLCRaoQ + SRLFRaoQ + PropFRaoQ, rao_forb)
step_forb <- stepAIC(model2, direction = "backward", trace = FALSE)
step_forb$anova
model3 <- lm(total ~ SLARaoQ + SRLFRaoQ, rao_forb)
summary(model3)

model4 <- lm(total ~ SLARaoQ + LDMCRaoQ + HtRaoQ + DensRaoQ + DiamCRaoQ + SRLCRaoQ + SRLFRaoQ + PropFRaoQ, rao_grass)
step_grass <- stepAIC(model4, direction = "backward", trace = FALSE)
step_grass$anova
model5 <- lm(total ~ SLARaoQ + LDMCRaoQ + DensRaoQ + SRLCRaoQ, rao_grass)
summary(model5)

model6 <- lm(total ~ SLARaoQ + LDMCRaoQ + HtRaoQ + DensRaoQ + DiamCRaoQ + SRLCRaoQ + SRLFRaoQ + PropFRaoQ, stand_rao)
step_all <- stepAIC(model6, direction = "backward", trace = FALSE)
step_all$anova
model7 <- lm(total ~ SLARaoQ, stand_rao)
summary(model7)

#plot ANPP BNPP by RaoQ 
library(ggplot2)
library(ggpubr)
p1 <- ggplot(joined_rao, aes(x = RaoQ, y = weight_g_m, color = subplot)) +
  geom_point() +
  theme_bw() +
  labs(y = "ANPP g/m2", x = "Rao's Q", color = "treatment") +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE) +
  annotate("text", x = 3.5, y = 750, label = "R2 = 0.50", size = 4, color = "#f8766d") +
  annotate("text", x = 1.5, y = 530, label = "R2 = 0.09", size = 4, color = "#619bff") +
  annotate("text", x = 7, y = 300, label = "R2 = 0.03", size = 4, color = "#00ba38") +
  scale_color_discrete(name = "treatment", labels = c("Mixed", "Forb", "Grass")) 

p2 <- ggplot(joined_rao, aes(x = RaoQ, y = agg_BNPP, color = subplot)) +
  geom_point() +
  theme_bw() +
  labs(y = "BNPP g/m2 depth 0-30 cm", x = "Rao's Q", color = "treatment") +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE) +
  annotate("text", x = 2, y = 400, label = "R2 = <0.001", size = 4, color = "#f8766d") +
  annotate("text", x = 3.8, y = 300, label = "R2 = 0.15", size = 4, color = "#619bff") +
  annotate("text", x = 7, y = 210, label = "R2 = 0.002", size = 4, color = "#00ba38") +
  scale_color_discrete(name = "treatment", labels = c("Mixed", "Forb", "Grass")) 

ggarrange(p1, p2, ncol = 2, nrow = 1, 
          common.legend = TRUE, legend = "right",
          align = "v",labels = c("a)", "b)"))

#community RaoQ with ANPP and BNPP
fitB <- lm(weight_g_m ~ RaoQ, joined_rao%>%filter(subplot == "B")) 
summary(fitB) 
fitF <- lm(weight_g_m ~ RaoQ, joined_rao%>%filter(subplot == "F"))
summary(fitF)
fitG <- lm(weight_g_m ~ RaoQ, joined_rao%>%filter(subplot == "G"))
summary(fitG)

fitB <- lm(agg_BNPP ~ RaoQ, joined_rao%>%filter(subplot == "B")) 
summary(fitB) 
fitF <- lm(agg_BNPP ~ RaoQ, joined_rao%>%filter(subplot == "F"))
summary(fitF)
fitG <- lm(agg_BNPP ~ RaoQ, joined_rao%>%filter(subplot == "G"))
summary(fitG)

#ind trait RaoQ with ANPP and BNPP
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
summary(FDLDMC_all) #not significant
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
summary(FDDens_grass ) #significant

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


