veg <- read.csv("./Dropbox/ClimVar/DATA/Plant_composition_data/Cover/Cover_CleanedData/ClimVar_2015_species-cover_wide.csv", header = TRUE)
traits <- read.csv("./Dropbox/ClimVar/DATA/Plant_composition_data/Traits/Traits_GHscreening_Brad/BradTraits_Cleaned.csv", header = TRUE)
BNPP <- read.csv("./Dropbox/ClimVar/DATA/Plant_composition_data/BNPP/BNPP_CleanedData/BNPP_MayHarvest_2015.csv", header = TRUE)

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
below_traits <- traits %>%
  filter(ID %in% c("ACHMIL", "AVESP", "BRADIS", "BRODIA", "BROHOR", "BROMAD", 
                   "CENSOL", "CYNDAC", "CYNECH", "EROBOT", "HORMUR", "LACSER",
                   "LOLMUL", "TAECAP", "TRIHIR","TRISP", "VULMYU")) %>%
  select(ID, Dens, DiamC, SRLC, SRLF, PropF)
dim(below_traits) #check the dimension of trait data
tr <- as.matrix(below_traits[1:nrow(below_traits), 2:ncol(below_traits)])
row.names(tr) <- below_traits$ID

#Now the data is ready to calculate trait diversity
results <- dbFD(tr, comp, corr="cailliez")
results <- as.data.frame(results) %>%
  tbl_df()
Div <- cbind(veg_keys, results) #combine with keys

###ANOVA Rao's Q ~ subplot
library(nlme)
library(ggplot2)
fit1 <- lme(RaoQ ~ subplot, random=~1|shelterBlock/treatment, Div) 
anova(fit1)
summary(fit1) #B significantly lower than F, B significantly higher than G
ggplot(data = Div, aes(x = subplot, y = RaoQ, fill = subplot)) +
  geom_boxplot() +
  theme_classic()

###Combine BNPP with trait diversity
#Replace entry "20-Oct" to "10-20" in depth column in the BNPP dataset
BNPP$depth <- as.character(BNPP$depth) #set depth as character
BNPP <- BNPP %>%
  mutate(depth = replace(depth, depth == "20-Oct", "10-20")) %>%
  mutate(treatment_code = fct_recode(treatment, consDry = "consistentDry", control = "controlRain")) #make a column with shorter treatment names

#Calculate the aggregated BNPP in three levels of soil depth
BNPP1 <- BNPP %>%
  group_by(plot, subplot, treatment, treatment_code, shelterBlock, shelter, fall, spring) %>% #group data
  summarise(agg_BNPP = sum(bmass_g_m2))#sum BNPP 

#Join the aggregated BNPP dataset with Div
joined <- Div %>%
  inner_join(BNPP1, by = c( "plot", "subplot", "treatment", "shelterBlock")) 

###Regression BNPP by Rao's Q
ggplot(data = joined, aes(x=RaoQ, y = agg_BNPP, col = subplot)) +
  geom_point() +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = TRUE)+
  theme_classic()
#subset by subplot
both <- joined %>%
  filter(subplot == "B")
forb <- joined %>%
  filter(subplot == "F")
grass <- joined %>%
  filter(subplot == "G")
#linear regression by subplot
fit2 <- lm(agg_BNPP ~ RaoQ, both) 
summary(fit2) 
fit3 <- lm(agg_BNPP ~ RaoQ, forb)
summary(fit3)
fit4 <- lm(agg_BNPP ~ RaoQ, grass)
summary(fit4)

###Regression BNPP by FDiv (functional divergence)
ggplot(data = joined, aes(x=FDiv, y = agg_BNPP, col = subplot)) +
  geom_point() +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = TRUE)+
  theme_classic()
#linear regression by subplot
fit5 <- lm(agg_BNPP ~ FDiv, both) 
summary(fit5) 
fit6 <- lm(agg_BNPP ~ FDiv, forb)
summary(fit6)
fit7 <- lm(agg_BNPP ~ FDiv, grass)
summary(fit7)

###Regression BNPP by FRic (functional richness)
ggplot(data = joined, aes(x=FRic, y = agg_BNPP, col = subplot)) +
  geom_point() +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = TRUE)+
  theme_classic()
#linear regression by subplot
fit8 <- lm(agg_BNPP ~ FRic, both) 
summary(fit5) 
fit9 <- lm(agg_BNPP ~ FRic, forb)
summary(fit6)
fit10 <- lm(agg_BNPP ~ FRic, grass)
summary(fit7)

###Relationships bw CWM of five root traits and BNPP
library(ggpubr)
p1 <- ggplot(data = joined, aes(x = CWM.Dens, y = agg_BNPP, col = subplot)) +
  geom_point() +
  theme_classic() +
  xlab("CWM Root Density") +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE)+
  labs(y= NULL) +
  theme(legend.position="none")

p2 <- ggplot(data = joined, aes(x = CWM.SRLF, y = agg_BNPP, col = subplot)) +
  geom_point() +
  theme_classic() +
  xlab("CWM SRLF") +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE)+
  labs(y= NULL) +
  theme(legend.position= "none")

p3 <- ggplot(data = joined, aes(x = CWM.SRLC, y = agg_BNPP, col = subplot)) +
  geom_point() +
  theme_classic() +
  xlab("CWM SRLC") +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE)+
  labs(y= NULL) +
  theme(legend.position = "none")

p4 <- ggplot(data = joined, aes(x = CWM.DiamC, y = agg_BNPP, col = subplot)) +
  geom_point() +
  theme_classic() +
  xlab("CWM Diameter Coarse") +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE)+
  labs(y= NULL) +
  theme(legend.position = "none")

p5 <- ggplot(data = joined, aes(x = CWM.PropF, y = agg_BNPP, col = subplot)) +
  geom_point() +
  theme_classic() +
  xlab("CWM Proportion of Fine") +
  labs(y= NULL) +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE)+
  theme(legend.position = "none")

figure <- ggarrange(p1, p2, p3, p4, p5,
          ncol =3, nrow =2, common.legend = TRUE,
          align = "v",labels = c("A", "B", "C", "D", "E"))
annotate_figure(figure, 
                left = text_grob("Belowground biomass production (soil 0-30 cm)", rot = 90))

#Regression BNPP ~ CWM.Dens
Dens_all <- lm(agg_BNPP ~ CWM.Dens, joined)
summary(Dens_all) #not significant
Dens_both<- lm(agg_BNPP ~ CWM.Dens, both) 
summary(Dens_both) #significant
Dens_forb <- lm(agg_BNPP ~ CWM.Dens, forb)
summary(Dens_forb) #not significant
Dens_grass <- lm(agg_BNPP ~ CWM.Dens, grass)
summary(Dens_grass) #not significant

#Regression BNPP ~ CWM.SRLF
SRLF_all <- lm(agg_BNPP ~ CWM.SRLF, joined)
summary(SRLF_all) #not significant
SRLF_both<- lm(agg_BNPP ~ CWM.SRLF, both) 
summary(SRLF_both) #not significant
SRLF_forb <- lm(agg_BNPP ~ CWM.SRLF, forb)
summary(SRLF_forb) #not significant
SRLF_grass <- lm(agg_BNPP ~ CWM.SRLF, grass)
summary(SRLF_grass) #not significant

#Regression BNPP ~ CWM.SRLC
SRLC_all <- lm(agg_BNPP ~ CWM.SRLC, joined)
summary(SRLC_all) #not significant
SRLC_both<- lm(agg_BNPP ~ CWM.SRLC, both) 
summary(SRLC_both) #not significant
SRLC_forb <- lm(agg_BNPP ~ CWM.SRLC, forb)
summary(SRLC_forb) #not significant
SRLC_grass <- lm(agg_BNPP ~ CWM.SRLC, grass)
summary(SRLC_grass) #not significant

#Regression BNPP ~ CWM.DiamC
DiamC_all <- lm(agg_BNPP ~ CWM.DiamC, joined)
summary(DiamC_all) #not significant
DiamC_both<- lm(agg_BNPP ~ CWM.DiamC, both) 
summary(DiamC_both) #not significant
DiamC_forb <- lm(agg_BNPP ~ CWM.DiamC, forb)
summary(DiamC_forb) #not significant
DiamC_grass <- lm(agg_BNPP ~ CWM.DiamC, grass)
summary(DiamC_grass) #significant

#Regression BNPP ~ CWM.PropF
PropF_all <- lm(agg_BNPP ~ CWM.PropF, joined)
summary(PropF_all) #not significant
PropF_both<- lm(agg_BNPP ~ CWM.PropF, both) 
summary(PropF_both) #not significant
PropF_forb <- lm(agg_BNPP ~ CWM.PropF, forb)
summary(PropF_forb) #not significant
PropF_grass <- lm(agg_BNPP ~ CWM.PropF, grass)
summary(PropF_grass) #not significant

#

###Calculate Rao's Q of each trait
Dens_results <- dbFD(tr[,1], comp, corr="cailliez")
Dens_results <- as.data.frame(Dens_results) %>%
  tbl_df()
Dens_Div <- cbind(veg_keys, Dens_results) #combine with keys
Dens_joined <- Dens_Div %>%
  inner_join(BNPP1, by = c( "plot", "subplot", "treatment", "shelterBlock")) 

SRLF_results <- dbFD(tr[,4], comp, corr="cailliez")
SRLF_results <- as.data.frame(SRLF_results) %>%
  tbl_df()
SRLF_Div <- cbind(veg_keys, SRLF_results) #combine with keys
SRLF_joined <- SRLF_Div %>%
  inner_join(BNPP1, by = c( "plot", "subplot", "treatment", "shelterBlock")) 

SRLC_results <- dbFD(tr[,3], comp, corr="cailliez")
SRLC_results <- as.data.frame(SRLC_results) %>%
  tbl_df()
SRLC_Div <- cbind(veg_keys, SRLC_results) #combine with keys
SRLC_joined <- SRLC_Div %>%
  inner_join(BNPP1, by = c( "plot", "subplot", "treatment", "shelterBlock"))

DiamC_results <- dbFD(tr[,2], comp, corr="cailliez")
DiamC_results <- as.data.frame(DiamC_results) %>%
  tbl_df()
DiamC_Div <- cbind(veg_keys, DiamC_results) #combine with keys
DiamC_joined <- DiamC_Div %>%
  inner_join(BNPP1, by = c( "plot", "subplot", "treatment", "shelterBlock"))

PropF_results <- dbFD(tr[,5], comp, corr="cailliez")
PropF_results <- as.data.frame(PropF_results) %>%
  tbl_df()
PropF_Div <- cbind(veg_keys, PropF_results) #combine with keys
PropF_joined <- PropF_Div %>%
  inner_join(BNPP1, by = c( "plot", "subplot", "treatment", "shelterBlock"))

#Plot BNPP vs Raos Q of each trait
f1 <- ggplot(data = Dens_joined, aes(x = RaoQ, y = agg_BNPP, col = subplot)) +
  geom_point() +
  theme_classic() +
  xlab("FD Root Density") +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE)+
  labs(y= NULL) +
  theme(legend.position="none")

f2 <- ggplot(data = SRLF_joined, aes(x = RaoQ, y = agg_BNPP, col = subplot)) +
  geom_point() +
  theme_classic() +
  xlab("FD SRLF") +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE)+
  labs(y= NULL) +
  theme(legend.position="none")

f3 <- ggplot(data = SRLC_joined, aes(x = RaoQ, y = agg_BNPP, col = subplot)) +
  geom_point() +
  theme_classic() +
  xlab("FD SRLC") +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE)+
  labs(y= NULL) +
  theme(legend.position="none")

f4 <- ggplot(data = DiamC_joined, aes(x = RaoQ, y = agg_BNPP, col = subplot)) +
  geom_point() +
  theme_classic() +
  xlab("FD Diameter Coarse") +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE)+
  labs(y= NULL) +
  theme(legend.position="none")

f5 <- ggplot(data = PropF_joined, aes(x = RaoQ, y = agg_BNPP, col = subplot)) +
  geom_point() +
  theme_classic() +
  xlab("FD Proportion Fine") +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE)+
  labs(y= NULL) +
  theme(legend.position="none")

figure2 <- ggarrange(f1, f2, f3, f4, f5,
                    ncol =3, nrow =2, common.legend = TRUE,
                    labels = c("A", "B", "C", "D", "E"))
annotate_figure(figure2, 
                left = text_grob("Belowground biomass production (soil 0-30 cm)", rot = 90))

#subset by subplot
Dens_both <- Dens_joined %>%
  filter(subplot == "B")
Dens_forb <- Dens_joined %>%
  filter(subplot == "F")
Dens_grass <- Dens_joined %>%
  filter(subplot == "G")
SRLF_both <- SRLF_joined %>%
  filter(subplot == "B")
SRLF_forb <- SRLF_joined %>%
  filter(subplot == "F")
SRLF_grass <- SRLF_joined %>%
  filter(subplot == "G")
SRLC_both <- SRLC_joined %>%
  filter(subplot == "B")
SRLC_forb <- SRLC_joined %>%
  filter(subplot == "F")
SRLC_grass <- SRLC_joined %>%
  filter(subplot == "G")
DiamC_both <- DiamC_joined %>%
  filter(subplot == "B")
DiamC_forb <- DiamC_joined %>%
  filter(subplot == "F")
DiamC_grass <- DiamC_joined %>%
  filter(subplot == "G")
PropF_both <- PropF_joined %>%
  filter(subplot == "B")
PropF_forb <- PropF_joined %>%
  filter(subplot == "F")
PropF_grass <- PropF_joined %>%
  filter(subplot == "G")

#Regression BNPP ~ Rao's Q for each trait
FDDens_all <- lm(agg_BNPP ~ RaoQ, Dens_joined)
summary(FDDens_all) #not significant
FDDens_both<- lm(agg_BNPP ~ RaoQ, Dens_both) 
summary(FDDens_both) #not significant
FDDens_forb <- lm(agg_BNPP ~ RaoQ, Dens_forb)
summary(FDDens_forb) #not significant
FDDens_grass <- lm(agg_BNPP ~ RaoQ, Dens_grass)
summary(FDDens_grass ) #significant

FDSRLF_all <- lm(agg_BNPP ~ RaoQ, SRLF_joined)
summary(FDSRLF_all) #not significant
FDSRLF_both<- lm(agg_BNPP ~ RaoQ, SRLF_both) 
summary(FDSRLF_both) #not significant
FDSRLF_forb <- lm(agg_BNPP ~ RaoQ, SRLF_forb)
summary(FDSRLF_forb) #not significant
FDSRLF_grass <- lm(agg_BNPP ~ RaoQ, SRLF_grass)
summary(FDSRLF_grass) #not significant

FDSRLC_all <- lm(agg_BNPP ~ RaoQ, SRLC_joined)
summary(FDSRLC_all) #not significant
FDSRLC_both<- lm(agg_BNPP ~ RaoQ, SRLC_both) 
summary(FDSRLC_both) #not significant
FDSRLC_forb <- lm(agg_BNPP ~ RaoQ, SRLC_forb)
summary(FDSRLC_forb) #not significant
FDSRLC_grass <- lm(agg_BNPP ~ RaoQ, SRLC_grass)
summary(FDSRLC_grass ) #not significant

FDDiamC_all <- lm(agg_BNPP ~ RaoQ, DiamC_joined)
summary(FDDiamC_all) #not significant
FDDiamC_both<- lm(agg_BNPP ~ RaoQ, DiamC_both) 
summary(FDDiamC_both) #not significant
FDDiamC_forb <- lm(agg_BNPP ~ RaoQ, DiamC_forb)
summary(FDDiamC_forb) #not significant
FDDiamC_grass <- lm(agg_BNPP ~ RaoQ, DiamC_grass)
summary(FDDiamC_grass ) #not significant

FDPropF_all <- lm(agg_BNPP ~ RaoQ, PropF_joined)
summary(FDPropF_all) #not significant
FDPropF_both<- lm(agg_BNPP ~ RaoQ, PropF_both) 
summary(FDPropF_both) #not significant
FDPropF_forb <- lm(agg_BNPP ~ RaoQ, PropF_forb)
summary(FDPropF_forb) #not significant
FDPropF_grass <- lm(agg_BNPP ~ RaoQ, PropF_grass)
summary(FDPropF_grass ) #not significant



