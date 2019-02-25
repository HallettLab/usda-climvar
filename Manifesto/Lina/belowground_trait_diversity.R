veg <- read.csv("./Dropbox/ClimVar/DATA/Plant_composition_data/Cover/Cover_CleanedData/ClimVar_2015_species-cover_wide.csv", header = TRUE)
traits <- read.csv("./Dropbox/ClimVar/DATA/Plant_composition_data/Traits/Traits_GHscreening_Brad/BradTraits_Cleaned.csv", header = TRUE)
BNPP <- read.csv("./Dropbox/ClimVar/DATA/Plant_composition_data/BNPP/BNPP_CleanedData/BNPP_MayHarvest_2015.csv", header = TRUE)

library(tidyverse)
library(FD)

##Species names must be in same order in both files
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

##calculate trait diversity
results <- dbFD(tr, comp, corr="cailliez")
results <- as.data.frame(results) %>%
  tbl_df()
Div <- cbind(veg_keys, results) ##combine with keys

##ANOVA Rao's Q with subplot
library(nlme)
library(ggplot2)
fit1 <- lme(RaoQ ~ subplot, random=~1|shelterBlock/treatment, Div) 
anova(fit1)
summary(fit1) #B significantly lower than F, B significantly higher than G
ggplot(data = Div, aes(x = subplot, y = RaoQ, fill = subplot)) +
  geom_boxplot() +
  theme_classic()

##combine BNPP
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

##Regression BNPP by Rao's Q
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

##Regression BNPP by FDiv (functional divergence)
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

##Regression BNPP by FRic (functional richness)
ggplot(data = joined, aes(x=FRic, y = agg_BNPP, col = subplot)) +
  geom_point() +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = TRUE)+
  theme_classic()
