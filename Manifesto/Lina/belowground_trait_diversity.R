veg <- read.csv("./Dropbox/ClimVar/DATA/Plant_composition_data/Cover/Cover_CleanedData/ClimVar_2015_species-cover_wide.csv", header = TRUE)
traits <- read.csv("./Dropbox/ClimVar/DATA/Plant_composition_data/Traits/Traits_GHscreening_Lina_Ashley/Combined_GHtraits.csv", header = TRUE)
BNPP <- read.csv("./Dropbox/ClimVar/DATA/Plant_composition_data/BNPP/BNPP_CleanedData/BNPP_MayHarvest_2015.csv", header = TRUE)

library(tidyverse)
library(FD)

###Set up data for FD package###
#Species names must be in same order in both files
veg_keys <- veg[,1:7]
vegdat <- veg[,8:64]
vegdat1 <- vegdat[, -c(6,10,12,15:16,20,22:25,28,31,32,34,36:39,41:44,46,53,57)] #remove all columns that we do not have trait data for and any species that has no cover in 2015
#group some species
vegdat1$BRODIA <- vegdat1$Bromus.diandrus + vegdat1$Bromus.sterilis
vegdat1$HORMUR <- vegdat1$Hordeum.marinum + vegdat1$Hordeum.murinum
vegdat1$TRISP <- vegdat1$Trifolium.dubium + vegdat1$Trifolium.glomeratum + vegdat1$Trifolium.sp. + vegdat1$Trifolium.wildenovii
vegdat2 <- vegdat1[ , -c(9,6,16,17,24,25,27,29)] #remove grouped columns 
colnames(vegdat2) <- c("ACHMIL","ANAARV", "AVEBAR", "AVEFAT", "BRADIS", "BROHOR", "BROMAD", "CENSOL","CERGLO", 
                       "CYNDAC", "CYNECH", "EROBOT", "EROMOS", "HYPGLA", "HYPRAD", 
                       "LACSER", "LOLMUL","RUMPUL", "TAECAP", "TRIHIR", "TRISUB", "VICSAT", "VULBRO", "VULMYU",
                       "BRODIA", "HORMUR", "TRISP") #rename columns
vegdat3 <- vegdat2[c("ACHMIL","ANAARV", "AVEBAR", "AVEFAT", "BRADIS", "BRODIA","BROHOR", "BROMAD", "CENSOL","CERGLO", 
                     "CYNDAC", "CYNECH", "EROBOT", "EROMOS",  "HORMUR","HYPGLA", "HYPRAD", 
                     "LACSER", "LOLMUL","RUMPUL", "TAECAP", "TRIHIR", "TRISP", "TRISUB", "VICSAT", "VULBRO", "VULMYU"
)] #reorder columns
dim(vegdat3) #check the dimension of veg data
comp <- as.matrix(vegdat3[1:nrow(vegdat3), 1:ncol(vegdat3)])

traits$ID <- as.character(traits$ID) #set ID column as character
traits$ID <- ifelse(traits$ID == "TRIREP", "TRISP", traits$ID) #use Trifolium repens traits for Trifolium sp. 
all_traits <- traits %>%
  filter(ID %in% c("ACHMIL","ANAARV", "AVEBAR", "AVEFAT", "BRADIS", "BRODIA","BROHOR", "BROMAD", "CENSOL","CERGLO", 
                   "CYNDAC", "CYNECH", "EROBOT", "EROMOS",  "HORMUR","HYPGLA", "HYPRAD", 
                   "LACSER", "LOLMUL","RUMPUL", "TAECAP", "TRIHIR", "TRISP", "TRISUB", "VICSAT", "VULBRO", "VULMYU")) %>%
  dplyr::select(ID, SLA, LDMC, Ht, Dens, DiamC, SRLC, SRLF, PropF)
dim(all_traits) #check the dimension of trait data
tr <- as.matrix(all_traits[1:nrow(all_traits), 2:ncol(all_traits)])
row.names(tr) <- all_traits$ID

#Now the data is ready to calculate trait diversity
results <- dbFD(tr, comp, corr="cailliez")
results <- as.data.frame(results) %>%
  tbl_df()
Div <- cbind(veg_keys, results) #combine with keys
#write.csv(Div, "Belowground_CWM_traits.csv")

#standardize Div dataset
library(vegan)
stand_Div_num <- decostand(Div[,10:20], "standardize")
stand_Div <- as.data.frame(append(Div[,1:9], stand_Div_num))

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

#Join the aggregated BNPP dataset with standardized Div
stand_BNPP <- as.data.frame(append(BNPP1[,1:8], decostand(BNPP1[,9], "standardize")))
stand_joined <- stand_Div %>%
  inner_join(stand_BNPP, by = c("plot", "subplot", "treatment", "shelterBlock"))

###Regression BNPP by Rao's Q
ggplot(data = joined, aes(x=RaoQ, y = agg_BNPP, col = subplot)) +
  geom_point() +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = TRUE)+
  theme_bw() +
  annotate("text", x = 1, y = 400, label = "R2 = 0.001", size = 4, color = "#F8766D") +
  annotate("text", x = 5, y = 450, label = "R2 = 0.122", size = 4, color = "#00BFC4") +
  annotate("text", x = 5, y = 250, label = "R2 = 0.011", size = 4, color = "#00ba38") +
  labs(y = "BNPP g/m2 depth 0-30 cm", col = "treatment")
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
lbl1 <- expression("y = 656.3 - 6142.1 x," ~ r^2 ~ "= 0.03")
lbl2 <- expression("y = 33.6 + 0.004 x," ~ r^2 ~ "= 0.02")
lbl3 <- expression("y = 129 + 0.03 x," ~ r^2 ~ "= 0.04")
lbl4 <- expression("y = 306.9 - 75.4 x," ~ r^2 ~ "= 0.03")
lbl5 <- expression("y = 10291 - 10181 x," ~ r^2 ~ "= 0.07")


p1 <- ggplot(data = joined, aes(x = CWM.Dens, y = agg_BNPP)) +
  geom_point(aes(color = as.factor(subplot))) +
  theme_bw() +
  xlab("CWM Root Density") +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE, color = "black")+
  labs(y= NULL) +
  theme(legend.position="none") +
  scale_color_discrete(name = "treatment", labels = c("Mixed", "Forb", "Grass")) +
  annotate("text", x = 0.065, y = 600, label = lbl1, size = 4, color = "black")

p2 <- ggplot(data = joined, aes(x = CWM.SRLF, y = agg_BNPP)) +
  geom_point(aes(color = as.factor(subplot))) +
  theme_bw() +
  xlab("CWM SRLF") +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE, color = "black")+
  labs(y= NULL) +
  theme(legend.position= "none")+
  scale_color_discrete(name = "treatment", labels = c("Mixed", "Forb", "Grass")) +
  annotate("text", x = 50000, y = 600, label = lbl2, size = 4, color = "black")

p3 <- ggplot(data = joined, aes(x = CWM.SRLC, y = agg_BNPP)) +
  geom_point(aes(color = as.factor(subplot))) +
  theme_bw() +
  xlab("CWM SRLC") +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE, color = "black")+
  labs(y= NULL) +
  theme(legend.position = "none")+
  scale_color_discrete(name = "treatment", labels = c("Mixed", "Forb", "Grass")) +
  annotate("text", x = 3000, y = 600, label = lbl3, size = 4, color = "black")

p4 <- ggplot(data = joined, aes(x = CWM.DiamC, y = agg_BNPP)) +
  geom_point(aes(color = as.factor(subplot))) +
  theme_bw() +
  xlab("CWM Diameter Coarse") +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE, color = "black")+
  labs(y= NULL) +
  theme(legend.position = "none")+
  scale_color_discrete(name = "treatment", labels = c("Mixed", "Forb", "Grass")) +
  annotate("text", x = 1, y = 600, label = lbl4, size = 4, color = "black")

p5 <- ggplot(data = joined, aes(x = CWM.PropF, y = agg_BNPP)) +
  geom_point(aes(color = as.factor(subplot))) +
  theme_bw() +
  xlab("CWM Proportion of Fine") +
  labs(y= NULL) +
  geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE, color = "black")+
  theme(legend.position = "none")+
  scale_color_discrete(name = "treatment", labels = c("Mixed", "Forb", "Grass")) +
  annotate("text", x = 0.982, y = 600, label = lbl5, size = 4, color = "black")

figure <- ggarrange(p1, p2, p3, p4, p5,
          ncol =3, nrow =2, common.legend = TRUE, legend = "bottom",
          align = "v",labels = c("a)", "b)", "c)", "d)", "e)"))
annotate_figure(figure, 
                left = text_grob("BNPP g/m2 depth 0-30 cm", rot = 90))

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

#Backward step-wise regression model - which trait contributes to BNPP?
stand_both <- stand_joined %>%
  filter(subplot == "B")
stand_forb <- stand_joined %>%
  filter(subplot == "F")
stand_grass <- stand_joined %>%
  filter(subplot == "G")

model0 <- lm(agg_BNPP ~ CWM.PropF, stand_joined)
summary(model0)
AIC(model0)

model1 <- lm(agg_BNPP ~ CWM.Dens + CWM.SRLC, stand_both)
summary(model1)
AIC(model1)

model2 <- lm(agg_BNPP ~ CWM.DiamC, stand_forb)
summary(model2)
AIC(model2)

model3 <- lm(agg_BNPP ~ CWM.SRLF + CWM.DiamC, stand_grass)
summary(model3)
AIC(model3)

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

#standardize ind RaoQ and BNPP
rao <- as.data.frame(cbind(Dens_joined[,14], SRLF_joined[,14], SRLC_joined[,14], DiamC_joined[,14], 
                                  PropF_joined[,14], Dens_joined[,20]))
colnames(rao) <- c("DensRao", "SRLFRao", "SRLCRao", "DiamCRao", "PropFRao", "BNPP")
stand_rao <- decostand(rao, "standardize")
stand_rao_plot <- as.data.frame(cbind(Dens_joined[,2:3],stand_rao))

#Plot BNPP vs Raos Q of each trait
lbl1 <- expression("y = 250.31 + 15.25 x," ~ r^2 ~ "= 0.001")
lbl2 <- expression("y = 291.31 - 89.81 x," ~ r^2 ~ "= 0.02")
lbl3 <- expression("y = 290.77 - 59.32 x," ~ r^2 ~ "= 0.03")
lbl4 <- expression("y = 274.82 - 30.08 x," ~ r^2 ~ "= 0.01")
lbl5 <- expression("y = 240.65 - 22.11 x," ~ r^2 ~ "= 0.004")

f1 <- ggplot(data = Dens_joined, aes(x = RaoQ, y = agg_BNPP)) +
        geom_point(aes(color = as.factor(subplot))) +
        theme_bw() +
        xlab("Rao's Q Root Density") +
        geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE, color = "black")+
        labs(y= NULL) +
        theme(legend.position="none")+
        scale_color_discrete(name = "treatment", labels = c("Mixed", "Forb", "Grass")) +
        annotate("text", x = 1, y = 600, label = lbl1, size = 4, color = "black")

f2 <- ggplot(data = SRLF_joined, aes(x = RaoQ, y = agg_BNPP)) +
        geom_point(aes(color = as.factor(subplot))) +
        theme_bw() +
        xlab("Rao's Q SRLF") +
        geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE, color = "black")+
        labs(y= NULL) +
        theme(legend.position="none") +
        annotate("text", x = 0.75, y = 600, label = lbl2 , size = 4, color = "black")

f3 <- ggplot(data = SRLC_joined, aes(x = RaoQ, y = agg_BNPP)) +
        geom_point(aes(color = as.factor(subplot))) +
        theme_bw() +
        xlab("Rao's Q SRLC") +
        geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE, color = "black")+
        labs(y= NULL) +
        theme(legend.position="none") +
        annotate("text", x = 0.9, y = 600, label = lbl3, size = 4, color = "black")

f4 <- ggplot(data = DiamC_joined, aes(x = RaoQ, y = agg_BNPP)) +
        geom_point(aes(color = as.factor(subplot))) +
        theme_bw() +
        xlab("Rao's Q Diameter Coarse") +
        geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE, color = "black")+
        labs(y= NULL) +
        theme(legend.position="none") +
        annotate("text", x = 1.5, y = 600, label = lbl4, size = 4, color = "black")

f5 <- ggplot(data = PropF_joined, aes(x = RaoQ, y = agg_BNPP)) +
        geom_point(aes(color = as.factor(subplot))) +
        theme_bw() +
        xlab("Rao's Q Proportion Fine") +
        geom_smooth(method = lm, size = 1, se = FALSE, fullrange = FALSE, color = "black")+
        labs(y= NULL) +
        theme(legend.position="none") +
        annotate("text", x = 0.85, y = 600, label = lbl5, size = 4, color = "black")

figure2 <- ggarrange(f1, f2, f3, f4, f5,
                    ncol =3, nrow =2, common.legend = TRUE, legend = "bottom",
                    labels = c("a)", "b)", "c)", "d)", "e)"))
annotate_figure(figure2, 
                left = text_grob("BNPP g/m2 depth 0-30 cm", rot = 90))

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
