library(vegan)
library(tidyverse)
library(FD)
library(gridExtra)
library(ggpubr)
library(MASS)
library(ggplot2)
library(nlme)
library(MuMln)
library(lsmeans)

#function for standard error
calcSE<-function(x){
  x <- x[!is.na(x)]
  sd(x)/sqrt(length(x))}


## PCA for trait groups by plant species; visualized by orgin and functional group ##
## Followed by CMW of PCA scores, regressed against ANPP and BNPP
## NOTE: Uses all_trait, abovetr15_fd2, cover15_fd2 from traitxANPP.R and BNPP1 from Lina's BNPP_exploratory_analysis.R


# Select out correlated traits(MD, actual_area, Total) and those I don't have as much faith in (RMF, RGR)
#all_trait2 <- all_trait %>%
#  dplyr::select( -Seed.mass.grams, -C.N.Ratio)
#all_trait2 <- all_trait2[-c(42:77), ]

# Remove legumes from analysis if desired
#trait.dat2 <- subset(trait.dat2, GF != "L")
#trait.dat2$GF <- as.character(trait.dat2$GF)
#trait.dat2$GF <- as.factor(trait.dat2$GF)

# # First a PCA for ALL TRAITS (ABOVE AND BELOW)
# # matrix for PCA
# traits <- as.matrix(abovetr15_fd3[,c(6:ncol(abovetr15_fd3))])
# #row.names(traits) <- all_trait2$ID
# all_trait_redo <- trait.dat1[-c(2,3,5,9,10,14,17,20,21,22,23,24,27,28,29,30,31,34,37,39,40,
#                                  41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,58), ]
# all_trait_redo <- all_trait_redo %>% dplyr::select( -RGR, -actual_area,-MD,-Total,-Area)
# traits <- as.matrix(all_trait_redo[,c(6:ncol(all_trait_redo))])

# Lina updated trait data 8/4/2022 
# Use all_traits and traits in Lina's "total_biomass_raoQ.R" 
traits_only <- all_traits %>% select(-ID)
all_trait_redo <- left_join(all_traits, traits) %>%
  select( -RGR, -actual_area,-MD,-Total,-Area, -RMF)

# run PCA
myrda <- rda(na.omit(traits_only), scale = T)

# extract values
siteout <- as.data.frame(scores(myrda, choices=c(1,2), display=c("sites")))
siteout$ID<-all_trait_redo$ID
siteout$name <- siteout$ID

enviroout<-as.data.frame(scores(myrda, choices=c(1,2), display=c("species")))
enviroout$type<-"traits"
enviroout$name<-rownames(enviroout)

# merge PC axes with trait data
tog <- left_join(all_trait_redo, siteout) %>%
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
# 
# siteout<- siteout %>% dplyr::select(-ID, -name) 
# siteout<-bind_cols(siteout,cover_names)
# row.names(siteout)<-siteout$Taxon 
# siteout <- siteout %>% dplyr::select(-Taxon)
# #siteout <- siteout[-c(14), ] #drop LUPBIC from PCA scores bc it has 0 % cover in any plot 
# 
# siteout_fd<-dbFD (siteout, comp, w.abun = T, stand.x = F,
#                  calc.FRic = TRUE, m = "max", stand.FRic = FALSE,
#                  scale.RaoQ = FALSE, calc.FGR = FALSE, clust.type = "ward",
#                  km.inf.gr = 2, km.sup.gr = nrow(x) - 1, km.iter = 100, calc.CWM = TRUE,
#                  calc.FDiv = TRUE, print.pco = FALSE, messages = TRUE)
# siteout_fd<-as.data.frame(siteout_fd)
# 
# traits_2015_2<- arrange(traits_2015, treatment, shelterBlock)
# siteout_fd$ANPPgm<-traits_2015_2$ANPPgm
# siteout_fd$shelterBlock<-traits_2015_2$shelterBlock
# siteout_fd$subplot<-traits_2015_2$subplot
# siteout_fd$treatment<-traits_2015_2$treatment
# 
# siteout_fd <- siteout_fd %>%
#   left_join(ANPP1, by = c("plot", "subplot", "treatment")) 
# 
# f_pca1_ANPP <- ggplot(data=siteout_fd, aes(x=CWM.PC1, y=weight_g_m, group=subplot, color=subplot))+
#   #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
#   geom_point()+
#   theme_bw()+
#   scale_color_manual(values=c("tomato", "green3", "dodgerblue"), guide = guide_legend(title = "Treatment"), #change legend title
#                     labels=c("Mixed", "Forb", "Grass"))+ #change labels in the legend)+
#   xlab("Community Weighted Means of PC1 Scores")+
#   ylab("ANPP g/m2")+
#   #xlim(40,100)+
#   #ylim(0,100)
#   geom_smooth(method="lm", formula= y ~ x, se=FALSE, color="black", aes(group=1))
# 
# ggplot(data=siteout_fd, aes(x=CWM.PC1, y=weight_g_m, group=treatment, color=treatment))+
#   #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=treatment)) +
#   geom_point()+
#   theme_bw()+
# #  scale_color_manual(values=c("tomato", "green3", "dodgerblue"), guide = guide_legend(title = "Treatment"), #change legend title
# #                     labels=c("Mixed", "Forb", "Grass"))+ #change labels in the legend)+
#   xlab("Community Weighted Means of PC1 Scores")+
#   ylab("ANPP g/m2")+
#   #xlim(40,100)+
#   #ylim(0,100)
#   geom_smooth(method="lm", formula= y ~ x, se=FALSE,  aes(group=treatment))
# 
# f_pca2_ANPP <- ggplot(data=siteout_fd, aes(x=CWM.PC2, y=weight_g_m, group=subplot, color=subplot))+
#   #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
#   geom_point()+
#   theme_bw()+
#   xlab("Community Weighted Means of PC2 Scores")+
#   ylab("ANPP g/m2")+
#   #xlim(40,100)+
#   #ylim(0,100)
#   geom_smooth(method="lm", formula= y ~ x, se=FALSE, color="black", aes(group=1))
# 
# ggplot(data=siteout_fd, aes(x=CWM.PC2, y=weight_g_m, group=treatment, color=treatment))+
#   #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=treatment)) +
#   geom_point()+
#   theme_bw()+
#   #  scale_color_manual(values=c("tomato", "green3", "dodgerblue"), guide = guide_legend(title = "Treatment"), #change legend title
#   #                     labels=c("Mixed", "Forb", "Grass"))+ #change labels in the legend)+
#   xlab("Community Weighted Means of PC1 Scores")+
#   ylab("ANPP g/m2")+
#   #xlim(40,100)+
#   #ylim(0,100)
#   geom_smooth(method="lm", formula= y ~ x, se=FALSE,  aes(group=treatment))
# 
# siteout_fd<-merge(siteout_fd, BNPP1)
# 
# f_pca1_BNPP <- ggplot(data=siteout_fd, aes(x=CWM.PC1, y=(agg_BNPP), group=subplot, color=subplot))+
#   #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
#   geom_point()+
#   theme_bw()+
#   scale_color_manual(values=c("tomato", "green3", "dodgerblue"), guide = guide_legend(title = "Treatment"), #change legend title
#                      labels=c("Mixed", "Forb", "Grass"))+ #change labels in the legend)+
#   xlab("Community Weighted Means of PC1 Scores")+
#   ylab("BNPP g/m2")+
#   #xlim(40,100)+
#   #ylim(0,100)
#   geom_smooth(method="lm", formula= y ~ x, se=FALSE, color="black", aes(group=1))
# f_pca2_BNPP <- ggplot(data=siteout_fd, aes(x=CWM.PC2, y=(agg_BNPP), group=subplot, color=subplot))+
#   #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
#   geom_point()+
#   theme_bw()+
#   scale_color_manual(values=c("tomato", "green3", "dodgerblue"), guide = guide_legend(title = "Treatment"), #change legend title
#                      labels=c("Mixed", "Forb", "Grass"))+ #change labels in the legend)+
#   xlab("Community Weighted Means of PC2 Scores")+
#   ylab("BNPP g/m2")+
#   #xlim(40,100)+
#   #ylim(0,100)
#   geom_smooth(method="lm", formula= y ~ x, se=FALSE, color="black", aes(group=1))
# 
# ggarrange(f_pca1_ANPP, f_pca2_ANPP, f_pca1_BNPP, f_pca2_BNPP, ncol = 2, nrow = 2, common.legend = TRUE)
# 
# ## Rao Q of PCA
# Rb<-ggplot(data=siteout_fd, aes(x=RaoQ, y=agg_BNPP, group=subplot, color=subplot))+
#   stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=subplot)) +
#   geom_point()+
#   theme_bw()+
#   ylab("BNPP g/m2")+
#   #xlim(40,100)+
#   #ylim(0,100)
#   geom_smooth(method="lm", formula= y ~ x, se=FALSE,  aes(group=subplot))
# Rb
# 
# Ra<-ggplot(data=siteout_fd, aes(x=RaoQ, y=ANPPgm, group=subplot, color=subplot))+
#   stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=subplot)) +
#   geom_point()+
#   theme_bw()+
#   ylab("ANPP g/m2")+
#   #xlim(40,100)+
#   #ylim(0,100)
#   geom_smooth(method="lm", formula= y ~ x, se=FALSE,  aes(group=subplot))
# Ra
# 
# 
# ggplot(data=siteout_fd, aes(x=RaoQ, y=agg_BNPP, group=treatment, color=treatment))+
#   stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=treatment)) +
#   geom_point()+
#   theme_bw()+
#   ylab("BNPP g/m2")+
#   #xlim(40,100)+
#   #ylim(0,100)
#   geom_smooth(method="lm", formula= y ~ x, se=FALSE,  aes(group=treatment))
# 
# 
# ggplot(data=siteout_fd, aes(x=RaoQ, y=ANPPgm, group=subplot, color=subplot))+
#   stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=subplot)) +
#   geom_point()+
#   theme_bw()+
#   ylab("ANPP g/m2")+
#   #xlim(40,100)+
#   #ylim(0,100)
#   geom_smooth(method="lm", formula= y ~ x, se=FALSE,  aes(group=subplot))
# 
# ggplot(data=siteout_fd, aes(x=RaoQ, y=ANPPgm, group=treatment, color=treatment))+
#   stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=treatment)) +
#   geom_point()+
#   theme_bw()+
#   ylab("ANPP g/m2")+
#   #xlim(40,100)+
#   #ylim(0,100)
#   geom_smooth(method="lm", formula= y ~ x, se=FALSE,  aes(group=treatment))
# 
# 
# figure2_pca_raoq <- ggarrange(Ra, Rb, 
#                      ncol =2, common.legend = TRUE, legend = "bottom",
#                      align = "v",labels = c("a)", "b)"))
# figure2_pca_raoq
# 
# ggplot(data=siteout_fd, aes(x=CWM.PC1, y=agg_BNPP, group=subplot, color=subplot))+
#   stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
#   geom_point()+
#   theme_bw()+
#   #xlim(40,100)+
#   #ylim(0,100)
#   geom_smooth(method="lm", formula= y ~ x, se=FALSE, color="black", aes(group=1))
# 
# ggplot(subset(siteout_fd), aes(x=CWM.PC1, y=agg_BNPP, group=treatment, color=treatment))+
#   stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=treatment)) +
#   geom_point()+
#   theme_bw()+
#   #xlim(40,100)+
#   #ylim(0,100)
#   geom_smooth(method="lm", formula= y ~ x, se=FALSE,  aes(group=treatment))
# 
# ggplot(data=siteout_fd, aes(x=CWM.PC2, y=agg_BNPP, group=subplot, color=subplot))+
#   stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
#   geom_point()+
#   theme_bw()+
#   #xlim(40,100)+
#   #ylim(0,100)
#   geom_smooth(method="lm", formula= y ~ x, se=FALSE, color="black", aes(group=1))
# 
# ggplot(subset(siteout_fd), aes(x=CWM.PC2, y=agg_BNPP, group=treatment, color=treatment))+
#   stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=treatment)) +
#   geom_point()+
#   theme_bw()+
#   #xlim(40,100)+
#   #ylim(0,100)
#   geom_smooth(method="lm", formula= y ~ x, se=FALSE,  aes(group=treatment))
# 
# ggplot(data=siteout_fd, aes(x=CWM.PC2, y=(ANPPgm+agg_BNPP), group=subplot, color=subplot))+
#   stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
#   geom_point()+
#   theme_bw()+
#   xlab("Community Weighted Means of PC2 Scores")+
#   ylab("Total biomass g/m2")+
#   #xlim(40,100)+
#   #ylim(0,100)
#   geom_smooth(method="lm", formula= y ~ x, se=FALSE, color="black", aes(group=1))
# 
# ggplot(data=siteout_fd, aes(x=CWM.PC1, y=(ANPPgm+agg_BNPP), group=subplot, color=subplot))+
#   stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
#   geom_point()+
#   theme_bw()+
#   xlab("Community Weighted Means of PC2 Scores")+
#   ylab("Total biomass g/m2")+
#   #xlim(40,100)+
#   #ylim(0,100)
#   geom_smooth(method="lm", formula= y ~ x, se=FALSE, color="black", aes(group=1))

#############################################
##run PCA again but for aboveground traits only
traits_above <- as.matrix(all_trait_redo[,c(2,3,4)])

myrda2 <- rda(na.omit(traits_above), scale = TRUE)

# extract values
siteout2 <- as.data.frame(scores(myrda2, choices=c(1,2), display=c("sites")))
siteout2$ID<-all_trait_redo$ID
siteout2$name <- siteout2$ID

enviroout2<-as.data.frame(scores(myrda2, choices=c(1,2), display=c("species")))
enviroout2$type<-"traits"
enviroout2$name<-rownames(enviroout2)

# merge PC axes with trait data
tog2 <- left_join(all_trait_redo, siteout2) %>%
  mutate(func = paste(Origin, GF, sep = "_"))

#write.csv(tog2, "~/Documents/Repositories/usda-climvar/Manifesto/Ashley/2015/figures/above_pcscores.csv")

a.pca<-ggplot(tog2, aes(x=PC1, y=PC2))+ 
  ggtitle("a) PCA on aboveground traits")+
  geom_hline(aes(yintercept=0), color="grey") + 
  geom_vline(aes(xintercept=0), color="grey") +
  #geom_text(aes(label = name, color = GF), size = 5) +
  geom_point(aes(color=GF, shape=GF), size=4)+
  scale_color_manual(values = c( "#66c2a5", "#8da0cb","#66c2a5"), labels=c("Forb","Grass", "Legume"), guide = guide_legend(title = "Functional Group:")) +
  scale_shape_manual(values = c(15, 19,17), labels=c("Forb","Grass", "Legume"), guide = guide_legend(title = "Functional Group:")) +
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
  xlim(-2.4, 1.8) +
  theme_bw() +theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                    text=element_text(size = 20), legend.position="none")+ 
  xlab(paste("Axis 1 (",sprintf("%.1f",myrda2$CA$eig["PC1"]/myrda2$tot.chi*100,3),"%)",sep="")) +
  ylab(paste("Axis 2 (",sprintf("%.1f",myrda2$CA$eig["PC2"]/myrda2$tot.chi*100,3),"%)",sep="")) 
a.pca

#keep species that we have trait data for
#cover15_fd4 <- cover15_fd2 %>% dplyr::select(ACHMIL, AVEBAR, AVEFAT, BRADIS, BRIMIN, BRODIA, BROHOR, BROMAD, CARPYC, CENSOL, CERGLO, CLAPUR, CYNDAC, CYNECH, EROBOT, EROMOS, FILGAL, GALPAR,HORMAR, HORMUR,HYPGLA, HYPRAD, LACSER, LOLMUL, SILGAL, TAECAP,TORARV, TRIDUB,TRIGLO, TRIHIR,TRISP, TRISUB,VICSAT,VULBRO, VULMYU)
#cover15_fd4<-data.matrix(cover15_fd4)

#remove excess rows from PCA scores
#siteout2 <- siteout2[-c(15), ]
siteout2<- siteout2 %>% dplyr::select(-ID, -name)
# siteout2<-bind_cols(siteout2,cover_names)
row.names(siteout2) <- all_traits$ID
# row.names(siteout2)<-siteout2$Taxon 
# siteout2 <- siteout2 %>% dplyr::select(-Taxon)


siteout2_fd<-dbFD (siteout2, comp, w.abun = T, stand.x = F,
                  calc.FRic = TRUE, m = "max", stand.FRic = FALSE,
                  scale.RaoQ = FALSE, calc.FGR = FALSE, clust.type = "ward",
                  km.inf.gr = 2, km.sup.gr = nrow(x) - 1, km.iter = 100, calc.CWM = TRUE,
                  calc.FDiv = TRUE, print.pco = FALSE, messages = TRUE)
siteout2_fd<-as.data.frame(siteout2_fd)
#write.csv(siteout2, "./Repositories/usda-climvar/Manifesto/Lina/above_pcscores.csv")

siteout2_fd_biomass <- cbind(siteout2_fd, ANPP1) %>% 
  mutate(subplot = recode(subplot, B = "Mixed", F = "Forb", G = "Grass"))
siteout2_fd_biomass$treatment <- factor(siteout2_fd_biomass$treatment, levels =c("controlRain",  "springDry", "fallDry","consistentDry"))
siteout2_fd_biomass$subplot <- factor(siteout2_fd_biomass$subplot, levels =c("Mixed",  "Forb", "Grass"))
# siteout2_fd$ANPPgm<-traits_2015_2$ANPPgm
# siteout2_fd$shelterBlock<-traits_2015_2$shelterBlock
# siteout2_fd$subplot<-traits_2015_2$subplot
# siteout2_fd$treatment<-traits_2015_2$treatment

# a.1<-ggplot(data=siteout2_fd_biomass, aes(x=CWM.PC1, y=weight_g_m, group=subplot, color=subplot))+
#   #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
#   #ggtitle("a)")+
#   geom_point()+
#   theme_classic()+
#   stat_cor(aes(group=1,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), show.legend = FALSE)+
#   theme(legend.position="none")+
#   scale_color_manual(values=c("#fc8d62", "#66c2a5", "#8da0cb"), guide = guide_legend(title = "Treatment"), #change legend title
#                      labels=c("Mixed", "Forb", "Grass"))+ #change labels in the legend)+
#   labs(x=expression(atop("High LDMC" %<->% "High SLA","(CWM Aboveground PC1 Scores)")), y = expression(paste("ANPP ", (g/m^{2}))))+
#   #xlim(40,100)+
#   #ylim(0,100)
#   geom_smooth(method="lm", formula= y ~ x, se=FALSE, color="black", aes(group=1))
# a.1

a.t1<-ggplot(data=siteout2_fd_biomass, aes(x=CWM.PC1, y=weight_g_m, group=treatment, color=treatment, shape = subplot))+
  #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=treatment)) +
  stat_cor(aes(group=1,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), show.legend = FALSE)+
  #ggtitle("a)")+
  geom_point()+
  theme_classic()+
#  theme(legend.position="none")+
  scale_color_manual(name = "Rainfall Treatment", labels = c("Control", "Spring Dry", "Fall Dry","Consistent Dry" ), values= c("#0070b8", "#b2c7e4", "#fccaaf", "#c85b23"))+
  labs(x=expression(atop("High LDMC" %<->% "High SLA","(CWM Aboveground PC1 Scores)")), y = expression(paste("ANPP ", (g/m^{2}))), shape = "Composition Treatment")+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method="lm", formula= y ~ x, se=FALSE, color="black", aes(group=1))
a.t1

# a.2<-ggplot(data=siteout2_fd_biomass, aes(x=CWM.PC2, y=weight_g_m, group=subplot, color=subplot))+
#   #ggtitle("b)")+
#   #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
#   geom_point()+
#   theme_classic()+
#   theme(legend.position="none")+
#   scale_color_manual(values=c("#fc8d62", "#66c2a5", "#8da0cb"), guide = guide_legend(title = "Treatment"), #change legend title
#                      labels=c("Mixed", "Forb", "Grass"))+ #change labels in the legend)+
#   labs(x=expression(atop("Tall" %<->% "Short","(CWM Aboveground PC2 Scores)")), y = expression(paste("ANPP ", (g/m^{2}))))+
#   stat_cor(aes(group=1,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), show.legend = FALSE)+
#   #xlim(40,100)+
#   #ylim(0,100)
#   geom_smooth(method="lm", formula= y ~ x, se=FALSE, color="black", aes(group=1))
# a.2

a.t2<-ggplot(data=siteout2_fd_biomass, aes(x=CWM.PC2, y=weight_g_m, group=treatment, color=treatment, shape = subplot))+
  #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=treatment)) +
  stat_cor(aes(group=1,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), show.legend = FALSE)+
  #ggtitle("a)")+
  labs(x=expression(atop("Tall" %<->% "Short","(CWM Aboveground PC2 Scores)")), y = expression(paste("ANPP ", (g/m^{2}))))+
  geom_point()+
  theme_classic()+
  #  theme(legend.position="none")+
  scale_color_manual(name = "Treatment", labels = c("Control", "Spring Dry", "Fall Dry","Consistent Dry" ), values= c("#0070b8", "#b2c7e4", "#fccaaf", "#c85b23"))+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method="lm", formula= y ~ x, se=FALSE, color="black", aes(group=1))
a.t2

# siteout2_fd_biomass_all<-merge(siteout2_fd_biomass, BNPP1)

# t.1<-ggplot(data=siteout2_fd, aes(x=CWM.PC1, y=(ANPPgm+agg_BNPP), group=subplot, color=subplot))+
#   #ggtitle("b)")+
#   #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
#   geom_point()+
#   theme_bw()+
#   theme(legend.position="none")+
#   scale_color_manual(values=c("tomato", "green3", "dodgerblue"), guide = guide_legend(title = "Treatment"), #change legend title
#                      labels=c("Mixed", "Forb", "Grass"))+ #change labels in the legend)+
#   xlab("Community Weighted Means of Aboveground PC2 Scores")+
#   ylab("")+
#   stat_cor(aes(group=1,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
#   #xlim(40,100)+
#   #ylim(0,100)
#   geom_smooth(method="lm", formula= y ~ x, se=FALSE, color="black", aes(group=1))
# t.1

# t.2<-ggplot(data=siteout2_fd_biomass_all, aes(x=CWM.PC2, y=(agg_BNPP), group=subplot, color=subplot))+
#   #ggtitle("b)")+
#   #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
#   geom_point()+
#   theme_bw()+
#   theme(legend.position="none")+
#   scale_color_manual(values=c("tomato", "green3", "dodgerblue"), guide = guide_legend(title = "Treatment"), #change legend title
#                      labels=c("Mixed", "Forb", "Grass"))+ #change labels in the legend)+
#   xlab("Community Weighted Means of Aboveground PC2 Scores")+
#   ylab("")+
#   stat_cor(aes(group=1,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
#   #xlim(40,100)+
#   #ylim(0,100)
#   geom_smooth(method="lm", formula= y ~ x, se=FALSE, color="black", aes(group=1))
# t.2

##do again but for BELOWGROUND traits only
# select belowground traits
traits.below<-as.data.frame(all_trait_redo) %>% dplyr::select(Dens, DiamC, SRLC,SRLF, PropF)

# run PCA
myrda.b <- rda(na.omit(traits.below), scale = TRUE)

# extract values
siteout.b <- as.data.frame(scores(myrda.b, choices=c(1,2), display=c("sites")))
siteout.b$ID<-all_trait_redo$ID
siteout.b$name <- siteout.b$ID

enviroout.b<-as.data.frame(scores(myrda.b, choices=c(1,2), display=c("species")))
enviroout.b$type<-"traits"
enviroout.b$name<-rownames(enviroout.b)

# merge PC axes with trait data
tog3 <- left_join(all_trait_redo, siteout.b) %>%
  mutate(func = paste(Origin, GF, sep = "_"))

#to create supplementary table of PC scores by species
#write.csv(tog3, "~/Documents/Repositories/usda-climvar/Manifesto/Ashley/2015/figures/below_pcscores.csv")

# Remove legumes from legend key (if running PCA without legumes)
#tog <- subset(tog, GF != "L")

#pdf("TraitPCA.pdf", width = 9, height = 7.5)
#pdf("TraitPCA_noLegumes.pdf", width = 9, height = 7.5)

#write.csv(siteout.b, "./Repositories/usda-climvar/Manifesto/Lina/below_pcscores.csv")


b.pca<-ggplot(tog3, aes(x=PC1, y=PC2))+ 
  ggtitle("b) PCA on belowground traits")+
  geom_hline(aes(yintercept=0), color="grey") + 
  geom_vline(aes(xintercept=0), color="grey") +
  geom_point(aes(color=GF, shape=GF), size=4)+
  #geom_text(aes(label = name, color = GF), size = 5) +
  scale_color_manual(values = c("#66c2a5", "#8da0cb","#66c2a5"), labels=c("Forb","Grass", "Legume"), guide = guide_legend(title = "Functional Group:")) +
  scale_shape_manual(values = c(15, 19,17), labels=c("Forb","Grass", "Legume"), guide = guide_legend(title = "Functional Group:")) +
  geom_segment(data = enviroout.b,
               aes(x = 0, xend =  PC1,
                   y = 0, yend =  PC2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "black") + #grid is required for arrow to work.
  geom_text(data = enviroout.b,
            aes(x=  PC1*1.2, y =  PC2*1.2, #we add 10% to the text to push it slightly out from arrows
                label = name), #otherwise you could use hjust and vjust. I prefer this option
            size = 6,
            hjust = 0.5, 
            color="black") + 
  xlim(-2.4, 2.2) +
  theme_bw() +theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                    text=element_text(size = 20))+ 
  xlab(paste("Axis 1 (",sprintf("%.1f",myrda.b$CA$eig["PC1"]/myrda.b$tot.chi*100,3),"%)",sep="")) +
  ylab(paste("Axis 2 (",sprintf("%.1f",myrda.b$CA$eig["PC2"]/myrda.b$tot.chi*100,3),"%)",sep="")) 
b.pca
#dev.off()

grid.arrange(a.pca, b.pca,  ncol = 2, widths = c(1,1.3))


figure3 <- ggarrange(a.pca, b.pca,
                      ncol =2, nrow =1, common.legend = TRUE, legend = "bottom",
                      align = "v")
figure3

#remove excess rows from PCA scores
#siteout.b <- siteout.b[-c(14), ] #remove LUPBIC
#siteout.b<- siteout.b %>% dplyr::select(-ID, -name)
#keep species that we have trait data for
#cover15_fd3 <- cover15_fd2 %>% dplyr::select(ACHMIL, AVEFAT, BRADIS, BRODIA, BROHOR, BROMAD, CENSOL, CYNDAC, CYNECH, EROBOT, HORMUR, LACSER, LOLMUL, TAECAP,TRIHIR,VULMYU)
#cover15_fd3<-data.matrix(cover15_fd3)
siteout.b<- siteout.b %>% dplyr::select(-ID, -name)
row.names(siteout.b) <- all_traits$ID
#siteout.b<-bind_cols(siteout.b,cover_names)
#row.names(siteout.b)<-siteout.b$Taxon 
#siteout.b <- siteout.b %>% dplyr::select(-Taxon)

siteout_fd.b<-dbFD (siteout.b, comp, w.abun = T, stand.x = F,
                  calc.FRic = TRUE, m = "max", stand.FRic = FALSE,
                  scale.RaoQ = FALSE, calc.FGR = FALSE, clust.type = "ward",
                  km.inf.gr = 2, km.sup.gr = nrow(x) - 1, km.iter = 100, calc.CWM = TRUE,
                  calc.FDiv = TRUE, print.pco = FALSE, messages = TRUE)
siteout_fd.b<-as.data.frame(siteout_fd.b)


# 
# traits_2015_2<- arrange(traits_2015, treatment, shelterBlock)
# siteout_fd.b$ANPPgm<-traits_2015_2$ANPPgm
# siteout_fd.b$shelterBlock<-traits_2015_2$shelterBlock
# siteout_fd.b$subplot<-traits_2015_2$subplot
# siteout_fd.b$treatment<-traits_2015_2$treatment
# 
# 
# ggplot(data=siteout_fd.b, aes(x=CWM.PC1, y=ANPPgm, group=subplot, color=subplot))+
#   #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
#   geom_point()+
#   theme_bw()+
#   #xlim(40,100)+
#   #ylim(0,100)
#   geom_smooth(method="lm", formula= y ~ x, se=FALSE, color="black", aes(group=1))
# 
# ggplot(data=siteout_fd.b, aes(x=CWM.PC2, y=ANPPgm, group=subplot, color=subplot))+
#   stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
#   geom_point()+
#   theme_bw()+
#   #xlim(40,100)+
#   #ylim(0,100)
#   geom_smooth(method="lm", formula= y ~ x, se=FALSE, color="black", aes(group=1))

#run Lina's BNPP_exploratory_analysis.R script to get BNPP1
siteout_fd.b_BNPP<-cbind(siteout_fd.b, BNPP1)
siteout_fd.b_BNPP$treatment <- factor(siteout_fd.b_BNPP$treatment, levels =c("controlRain",  "springDry", "fallDry","consistentDry"))
# b.1<-ggplot(data=siteout_fd.b_BNPP, aes(x=CWM.PC1, y=agg_BNPP, group=treatment, color=treatment))+
#   #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
#   scale_color_manual(name = "Treatment", labels = c("Control", "Spring Dry", "Fall Dry","Consistent Dry" ), values= c("#0070b8", "#b2c7e4", "#fccaaf", "#c85b23"))+
#   labs(x=expression(atop("Coarse" %<->% "Dense, More fine roots","(CWM Belowground PC1 Scores)")), y = expression(paste("BNPP ", (g/m^{2}))))+
#     #ggtitle("c)")+
#   geom_point()+
#   theme_classic()+
#   #xlim(40,100)+
#   #ylim(0,100)
#   stat_cor(aes(group=1,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),  show.legend = FALSE)+
#   geom_smooth(method="lm", formula= y ~ x, se=FALSE, color="black", aes(group=1))+
#   theme(legend.position="none")
# b.1

bt.1<-ggplot(data=siteout_fd.b_BNPP, aes(x=CWM.PC1, y=agg_BNPP, group=treatment, color=treatment, shape = subplot))+
  #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
  scale_color_manual(name = "Treatment", labels = c("Control", "Spring Dry", "Fall Dry","Consistent Dry" ), values= c("#0070b8", "#b2c7e4", "#fccaaf", "#c85b23"))+
  labs(x=expression(atop("Coarse roots" %<->% "Dense, More fine roots","(CWM Belowground PC1 Scores)")), y = expression(paste("BNPP ", (g/m^{2}))))+
  #ggtitle("c)")+
  geom_point()+
  theme_classic()+
  #xlim(40,100)+
  #ylim(0,100)
  stat_cor(aes(group=1,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), show.legend = FALSE)+
  geom_smooth(method="lm", formula= y ~ x, se=FALSE,  color = "black", aes(group=1))+
  theme(legend.position="none")
bt.1

# below.1<-ggplot(data=siteout_fd.b_BNPP, aes(y=CWM.PC1, x=as.factor(shelter), color=as.factor(shelter)))+
#   #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
#   #scale_color_manual(values=c("tomato", "green3", "dodgerblue"), guide = guide_legend(title = "Treatment"), #change legend title
#   #                   labels=c("Mixed", "Forb", "Grass"))+ #change labels in the legend)+
#   ylab("Community Weighted Means of Belowground PC1 Scores")+
#   xlab("Treatment")+
#   #ggtitle("c)")+
#   geom_boxplot()+
#   theme_bw()+
#   #xlim(40,100)+
#   #ylim(0,100)
#   #stat_cor(aes(group=treatment,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
#   #geom_smooth(method="lm", formula= y ~ x, se=FALSE,  aes(group=treatment))+
#   theme(legend.position="none")
# below.1
# 
# below.2<-ggplot(data=siteout_fd, aes(y=CWM.PC2, x=as.factor(shelter), color=as.factor(shelter)))+
#   #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
#   #scale_color_manual(values=c("tomato", "green3", "dodgerblue"), guide = guide_legend(title = "Treatment"), #change legend title
#   #                   labels=c("Mixed", "Forb", "Grass"))+ #change labels in the legend)+
#   ylab("Community Weighted Means of Belowground PC2 Scores")+
#   xlab("Treatment")+
#   #ggtitle("c)")+
#   geom_boxplot()+
#   theme_bw()+
#   #xlim(40,100)+
#   #ylim(0,100)
#   #stat_cor(aes(group=treatment,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
#   #geom_smooth(method="lm", formula= y ~ x, se=FALSE,  aes(group=treatment))+
#   theme(legend.position="none")
# below.2
# 
# b.2<-ggplot(data=siteout_fd, aes(x=CWM.PC2, y=agg_BNPP, group=subplot, color=subplot))+
#   #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
#   scale_color_manual(values=c("#fc8d62", "#66c2a5", "#8da0cb"), guide = guide_legend(title = "Treatment"), #change legend title
#                      labels=c("Mixed", "Forb", "Grass"))+ #change labels in the legend)+
#   labs(x=expression(atop("Fine" %<->% "Coarse","(CWM Belowground PC2 Scores)")), y = expression(paste("BNPP ", (g/m^{2}))))+
#   #ggtitle("d)")+
#   geom_point()+
#   theme_classic()+
#   #xlim(40,100)+
#   #ylim(0,100)
#   stat_cor(aes(group=1,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), show.legend = FALSE)+
#   #stat_cor(aes())+ #for p-values and R2
#   geom_smooth(method="lm", formula= y ~ x, se=FALSE, color="black", aes(group=1))+
#   theme(legend.position="none")
# b.2

bt.2<-ggplot(data=siteout_fd.b_BNPP, aes(x=CWM.PC2, y=agg_BNPP, group=treatment, color=treatment, shape = subplot))+
  #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
  scale_color_manual(name = "Treatment", labels = c("Control", "Spring Dry", "Fall Dry","Consistent Dry" ), values= c("#0070b8", "#b2c7e4", "#fccaaf", "#c85b23"))+
  labs(x=expression(atop("Long roots" %<->% "Short roots","(CWM Belowground PC2 Scores)")), y = expression(paste("BNPP ", (g/m^{2}))))+
  #ggtitle("c)")+
  geom_point()+
  theme_classic()+
  #xlim(40,100)+
  #ylim(0,100)
  stat_cor(aes(group=1,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), show.legend = FALSE)+
  geom_smooth(method="lm", formula= y ~ x, se=FALSE, color = "black", aes(group=1))+
  theme(legend.position="none")
bt.2

Figure4 <- ggarrange(a.t1,bt.1, a.t2,  bt.2, ncol = 2, nrow =2, common.legend = TRUE, labels = c("a)", "b)", "c)", "d)"), legend = "right")
# b.3<-ggplot(data=siteout_fd, aes(x=CWM.PC2, y=agg_BNPP, group=subplot, color=subplot))+
#   #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
#   scale_color_manual(values=c("tomato", "green3", "dodgerblue"), guide = guide_legend(title = "Treatment:"), #change legend title
#                      labels=c("Mixed", "Forb", "Grass"))+ #change labels in the legend)+
#   xlab("Community Weighted Means of Belowground PC2 Scores")+
#   ylab("BNPP g/m2")+
#   #ggtitle("d)")+
#   geom_point()+
#   theme_bw()+
#   #xlim(40,100)+
#   #ylim(0,100)
#   #geom_smooth(method="lm", formula= y ~ x, se=FALSE, color="black", aes(group=1))+
#   theme(legend.position="bottom")
# b.3
# 
# t.3<-ggplot(data=siteout_fd, aes(x=CWM.PC1, y=(agg_BNPP+ANPPgm), group=subplot, color=subplot))+
#   #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
#   scale_color_manual(values=c("tomato", "green3", "dodgerblue"), guide = guide_legend(title = "Treatment"), #change legend title
#                      labels=c("Mixed", "Forb", "Grass"))+ #change labels in the legend)+
#   xlab("Community Weighted Means of Belowground PC2 Scores")+
#   ylab("")+
#   #ggtitle("d)")+
#   geom_point()+
#   theme_bw()+
#   #xlim(40,100)+
#   #ylim(0,100)
#   stat_cor(aes(group=1,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
#   #stat_cor(aes())+ #for p-values and R2
#   geom_smooth(method="lm", formula= y ~ x, se=FALSE, color="black", aes(group=1))+
#   theme(legend.position="none")
# t.3
# 
# t.4<-ggplot(data=siteout_fd, aes(x=CWM.PC2, y=(agg_BNPP+ANPPgm), group=subplot, color=subplot))+
#   #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
#   scale_color_manual(values=c("tomato", "green3", "dodgerblue"), guide = guide_legend(title = "Treatment"), #change legend title
#                      labels=c("Mixed", "Forb", "Grass"))+ #change labels in the legend)+
#   xlab("Community Weighted Means of Belowground PC2 Scores")+
#   ylab("")+
#   #ggtitle("d)")+
#   geom_point()+
#   theme_bw()+
#   #xlim(40,100)+
#   #ylim(0,100)
#   stat_cor(aes(group=1,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
#   #stat_cor(aes())+ #for p-values and R2
#   geom_smooth(method="lm", formula= y ~ x, se=FALSE, color="black", aes(group=1))+
#   theme(legend.position="none")
# t.4
# 
# g_legend<-function(a.gplot){
#   tmp <- ggplot_gtable(ggplot_build(a.gplot))
#   leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
#   legend <- tmp$grobs[[leg]]
#   return(legend)}
# 
# mylegend<-g_legend(b.3)
# 
# lay <- rbind(c(1,1,2,2),
#              c(1,1,2,2),
#              c(3,3,4,4),
#              c(3,3,4,4),
#              c(NA,5,5,NA))
# 
# grid.arrange(a.1, a.2, b.1, b.2,mylegend,layout_matrix = lay)
# 
# figure4 <- ggarrange(a.1, a.2, b.1, b.2,
#                       ncol =2, nrow =2, common.legend = TRUE, legend = "bottom",
#                       align = "v",labels = c("a)", "b)", "c)", "d)"))
# figure4
# 
# figures7 <- ggarrange(a.t1, a.t2, bt.1, bt.2,
#                      ncol =2, nrow =2, common.legend = TRUE, legend = "bottom",
#                      align = "v",labels = c("a)", "b)", "c)", "d)"))
# figures7

#Predictive models
unique_sp_cover <- cover15 %>%
  dplyr::select(plot, subplot, treatment, Centaurea_solstitialis, Taeniatherum_caput_medusae)
siteout_fd.b_BNPP_unique <- left_join(siteout_fd.b_BNPP, unique_sp_cover)
siteout2_fd_biomass_unique <- left_join(siteout2_fd_biomass, unique_sp_cover)

#Following Diaz et al. method
summary(lm(weight_g_m ~ CWM.PC1, siteout2_fd_biomass_unique)) #sig
summary(lm(weight_g_m ~ CWM.PC2, siteout2_fd_biomass_unique))
summary(lm(weight_g_m ~ FEve, siteout2_fd_biomass_unique))
summary(lm(weight_g_m ~ RaoQ, siteout2_fd_biomass_unique))
summary(lm(weight_g_m ~ Centaurea_solstitialis, siteout2_fd_biomass_unique))
summary(lm(weight_g_m ~ Taeniatherum_caput_medusae, siteout2_fd_biomass_unique))

summary(lm(agg_BNPP ~ CWM.PC1, siteout_fd.b_BNPP_unique))
summary(lm(agg_BNPP ~ CWM.PC2, siteout_fd.b_BNPP_unique))
summary(lm(agg_BNPP ~ FEve, siteout_fd.b_BNPP_unique))
summary(lm(agg_BNPP ~ RaoQ, siteout_fd.b_BNPP_unique))
summary(lm(agg_BNPP ~ Centaurea_solstitialis, siteout_fd.b_BNPP_unique))
summary(lm(agg_BNPP ~ Taeniatherum_caput_medusae, siteout_fd.b_BNPP_unique))

##Table Backwards stepwise regression
fullmodel_ANPP <- lm(weight_g_m ~ FEve + RaoQ + CWM.PC1 + CWM.PC2 +Centaurea_solstitialis+Taeniatherum_caput_medusae, siteout2_fd_biomass_unique)
step_fullmodel_ANPP  <- stepAIC(fullmodel_ANPP, direction = "backward", trace = FALSE)
step_fullmodel_ANPP$anova
best_ANPP <- lm(weight_g_m ~ CWM.PC1, siteout2_fd_biomass_unique)
summary(best_ANPP)

fullmodel_ANPP_CD <- lm(weight_g_m ~ FEve + RaoQ + CWM.PC1 + CWM.PC2 +Centaurea_solstitialis+Taeniatherum_caput_medusae, siteout2_fd_biomass_unique %>% filter(treatment == "consistentDry"))
step_fullmodel_ANPP_CD  <- stepAIC(fullmodel_ANPP_CD, direction = "backward", trace = FALSE)
step_fullmodel_ANPP_CD$anova
model_ANPP_CD <- lm(weight_g_m ~ CWM.PC2 , siteout2_fd_biomass_unique %>% filter(treatment == "consistentDry"))
summary(model_ANPP_CD)

fullmodel_ANPP_control <- lm(weight_g_m ~ FEve + RaoQ + CWM.PC1 + CWM.PC2 +Centaurea_solstitialis+Taeniatherum_caput_medusae, siteout2_fd_biomass_unique %>% filter(treatment == "controlRain"))
step_fullmodel_ANPP_control  <- stepAIC(fullmodel_ANPP_control, direction = "backward", trace = FALSE)
step_fullmodel_ANPP_control$anova
model_ANPP_control <- lm(weight_g_m ~ CWM.PC1 , siteout2_fd_biomass_unique %>% filter(treatment == "controlRain"))
summary(model_ANPP_control)

fullmodel_ANPP_SD <- lm(weight_g_m ~ FEve + RaoQ + CWM.PC1 + CWM.PC2 +Centaurea_solstitialis+Taeniatherum_caput_medusae, siteout2_fd_biomass_unique %>% filter(treatment == "springDry"))
step_fullmodel_ANPP_SD  <- stepAIC(fullmodel_ANPP_SD, direction = "backward", trace = FALSE)
step_fullmodel_ANPP_SD$anova
model_ANPP_SD <- lm(weight_g_m ~ FEve , siteout2_fd_biomass_unique %>% filter(treatment == "springDry"))
summary(model_ANPP_SD)

fullmodel_ANPP_FD <- lm(weight_g_m ~ FEve + RaoQ + CWM.PC1 + CWM.PC2 +Centaurea_solstitialis+Taeniatherum_caput_medusae, siteout2_fd_biomass_unique %>% filter(treatment == "fallDry"))
step_fullmodel_ANPP_FD  <- stepAIC(fullmodel_ANPP_FD, direction = "backward", trace = FALSE)
step_fullmodel_ANPP_FD$anova
model_ANPP_FD <- lm(weight_g_m ~ FEve + CWM.PC1 + CWM.PC2 , siteout2_fd_biomass_unique %>% filter(treatment == "fallDry"))
summary(model_ANPP_FD)

fullmodel_BNPP <- lm(agg_BNPP ~ FEve + RaoQ + CWM.PC1 + CWM.PC2 +Centaurea_solstitialis+Taeniatherum_caput_medusae, siteout_fd.b_BNPP_unique)
step_fullmodel_BNPP  <- stepAIC(fullmodel_BNPP, direction = "backward", trace = FALSE)
step_fullmodel_BNPP$anova

fullmodel_BNPP_CD <- lm(agg_BNPP ~ FEve + RaoQ + CWM.PC1 + CWM.PC2 +Centaurea_solstitialis+Taeniatherum_caput_medusae, siteout_fd.b_BNPP_unique %>% filter(treatment == "consistentDry"))
step_fullmodel_BNPP_CD  <- stepAIC(fullmodel_BNPP_CD, direction = "backward", trace = FALSE)
step_fullmodel_BNPP_CD$anova
model_BNPP_CD <- lm(agg_BNPP ~ RaoQ + Centaurea_solstitialis, siteout_fd.b_BNPP_unique %>% filter(treatment == "consistentDry"))
summary(model_BNPP_CD)

fullmodel_BNPP_control <- lm(agg_BNPP ~ FEve + RaoQ + CWM.PC1 + CWM.PC2 +Centaurea_solstitialis+Taeniatherum_caput_medusae, siteout_fd.b_BNPP_unique %>% filter(treatment == "controlRain"))
step_fullmodel_BNPP_control  <- stepAIC(fullmodel_BNPP_control, direction = "backward", trace = FALSE)
step_fullmodel_BNPP_control$anova
model_BNPP_control <- lm(agg_BNPP ~ FEve + CWM.PC1 + Taeniatherum_caput_medusae, siteout_fd.b_BNPP_unique %>% filter(treatment == "controlRain"))
summary(model_BNPP_control)

fullmodel_BNPP_SD <- lm(agg_BNPP ~ FEve + RaoQ + CWM.PC1 + CWM.PC2 +Centaurea_solstitialis+Taeniatherum_caput_medusae, siteout_fd.b_BNPP_unique %>% filter(treatment == "springDry"))
step_fullmodel_BNPP_SD  <- stepAIC(fullmodel_BNPP_SD, direction = "backward", trace = FALSE)
step_fullmodel_BNPP_SD$anova
model_BNPP_SD <- lm(agg_BNPP ~ 1, siteout_fd.b_BNPP_unique %>% filter(treatment == "springDry"))
summary(model_BNPP_SD)

fullmodel_BNPP_FD <- lm(agg_BNPP ~ FEve + RaoQ + CWM.PC1 + CWM.PC2 +Centaurea_solstitialis+Taeniatherum_caput_medusae, siteout_fd.b_BNPP_unique %>% filter(treatment == "fallDry"))
step_fullmodel_BNPP_FD  <- stepAIC(fullmodel_BNPP_FD, direction = "backward", trace = FALSE)
step_fullmodel_BNPP_FD$anova

#### Fig 5 trait shift with rainfall treatment
above_rainfall_tr <- siteout2_fd_biomass %>% mutate(treatment=ordered(treatment, levels = c(controlRain="controlRain", springDry="springDry" , fallDry="fallDry", consistentDry="consistentDry"))) %>%
  group_by(treatment) %>% 
  summarize(CWM.PC1.m=mean(CWM.PC1), CWM.PC2.m=mean(CWM.PC2), se.PC1=calcSE(CWM.PC1),se.PC2=calcSE(CWM.PC2))

below_rainfall_tr <- siteout_fd.b_BNPP %>% mutate(treatment=ordered(treatment, levels = c(controlRain="controlRain", springDry="springDry" , fallDry="fallDry", consistentDry="consistentDry"))) %>%
  group_by(treatment) %>% 
  summarize(CWM.PC1.m=mean(CWM.PC1), CWM.PC2.m=mean(CWM.PC2), se.PC1=calcSE(CWM.PC1),se.PC2=calcSE(CWM.PC2))

fit_above_a <- aov(CWM.PC1~treatment, siteout2_fd_biomass)
TukeyHSD(aov(fit_above_a))
fit_above_c <- aov(CWM.PC2~treatment, siteout2_fd_biomass)
TukeyHSD(aov(fit_above_c))
fit_below_b <- aov(CWM.PC1~treatment, siteout_fd.b_BNPP)
TukeyHSD(aov(fit_below_b))
fit_below_d <- aov(CWM.PC2~treatment, siteout_fd.b_BNPP)
TukeyHSD(aov(fit_below_d))

#no interaction effects
summary(aov(CWM.PC1~subplot*treatment, siteout2_fd_biomass))
summary(aov(CWM.PC2~subplot*treatment, siteout2_fd_biomass))
summary(aov(CWM.PC1~subplot*treatment, siteout_fd.b_BNPP))
summary(aov(CWM.PC2~subplot*treatment, siteout_fd.b_BNPP))

F5a<-ggplot(data=subset(above_rainfall_tr), aes(x=treatment, y=CWM.PC1.m, col = treatment))+
  geom_point(position=position_dodge(0.9),size=4)+
  labs(y=expression(atop("High LDMC" %<->% "High SLA","(CWM Aboveground PC1 Scores)")), color = "Rain treatment")+
  theme(axis.title.x=element_blank(),
        text = element_text(size=14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2))+
  scale_x_discrete(labels = c("Control\n ","Spring\n Dry\n", "Fall\n Dry\n", "Consistent\n Dry\n"))+
  scale_color_manual(name = "Treatment", labels = c("Control", "Spring Dry", "Fall Dry","Consistent Dry" ), values= c("#0070b8", "#b2c7e4", "#fccaaf", "#c85b23"))+
  geom_errorbar(aes(ymin = CWM.PC1.m-se.PC1, ymax = CWM.PC1.m+se.PC1), size = 0.5, width = 0,position=position_dodge(0.9))+
  annotate("text", x = 1, y = -0.1, label = "ab", size = 6) +
  annotate("text", x = 2, y = -0.1, label = "a", size = 6) +
  annotate("text", x = 3, y = -0.1, label = "ab", size = 6) +
  annotate("text", x = 4, y = -0.1, label = "b", size = 6)
F5a

F5c<-ggplot(data=subset(above_rainfall_tr), aes(x=treatment, y=CWM.PC2.m, col = treatment))+
  geom_point(position=position_dodge(0.9), size=4)+
  labs(y=expression(atop("Tall" %<->% "Short","(CWM Aboveground PC2 Scores)")))+
  scale_x_discrete(labels = c("Control\n ","Spring\n Dry\n", "Fall\n Dry\n", "Consistent\n Dry\n"))+
  scale_color_manual(name = "Treatment", labels = c("Control", "Spring Dry", "Fall Dry","Consistent Dry" ), values= c("#0070b8", "#b2c7e4", "#fccaaf", "#c85b23"))+
  geom_errorbar(aes(ymin = CWM.PC2.m-se.PC2, ymax = CWM.PC2.m+se.PC2), size = 0.5, width = 0,position=position_dodge(0.9))+
  theme(axis.title.x=element_blank(),
        text = element_text(size=14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2))+
  annotate("text", x = 1, y = -0.01, label = "a", size = 6) +
  annotate("text", x = 2, y = -0.01, label = "a", size = 6) +
  annotate("text", x = 3, y = -0.01, label = "a", size = 6) +
  annotate("text", x = 4, y = -0.01, label = "a", size = 6)
F5c

F5b<-ggplot(data=subset(below_rainfall_tr), aes(x=treatment, y=CWM.PC1.m, col = treatment))+
  geom_point(position=position_dodge(0.9),size=4)+
  labs(y=expression(atop("Coarse roots" %<->% "Dense, fine roots","(CWM Belowground PC1 Scores)")))+
  scale_x_discrete(labels = c("Control\n ","Spring\n Dry\n", "Fall\n Dry\n", "Consistent\n Dry\n"))+
  scale_color_manual(name = "Treatment", labels = c("Control", "Spring Dry", "Fall Dry","Consistent Dry" ), values= c("#0070b8", "#b2c7e4", "#fccaaf", "#c85b23"))+
  geom_errorbar(aes(ymin = CWM.PC1.m-se.PC1, ymax = CWM.PC1.m+se.PC1), size = 0.5, width = 0,position=position_dodge(0.9))+
  theme(axis.title.x=element_blank(),
        text = element_text(size=14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2))+
  annotate("text", x = 1, y = 0.4, label = "ab", size = 6) +
  annotate("text", x = 2, y = 0.4, label = "a", size = 6) +
  annotate("text", x = 3, y = 0.4, label = "b", size = 6) +
  annotate("text", x = 4, y = 0.4, label = "ab", size = 6)
F5b

F5d<-ggplot(data=subset(below_rainfall_tr), aes(x=treatment, y=CWM.PC2.m, col = treatment))+
  geom_point(position=position_dodge(0.9), size=4)+
  labs(y=expression(atop("Long roots" %<->% "Short roots","(CWM Belowground PC2 Scores)")))+
  scale_x_discrete(labels = c("Control\n ","Spring\n Dry\n", "Fall\n Dry\n", "Consistent\n Dry\n"))+
  scale_color_manual(name = "Treatment", labels = c("Control", "Spring Dry", "Fall Dry","Consistent Dry" ), values= c("#0070b8", "#b2c7e4", "#fccaaf", "#c85b23"))+
  geom_errorbar(aes(ymin = CWM.PC2.m-se.PC2, ymax = CWM.PC2.m+se.PC2), size = 0.5, width = 0,position=position_dodge(0.9))+
  theme(axis.title.x=element_blank(),
        text = element_text(size=14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2))+
  annotate("text", x = 1, y = 0.5, label = "ab", size = 6) +
  annotate("text", x = 2, y = 0.5, label = "a", size = 6) +
  annotate("text", x = 3, y = 0.5, label = "b", size = 6) +
  annotate("text", x = 4, y = 0.5, label = "b", size = 6)
F5d

figure5 <- ggarrange(F5a, F5b, F5c, F5d,
                      ncol =2, nrow =2, 
                      align = "v",labels = c("e)", "f)", "g)", "h)"), common.legend = TRUE, legend = "top")
figure5
#annotate_figure(figure5, bottom = text_grob("Rainfall treatment", color = "black", size = 14))

above_rain_keys <- left_join(veg_keys, siteout2_fd_biomass)
below_rain_keys <- left_join(veg_keys, siteout_fd.b_BNPP)

mf6a<-lme(CWM.PC1 ~treatment, random=~1|shelterBlock, above_rain_keys , na.action=na.exclude)
summary(mf6a)
anova(mf6a) #treatment is significant
r.squaredGLMM(mf6a) 
qqnorm(residuals(mf6a))
qqline(residuals(mf6a))
shapiro.test(residuals(mf6a))
LS1b<-lsmeans(mf6a, ~treatment)
contrast(LS1b, "pairwise")

mf6a<-lme(CWM.PC2 ~treatment, random=~1|shelterBlock, above_rain_keys , na.action=na.exclude)
summary(mf6a)
anova(mf6a) #treatment is not significant
r.squaredGLMM(mf6a) 
qqnorm(residuals(mf6a))
qqline(residuals(mf6a))
shapiro.test(residuals(mf6a))
LS1b<-lsmeans(mf6a, ~treatment)
contrast(LS1b, "pairwise")

mf6a<-lme(CWM.PC1 ~treatment, random=~1|shelterBlock, below_rain_keys , na.action=na.exclude)
summary(mf6a)
anova(mf6a) #treatment is significant
r.squaredGLMM(mf6a) 
qqnorm(residuals(mf6a))
qqline(residuals(mf6a))
shapiro.test(residuals(mf6a))
LS1b<-lsmeans(mf6a, ~treatment)
contrast(LS1b, "pairwise")

mf6a<-lme(CWM.PC2 ~treatment, random=~1|shelterBlock, below_rain_keys , na.action=na.exclude)
summary(mf6a)
anova(mf6a) #treatment is significant
r.squaredGLMM(mf6a) 
qqnorm(residuals(mf6a))
qqline(residuals(mf6a))
shapiro.test(residuals(mf6a))
LS1b<-lsmeans(mf6a, ~treatment)
contrast(LS1b, "pairwise")

#### Supp Fig  trait shift with seeding treatment
above_seeding_tr <- siteout2_fd_biomass %>% 
  group_by(subplot) %>% 
  summarize(CWM.PC1.m=mean(CWM.PC1), CWM.PC2.m=mean(CWM.PC2), se.PC1=calcSE(CWM.PC1),se.PC2=calcSE(CWM.PC2))

below_seeding_tr <- siteout_fd.b_BNPP %>% 
  group_by(subplot) %>% 
  summarize(CWM.PC1.m=mean(CWM.PC1), CWM.PC2.m=mean(CWM.PC2), se.PC1=calcSE(CWM.PC1),se.PC2=calcSE(CWM.PC2))

fit_above_seed_a <- aov(CWM.PC1~subplot, siteout2_fd_biomass)
TukeyHSD(aov(fit_above_seed_a))
fit_above_seed_c <- aov(CWM.PC2~subplot, siteout2_fd_biomass)
TukeyHSD(aov(fit_above_seed_c))
fit_below_seed_b <- aov(CWM.PC1~subplot, siteout_fd.b_BNPP)
TukeyHSD(aov(fit_below_seed_b))
fit_below_seed_d <- aov(CWM.PC2~subplot, siteout_fd.b_BNPP)
TukeyHSD(aov(fit_below_seed_d))

S9a<-ggplot(data=subset(above_seeding_tr), aes(x=subplot, y=CWM.PC1.m, col = as.factor(subplot)))+
  geom_point(position=position_dodge(0.9),size=4)+
  labs(y=expression(atop("High LDMC" %<->% "High SLA","(CWM Aboveground PC1 Scores)")), color = "Treatment")+
  theme(axis.title.x=element_blank(),
        text = element_text(size=14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2))+
  scale_x_discrete(labels = c("Mixed","Forb", "Grass"))+
  scale_color_manual(labels = c("Mixed","Forb", "Grass"),  values = c("#fc8d62", "#66c2a5", "#8da0cb"))+
  geom_errorbar(aes(ymin = CWM.PC1.m-se.PC1, ymax = CWM.PC1.m+se.PC1), size = 0.5, width = 0,position=position_dodge(0.9))+
  annotate("text", x = 1, y = -0.1, label = "a", size = 6) +
  annotate("text", x = 2, y = -0.1, label = "b", size = 6) +
  annotate("text", x = 3, y = -0.1, label = "ab", size = 6) 
S9a

S9c <-ggplot(data=subset(above_seeding_tr), aes(x=subplot, y=CWM.PC2.m, col = as.factor(subplot)))+
  geom_point(position=position_dodge(0.9), size=4)+
  labs(y=expression(atop("Tall" %<->% "Short","(CWM Aboveground PC2 Scores)")))+
  scale_x_discrete(labels = c("Mixed","Forb", "Grass"))+
  scale_color_manual(labels = c("Mixed","Forb", "Grass"),  values = c("#fc8d62", "#66c2a5", "#8da0cb"))+
  geom_errorbar(aes(ymin = CWM.PC2.m-se.PC2, ymax = CWM.PC2.m+se.PC2), size = 0.5, width = 0,position=position_dodge(0.9))+
  theme(axis.title.x=element_blank(),
        text = element_text(size=14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2))+
  annotate("text", x = 1, y = 0.3, label = "a", size = 6) +
  annotate("text", x = 2, y = 0.3, label = "b", size = 6) +
  annotate("text", x = 3, y = 0.3, label = "a", size = 6) 
S9c

S9b<-ggplot(data=subset(below_seeding_tr), aes(x=subplot, y=CWM.PC1.m, col = as.factor(subplot)))+
  geom_point(position=position_dodge(0.9),size=4)+
  labs(y=expression(atop("Coarse roots" %<->% "Dense, fine roots","(CWM Belowground PC1 Scores)")))+
  scale_x_discrete(labels = c("Mixed","Forb", "Grass"))+
  scale_color_manual(labels = c("Mixed","Forb", "Grass"),  values = c("#fc8d62", "#66c2a5", "#8da0cb"))+
  geom_errorbar(aes(ymin = CWM.PC1.m-se.PC1, ymax = CWM.PC1.m+se.PC1), size = 0.5, width = 0,position=position_dodge(0.9))+
  theme(axis.title.x=element_blank(),
        text = element_text(size=14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2))+
  annotate("text", x = 1, y = 0.3, label = "a", size = 6) +
  annotate("text", x = 2, y = 0.3, label = "a", size = 6) +
  annotate("text", x = 3, y = 0.3, label = "a", size = 6) 
S9b

S9d<-ggplot(data=subset(below_seeding_tr), aes(x=subplot, y=CWM.PC2.m, col = as.factor(subplot)))+
  geom_point(position=position_dodge(0.9), size=4)+
  labs(y=expression(atop("Long roots" %<->% "Short roots","(CWM Belowground PC2 Scores)")))+
  scale_x_discrete(labels = c("Mixed","Forb", "Grass"))+
  scale_color_manual(labels = c("Mixed","Forb", "Grass"),  values = c("#fc8d62", "#66c2a5", "#8da0cb"))+
  geom_errorbar(aes(ymin = CWM.PC2.m-se.PC2, ymax = CWM.PC2.m+se.PC2), size = 0.5, width = 0,position=position_dodge(0.9))+
  theme(axis.title.x=element_blank(),
        text = element_text(size=14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2))+
  annotate("text", x = 1, y = 0.25, label = "a", size = 6) +
  annotate("text", x = 2, y = 0.25, label = "a", size = 6) +
  annotate("text", x = 3, y = 0.25, label = "a", size = 6) 
S9d

figureS9 <- ggarrange(S9a, S9b, S9c, S9d,
                     ncol =2, nrow =2, 
                     align = "v",labels = c("a)", "b)", "c)", "d)"), common.legend = TRUE, legend = "top")
figureS9

#put Fig 5 and fig S9 together
ggarrange(figureS9, figure5, ncol = 1, nrow = 2)
##soil moisture plot
source("./Manifesto/Ashley/ClimVar_Decagon_cleaning.R")
library(lubridate)

dat7 <- dat2 %>%
  filter(year == 2014 | year == 2015) %>%
  mutate(date = ymd_hms(date)) %>%
  mutate(season = ifelse(month%in%c(9:12, 1), "early", "late")) %>%
  mutate(complete_date = (paste(year, month, day, sep="-")),
         complete_date = ymd(complete_date)) %>%
  filter(month%in%c(9:12, 1:5)) %>%
  #mutate(toremove = ifelse(plot == 6 & subplot == "G", 1, 0)) %>%
  filter(date < '2015-05-10 18:00:00',
         subplot != 'C',
         date > '2014-10-10 22:00:00') %>%
  group_by(treatment, subplot, shelterBlock, season, complete_date) %>%
  summarize(sm = mean(sm))

dat8 <- dat7 %>%
  tbl_df() %>%
  #mutate(falltreatment = ifelse(treatment =="consistentDry" | treatment == "fallDry", "dry", "wet")) %>%
  mutate(treatment=ordered(treatment, levels = c( consistentDry="consistentDry", fallDry="fallDry",
                                                  springDry="springDry", controlRain="controlRain"))) %>%
  group_by(treatment, complete_date) %>%
  summarize(sm = mean(sm, na.rm=T)) 

S1<-ggplot(dat8, aes(x=complete_date, y=sm, color = treatment))  + geom_line(lwd = 1) +
  theme_bw() + 
  ylim(0,0.57)+
  scale_x_date(limits = as.Date(c("2014-10-11","2015-05-10")))+
  scale_color_manual(values = c("sienna", "peachpuff", "lightsteelblue1","royalblue2"), 
                     guide = guide_legend(title = "Treatment"),
                     labels=c("Consistent Drought", "Fall Drought", "Spring Drought","Control" )) +
  theme(strip.background = element_blank(), 
        text = element_text(size = 12), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  geom_segment(aes(x = as.Date("2014-10-11"), y = 0.55, xend = as.Date("2015-05-10"), yend = 0.55), color="sienna", size=5.5)+
  geom_segment(aes(x = as.Date("2014-10-11"), y = 0.53, xend = as.Date("2015-02-01"), yend = 0.53), color="peachpuff", size=5)+
  geom_segment(aes(x = as.Date("2015-02-01"), y = 0.53, xend = as.Date("2015-05-10"), yend = 0.53), color="lightsteelblue1", size=5)+
  labs(y = expression(paste("Soil moisture (m"^3, " m"^-3,")")), x = "Day of 2014-2015 growing season") 
#scale_color_manual(values = c( "indianred3", "dodgerblue2"))
#ggsave("Competition/Figures/soil-moisture_201617.pdf", width = 8, height = 5)
S1

####Extra
ggplot(gf_prop_forb, aes(x=percentForb, y=agg_BNPP, color=subplot))+
         geom_point()

ggplot(gf_prop_forb, aes(x=totcover, y=agg_BNPP, color=subplot))+
  geom_point()

m_pc1<-lme(CWM.PC2~shelter*subplot, random=~1|shelterBlock, subset(siteout_fd), na.action=na.exclude)
summary(m_pc1)
anova(m_pc1) #significant
r.squaredGLMM(m_pc1) #24% of variation explained by fixed effects, 42% by whole model (spatial variation?)
qqnorm(residuals(m_pc1))
qqline(residuals(m_pc1))
shapiro.test(residuals(m_pc1))
#normal
LSpc1<-lsmeans(m_pc1, ~subplot, by="shelter")
contrast(LSpc1, "pairwise") #forb is more diverse

       
       