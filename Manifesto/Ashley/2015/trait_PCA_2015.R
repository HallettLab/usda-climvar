library(vegan)
library(tidyverse)
library(FD)
library(gridExtra)

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

# First a PCA for ALL TRAITS (ABOVE AND BELOW)
# matrix for PCA
traits <- as.matrix(abovetr15_fd2[,c(6:ncol(all_trait2))])
row.names(traits) <- all_trait2$ID
             
# run PCA
myrda <- rda(na.omit(abovetr15_fd2), scale = TRUE)

# extract values
siteout <- as.data.frame(scores(myrda, choices=c(1,2), display=c("sites")))
siteout$ID<-rownames(siteout)
siteout$name <- siteout$ID

enviroout<-as.data.frame(scores(myrda, choices=c(1,2), display=c("species")))
enviroout$type<-"traits"
enviroout$name<-rownames(enviroout)

# merge PC axes with trait data
tog <- left_join(abovetr15_FD2, siteout) %>%
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
cover15_fd3 <- cover15_fd2 %>% dplyr::select(ACHMIL, AVEFAT, BRADIS, BRODIA, BROHOR, BROMAD, CENSOL, CYNDAC, CYNECH, EROBOT, HORMUR, LACSER, LOLMUL, TAECAP,TRIHIR,VULMYU)
cover15_fd3<-data.matrix(cover15_fd3)

siteout<- siteout %>% dplyr::select(-ID, -name) 
siteout <- siteout[-c(14), ] #drop LUPBIC from PCA scores bc it has 0 % cover in any plot 

siteout_fd<-dbFD (siteout, cover15_fd3, w.abun = T, stand.x = F,
                 calc.FRic = TRUE, m = "max", stand.FRic = FALSE,
                 scale.RaoQ = FALSE, calc.FGR = FALSE, clust.type = "ward",
                 km.inf.gr = 2, km.sup.gr = nrow(x) - 1, km.iter = 100, calc.CWM = TRUE,
                 calc.FDiv = TRUE, print.pco = FALSE, messages = TRUE)
siteout_fd<-as.data.frame(siteout_fd)

traits_2015_2<- arrange(traits_2015, treatment, shelterBlock)
siteout_fd$ANPPgm<-traits_2015_2$ANPPgm
siteout_fd$shelterBlock<-traits_2015_2$shelterBlock
siteout_fd$subplot<-traits_2015_2$subplot
siteout_fd$treatment<-traits_2015_2$treatment

ggplot(data=siteout_fd, aes(x=CWM.PC1, y=ANPPgm, group=subplot, color=subplot))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
  geom_point()+
  theme_bw()+
  scale_color_manual(values=c("tomato", "green3", "dodgerblue"), guide = guide_legend(title = "Treatment"), #change legend title
                    labels=c("Mixed", "Forb", "Grass"))+ #change labels in the legend)+
  xlab("Community Weighted Means of PC1 Scores")+
  ylab("ANPP g/m2")+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method="lm", formula= y ~ x, se=FALSE, color="black", aes(group=1))

ggplot(data=siteout_fd, aes(x=CWM.PC2, y=ANPPgm, group=subplot, color=subplot))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
  geom_point()+
  theme_bw()+
  xlab("Community Weighted Means of PC2 Scores")+
  ylab("ANPP g/m2")+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method="lm", formula= y ~ x, se=FALSE, color="black", aes(group=1))

siteout_fd<-merge(siteout_fd, BNPP1)

ggplot(data=siteout_fd, aes(x=CWM.PC1, y=agg_BNPP, group=subplot, color=subplot))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
  geom_point()+
  theme_bw()+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method="lm", formula= y ~ x, se=FALSE, color="black", aes(group=1))

ggplot(data=siteout_fd, aes(x=CWM.PC2, y=agg_BNPP, group=subplot, color=subplot))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
  geom_point()+
  theme_bw()+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method="lm", formula= y ~ x, se=FALSE, color="black", aes(group=1))

#############################################
##run PCA again but for aboveground traits only

myrda2 <- rda(na.omit(abovetr15_fd2), scale = TRUE)

# extract values
siteout2 <- as.data.frame(scores(myrda2, choices=c(1,2), display=c("sites")))
siteout2$ID<-rownames(siteout2)
siteout2$name <- siteout2$ID

enviroout2<-as.data.frame(scores(myrda2, choices=c(1,2), display=c("species")))
enviroout2$type<-"traits"
enviroout2$name<-rownames(enviroout2)

# merge PC axes with trait data
tog <- left_join(abovetr15_fd2a, siteout2) %>%
  mutate(func = paste(Origin, GF, sep = "_"))

a.pca<-ggplot(tog, aes(x=PC1, y=PC2))+ 
  ggtitle("a) PCA on aboveground traits")+
  geom_hline(aes(yintercept=0), color="grey") + 
  geom_vline(aes(xintercept=0), color="grey") +
  geom_text(aes(label = name, color = GF), size = 5) +
  scale_color_manual(values = c("mediumpurple", "seagreen3","sienna1"), labels=c("Forb","Grass", "Legume"), guide = guide_legend(title = "Functional Group")) +
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
  theme_bw() +theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                    text=element_text(size = 20), legend.position="none")+ 
  xlab(paste("Axis 1 (",sprintf("%.1f",myrda2$CA$eig["PC1"]/myrda2$tot.chi*100,3),"%)",sep="")) +
  ylab(paste("Axis 2 (",sprintf("%.1f",myrda2$CA$eig["PC2"]/myrda2$tot.chi*100,3),"%)",sep="")) 
a.pca

#keep species that we have trait data for
cover15_fd4 <- cover15_fd2 %>% dplyr::select(ACHMIL, AVEBAR, AVEFAT, BRADIS, BRIMIN, BRODIA, BROHOR, BROMAD, CARPYC, CENSOL, CERGLO, CLAPUR, CYNDAC, CYNECH, EROBOT, EROMOS, FILGAL, GALPAR,HORMAR, HORMUR,HYPGLA, HYPRAD, LACSER, LOLMUL, SILGAL, TAECAP,TORARV, TRIDUB,TRIGLO, TRIHIR,TRISP, TRISUB,VICSAT,VULBRO, VULMYU)
cover15_fd4<-data.matrix(cover15_fd4)

#remove excess rows from PCA scores
#siteout2 <- siteout2[-c(15), ]
siteout2<- siteout2 %>% dplyr::select(-ID, -name)

siteout2_fd<-dbFD (siteout2, cover15_fd4, w.abun = T, stand.x = F,
                  calc.FRic = TRUE, m = "max", stand.FRic = FALSE,
                  scale.RaoQ = FALSE, calc.FGR = FALSE, clust.type = "ward",
                  km.inf.gr = 2, km.sup.gr = nrow(x) - 1, km.iter = 100, calc.CWM = TRUE,
                  calc.FDiv = TRUE, print.pco = FALSE, messages = TRUE)
siteout2_fd<-as.data.frame(siteout2_fd)


siteout2_fd$ANPPgm<-traits_2015_2$ANPPgm
siteout2_fd$shelterBlock<-traits_2015_2$shelterBlock
siteout2_fd$subplot<-traits_2015_2$subplot
siteout2_fd$treatment<-traits_2015_2$treatment

a.1<-ggplot(data=siteout2_fd, aes(x=CWM.PC1, y=ANPPgm, group=subplot, color=subplot))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
  #ggtitle("a)")+
  geom_point()+
  theme_bw()+
  theme(legend.position="none")+
  scale_color_manual(values=c("tomato", "green3", "dodgerblue"), guide = guide_legend(title = "Treatment"), #change legend title
                     labels=c("Mixed", "Forb", "Grass"))+ #change labels in the legend)+
  xlab("Community Weighted Means of Aboveground PC1 Scores")+
  ylab("ANPP g/m2")+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method="lm", formula= y ~ x, se=FALSE, color="black", aes(group=1))
a.1

a.2<-ggplot(data=siteout2_fd, aes(x=CWM.PC2, y=ANPPgm, group=subplot, color=subplot))+
  #ggtitle("b)")+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
  geom_point()+
  theme_bw()+
  theme(legend.position="none")+
  scale_color_manual(values=c("tomato", "green3", "dodgerblue"), guide = guide_legend(title = "Treatment"), #change legend title
                     labels=c("Mixed", "Forb", "Grass"))+ #change labels in the legend)+
  xlab("Community Weighted Means of Aboveground PC2 Scores")+
  ylab("")+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method="lm", formula= y ~ x, se=FALSE, color="black", aes(group=1))
a.2

##do again but for BELOWGROUND traits only
# select belowground traits
traits.below<-as.data.frame(traits) %>% dplyr::select(Dens, DiamC, SRLC,SRLF, PropF)

# run PCA
myrda.b <- rda(na.omit(traits.below), scale = TRUE)

# extract values
siteout.b <- as.data.frame(scores(myrda.b, choices=c(1,2), display=c("sites")))
siteout.b$ID<-row.names(siteout.b)
siteout.b$name <- siteout.b$ID

enviroout.b<-as.data.frame(scores(myrda.b, choices=c(1,2), display=c("species")))
enviroout.b$type<-"traits"
enviroout.b$name<-rownames(enviroout.b)

# merge PC axes with trait data
tog <- left_join(trait.dat1, siteout.b) %>%
  mutate(func = paste(Origin, GF, sep = "_"))

# Remove legumes from legend key (if running PCA without legumes)
#tog <- subset(tog, GF != "L")

#pdf("TraitPCA.pdf", width = 9, height = 7.5)
#pdf("TraitPCA_noLegumes.pdf", width = 9, height = 7.5)

b.pca<-ggplot(tog, aes(x=PC1, y=PC2))+ 
  ggtitle("b) PCA on belowground traits")+
  geom_hline(aes(yintercept=0), color="grey") + 
  geom_vline(aes(xintercept=0), color="grey") +
  geom_text(aes(label = name, color = GF), size = 5) +
  scale_color_manual(values = c("mediumpurple","seagreen3", "sienna1"), labels=c("Forb","Grass", "Legume"), guide = guide_legend(title = "Functional Group:")) +
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
siteout.b <- siteout.b[-c(14), ] #remove LUPBIC
siteout.b<- siteout.b %>% dplyr::select(-ID, -name)
#keep species that we have trait data for
cover15_fd3 <- cover15_fd2 %>% dplyr::select(ACHMIL, AVEFAT, BRADIS, BRODIA, BROHOR, BROMAD, CENSOL, CYNDAC, CYNECH, EROBOT, HORMUR, LACSER, LOLMUL, TAECAP,TRIHIR,VULMYU)
cover15_fd3<-data.matrix(cover15_fd3)

siteout_fd.b<-dbFD (siteout.b, cover15_fd3, w.abun = T, stand.x = F,
                  calc.FRic = TRUE, m = "max", stand.FRic = FALSE,
                  scale.RaoQ = FALSE, calc.FGR = FALSE, clust.type = "ward",
                  km.inf.gr = 2, km.sup.gr = nrow(x) - 1, km.iter = 100, calc.CWM = TRUE,
                  calc.FDiv = TRUE, print.pco = FALSE, messages = TRUE)
siteout_fd.b<-as.data.frame(siteout_fd.b)

traits_2015_2<- arrange(traits_2015, treatment, shelterBlock)
siteout_fd.b$ANPPgm<-traits_2015_2$ANPPgm
siteout_fd.b$shelterBlock<-traits_2015_2$shelterBlock
siteout_fd.b$subplot<-traits_2015_2$subplot
siteout_fd.b$treatment<-traits_2015_2$treatment

ggplot(data=siteout_fd.b, aes(x=CWM.PC1, y=ANPPgm, group=subplot, color=subplot))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
  geom_point()+
  theme_bw()+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method="lm", formula= y ~ x, se=FALSE, color="black", aes(group=1))

ggplot(data=siteout_fd.b, aes(x=CWM.PC2, y=ANPPgm, group=subplot, color=subplot))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
  geom_point()+
  theme_bw()+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method="lm", formula= y ~ x, se=FALSE, color="black", aes(group=1))

#run Lina's BNPP_exploratory_analysis.R script to get BNPP1
siteout_fd<-merge(siteout_fd.b, BNPP1)

b.1<-ggplot(data=siteout_fd, aes(x=CWM.PC1, y=agg_BNPP, group=subplot, color=subplot))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
  scale_color_manual(values=c("tomato", "green3", "dodgerblue"), guide = guide_legend(title = "Treatment"), #change legend title
                     labels=c("Mixed", "Forb", "Grass"))+ #change labels in the legend)+
  xlab("Community Weighted Means of Belowground PC1 Scores")+
  ylab("BNPP g/m2")+
    #ggtitle("c)")+
  geom_point()+
  theme_bw()+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method="lm", formula= y ~ x, se=FALSE, color="black", aes(group=1))+
  theme(legend.position="none")
b.1

b.2<-ggplot(data=siteout_fd, aes(x=CWM.PC2, y=agg_BNPP, group=subplot, color=subplot))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
  scale_color_manual(values=c("tomato", "green3", "dodgerblue"), guide = guide_legend(title = "Treatment"), #change legend title
                     labels=c("Mixed", "Forb", "Grass"))+ #change labels in the legend)+
  xlab("Community Weighted Means of Belowground PC2 Scores")+
  ylab("")+
  #ggtitle("d)")+
  geom_point()+
  theme_bw()+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method="lm", formula= y ~ x, se=FALSE, color="black", aes(group=1))+
  theme(legend.position="none")
b.2

b.3<-ggplot(data=siteout_fd, aes(x=CWM.PC2, y=agg_BNPP, group=subplot, color=subplot))+
  #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
  scale_color_manual(values=c("tomato", "green3", "dodgerblue"), guide = guide_legend(title = "Treatment:"), #change legend title
                     labels=c("Mixed", "Forb", "Grass"))+ #change labels in the legend)+
  xlab("Community Weighted Means of Belowground PC2 Scores")+
  ylab("BNPP g/m2")+
  #ggtitle("d)")+
  geom_point()+
  theme_bw()+
  #xlim(40,100)+
  #ylim(0,100)
  #geom_smooth(method="lm", formula= y ~ x, se=FALSE, color="black", aes(group=1))+
  theme(legend.position="bottom")
b.3

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(b.3)

lay <- rbind(c(1,1,2,2),
             c(1,1,2,2),
             c(3,3,4,4),
             c(3,3,4,4),
             c(NA,5,5,NA))

grid.arrange(a.1, a.2, b.1, b.2,mylegend,layout_matrix = lay)

figure4 <- ggarrange(a.1, a.2, b.1, b.2,
                      ncol =2, nrow =2, common.legend = TRUE, legend = "bottom",
                      align = "v",labels = c("a)", "b)", "c)", "d)"))
figure4

ggplot(gf_prop_forb, aes(x=percentForb, y=agg_BNPP, color=subplot))+
         geom_point()

ggplot(gf_prop_forb, aes(x=totcover, y=agg_BNPP, color=subplot))+
  geom_point()

       
       