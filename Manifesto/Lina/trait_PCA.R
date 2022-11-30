library(vegan)
library(tidyverse)

#updated 11/30/2022

traits <- read.csv("./Dropbox/ClimVar/DATA/Plant_composition_data/Traits/Traits_GHscreening_Lina_Ashley/Combined_GHtraits.csv", header = TRUE)
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
cover15_fd3 <- cover15_fd2 %>% dplyr::select(ACHMIL, AVEFAT, BRADIS, BRODIA, BROHOR, BROMAD, CENSOL, CYNDAC, CYNECH, EROBOT, HORMUR, LACSER, LOLMUL, TAECAP,TRIHIR,VULMYU)
cover15_fd3<-data.matrix(cover15_fd3)

siteout<- siteout %>% dplyr::select(-ID, -name) 
siteout<-bind_cols(siteout,cover_names)
row.names(siteout)<-siteout$Taxon 
siteout <- siteout %>% dplyr::select(-Taxon)
#siteout <- siteout[-c(14), ] #drop LUPBIC from PCA scores bc it has 0 % cover in any plot 

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

siteout_fd <- siteout_fd %>%
  left_join(ANPP1, by = c("plot", "subplot", "treatment")) 

f_pca1_ANPP <- ggplot(data=siteout_fd, aes(x=CWM.PC1, y=weight_g_m, group=subplot, color=subplot))+
  #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
  geom_point()+
  theme_bw()+
  scale_color_manual(values=c("tomato", "green3", "dodgerblue"), guide = guide_legend(title = "Treatment"), #change legend title
                     labels=c("Mixed", "Forb", "Grass"))+ #change labels in the legend)+
  xlab("Community Weighted Means of PC1 Scores")+
  ylab("ANPP g/m2")+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method="lm", formula= y ~ x, se=FALSE, color="black", aes(group=1))

ggplot(data=siteout_fd, aes(x=CWM.PC1, y=weight_g_m, group=treatment, color=treatment))+
  #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=treatment)) +
  geom_point()+
  theme_bw()+
  #  scale_color_manual(values=c("tomato", "green3", "dodgerblue"), guide = guide_legend(title = "Treatment"), #change legend title
  #                     labels=c("Mixed", "Forb", "Grass"))+ #change labels in the legend)+
  xlab("Community Weighted Means of PC1 Scores")+
  ylab("ANPP g/m2")+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method="lm", formula= y ~ x, se=FALSE,  aes(group=treatment))

f_pca2_ANPP <- ggplot(data=siteout_fd, aes(x=CWM.PC2, y=weight_g_m, group=subplot, color=subplot))+
  #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
  geom_point()+
  theme_bw()+
  xlab("Community Weighted Means of PC2 Scores")+
  ylab("ANPP g/m2")+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method="lm", formula= y ~ x, se=FALSE, color="black", aes(group=1))

ggplot(data=siteout_fd, aes(x=CWM.PC2, y=weight_g_m, group=treatment, color=treatment))+
  #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=treatment)) +
  geom_point()+
  theme_bw()+
  #  scale_color_manual(values=c("tomato", "green3", "dodgerblue"), guide = guide_legend(title = "Treatment"), #change legend title
  #                     labels=c("Mixed", "Forb", "Grass"))+ #change labels in the legend)+
  xlab("Community Weighted Means of PC1 Scores")+
  ylab("ANPP g/m2")+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method="lm", formula= y ~ x, se=FALSE,  aes(group=treatment))

siteout_fd<-merge(siteout_fd, BNPP1)

f_pca1_BNPP <- ggplot(data=siteout_fd, aes(x=CWM.PC1, y=(agg_BNPP), group=subplot, color=subplot))+
  #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
  geom_point()+
  theme_bw()+
  scale_color_manual(values=c("tomato", "green3", "dodgerblue"), guide = guide_legend(title = "Treatment"), #change legend title
                     labels=c("Mixed", "Forb", "Grass"))+ #change labels in the legend)+
  xlab("Community Weighted Means of PC1 Scores")+
  ylab("BNPP g/m2")+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method="lm", formula= y ~ x, se=FALSE, color="black", aes(group=1))
f_pca2_BNPP <- ggplot(data=siteout_fd, aes(x=CWM.PC2, y=(agg_BNPP), group=subplot, color=subplot))+
  #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(group=1)) +
  geom_point()+
  theme_bw()+
  scale_color_manual(values=c("tomato", "green3", "dodgerblue"), guide = guide_legend(title = "Treatment"), #change legend title
                     labels=c("Mixed", "Forb", "Grass"))+ #change labels in the legend)+
  xlab("Community Weighted Means of PC2 Scores")+
  ylab("BNPP g/m2")+
  #xlim(40,100)+
  #ylim(0,100)
  geom_smooth(method="lm", formula= y ~ x, se=FALSE, color="black", aes(group=1))

ggarrange(f_pca1_ANPP, f_pca2_ANPP, f_pca1_BNPP, f_pca2_BNPP, ncol = 2, nrow = 2, common.legend = TRUE)
