library(vegan)
library(tidyverse)
setwd("./Dropbox/ClimVar/DATA/Plant_composition_data/Traits/Traits_GHscreening_Brad")
trait.dat <- read.csv("BradTraits_Cleaned.csv", header = TRUE)
trait.dat2 <- trait.dat %>%
  select(-MD, -PropF, - Total, - actual_area, -RMF, - RGR)

# matrix for PCA
traits <- as.matrix(trait.dat2[,c(6:ncol(trait.dat2))])
row.names(traits) <- trait.dat2$ID

# run PCA
myrda <- rda(traits, scale = TRUE)

# extract values
siteout <- as.data.frame(scores(myrda, choices=c(1,2), display=c("sites")))
siteout$ID<-rownames(siteout)
siteout$name <- siteout$ID

enviroout<-as.data.frame(scores(myrda, choices=c(1,2), display=c("species")))
enviroout$type<-"traits"
enviroout$name<-rownames(enviroout)

# merge PC axes with trait data
tog <- left_join(trait.dat, siteout) #%>%
  #mutate(func = paste(Origin, GF, sep = "_"))

pdf("TraitPCA_by_FG.pdf", width = 9, height = 7.5)

ggplot(tog, aes(x=PC1, y=PC2))+ 
  geom_hline(aes(yintercept=0), color="grey") + 
  geom_vline(aes(xintercept=0), color="grey") +
  geom_text(aes(label = name, color = GF), size = 5) +
  # scale_color_manual(values = c("grey20", "grey70")) +
  geom_segment(data = enviroout,
               aes(x = 0, xend =  PC1,
                   y = 0, yend =  PC2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "black") + #grid is required for arrow to work.
  geom_text(data = enviroout,
            aes(x=  PC1*1.2, y =  PC2*1.2, #we add 10% to the text to push it slightly out from arrows
                label = name), #otherwise you could use hjust and vjust. I prefer this option
            size = 6,
            hjust = 0.5, 
            color="black") + 
  theme_bw() +theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                    text=element_text(size = 20))+ 
  xlab(paste("Axis 1 (",sprintf("%.1f",myrda$CA$eig["PC1"]/myrda$tot.chi*100,3),"%)",sep="")) +
  ylab(paste("Axis 2 (",sprintf("%.1f",myrda$CA$eig["PC2"]/myrda$tot.chi*100,3),"%)",sep="")) 

dev.off()
