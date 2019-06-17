library(tidyverse)
library(nlme)
library(ggplot2)
library(dplyr)
library(multcomp)
library(lsmeans)
library(MuMIn)
library(vegan)
library(RColorBrewer)
library(gridExtra)

setwd("~/Dropbox/ClimVar/DATA/")
soil<-read.csv ("~/Dropbox/ClimVar/DATA/Soil Extracts/Soil_CleanedData/ClimVar_soil_2015.csv")
soil_2015<- soil %>% filter (subplot=="F" | subplot =="G" | subplot == "B", Depth == "1", DOY == "168")  %>% 
  arrange(treatment, shelterBlock) #arrange by treatment, shelterBlock so all files used in PCA match

SOC_graphic_2015 <- soil_2015 %>%
  group_by(shelterBlock, treatment) %>%
  summarize(meanSOC=mean(SOC), seSOC=sd(SOC)/sqrt(length(SOC)))

ggplot(SOC_graphic_2015, aes(x=treatment, y=meanSOC, fill=treatment)) + 
  geom_bar(stat="identity", position="dodge") + 
  facet_wrap(~shelterBlock) +
  theme_bw() + 
  geom_errorbar(aes(ymax = meanSOC+seSOC, ymin = meanSOC-seSOC), width=.25, position=position_dodge(width=0.9)) + 
  labs(x="Treatment", y="SOC", fill="Composition Treatment") +
  theme(text = element_text(size=20))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

MBC_graphic_2015 <- soil_2015 %>%
  group_by(shelterBlock, treatment) %>%
  summarize(meanSOC=mean(MBC), seSOC=sd(MBC)/sqrt(length(MBC)))

ggplot(MBC_graphic_2015, aes(x=treatment, y=meanSOC, fill=treatment)) + 
  geom_bar(stat="identity", position="dodge") + 
  facet_wrap(~shelterBlock) +
  theme_bw() + 
  geom_errorbar(aes(ymax = meanSOC+seSOC, ymin = meanSOC-seSOC), width=.25, position=position_dodge(width=0.9)) + 
  labs(x="Treatment", y="MBC", fill="Composition Treatment") +
  theme(text = element_text(size=20))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

NH4_graphic_2015 <- soil_2015 %>%
  group_by(shelterBlock, treatment) %>%
  summarize(meanSOC=mean(NH4), seSOC=sd(NH4)/sqrt(length(NH4)))

ggplot(NH4_graphic_2015, aes(x=treatment, y=meanSOC, fill=treatment)) + 
  geom_bar(stat="identity", position="dodge") + 
  facet_wrap(~shelterBlock) +
  theme_bw() + 
  geom_errorbar(aes(ymax = meanSOC+seSOC, ymin = meanSOC-seSOC), width=.25, position=position_dodge(width=0.9)) + 
  labs(x="Treatment", y="NH4", fill="Composition Treatment") +
  theme(text = element_text(size=20))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

sm2015<-read.csv("~/Dropbox/ClimVar/DATA/Decagon data/ClimVar_sm_2015.csv")
smdat2 <- sm2015 %>%
  # filter(subplot == "B" | subplot == "G" | subplot == "F" ) %>%
  # filter(subplot == "XC" | subplot == "C") %>%
  #  filter(subplot == "XC" | subplot == "B") %>%
  #filter( subplot == "G" | subplot == "F" ) %>%
  filter(subplot !="C", subplot != "XC") %>%
  tbl_df() %>%
  mutate(doy3= julian/365,
         doy4 = year + doy3) %>%
  filter(doy4 < 2015.75)

smdat3<-smdat2%>%
  group_by(treatment, shelterBlock, subplot, month, doy4)%>%
  summarise(meansm=mean(sm, na.rm=T))

smdat_mean1<-smdat2%>%
  group_by(treatment, subplot)%>%
  summarise(avg_sm=mean(sm, na.rm=T))

smdat_mean<-smdat2%>%
  group_by(treatment, shelterBlock, subplot)%>%
  summarise(avg_sm=mean(sm, na.rm=T))

soil_2015<-merge(soil_2015, smdat_mean)
soil_2015<-merge(soil_2015,lit) 
soil_2015<- soil_2015 %>% select(-year)

#create a new variable for growing season?
#smdat4<-smdat2 %>% mutate(season=ifelse(doy4 %in% 2014:2015.085, "fall", ifelse(doy4 %in% 2015.086:2015.33, "spring", ifelse(doy4 %in% 2015.34:2015.75, "summer", "2016"))))

CV <- function(x){(sd(x)/mean(x))*100}
moistCV<-aggregate(sm ~ treatment*shelterBlock*subplot, data= smdat2, FUN = CV)
colnames(moistCV)[colnames(moistCV)=="sm"] <- "sm_cv"
soil_2015<-merge(soil_2015, moistCV, all=T)


#create objects from community analysis
data<- cover %>% dplyr::select(-X, -species, -genus, -status, -func, -func2) %>% spread(species_name, cover)
data<-data %>% dplyr::select(-Unknown)
levels(cover$plot)
str(data)
levels(data$treatment)
levels(data$year)
levels(data$subplot)

data <- tibble::rowid_to_column(data, "subplot")
data2 <- filter(data, subplot=="XC") %>% 
  arrange(treatment, shelterBlock) #arrange to match soil datasheet for PCA

Treatment<-data2[,4]
Year<-data2[,3]
data2$ID <- as.character(seq.int(nrow(data2)))
plotnames<-data2[,64]
cover.Bio<- data2 %>% dplyr::select(-c(1:6), -ID)
rownames(cover.Bio)<-as.character(plotnames)
#check for empty rows
cover.Biodrop<-cover.Bio[rowSums(cover.Bio[, (1:57)]) ==0, ] #no empty rows, next step not needed
# remove empty rows cover.Biodrop<-cover.Bio[rowSums(cover.Bio[, (1:58)])  >0, ]

#look at how many plots each species occurs in
colSums(cover.Bio != 0) #note there are 17 species that occur 0 times, remove those
names(cover.Bio)<-str_replace_all(names(cover.Bio), c(" " = "_"))
cover.Bio2<- cover.Bio %>% dplyr::select(-Achillea_millefolium, -Bromus_sp.,-Bromus_sterilis,-Cynosaurus_echinatus,
                                         - Erodium_cicutarium, -Fillago_gallica, -Gastridium_phleoides, -Geranium_sp., -Juncus_bufonius,
                                         -Kickxia_spuria, -Linum_bienne, -Lupinus_bicolor,-Lythrum_hyssopifolia, -Medicago_arabica,
                                         -Sonchus_oleraceus, -Torilis_arvensis, -Triteleia_hyacintha)

cover.rowsums <- rowSums(cover.Bio2 [1:40])

cover.relrow <- data.frame(cover.Bio2 /cover.rowsums)
cover.colmax<-sapply(cover.Bio2 ,max)
cover.relcolmax <- data.frame(sweep(cover.Bio2 ,2,cover.colmax,'/'))
cover.pa <- cover.Bio2 %>% mutate_each(funs(ifelse(.>0,1,0)), 1:40)
#calculate shannon-weiner
data2$H <- diversity(cover.Bio2)
#calculate pielou's J
data2$J <- May_2015$H/log(specnumber(cover.Bio2))

# run PCA
myrda2 <- rda(cover.Bio2, scale = TRUE)

# extract values
siteout <- as.data.frame(scores(myrda2, choices=c(1,2), display=c("sites")))
siteout$ID<-rownames(siteout)
siteout$name <- siteout$ID

enviroout<-as.data.frame(scores(myrda2, choices=c(1,2), display=c("species")))
enviroout$type<-"species"
enviroout$name<-rownames(enviroout)

env<- soil_2015 %>% dplyr::select(6:11)
envvec<-envfit(myrda2,cover.Bio2, na.rm=TRUE)
envvec
plot(myrda2, display = "sites")
plot(envvec, p.max=0.05) #add vectors to previous ordination
trait.dat2<-trait.dat2 %>% arrange(treatment, shelterBlock)
tr2 <- as.matrix(trait.dat2[,c(15:24)])
envvec2<-envfit(myrda2,tr2, na.rm=TRUE)
envvec2
plot(envvec2, p.max=0.05, col="skyblue")
envvec3<-envfit(myrda2,env, na.rm=TRUE)
envvec3
plot(envvec3, p.max=0.10, col="green")

#plots colored based on treatment
plot <- ordiplot(myrda2,choices=c(1,2), type = "none")   #Set up the plot
colsT2 <- rep(c("Red","Blue","Orange", "Pink"), each = 12) #color based on drought treatment
cols1 <- rep(c("Red","Blue","Orange", "Pink"))
shapes <- rep(c(15, 8, 17 ), each=1) #shapes on year
points(siteout$PC1,siteout$PC2,col=colsT2,pch=shapes) 
text(enviroout, display = "species", cex=0.5, col="grey40") #label species
# add legend for treatment
legend("topleft",legend=levels(data2$treatment), col=cols1, pch=19, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)
# add legend for year
legend("topright",legend=levels(data2$subplot), col="black", pch=shapes, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)


# merge PC axes with trait data
tog <- left_join(data2, siteout) 

ggplot(tog, aes(x=PC1, y=PC2))+ 
  geom_hline(aes(yintercept=0), color="grey") + 
  geom_vline(aes(xintercept=0), color="grey") +
  geom_point(aes(shape = shelterBlock, color = treatment), size = 3) +
  # scale_color_manual(values = c("grey20", "grey70")) +
  #geom_segment(data = envvec,
  #aes(x = 0, xend =  PC1,
  #y = 0, yend =  PC2),
  #arrow = arrow(length = unit(0.25, "cm")), colour = "black") + #grid is required for arrow to work.
  geom_text(data = enviroout,
            aes(x=  PC1*1.1, y =  PC2*1.1, #we add 10% to the text to push it slightly out from arrows
                label = name), #otherwise you could use hjust and vjust. I prefer this option
            size = 2,
            hjust = 0.5, 
            color="black") + 
  theme_bw() +theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                    text=element_text(size = 15))+ 
  xlab(paste("Axis 1 (",sprintf("%.1f",myrda$CA$eig["PC1"]/myrda$tot.chi*100,3),"%)",sep="")) +
  ylab(paste("Axis 2 (",sprintf("%.1f",myrda$CA$eig["PC2"]/myrda$tot.chi*100,3),"%)",sep=""))


#traits PCA with correlated species
# run PCA
trrda <- rda(tr2, scale = TRUE)

# extract values
siteout <- as.data.frame(scores(trrda, choices=c(1,2), display=c("sites")))
siteout$ID<-rownames(siteout)
siteout$name <- siteout$ID

enviroout<-as.data.frame(scores(trrda, choices=c(1,2), display=c("species")))
enviroout$type<-"traits"
enviroout$name<-rownames(enviroout)

soil_2015 <- soil_2015 %>% arrange(treatment, shelterBlock)
env<- soil_2015 %>% dplyr::select(-X) %>% dplyr::select(8:20) 
envvec<-envfit(trrda,cover.Bio2, na.rm=TRUE)
envvec
sp<-as.factor(data2$subplot)
trt<-as.factor(data2$treatment)
plot(trrda, display="species")
plot(envvec, p.max=0.02) #add vectors to previous ordination
envvec3<-envfit(trrda,env, na.rm=TRUE)
envvec3
plot(envvec3, p.max=0.10, col="green3")

#plots colored based on treatment
plot <- ordiplot(trrda,choices=c(1,2), type = "none")   #Set up the plot
colsT2 <- rep(c("sienna","royalblue2","peachpuff2", "lightsteelblue1"), each = 12) #color based on drought treatment
cols1 <- rep(c("sienna","royalblue2","peachpuff2", "lightsteelblue1"))
shapes <- rep(c(15, 8, 21,17 ), each=3) #shapes on shelterBlock
shapes1 <- rep(c(15,8,21, 17), each = 1)
points(siteout$PC1,siteout$PC2,col=colsT2,pch=shapes) 
text(trrda, display = "species", cex=0.9, col="grey40") #label species
plot(envvec3, p.max=0.10, col="red", cex=0.9)
plot(envvec, p.max=0.02, col= "seagreen", cex=0.6)
# add legend for treatment
legend("topleft",legend=levels(data2_B$treatment), col=cols1, pch=19, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)
# add legend for year
legend("bottomleft",legend=levels(data2$shelterBlock), col="black", pch=shapes1, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)
help(ordiplot)

#note: this section of code from http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/
library("FactoMineR")
library(factoextra)
library("corrplot")
trPCA<- PCA(tr, scale.unit = TRUE, ncp = 4, graph = TRUE)
eig.val <- get_eigenvalue(trPCA)
eig.val #note that 95% of variance is explained by the 1st 4 dimensions
fviz_eig(trPCA, addlabels = TRUE, ylim = c(0, 50)) #visualize variance explained in scree plot
var <- get_pca_var(trPCA)
var
head(var$coord)
corrplot(var$cos2, is.corr=FALSE)
fviz_contrib(trPCA, choice = "var", axes = 1, top = 10)
fviz_contrib(trPCA, choice = "var", axes = 2, top = 10)
fviz_contrib(trPCA, choice = "var", axes = 1:2, top = 10)
res.km <- kmeans(var$coord, centers = 4, nstart = 25)
grp <- as.factor(res.km$cluster)
# Color variables by groups
fviz_pca_var(trPCA, col.var = grp, 
             palette = c("#0073C2FF", "#EFC000FF", "#868686FF", "red"),
             legend.title = "Cluster")
library(devtools)
library(ggbiplot)
tr.pca<-prcomp(tr,scale=TRUE)
spp.scrs <- as.data.frame(scores(envvec3, display = "vectors"))
spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs))
g <- ggbiplot(tr.pca, obs.scale = 1, var.scale = 1, 
              groups = sp, ellipse = F, circle = TRUE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
g
g2 <- ggbiplot(tr.pca, obs.scale = 1, var.scale = 1, alpha=0)
g2+geom_point(aes(shape=factor(sp), color=factor(trt)))+
  geom_segment(data = spp.scrs,
               aes(x = 0, xend = PC1, y = 0, yend = PC2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
  geom_text(data = spp.scrs, aes(x = PC1, y = PC2, label = Species),
            size = 3)
g2


bioPCA<- PCA(cover.relrow, scale.unit = TRUE, ncp = 10, graph = TRUE)
eig.val <- get_eigenvalue(bioPCA)
eig.val #note that 95% of variance is explained by the 1st 4 dimensions
fviz_eig(bioPCA, addlabels = TRUE, ylim = c(0, 50)) #visualize variance explained in scree plot
var <- get_pca_var(bioPCA)
var
head(var$coord)
corrplot(var$cos2, is.corr=FALSE)
fviz_contrib(bioPCA, choice = "var", axes = 1, top = 10)
fviz_contrib(bioPCA, choice = "var", axes = 2, top = 10)
fviz_contrib(bioPCA, choice = "var", axes = 1:2, top = 10)
res.km <- kmeans(var$coord, centers = 4, nstart = 25)
grp <- as.factor(res.km$cluster)
# Color variables by groups
fviz_pca_var(bioPCA, col.var = grp, 
             palette = c("#0073C2FF", "#EFC000FF", "#868686FF", "red"),
             legend.title = "Cluster")

##Let's do this again but just for FORB treatment
data2_F <- filter(data, year =="2015", subplot=="F") %>% arrange(treatment, shelterBlock)

Treatment<-data2_F[,4]
Year<-data2_F[,3]
data2_F$ID <- as.character(seq.int(nrow(data2_F)))
plotnames<-data2_F[,64]
cover.Bio_F<- data2_F %>% dplyr::select(-c(1:6), -ID)
rownames(cover.Bio_F)<-as.character(plotnames)
#check for empty rows
cover.Biodrop<-cover.Bio_F[rowSums(cover.Bio_F[, (1:57)]) ==0, ] #no empty rows, next step not needed
# remove empty rows cover.Biodrop<-cover.Bio[rowSums(cover.Bio[, (1:58)])  >0, ]

#look at how many plots each species occurs in
#colSums(cover.Bio_B != 0) #note there are lots of species that occur 1 or 0 times, remove those
#names(cover.Bio)<-str_replace_all(names(cover.Bio), c(" " = "_"))
#cover.Bio2<- cover.Bio %>% dplyr::select(-Achillea_millefolium, -Brachypodium_distachyon, -Bromus_sterilis, -Clarkia_amoena, -Erodium moschatum, -Triteleia_hyacintha, -Trifolium_wildenovii,-Sonchus_oleraceus,-Sherardia_arvensis,
#-Fillago gallica, -Galium parisiense, -Gastridium phleoides,-Kickxia spuria,-Linum bienne, -Lythrum hyssopifolia, -Medicago polymorpha, -Torilis arvensis, -Trifolium dubium, -Trifolium glomeratum,
#-Trifolium sp.,-Senecio_vulgaris, -Rumex_pulcher, -Medicago_arabica, -Lupinus_bicolor, -Juncus_bufonius,
#-Hordeum_sp., -Geranium_sp., -Gastridium_phleoides,- Erodium_cicutarium, -Cynosaurus_echinatus,
#-Convolvulus_arvensis, -Carduus_pycnocephalus, -Bromus_sp.)

cover.rowsums <- rowSums(cover.Bio_F [1:57])

cover.relrow <- data.frame(cover.Bio_F /cover.rowsums)
cover.colmax<-sapply(cover.Bio_F ,max)
cover.relcolmax <- data.frame(sweep(cover.Bio_F ,2,cover.colmax,'/'))
cover.pa <- cover.Bio_F %>% mutate_each(funs(ifelse(.>0,1,0)), 1:57)
#calculate shannon-weiner
data2_F$H <- diversity(cover.Bio_F)
#calculate pielou's J
data2_F$J <- May_2015$H/log(specnumber(cover.Bio_F))

trait.dat2_F<- filter(trait.dat2, subplot=="F")
trF <- as.matrix(trait.dat2_F[,c(14:21)])
row.names(trF) <- trait.dat2_F$ID

# run PCA
TF.rda <- rda(trF, scale = TRUE)

# extract values
siteout <- as.data.frame(scores(TF.rda, choices=c(1,2), display=c("sites")))
siteout$ID<-rownames(siteout)
siteout$name <- siteout$ID

enviroout<-as.data.frame(scores(TF.rda, choices=c(1,2), display=c("species")))
enviroout$type<-"species"
enviroout$name<-rownames(enviroout)

env<- soil_2015 %>% filter(subplot=="F")%>% dplyr::select(6:11) 
envvec<-envfit(TF.rda,cover.Bio_F, na.rm=TRUE)
envvec
plot(TF.rda)
plot(envvec, p.max=0.05) #add vectors to previous ordination
envvec3<-envfit(TF.rda,env, na.rm=TRUE)
envvec3
plot(envvec3, p.max=0.10, col="green")

#plots colored based on treatment
plot <- ordiplot(myrda3,choices=c(1,2), type = "none")   #Set up the plot
colsT2 <- rep(c("Red","Blue","Orange", "Pink"), each = 4) #color based on drought treatment
cols1 <- rep(c("Red","Blue","Orange", "Pink"))
shapes <- rep(c(15, 8, 17, 22 ), each=1) #shapes on shelterBlock
points(siteout$PC1,siteout$PC2,col=colsT2,pch=shapes) 
text(enviroout, display = "species", cex=0.5, col="grey40") #label species
# add legend for treatment
legend("topleft",legend=levels(data2_B$treatment), col=cols1, pch=19, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)
# add legend for year
legend("topright",legend=levels(data2_B$shelterBlock), col="black", pch=shapes, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)
help(ordiplot)

B.PCA<- PCA(trB, scale.unit = TRUE, ncp = 10, graph = TRUE)
eig.val <- get_eigenvalue(B.PCA)
eig.val #note that 95% of variance is explained by the 1st 4 dimensions
fviz_eig(B.PCA, addlabels = TRUE, ylim = c(0, 50)) #visualize variance explained in scree plot
var <- get_pca_var(B.PCA)
var
head(var$coord)
corrplot(var$cos2, is.corr=FALSE)
fviz_contrib(B.PCA, choice = "var", axes = 1, top = 10)
fviz_contrib(B.PCA, choice = "var", axes = 2, top = 10)
fviz_contrib(B.PCA, choice = "var", axes = 1:2, top = 10)
res.km <- kmeans(var$coord, centers = 3, nstart = 25)
grp <- as.factor(res.km$cluster)
# Color variables by groups
fviz_pca_var(B.PCA, col.var = grp, 
             palette = c("#0073C2FF", "#EFC000FF", "#868686FF"),
             legend.title = "Cluster")

data3<- data2 %>% dplyr::select (-ID, -H, -J)
write.csv(data3, "~/Dropbox/ClimVar/DATA/Plant_composition_data/Cover/Cover_CleanedData/ClimVar_species-cover_wide.csv")
