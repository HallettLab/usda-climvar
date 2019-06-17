library(vegan)
library(MASS)
library(tidyverse)
library(dplyr)
library (vegan3d)
library(RVAideMemoire) #for posthoc tests on permanova
# Import csv file, transform to wide data, call it data
setwd("~/Dropbox/ClimVar/DATA/Plant_composition_data")
cover<-read.csv("Cover/Cover_CleanedData/ClimVar_species-cover.csv")
data<- cover %>% dplyr::select(-X, -species, -genus, -status, -func, -func2) %>% spread(species_name, cover)
data<-data %>% dplyr::select(-Unknown)
levels(cover$plot)
str(data)
levels(data$treatment)
levels(data$year)
levels(data$subplot)

data <- tibble::rowid_to_column(data, "subplot")

######################
#0. Create objects needed for ordination
#This first round is for full dataset (all subplots)... 
#I will do XC only next
#And then composition treatments
######################

###some code to load and manipulate env variables (for later)
#LTM_env<-data %>% select(-c(18:51))
#plotnames<-LTM_env[,1]
#LTM.env<-LTM_env
#rownames(LTM.env)<-plotnames 
#Treatment=as.factor(LTM_env[,9])
#Year=as.factor(LTM_env[,5])
#LTM.env2 <-LTM.env[,-c(1:16)]
###standardize the environmental variables (scale function converts to z-scores)
#LTM.env.z <- data.frame(scale(LTM.env2))
soil<-read.csv ("~/Dropbox/ClimVar/DATA/Soil Extracts/Soil_CleanedData/ClimVar_soil_2015.csv")
soil_all<- soil %>% filter (subplot=="XC" | subplot =="G" | subplot == "B" | subplot =="F", Depth == "1")  %>% 
  arrange(treatment, shelterBlock, subplot) #arrange so all files used in PCA match

Treatment<-data[,4]
Year<-data[,3]
data$ID <- seq.int(nrow(data))
plotnames<-data[,64]
cover.Bio<- data %>% dplyr::select(-c(1:6), -ID)
rownames(cover.Bio)<-plotnames
#check for empty rows
cover.Biodrop<-cover.Bio[rowSums(cover.Bio[, (1:57)]) ==0, ] #no empty rows, next step not needed
# remove empty rows cover.Biodrop<-cover.Bio[rowSums(cover.Bio[, (1:58)])  >0, ]

cover.rowsums <- rowSums(cover.Bio [1:57])

cover.relrow <- data.frame(cover.Bio /cover.rowsums)
cover.colmax<-sapply(cover.Bio ,max)
cover.relcolmax <- data.frame(sweep(cover.Bio ,2,cover.colmax,'/'))
cover.pa <- cover.Bio %>% mutate_each(funs(ifelse(.>0,1,0)), 1:57)

######################
#2. NMS
#see also section 2.1 of vegan tutorial: 
#http://cc.oulu.fi/~jarioksa/opetus/metodi/vegantutor.pdf
######################

#make bray-curtis dissimilarity matrix
spp.bcd <- vegdist(cover.relrow)
#or
spp.pa.bcd<-vegdist(cover.pa, binary=T) #this calcs distances based on presence/absence, just in case

#run NMS ordination
spp.mds0 <-isoMDS(spp.bcd) #runs nms only once
spp.mds0  #by default 2 dimensions returned, stress is 18.995
ordiplot(spp.mds0)

#prefer to run multiple NMS ordinations
#with different starting configurations, and choose the best
#this function does that, and in addition does several other steps too
#including: 
  #standardizing the data (though fuction call below turns this off with autotransform=F)
  #calculating distance matrix (default bray-curtis)
  #running NMDS with random starts
  #rotation of axes to maximize variance of site scores on axis 1
  #calculate species scores based on weighted averaging
  
help(metaMDS)
spp.mds<-metaMDS(cover.relrow, trace = FALSE, autotransform=T, trymax=100, k=4) #runs several with different starting configurations
#trace= TRUE will give output for step by step what its doing
#default is 2 dimensions, can put k=4 for 4 dimensions
spp.mds #solution converged after 20 tries, stress = 9.04
summary(spp.mds)

#plot results
stressplot(spp.mds, spp.bcd) #stressplot to show fit, fit is decent
ordiplot(spp.mds)
spscores1<-scores(spp.mds,display="sites",choices=1)
spscores2<-scores(spp.mds,display="sites",choices=2)
tplots<-data[,2]
tplot_levels<-levels(tplots)
spscoresall<-data.frame(tplots,spscores1,spscores2)

#plots colored based on treatment
bio.plot <- ordiplot(spp.mds,choices=c(1,2), type = "none")   #Set up the plot
colsT1 <- rep(c("Red","Orange","purple", "black"), each = 15) #color based on drought treatment
Lcols <- rep(c("Red", "Orange", "purple", "black"))
shapes <- rep(c(15, 3, 17, 19, 5), each=3) #shapes on subplot
Lshapes <- rep(c(15,3,17,19,5))
points(spscoresall$NMDS1,spscoresall$NMDS2,col=colsT1,pch=shapes) 
text(spp.mds, display = "species", cex=0.5, col="grey30") #label species
legend("bottomright",legend=levels(Treatment), col=Lcols, pch=15, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)
# add legend for treatment
legend("topright",legend=levels(data$subplot), col="black", pch=Lshapes, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)
#help(ordiplot)

#plots colored based on year
bio.plot2 <- ordiplot(spp.mds,choices=c(1,2), type = "none")   #Set up the plot
colyr1 <- rep(c("black","cyan","orange"), each = 1) 
points(spscoresall$NMDS1,spscoresall$NMDS2,col=colyr1,pch=20) 
#ordiellipse(spp.mds, groups=Year, fill=colyr)
text(spp.mds, display = "species", cex=0.6, col="grey30") #label species
legend("bottomright",legend=levels(as.factor(Year)), col=colyr1, pch=15, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)
# add legend for year

#plots colored based on block
bio.plot3 <- ordiplot(spp.mds,choices=c(1,2), type = "none")   #Set up the plot
colb1 <- rep(c("grey", "magenta", "yellow", "navy"), each = 60) 
colb1L<- rep(c("grey", "magenta", "yellow", "navy"))
points(spscoresall$NMDS1,spscoresall$NMDS2,col=colb1,pch=20) 
text(spp.mds, display = "species", cex=0.4, col="black") #label species 
legend("bottomright",legend=levels(as.factor(data$shelterBlock)), col=colb1L, pch=15, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)
# add legend for year
#seems to be different based on block?

#plots with plant species colored on functional group: red= forb, green = grass, blue = medusahead, and pink=n-fixer
colspec<- rep(c("red", "red", "green", "green", "green", "green", "green", "green", "green", "green", "green", "red", "red", "red", "red","red", "green", "green", "red", "red", "red", "red", "red", "red", "red", "green", "green", "green", "red", "red", "green", "red", "red", "red", "green", "red", "red", "red", "red", "red", "red", "red", "red", "red", "blue", "red", "pink", "pink", "pink", "pink", "pink", "pink", "red", "pink", "green","green", "red"))
bio.plot4 <- ordiplot(spp.mds,choices=c(1,2), type = "none")   #Set up the plot
points(spscoresall$NMDS1,spscoresall$NMDS2,col="grey30",pch=20) 
text(spp.mds, display = "species", cex=0.8, col=colspec) #label species 

permanova1 <- adonis(spp.bcd~data$shelterBlock, perm=100, method="bray")
permanova1

permanova2 <- adonis(spp.bcd~Year, perm=100, method="bray")
permanova2

permanova3 <- adonis(spp.bcd~Treatment, perm=999, method="bray")
permanova3

permanova4 <- adonis(spp.bcd~Treatment*Year*data$subplot, perm=999, method="bray")
permanova4 #everything is significant, need model selection


#IF NEEDED LATER: code here is to overlay environmental data
#first remove rows that were removed from bio (if any)
#LTM.env3<-LTM.env2[ -c(642, 683, 691, 695), ]
#envvec.nms<-envfit(spp.mds,LTM.env3, na.rm=TRUE)
#envvec.nms
#plot(envvec.nms) #add vectors to previous ordination

#if desired, rotate ordination so that first axis is parallel to an environmental variable
#spp.mds.rotelev<-MDSrotate(spp.mds,LTM.env3$moisture)
#envvec.nms.rotelev<-envfit(spp.mds.rotelev,LTM.env3, na.rm=T)
#ordiplot(spp.mds.rotelev)
#text(spp.mds.rotelev, display = "species", cex=0.7, col="red") 
#plot(envvec.nms.rotelev) 

##Make a 3 dimensional plot of this ordination 
#color based on treatment
ordirgl (spp.mds, col=colsT1)
#color based on year
ordirgl (spp.mds, col=colyr1)
#color based on block
ordirgl (spp.mds, col=colb1)

####################
#repeat ordination for control (XC) subplots only
##################

data2<- cover %>% dplyr::select(-X, -species, -genus, -status, -func, -func2) %>% 
  filter(subplot=="XC") %>% spread(species_name, cover) %>% arrange(treatment, shelterBlock) 
data2$ID <- seq.int(nrow(data2))
data2$TB<-paste(data2$treatment,data2$shelterBlock,sep=".")
data2<-data2 %>% dplyr::select(-Unknown)

###some code to load and manipulate env variables (for later)
#LTM_env<-data %>% select(-c(18:51))
#plotnames<-LTM_env[,1]
#LTM.env<-LTM_env
#rownames(LTM.env)<-plotnames 
#Treatment=as.factor(LTM_env[,9])
#Year=as.factor(LTM_env[,5])
#LTM.env2 <-LTM.env[,-c(1:16)]
###standardize the environmental variables (scale function converts to z-scores)
#LTM.env.z <- data.frame(scale(LTM.env2))

TB<-as.factor(data2[,65])
Treatment2<-data2[,4]
Year2<-data2[,3]
year2<-as.factor(Year)
plotnames<-data2[,64]
cover.Bio2<- data2 %>% dplyr::select(-c(1:6), -ID, -TB)
rownames(cover.Bio2)<-plotnames

#check for empty rows
cover.Biodrop2<-cover.Bio2[rowSums(cover.Bio2[, (1:57)]) ==0, ] #no empty rows, next step not needed
# remove empty rows: cover.Biodrop<-cover.Bio[rowSums(cover.Bio[, (1:58)])  >0, ]

cover.rowsums2 <- rowSums(cover.Bio2 [1:57])
cover.relrow2 <- data.frame(cover.Bio2 /cover.rowsums2)
cover.pa2 <- cover.Bio2 %>% mutate_each(funs(ifelse(.>0,1,0)), 1:57)

#make bray-curtis dissimilarity matrix for controls
spp.bcd2 <- vegdist(cover.relrow2)
#or
spp.pa.bcd2<-vegdist(cover.pa2, binary=T) #this calcs distances based on presence/absence

#run NMS ordination
spp.mds0_2 <-isoMDS(spp.bcd2) #runs nms only once
spp.mds0_2  #by default 2 dimensions returned, stress 13.86
ordiplot(spp.mds0_2)

#prefer to run multiple NMS ordinations
#with different starting configurations, and choose the best
#this function does that, and in addition does several other steps too
#including: 
#standardizing the data (though fuction call below turns this off with autotransform=F)
#calculating distance matrix (default bray-curtis)
#running NMDS with random starts
#rotation of axes to maximize variance of site scores on axis 1
#calculate species scores based on weighted averaging

spp.mds2<-metaMDS(cover.relrow2, trace = T, autotransform=F, trymax=999, k=4) #runs several with different starting configurations

#trace= TRUE will give output for step by step what its doing
#default is 2 dimensions, can put k=4 for 4 dimensions
spp.mds2 #solution converged after 29 tries, stress is 8.70
summary(spp.mds2)

#plot results
stressplot(spp.mds2, spp.bcd2) #stressplot
ordiplot(spp.mds2)
spscores1_2<-scores(spp.mds2,display="sites",choices=1)
spscores2_2<-scores(spp.mds2,display="sites",choices=2)
tplots2<-data2[,65]
tplots2<-as.factor(tplots2)
tplot_levels2<-as.factor(levels(tplots2))
spscoresall_2<-data.frame(tplots2,spscores1_2,spscores2_2)

library(goeveg)
## Select the 30% most frequent species with 50% best axis fit
limited <- ordiselect(cover.relrow2, spp.mds, ablim = 0.3, fitlim=0.5, method="axes")

#plots colored based on treatment
xc.plot <- ordiplot(spp.mds2,choices=c(1,2), type = "none")   #Set up the plot
colsT2 <- rep(c("Red","Blue","Orange", "Pink"), each = 12) #color based on drought treatment
cols1 <- rep(c("Red","Blue","Orange", "Pink"))
shapes <- rep(c(15, 8, 17 ), each=1) #shapes on year
points(spscoresall_2$NMDS1,spscoresall_2$NMDS2,col=colsT2,pch=shapes) 
text(spp.mds2, display = "species", cex=0.5, col=colspec) #label species
# add legend for treatment
legend("bottomright",legend=levels(Treatment2), col=cols1, pch=19, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)
# add legend for year
legend("topright",legend=levels(as.factor(as.character(Year2))), col="black", pch=shapes, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)

#plots colored based on year
xc.plot2 <- ordiplot(spp.mds2,choices=c(1,2), type = "none")   #Set up the plot
cols <- rep(c("black","cyan","orange"), each = 1) 
points(spscoresall_2$NMDS1,spscoresall_2$NMDS2,col=cols,pch=20, cex=2) 
#ordiellipse(spp.mds, groups=Year, fill=cols)
text(spp.mds2, display = "species", cex=0.6, col=colspec) #label species 
legend("topright",legend=levels(as.factor(as.character(Year2))), col=cols, pch=19, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)


#plots colored based on block
xc.plot3 <- ordiplot(spp.mds2,choices=c(1,2), type = "none")   #Set up the plot
cols <- rep(c("grey", "magenta", "yellow", "navy"), each = 3) 
cols1 <- c("grey", "magenta", "yellow", "navy") 
points(spscoresall_2$NMDS1,spscoresall_2$NMDS2,col=cols,pch=20, cex=2) 
text(spp.mds2, display = "species", cex=0.6, col=colspec) #label species 
legend("bottomright",legend=levels(data2$shelterBlock), col=cols1, pch=19, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)

permanova1_2 <- adonis(spp.bcd2~tplots2, perm=100, method="bray")
permanova1_2
pairwise.perm.manova(spp.bcd2, tplots2, nperm=999)

permanova2_2 <- adonis(spp.bcd2~Year2, perm=100, method="bray")
permanova2_2 #Year not significant

permanova3_2 <- adonis(spp.bcd2~Treatment2, perm=999, method="bray")
permanova3_2 #not significant

permanova4_2 <- adonis(spp.bcd2~Year2*Treatment2*data2$shelterBlock, perm=999, method="bray")
permanova4_2 
pairwise.perm.manova(spp.bcd2, Treatment2, nperm=999)

permanova5_2 <- adonis(spp.bcd2~data2$shelterBlock, perm=999, method="bray")
permanova5_2 #block is significant
#which blocks differ?
pairwise.perm.manova(spp.bcd2, data2$shelterBlock, nperm=999)
#A is different from all other blocks, B different from D, C & D are not different

##Creating an ordination plot with succession vectors
xc.plot4 <- ordiplot(spp.mds2,choices=c(1,2), type = "none", xlim=c(-1.5,1.5), ylim=c(-1.5,1.5))   #Set up the plot
cols <- rep(c("sienna","royalblue2","lightsteelblue3", "peachpuff2"), each = 12) 
cols1 <- rep(c("sienna","royalblue2","lightsteelblue3", "peachpuff2"), each = 1) 
cols2 <- rep(c("sienna","royalblue2","lightsteelblue3", "peachpuff2"), each = 4) 
shapes <- rep(c(15, 3, 17, 19), each=3) #shapes on block
shapes1 <- rep(c(15, 3, 17, 19), each=1) #shapes on block
points(spscoresall_2$NMDS1,spscoresall_2$NMDS2,col=cols,pch=shapes, cex=1.2)#Plot the ordination points 
text(spp.mds2, display = "species", cex=0.4, col="grey30") #label species
ordiarrows(spp.mds2, groups=TB, order.by=Year2, label=F, col=cols2)
legend("topright",legend=levels(Treatment2), col=cols1, pch=19, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)
legend("bottomright",legend=levels(data2$shelterBlock), col="black", pch=shapes1, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)

#to color grasses and forbs labels:
colspec<- rep(c("plum1", "plum1", "darkgreen", "palegreen", "palegreen", "palegreen", "palegreen", "palegreen", "palegreen", "palegreen", "palegreen", "plum1", "plum1", "plum1", "plum1","plum1", "palegreen", "palegreen", "plum1", "plum1", "plum1", "plum1", "plum1", "plum1", "plum1", "palegreen", "palegreen", "palegreen", "plum1", "plum1", "palegreen", "plum1", "plum1", "plum1", "palegreen", "plum1", "plum1", "plum1", "plum1", "plum1", "plum1", "plum1", "plum1", "plum1", "palegreen", "plum1", "plum4", "plum4", "plum4", "plum4", "plum4", "plum4", "plum1", "plum4", "palegreen","palegreen", "plum1"))
##Creating an ordination plot with succession vectors
xc.plot4 <- ordiplot(spp.mds2,choices=c(1,2), type = "none", xlim=c(-1.5,1.5), ylim=c(-1.5,1.5))   #Set up the plot
cols <- rep(c("sienna","royalblue2","lightsteelblue3", "peachpuff2"), each = 12) 
cols1 <- rep(c("sienna","royalblue2","lightsteelblue3", "peachpuff2"), each = 1) 
cols2 <- rep(c("sienna","royalblue2","lightsteelblue3", "peachpuff2"), each = 4) 
shapes <- rep(c(15, 3, 17, 19), each=3) #shapes on block
shapes1 <- rep(c(15, 3, 17, 19), each=1) #shapes on block
points(spscoresall_2$NMDS1,spscoresall_2$NMDS2,col=cols,pch=shapes, cex=1.2)#Plot the ordination points 
text(spp.mds2, display = "species", cex=0.7, col=colspec) #label species
ordiarrows(spp.mds2, groups=TB, order.by=Year2, label=F, col=cols2)
legend("topright",legend=levels(Treatment2), col=cols1, pch=19, cex=0.9,inset=0.01,bty="n",y.intersp=0.5,x.intersp=0.4,pt.cex=1.1)
legend("bottomright",legend=levels(data2$shelterBlock), col="black", pch=shapes1, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)


##Make a 3 dimensional plot of this ordination
##color based on treatment
ordirgl (spp.mds2, col=colsT2, pch=shapes)
help(vegan3d)

####
#test vector length
###

spscoresall_2a<-cbind(spscoresall_2, Year2)
spscoresall_2a<-spscoresall_2a %>% filter(Year2 != "2016") %>% dplyr::select(-Year2) #so endpoints only in dataset
tplot_levels2a<-levels(spscoresall_2a$tplots2)
options(scipen=999)
veclengths<-c()
vecindices<-c()

for (i in 1:length(tplot_levels2a)){
  index<-spscoresall_2a$tplots2[i]
  tempscores<-spscoresall_2a[spscoresall_2a$tplots2==index,]
  if(nrow(tempscores) <= 1) next
  xshift<- tempscores$NMDS1[2]-tempscores$NMDS1[1]
  yshift<- tempscores$NMDS2[2]-tempscores$NMDS2[1]
  veclength<- (xshift^2+yshift^2)^0.5
  veclengths<-c(veclengths, veclength)
  vecindex<-index
  vecindices<-c(vecindices,vecindex)
  lengths<-data.frame(tplot=vecindices, length=veclengths)
  lengths
}

veclength<-read.csv("~/Desktop/veclength_xc.csv")

ggplot(data=veclength, aes(x=treatment, y=veclength, fill=treatment))+
  geom_boxplot()+
  scale_fill_manual(values = c("sienna","royalblue2","peachpuff", "lightsteelblue1"))

ggplot(data=veclength, aes(x=shelterBlock, y=veclength, fill=shelterBlock))+
  geom_boxplot()

library(nlme)
v<-lme(veclength ~ treatment, random=~1|shelterBlock, veclength, na.action=na.exclude)
summary(v)
anova(v)
r.squaredGLMM(v) 
qqnorm(residuals(v))
qqline(residuals(v))
shapiro.test(residuals(v))
#normal
jLS<-lsmeans(v, ~treatment)
contrast(jLS, "pairwise") 

 ################################
###NMDS to compare grass-forb-both treatments
################################
data3<- cover %>% dplyr::select(-X, -species, -genus, -status, -func, -func2) %>% 
  filter(subplot!="XC"& subplot!="C") %>% spread(species_name, cover) %>% arrange(treatment) 
data3$ID <- seq.int(nrow(data3))
data3$TB<-paste(data3$treatment,data2$shelterBlock,sep=".")
data3$SB<-paste(data3$TB,data3$subplot,sep=".")


###some code to load and manipulate env variables (for later)
#LTM_env<-data %>% select(-c(18:51))
#plotnames<-LTM_env[,1]
#LTM.env<-LTM_env
#rownames(LTM.env)<-plotnames 
#Treatment=as.factor(LTM_env[,9])
#Year=as.factor(LTM_env[,5])
#LTM.env2 <-LTM.env[,-c(1:16)]
###standardize the environmental variables (scale function converts to z-scores)
#LTM.env.z <- data.frame(scale(LTM.env2))
TB3<-as.factor(data3[,66])
SB3<-as.factor(data3[,67])
Treatment3<-data3[,4]
Year3<-data3[,3]
subplot3<-as.factor(c("B", "F", "G"))
plotnames<-data3[,65]
cover.Bio3<- data3 %>% dplyr::select(-c(1:6), -ID, -TB, -Unknown, -SB)
rownames(cover.Bio3)<-plotnames
#check for empty rows
cover.Biodrop3<-cover.Bio3[rowSums(cover.Bio3[, (1:57)]) ==0, ] #no empty rows, next step not needed
# remove empty rows: cover.Biodrop<-cover.Bio[rowSums(cover.Bio[, (1:58)])  >0, ]

cover.rowsums3 <- rowSums(cover.Bio3 [1:57])
cover.relrow3 <- data.frame(cover.Bio3 /cover.rowsums3)
cover.pa3 <- cover.Bio3 %>% mutate_each(funs(ifelse(.>0,1,0)), 1:57)

#make bray-curtis dissimilarity matrix for controls
spp.bcd3 <- vegdist(cover.relrow3)
#or
spp.pa.bcd3<-vegdist(cover.pa3, binary=T) #this calcs distances based on presence/absence

#run NMS ordination
spp.mds0_3 <-isoMDS(spp.bcd3) #runs nms only once
spp.mds0_3  #by default 2 dimensions returned, stress = 20.04
ordiplot(spp.mds0_3)

#prefer to run multiple NMS ordinations
#with different starting configurations, and choose the best
#this function does that, and in addition does several other steps too
#including: 
#standardizing the data (though fuction call below turns this off with autotransform=F)
#calculating distance matrix (default bray-curtis)
#running NMDS with random starts
#rotation of axes to maximize variance of site scores on axis 1
#calculate species scores based on weighted averaging

spp.mds3<-metaMDS(cover.relrow3, trace = T, autotransform=F, trymax=999, k=3) #runs several with different starting configurations
#trace= TRUE will give output for step by step what its doing
#default is 2 dimensions, can put k=4 for 4 dimensions
spp.mds3 #converged after 20 tries, stress=13.2
summary(spp.mds3)

#plot results
stressplot(spp.mds3, spp.bcd3) #stressplot
ordiplot(spp.mds3)
spscores1_3<-scores(spp.mds3,display="sites",choices=1)
spscores2_3<-scores(spp.mds3,display="sites",choices=2)
tplots3<-data3[,66]
tplots3<-as.factor(tplots3)
tplot_levels3<-levels(tplots3)
spscoresall_3<-data.frame(tplots3,spscores1_3,spscores2_3)

#plots colored based on subplot
comp.plot <- ordiplot(spp.mds3,choices=c(1,2), type = "none")   #Set up the plot
cols <- rep(c("Purple","Blue","Red"), each = 3) #color based on composition treatment
cols1 <- rep(c("Purple","Blue","Red"))
shapes <- rep(c(15, 8, 17, 1 ), each=36) #shapes on drought treatment
shapes1 <- rep(c(15, 8, 17, 1 ), each=1)
points(spscoresall_3$NMDS1,spscoresall_3$NMDS2,col=cols,pch=shapes) 
text(spp.mds3, display = "species", cex=0.6, col="grey30") #label species
# add legend for subplot
legend("topleft",legend=levels(subplot3), col=cols1, pch=19, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)
# add legend for drought treatment
legend("topright",legend=levels(as.factor(as.character(Treatment3))), col="black", pch=shapes1, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)
help(ordiplot)

#plots colored based on year
comp.plot2 <- ordiplot(spp.mds3,choices=c(1,2), type = "none")   #Set up the plot
cols <- rep(c("black","cyan","orange"), each = 1) 
points(spscoresall_3$NMDS1,spscoresall_3$NMDS2,col=cols,pch=20, cex=2) 
#ordiellipse(spp.mds, groups=Year, fill=cols)
text(spp.mds3, display = "species", cex=0.8, col="black") #label species 
legend("topright",legend=levels(as.factor(as.character(Year3))), col=cols, pch=19, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)

#plots colored based on block
comp.plot3 <- ordiplot(spp.mds3,choices=c(1,2), type = "none")   #Set up the plot
cols <- rep(c("grey", "magenta", "yellow", "navy"), each = 9) 
cols1 <- c("grey", "magenta", "yellow", "navy") 
points(spscoresall_3$NMDS1,spscoresall_3$NMDS2,col=cols,pch=20, cex=2) 
text(spp.mds3, display = "species", cex=0.8, col="black") #label species 
legend("topright",legend=levels(data3$shelterBlock), col=cols1, pch=19, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)

permanova1_3 <- adonis(spp.bcd3~tplots3, perm=100, method="bray")
permanova1_3

permanova2_3 <- adonis(spp.bcd3~Year3, perm=100, method="bray")
permanova2_3 

permanova3_3 <- adonis(spp.bcd3~Treatment3, perm=999, method="bray")
permanova3_3

permanova4_3 <- adonis(spp.bcd3~data3$subplot*Year3*Treatment3*data3$shelterBlock, perm=999, method="bray")
permanova4_3 #significant treatment * block interaction 
pairwise.perm.manova(spp.bcd3, Treatment3, nperm=999)


##Creating an ordination plot with succession vectors
comp.plot4 <- ordiplot(spp.mds3,choices=c(1,2), type = "none")   #Set up the plot
cols <- rep(c("Red","Black","Orange", "Pink"), each = 36) #drought treatment
cols1 <- rep(c("Red","Black","Orange", "Pink"), each = 1) 
cols2 <- rep(c("Red","Black","Orange", "Pink"), each = 12) #36 divided by 3 years
shapes <- rep(c("B", "F", "G"), each=3) #shapes on subplot
shapes1 <- rep(c("B", "F", "G"), each=1) #shapes on subplot
points(spscoresall_3$NMDS1,spscoresall_3$NMDS2,col=cols,pch=shapes, cex=1)#Plot the ordination points 
text(spp.mds3, display = "species", cex=0.5, col="black") #label species
ordiarrows(spp.mds3, groups=SB3, order.by=Year3, label=F, col=cols2)
legend("topleft",legend=levels(Treatment3), col=cols1, pch=19, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)
legend("topright",legend=levels(subplot3), col="black", pch=shapes1, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)

##Creating an ordination plot with succession vectors
comp.plot5 <- ordiplot(spp.mds3,choices=c(1,2), type = "none")   #Set up the plot
cols <- rep(c("Red","Black","Orange", "Pink"), each = 36) #drought treatment
cols1 <- rep(c("Red","Black","Orange", "Pink"), each = 1) 
cols2 <- rep(c("Red","Black","Orange", "Pink"), each = 12) #36 divided by 3 years
shapes <- rep(c(16, 18, 2, 22), each=9) #shapes on block
shapes1 <- rep(c(16, 18, 2, 22), each=1) #shapes on block for legend
points(spscoresall_3$NMDS1,spscoresall_3$NMDS2,col=cols,pch=shapes, cex=1)#Plot the ordination points 
text(spp.mds3, display = "species", cex=0.5, col="black") #label species
ordiarrows(spp.mds3, groups=SB3, order.by=Year3, label=F, col=cols2)
legend("topleft",legend=levels(Treatment3), col=cols1, pch=19, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)
legend("topright",legend=levels(data3$shelterBlock), col="black", pch=shapes1, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)


##Make a 3 dimensional plot of the ordination
ordirgl (spp.mds3)


#############
#another look at composition but group by block first
#everything seems to be converging towards Avena dominated communities?
###############
################################
###NMDS to compare grass-forb-both treatments
################################
data4<- cover %>% dplyr::select(-X, -species, -genus, -status, -func, -func2) %>% 
  filter(subplot!="XC"& subplot!="C") %>% group_by(subplot, treatment, year, species_name) %>% 
  summarize(cover=mean(cover)) %>%
  spread(species_name, cover) %>% arrange(treatment) 
   
data4$ID <- seq.int(nrow(data4))
data4$TB<-paste(data4$treatment,data4$year,sep=".")
data4$SB<-paste(data4$treatment,data4$subplot,sep=".")
data4<-as.data.frame(data4)

###some code to load and manipulate env variables (for later)
#LTM_env<-data %>% select(-c(18:51))
#plotnames<-LTM_env[,1]
#LTM.env<-LTM_env
#rownames(LTM.env)<-plotnames 
#Treatment=as.factor(LTM_env[,9])
#Year=as.factor(LTM_env[,5])
#LTM.env2 <-LTM.env[,-c(1:16)]
###standardize the environmental variables (scale function converts to z-scores)
#LTM.env.z <- data.frame(scale(LTM.env2))

TB4<-as.factor(data4[,63])
SB4<-as.factor(data4[,64])
Treatment4<-data4[,2]
Year4<-data4[,3]
Subplot4<-as.factor(c("B", "F", "G"))
plotnames<-data4[,62]
data4<-as.data.frame(data4)
cover.Bio4<- data4 %>% dplyr::select(-subplot, -treatment, -year, -ID, -TB, -Unknown, -SB)
rownames(cover.Bio4)<-plotnames
#check for empty rows
cover.Biodrop4<-cover.Bio4[rowSums(cover.Bio4[, (1:57)]) ==0, ] #no empty rows, next step not needed
# remove empty rows: cover.Biodrop<-cover.Bio[rowSums(cover.Bio[, (1:58)])  >0, ]

cover.rowsums4 <- rowSums(cover.Bio4 [1:57])
cover.relrow4 <- data.frame(cover.Bio4 /cover.rowsums4)
cover.pa4 <- cover.Bio4 %>% mutate_each(funs(ifelse(.>0,1,0)), 1:57)

#make bray-curtis dissimilarity matrix for controls
spp.bcd4 <- vegdist(cover.relrow4)
#or
spp.pa.bcd4<-vegdist(cover.pa4, binary=T) #this calcs distances based on presence/absence

#run NMS ordination
spp.mds0_4 <-isoMDS(spp.bcd4) #runs nms only once
spp.mds0_4  #by default 2 dimensions returned
ordiplot(spp.mds0_4)

#prefer to run multiple NMS ordinations
#with different starting configurations, and choose the best
#this function does that, and in addition does several other steps too
#including: 
#standardizing the data (though fuction call below turns this off with autotransform=F)
#calculating distance matrix (default bray-curtis)
#running NMDS with random starts
#rotation of axes to maximize variance of site scores on axis 1
#calculate species scores based on weighted averaging

spp.mds4<-metaMDS(cover.relrow4, trace = T, autotransform=F, trymax=999, k=3) #runs several with different starting configurations
#trace= TRUE will give output for step by step what its doing
#default is 2 dimensions, can put k=4 for 4 dimensions
spp.mds4 #solution converges after 20 tries, stress 9.32
summary(spp.mds4)

#plot results
stressplot(spp.mds4, spp.bcd4) #stressplot
ordiplot(spp.mds4)
spscores1_4<-scores(spp.mds4,display="sites",choices=1)
spscores2_4<-scores(spp.mds4,display="sites",choices=2)
tplots4<-data4[,63]
tplots4<-as.factor(tplots4)
tplot_levels4<-levels(tplots4)
spscoresall_4<-data.frame(tplots4,spscores1_4,spscores2_4)

#plots colored based on subplot


compM.plot <- ordiplot(spp.mds4,choices=c(1,2), type = "none", xlim=c(-2.5,1.5),ylim=c(-1.5,1.5))   #Set up the plot
cols <- rep(c("sienna","royalblue2","lightsteelblue3", "peachpuff2"), each = 9) #color based on drought treatment
cols1 <- rep(c("sienna","royalblue2","lightsteelblue3", "peachpuff2"), each=3) #for arrows
cols2 <- rep(c("sienna","royalblue2","lightsteelblue3", "peachpuff2")) #for legend
shapes <- rep(c(15, 8, 21 ), each=3) #shapes on composition treatment
shapes1 <- rep(c(15, 8, 21), each=1)
points(spscoresall_4$NMDS1,spscoresall_4$NMDS2,col=cols,pch=shapes) 
text(spp.mds4, display = "species", cex=0.4, col=colspec) #label species
ordiarrows(spp.mds4, groups=SB4, order.by=Year4, label=F, col=cols1)
# add legend for subplot
legend("topleft",legend=levels(as.factor(as.character(Treatment4))), col=cols2, pch=19, cex=0.8,inset=0.05,bty="n",y.intersp=0.7,x.intersp=1.5,pt.cex=1)
# add legend for drought treatment
legend("topright",legend=levels(Subplot4), col="black", pch=shapes1, cex=0.8,inset=0.05,bty="n",y.intersp=0.7,x.intersp=0.8,pt.cex=1)
help(ordiplot)

#plots colored based on year
compM.plot2 <- ordiplot(spp.mds4,choices=c(1,2), type = "none")   #Set up the plot
cols <- rep(c("black","cyan","orange"), each = 1) 
points(spscoresall_4$NMDS1,spscoresall_4$NMDS2,col=cols,pch=20, cex=2) 
#ordiellipse(spp.mds, groups=Year, fill=cols)
text(spp.mds4, display = "species", cex=0.6, col="grey30") #label species 
legend("topright",legend=levels(as.factor(as.character(Year4))), col=cols, pch=19, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)


permanova1_4 <- adonis(spp.bcd4~tplots4, perm=100, method="bray")
permanova1_4

permanova2_4 <- adonis(spp.bcd4~Year4, perm=100, method="bray")
permanova2_4 #year significant R2=0.197

permanova3_4 <- adonis(spp.bcd4~Treatment4, perm=999, method="bray")
permanova3_4 #treatment significant R2=0.19

permanova4_4 <- adonis(spp.bcd4~Treatment4*Year4, perm=999, method="bray")
permanova4_4 #interaction term is not significant


##Creating an ordination plot with succession vectors
compM.plot3 <- ordiplot(spp.mds4,choices=c(1,2), type = "none")   #Set up the plot
cols <- rep(c("Red","Black","Orange", "Pink"), each = 9) #drought treatment
cols1 <- rep(c("Red","Black","Orange", "Pink"), each = 1) 
cols2 <- rep(c("Red","Black","Orange", "Pink"), each = 3) #9 divided by 3 years
cols3 <-rep(c("Purple","Blue","Red"), each=1)#to color by composition treatment
shapes <- rep(c("B", "F", "G"), each=3) #shapes on subplot
shapes1 <- rep(c("B", "F", "G"), each=1) #shapes on block
points(spscoresall_4$NMDS1,spscoresall_4$NMDS2,col=cols,pch=shapes, cex=1)#Plot the ordination points 
text(spp.mds4, display = "species", cex=0.5, col="grey30") #label species
ordiarrows(spp.mds4, groups=(as.factor(data4$SB)), order.by=(data4$year), label=F, col=cols2)
legend("topleft",legend=levels(Treatment4), col=cols1, pch=19, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)
legend("topright",legend=levels(Subplot4), col="black", pch=shapes1, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)

##Make a 3 dimensional plot of the ordination
compM3d<-ordirgl (spp.mds4, col=cols, ax.col= "lightblue")
orgltext(spp.mds4, text, display = "species", adj = 0.5, col = "black")
help(vegan3d)

######################
#Principle coordinates analysis
# ="Metric multidimensional scaling"
# = "Classic multidimensional scaling"
######################
help(cmdscale)
spp.cmd4<-cmdscale(spp.bcd4,eig=T)
spp.cmd4
#cmdscale does not calculate species scores automatically
#so we need to do this manually using wascores()
#will not handle NA values
spp.cmd.spscores4<-wascores(spp.cmd4$points,cover.Bio4)

#create plot
plot(spp.cmd4$points[,1],spp.cmd4$points[,2],pch=1,col="blue",xlab="Axis1",ylab="Axis2")
points(spp.cmd.spscores4,cex=0.7, pch=3, col="red")
text(spp.cmd.spscores4,rownames(spp.cmd.spscores4), cex=0.7, col="red") 
