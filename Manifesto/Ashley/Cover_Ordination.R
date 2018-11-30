library(vegan)
library(MASS)
library(tidyverse)
library(dplyr)

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
spp.pa.bcd<-vegdist(cover.pa, binary=T) #this calcs distances based on presence/absence

#run NMS ordination
spp.mds0 <-isoMDS(spp.bcd) #runs nms only once
spp.mds0  #by default 2 dimensions returned
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
spp.mds #note if solution converges or not
summary(spp.mds)

#plot results
stressplot(spp.mds, spp.bcd) #stressplot
ordiplot(spp.mds)
spscores1<-scores(spp.mds,display="sites",choices=1)
spscores2<-scores(spp.mds,display="sites",choices=2)
tplots<-data[,2]
tplot_levels<-levels(tplots)
spscoresall<-data.frame(tplots,spscores1,spscores2)

#plots colored based on treatment
bio.plot <- ordiplot(spp.mds,choices=c(1,2), type = "none")   #Set up the plot
cols <- rep(c("Red","Orange","purple", "black"), each = 15) #color based on drought treatment
shapes <- rep(c(15, 3, 17, 19, 5), each=3) #shapes on subplot
points(spscoresall$NMDS1,spscoresall$NMDS2,col=cols,pch=shapes) 
text(spp.mds, display = "species", cex=0.5, col="grey30") #label species 
#help(ordiplot)

#plots colored based on year
bio.plot2 <- ordiplot(spp.mds,choices=c(1,2), type = "none")   #Set up the plot
cols <- rep(c("black","cyan","orange"), each = 1) 
points(spscoresall$NMDS1,spscoresall$NMDS2,col=cols,pch=20) 
#ordiellipse(spp.mds, groups=Year, fill=cols)
text(spp.mds, display = "species", cex=0.6, col="grey30") #label species 

#plots colored based on block
bio.plot3 <- ordiplot(spp.mds,choices=c(1,2), type = "none")   #Set up the plot
cols <- rep(c("grey", "magenta", "yellow", "navy"), each = 60) 
points(spscoresall$NMDS1,spscoresall$NMDS2,col=cols,pch=20) 
text(spp.mds, display = "species", cex=0.8, col="black") #label species 
#seems to be different based on block?

permanova1 <- adonis(spp.bcd~data$shelterBlock, perm=100, method="bray")
permanova1

permanova2 <- adonis(spp.bcd~Year, perm=100, method="bray")
permanova2

permanova3 <- adonis(spp.bcd~Treatment, perm=999, method="bray")
permanova3

#permanova4 <- adonis(spp.bcd~factor(LTM.BioT$temp), perm=999, method="bray")
#permanova4


#overlay environmental data
#first remove rows that were removed from bio
#LTM.env3<-LTM.env2[ -c(642, 683, 691, 695), ]
#envvec.nms<-envfit(spp.mds,LTM.env3, na.rm=TRUE)
#envvec.nms
#plot(envvec.nms) #add vectors to previous ordination, only showing p<0.05

#if desired, rotate ordination so that first axis is parallel to an environmental variable
#spp.mds.rotelev<-MDSrotate(spp.mds,LTM.env3$moisture)
#envvec.nms.rotelev<-envfit(spp.mds.rotelev,LTM.env3, na.rm=T)
#ordiplot(spp.mds.rotelev)
#text(spp.mds.rotelev, display = "species", cex=0.7, col="red") 
#plot(envvec.nms.rotelev) 

####################
#repeat ordination for control subplots only
##################

data2<- cover %>% dplyr::select(-X, -species, -genus, -status, -func, -func2) %>% 
  filter(subplot=="XC") %>% spread(species_name, cover) %>% arrange(treatment) 
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
Treatment<-data2[,4]
Year<-data2[,3]
year<-as.factor(Year)
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
spp.mds0 <-isoMDS(spp.bcd2) #runs nms only once
spp.mds0  #by default 2 dimensions returned
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

spp.mds<-metaMDS(cover.relrow2, trace = T, autotransform=F, trymax=999, k=3, wascores=T) #runs several with different starting configurations
#trace= TRUE will give output for step by step what its doing
#default is 2 dimensions, can put k=4 for 4 dimensions
spp.mds #note if solution converges or not
summary(spp.mds)

#plot results
stressplot(spp.mds, spp.bcd) #stressplot
ordiplot(spp.mds)
spscores1<-scores(spp.mds,display="sites",choices=1)
spscores2<-scores(spp.mds,display="sites",choices=2)
tplots<-data2[,65]
tplots<-as.factor(tplots)
tplot_levels<-levels(tplots)
spscoresall<-data.frame(tplots,spscores1,spscores2)

library(goeveg)
## Select the 30% most frequent species with 50% best axis fit
limited <- ordiselect(cover.relrow2, spp.mds, ablim = 0.3, fitlim=0.5, method="axes")

#plots colored based on treatment
bio.plot <- ordiplot(spp.mds,choices=c(1,2), type = "none")   #Set up the plot
cols <- rep(c("Red","Orange","purple", "black"), each = 12) #color based on drought treatment
cols1 <- rep(c("Red","Orange","purple", "black"))
shapes <- rep(c(15, 8, 17 ), each=1) #shapes on year
points(spscoresall$NMDS1,spscoresall$NMDS2,col=cols,pch=shapes) 
text(spp.mds, display = "species", cex=0.7, col="black") #label species
# add legend for treatment
legend("bottomright",legend=levels(Treatment), col=cols1, pch=19, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)
# add legend for year
legend("topright",legend=levels(year), col="black", pch=shapes, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)
help(ordiplot)

#plots colored based on year
bio.plot2 <- ordiplot(spp.mds,choices=c(1,2), type = "none")   #Set up the plot
cols <- rep(c("black","cyan","orange"), each = 1) 
points(spscoresall$NMDS1,spscoresall$NMDS2,col=cols,pch=20) 
#ordiellipse(spp.mds, groups=Year, fill=cols)
text(spp.mds, display = "species", cex=0.8, col="black") #label species 
legend("bottomright",legend=levels(year), col=cols, pch=19, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)

#plots colored based on block
bio.plot3 <- ordiplot(spp.mds,choices=c(1,2), type = "none")   #Set up the plot
cols <- rep(c("grey", "magenta", "yellow", "navy"), each = 3) 
cols1 <- c("grey", "magenta", "yellow", "navy") 
points(spscoresall$NMDS1,spscoresall$NMDS2,col=cols,pch=20) 
text(spp.mds, display = "species", cex=0.8, col="black") #label species 
legend("bottomright",legend=levels(data2$shelterBlock), col=cols1, pch=19, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)

permanova1 <- adonis(spp.bcd2~tplots, perm=100, method="bray")
permanova1

permanova2 <- adonis(spp.bcd2~Year, perm=100, method="bray")
permanova2

permanova3 <- adonis(spp.bcd2~Treatment*Year*data2$shelterBlock, perm=999, method="bray")
permanova3


##Create an ordination plot with succession vectors over years for treatment*block 
cover.plot <- ordiplot(spp.mds,choices=c(1,2), type = "none")   #Set up the plot
cols <- rep(c("Red","Orange","purple", "black"), each = 12) 
cols1 <- rep(c("Red","Orange","purple", "black"), each = 1) 
cols2 <- rep(c("Red","Orange","purple", "black"), each = 4) 
shapes <- rep(c(15, 3, 17, 1), each=3) #shapes on block
shapes1 <- rep(c(15, 3, 17, 1), each=1) #shapes on block
points(spscoresall$NMDS1,spscoresall$NMDS2,col=cols,pch=shapes)#Plot the ordination points 
text(spp.mds, display = "species", cex=0.8, col="grey30") #label species
ordiarrows(spp.mds, groups=TB, order.by=Year, label=F, col=cols2)
legend("bottomright",legend=levels(Treatment), col=cols1, pch=19, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)
legend("topright",legend=levels(data2$shelterBlock), col="black", pch=shapes1, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)
library (vegan3d)
ordirgl (spp.mds)


######################
#4. Principle coordinates analysis
# ="Metric multidimensional scaling"
# = "Classic multidimensional scaling"
######################
help(cmdscale)
spp.cmd<-cmdscale(spp.bcd2,eig=T)
spp.cmd
#cmdscale does not calculate species scores automatically
#so we need to do this manually using wascores()
#will not handle NA values
spp.cmd.spscores<-wascores(spp.cmd$points,cover.Bio2)

#create plot
#consider using package LabDSV if you plan to work a lot with PCoA
#otherwise stick with brute force as here
plot(spp.cmd$points[,1],spp.cmd$points[,2],pch=1,col="blue",xlab="Axis1",ylab="Axis2")
points(spp.cmd.spscores,cex=0.7, pch=3, col="red")
text(spp.cmd.spscores,rownames(spp.cmd.spscores), cex=0.7, col="red") 
