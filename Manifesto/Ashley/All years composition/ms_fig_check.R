library(vegan)
library(MASS)
library(tidyverse)
library(dplyr)
library(RVAideMemoire) #for posthoc tests on permanova
library(goeveg)#for scree plot of NMDS to test number of dimensions
library(gridExtra)
library(cluster)
library(nlme)
library(indicspecies)

theme_set(theme_bw())
theme_update( panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
              strip.background = element_blank(),
              text = element_text(size = 16),
              strip.text= element_text(size = 14), 
              axis.text = element_text(size = 12))

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

data2<- cover %>% dplyr::select(-X, -species, -genus, -status, -func, -func2) %>% 
  filter(subplot=="XC", year == "2015") %>% spread(species_name, cover) %>% arrange(treatment, shelterBlock) 
data2$ID <- seq.int(nrow(data2))
data2$TB<-paste(data2$treatment,data2$shelterBlock,sep=".")
data2<-data2 %>% dplyr::select(-Unknown)
data2[is.na(data2)] <- 0

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
#spp.pa.bcd2<-vegdist(cover.pa2, binary=T) #this calcs distances based on presence/absence

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

dimcheckMDS(cover.Bio2, distance = "bray", k = 8, trymax = 20, autotransform = F) 
spp.mds2<-metaMDS(cover.Bio2, trace = T, autotransform=F, trymax=999, k=5) #runs several with different starting configurations
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

#test for treatment differences
permanova1<-adonis(cover.Bio2~data2$treatment, strata=data2$shelterBlock, perm=1000, method="bray")
permanova1


#first put treatments in order
data2$treatment2 <- factor(data2$treatment, levels = c("controlRain", "fallDry","springDry","consistentDry"))
Fig1a<-ggplot(spscoresall_2, aes(x=NMDS1, y=NMDS2))+
  xlim(-1,1)+
  geom_path(data=spscoresall_2, arrow=arrow(), mapping=aes(x=NMDS1, y=NMDS2, group=TB))+
  geom_point(cex=4, aes(shape=as.factor(data2$treatment2), fill=data2$treatment2))+
  ggtitle("a) Community ordination showing rainfall treatments")+
  scale_fill_manual(values=c("Black", "gray40","lightgrey", "white"), guide = guide_legend(title = "Treatment:"), #change legend title
                    labels=c("Control", "Early\nDrought", "Late\nDrought", "Consistent\nDrought"))+ #change labels in the legend)+
  scale_shape_manual(values=c(21,24,22,23),guide = guide_legend(title = "Treatment:"), #change legend title
                     labels=c("Control", "Early\nDrought", "Late\nDrought", "Consistent\nDrought"))+
  #theme(legend.position=c(0.7, 0.1), legend.direction="horizontal",legend.key.size = unit(0.1,"cm"),legend.title=element_text(size=11), legend.text=element_text(size=9))+
  theme(legend.position="top", legend.direction="horizontal",legend.key.size = unit(0.1,"cm"),legend.title=element_text(size=11), legend.text=element_text(size=9))+
  theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))
Fig1a

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(Fig1a)
mylegend

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

Fig1b<- ggplot(data=veclength, aes(x=treatment, y=veclength, fill=treatment))+
  ylim(0,1.5)+
  ggtitle("b) Community change by treatment")+
  geom_boxplot()+
  theme(legend.position = "none")+
  labs(x="Rainfall Treatment", y = "Vector Length")+
  scale_x_discrete(labels = c("Consistent\n Drought","Control", "Early\n Drought", "Late\n Drought"))+
  #scale_fill_manual(values = c("sienna","royalblue2","peachpuff", "lightsteelblue1"))
  scale_fill_manual(values=c("white", "Black", "gray40","lightgrey" ), guide = guide_legend(title = "Treatment"), #change legend title
                    labels=c("Consistent Drought","Control", "Early Drought", "Late Drought" )) #change labels in the legend)+
Fig1b

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

#############
#Cluster Analysis
#############
spp.clust.sing <- hclust(spp.bcd2, "single") #using single (nearest neighbor) linkage
plot(spp.clust.sing, cex=0.8) #maybe 3 or 4 groups

spp.clust.comp <- hclust(spp.bcd2, "complete") #using complete (farthest neighbor) linkage
plot(spp.clust.comp, cex=0.8) #3 clusters

#plot dendrograms side by side
par(mfrow=c(1,2))
plot(spp.clust.sing, cex=0.5)
plot(spp.clust.comp, cex=0.5)
par(mfrow=c(1,1))

##flexible beta method
spp.clust.beta<-agnes(spp.bcd2, method="flexible", par.method=0.625)#par.method as 0.625 is equivalent to beta -0.25
spp.clust.beta<-as.hclust(spp.clust.beta)#match format of hclust
plot(spp.clust.beta, cex=0.8) #4 clusters

##from dendrogram, define groups, either by specifying k (number of groups desired) or
#by specifying h=height where to slice
spp.clust.beta.h1 <- cutree(spp.clust.beta, h=1)
spp.clust.beta.h1 #5 groups

spp.clust.beta.k4<-cutree(spp.clust.beta, k=4)
spp.clust.beta.k4

#loop calculates a df with cluster membership from k=2 to k=20
clustmem<-c()
for (i in 2:10){
  temp.clust<-cutree(spp.clust.beta,k=i)
  clustmem<-cbind(clustmem,temp.clust)
}
clustmem<-data.frame(clustmem)
colnames(clustmem)<-paste("clust",c(2:10),sep="")
clustmem

#cluster membership can be treated like a categorical variable
data2$clust15<-clustmem$clust2
d<-data2 %>% dplyr::select(treatment,shelterBlock, clust15)

data3<- cover %>% dplyr::select(-X, -species, -genus, -status, -func, -func2) %>% 
  filter(subplot=="XC") %>% spread(species_name, cover) %>% arrange(treatment, shelterBlock) 

data4<-right_join(d,data3, by=c("treatment", "shelterBlock"))
data4<-merge(data4, May_all_XC)

ggplot(data=data4, aes(y=ANPPgm, x=treatment, color=as.factor(clust15)))+
  geom_boxplot()

c<-lme(weight_g_m ~treatment, random=~1|year/shelterBlock, subset(data4, clust15=="1"), na.action=na.exclude)
summary(c)
anova(c) #treatment is significant
r.squaredGLMM(c) #12% of variation explained by fixed effects, 57% by whole model (interannual variation?)
qqnorm(residuals(c))
qqline(residuals(c))
shapiro.test(residuals(c))
#normally distributed, continue
LS.all2<-lsmeans(c, ~treatment)
contrast(LS.all2, "pairwise")

#overlay cluster membership on nmds ordination
plot(spp.mds2, choices=c(1,2), type="n")
points(spp.mds2, display="sites", cex=0.8, pch=16, col=clustmem$clust2)

