library(cluster)
library(indicspecies)

####################
#repeat ordination for control (XC) subplots only
##################

data2<- cover %>% dplyr::select(-X, -species, -genus, -status, -func, -func2) %>% 
  filter(subplot=="XC") %>% spread(species_name, cover) %>% arrange(treatment, shelterBlock) 
data2$ID <- seq.int(nrow(data2))
data2$TB<-paste(data2$treatment,data2$shelterBlock,sep=".")
data2<-data2 %>% dplyr::select(-Unknown)


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

soil<-read.csv ("~/Dropbox/ClimVar/DATA/Soil Extracts/Soil_CleanedData/ClimVar_soil_2015.csv")
soil_xc<- soil %>% filter (subplot=="XC", Depth == "1", DOY=="69")  %>% 
  arrange(treatment, shelterBlock, subplot) #arrange so all files used in PCA match
env<- soil_xc %>% dplyr::select(SOC,SIC,MBC,MBIC,NH4,NO3.NO2)
env.z <- data.frame(scale(env)) #standardize env variables using z-scores

env.euc<-vegdist(env.z, "euclidian")
mrpp(spp.bcd2,clustmem$clust4)

#plot results
stressplot(spp.mds2, spp.bcd2) #stressplot
ordiplot(spp.mds2)
spscores1_2<-scores(spp.mds2,display="sites",choices=1)
spscores2_2<-scores(spp.mds2,display="sites",choices=2)
tplots2<-data2[,65]
tplots2<-as.factor(tplots2)
tplot_levels2<-as.factor(levels(tplots2))
spscoresall_2<-data.frame(tplots2,spscores1_2,spscores2_2)

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
for (i in 2:20){
  temp.clust<-cutree(spp.clust.beta,k=i)
  clustmem<-cbind(clustmem,temp.clust)
}
clustmem<-data.frame(clustmem)
colnames(clustmem)<-paste("clust",c(2:20),sep="")
clustmem

#cluster membership can be treated like a categorical variable
May_all_XC<-arrange(May_all_XC, treatment, shelterBlock)
boxplot(May_all_XC$ANPPgm~clustmem$clust4)
boxplot(May_all_XC$RaoQ~clustmem$clust4)
boxplot(May_all_XC$CWM.Ht~clustmem$clust4)

May_all_XC$clust4<-clustmem$clust4
May_all_XC$clust3<-clustmem$clust3
ggplot(data=May_all_XC, aes(y=ANPPgm, x=treatment, color=as.factor(clust4)))+
  geom_boxplot()

#overlay cluster membership on nmds ordination
plot(spp.mds2, choices=c(1,2), type="n")
points(spp.mds2, display="sites", cex=0.8, pch=16, col=clustmem$clust4)

##Creating an ordination plot with succession vectors
xc.plot4 <- ordiplot(spp.mds2,choices=c(1,2), type = "none", xlim=c(-1.5,1.5), ylim=c(-1.5,1.5))   #Set up the plot
cols <- rep(c("sienna","royalblue2","lightsteelblue3", "peachpuff2"), each = 12) 
cols1 <- rep(c("sienna","royalblue2","lightsteelblue3", "peachpuff2"), each = 1) 
cols2 <- rep(c("sienna","royalblue2","lightsteelblue3", "peachpuff2"), each = 4) 
shapes <- rep(c(15, 3, 17), each=1) #shapes on year
shapes1 <- rep(c(15, 3, 17), each=1) #shapes on year
points(spscoresall_2$NMDS1,spscoresall_2$NMDS2,col=clustmem$clust4,pch=shapes, cex=1.2)#Plot the ordination points 
text(spp.mds2, display = "species", cex=0.4, col="grey30") #label species
ordiarrows(spp.mds2, groups=TB, order.by=Year2, label=F, col=cols2)
legend("topright",legend=levels(data2$treatment), col=cols1, pch=19, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)
legend("bottomright",legend=levels(as.factor(data2$year)), col="black", pch=shapes1, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)

#clustmem is a dataframe with columns indicating group membership at different numbers of clusters
#now create a dataframe showing how many significant indicator species at each clustering level
alpha<-0.05 #can change alpha here
c1<-colnames(clustmem)
c2<-c()
for (i in 1:length(clustmem)){
  temp.isa<-multipatt(cover.relrow2,clustmem[,i],duleg=T)
  nu.sig<-sum(temp.isa$sign$p.value<=alpha, na.rm=T)
  c2<-c(c2,nu.sig)
}
data.frame(c1,c2)

#indicator species by cluster membership
xc_isa = multipatt(cover.relrow2, clustmem$clust4, control=how(nperm=999))
summary(xc_isa)

#indicator species by treatment
trt_isa = multipatt(cover.relrow2, May_all_XC$treatment, control=how(nperm=999))
summary(trt_isa)
