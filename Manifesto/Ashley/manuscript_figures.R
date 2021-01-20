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
  filter(subplot=="XC") %>% spread(species_name, cover) %>% arrange(treatment, shelterBlock) 
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
permanova1<-adonis(cover.Bio2~data2$treatment*data2$year, perm=1000, method="bray")
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
for (i in 2:20){
  temp.clust<-cutree(spp.clust.beta,k=i)
  clustmem<-cbind(clustmem,temp.clust)
}
clustmem<-data.frame(clustmem)
colnames(clustmem)<-paste("clust",c(2:20),sep="")
clustmem

#cluster membership can be treated like a categorical variable
data2$clust<-clustmem$clust2
data2<-merge(May_all_XC, data2)


ggplot(data=data2, aes(y=ANPPgm, x=treatment, color=as.factor(clust)))+
  geom_boxplot()

c<-lme(weight_g_m ~treatment*as.factor(clust), random=~1|year/shelterBlock, subset(data2, clust2=="1"), na.action=na.exclude)
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

##Creating an ordination plot with succession vectors
xc.plot4 <- ordiplot(spp.mds2,choices=c(1,2), type = "none", xlim=c(-1.5,1.5), ylim=c(-1.5,1.5))   #Set up the plot
cols <- rep(c("sienna","royalblue2","lightsteelblue3", "peachpuff2"), each = 12) 
cols1 <- rep(c("sienna","royalblue2","lightsteelblue3", "peachpuff2"), each = 1) 
cols2 <- rep(c("sienna","royalblue2","lightsteelblue3", "peachpuff2"), each = 4) 
colspec<- rep(c("lightgrey", "lightgrey", "blue", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "blue", "lightgrey", "lightgrey","lightgrey", "blue", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "blue", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "blue", "blue","lightgrey", "lightgrey"))
shapes <- rep(c(15, 18, 17, 19), each=12) #shapes on trt
shapes1 <- rep(c(15, 18, 17, 19), each=1) #shapes on trt
points(spscoresall_2$NMDS1,spscoresall_2$NMDS2,col=clustmem$clust2,pch=shapes, cex=1.2)#Plot the ordination points 
text(spp.mds2, display = "species", cex=0.6, col=colspec) #label species
ordiarrows(spp.mds2, groups=TB, order.by=Year2, label=F, col="grey")
legend("topright",legend=levels(as.factor(clustmem$clust2)), col=c("black", "red"), pch=19, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)
legend("bottomright",legend=levels(as.factor(data2$treatment)), col="black", pch=shapes1, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)

#nice plot for pub
species<-as.data.frame(spp.mds2$species)
species$name<-row.names(species)
spc<- species %>% filter(name == "Avena.barbata"| name=="Lolium.multiflorum"| name=="Vicia.sativa"| name == "Cynodon.dactylon" |name=="Centaurea.solstitialis"| name=="Vulpia.bromoides"|name=="Trifolium.glomeratum")
spc<- species %>% filter(name == "Avena barbata"| name=="Lolium multiflorum"| name=="Vicia sativa"| name == "Cynodon dactylon" |name=="Centaurea solstitialis"| name=="Vulpia bromoides"|name=="Trifolium glomeratum")

Fig1c<-ggplot(spscoresall_2, aes(x=NMDS1, y=NMDS2))+
  xlim(-1,1)+
  geom_path(data=spscoresall_2, arrow=arrow(), mapping=aes(x=NMDS1, y=NMDS2, group=TB))+
  geom_point(cex=5, aes(shape=as.factor(data2$treatment2), fill=as.factor(clustmem$clust2)))+
  ggtitle("c) Community ordination showing clustering ")+
  scale_fill_manual(values=c("white", "black"), guide = guide_legend(title = "Cluster:", override.aes = list(shape = c(1, 16))))+ #change labels in the legend)+
  scale_shape_manual(values=c(21,24,22,23),guide = guide_legend(title = "Treatment:"), #change legend title
                     labels=c("Control", "Early\nDrought", "Late\nDrought", "Consistent\nDrought"))+
  theme(legend.position=c(0.7, 0.1), legend.direction="horizontal",legend.box.just="right", legend.key.size = unit(0.1,"cm"), legend.title=element_text(size=11), legend.text=element_text(size=9))+
  theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))
Fig1c
Fig1c<-Fig1c+geom_label(data=spc, mapping=aes(x=MDS1, y=MDS2, label=name, nudge_y=0.1), cex=3)
Fig1c

#redo fig1a with hull around clusters
spscoresall_2$clust2<-clustmem$clust2
clust1 <- spscoresall_2[spscoresall_2$clust2 == "1", ][chull(spscoresall_2[spscoresall_2$clust2 == 
                                                                   "1", c("NMDS1", "NMDS2")]), ]  # hull values for cluster 1
clust2 <- spscoresall_2[spscoresall_2$clust2 == "2", ][chull(spscoresall_2[spscoresall_2$clust2 == 
                                                                             "2", c("NMDS1", "NMDS2")]), ]  # hull values for cluster 1

hull.data <- rbind(clust1 , clust2)  #combine clust1 and clust2
hull.data

Fig1a<-ggplot(spscoresall_2, aes(x=NMDS1, y=NMDS2))+
  xlim(-1,1.3)+
  geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,fill=as.factor(clust2),group=as.factor(clust2)),alpha=0.6) + # add the convex hulls
  #scale_fill_manual(values=c("lightgrey", "black"), guide = guide_legend(title = "Cluster:"), #change legend title
                    #labels=c("1", "2")) #change labels in the legend)
geom_path(data=spscoresall_2, arrow=arrow(), mapping=aes(x=NMDS1, y=NMDS2, group=TB))+
  geom_point(cex=4, aes(shape=as.factor(data2$treatment2), fill=data2$treatment2))+
  ggtitle("a) Community composition over time from 2015 to 2017")+
  scale_fill_manual(values=c("gray90", "black","Black", "gray40","lightgrey", "white"), guide = guide_legend(title = "Treatment:"), #change legend title
                    labels=c("Cluster1","Cluster2","Control", "Early\nDrought", "Late\nDrought", "Consistent\nDrought"))+ #change labels in the legend)+
  scale_shape_manual(values=c(21,24,22,23),guide = guide_legend(title = "Treatment:"), #change legend title
                     labels=c("Control", "Early\nDrought", "Late\nDrought", "Consistent\nDrought"))+
  theme(legend.position="none")
  #theme(legend.position=c(0.7, 0.1), legend.direction="horizontal",legend.key.size = unit(0.1,"cm"),legend.title=element_text(size=11), legend.text=element_text(size=9))+
  #theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))
Fig1a
Fig1a<-Fig1a+geom_label(data=spc, mapping=aes(x=MDS1, y=MDS2, label=name), position=position_nudge(x=0.2), cex=3)
Fig1a<-Fig1a+annotate("text", y=c(1.2, 1.2), x=c(0.5,-0.7), label=c("Cluster 2","Cluster 1"), cex=6)
Fig1a

clust_xc<-data2 %>% filter(year==2015) #%>% dplyr::select(treatment, shelterBlock, clust2)
veclength<-veclength %>% arrange(treatment,shelterBlock)
veclength2<-merge(clust_xc,veclength) %>%filter(clust==1|clust==2)

v2<-lme(veclength ~ clust, random=~1|treatment/shelterBlock, veclength2, na.action=na.exclude)
summary(v2)
anova(v2)
r.squaredGLMM(v2) 
qqnorm(residuals(v2))
qqline(residuals(v2))
shapiro.test(residuals(v2))
#normal

Fig1d<- ggplot(data=veclength2, aes(x=as.factor(clust), y=veclength, fill=as.factor(clust), alpha=0.8))+
  ylim(0,1.5)+
  geom_boxplot()+
  ggtitle("c) Community change by cluster")+
  #theme(legend.position="none")+
  labs(x="Cluster", y="Vector Length")+
  annotate("text", x= c("1", "2"), y = c(1.5, 1.5), label = c("a", "b"), color = "black")+
  #facet_wrap(~clust)+
  theme(legend.position="none")+
  scale_x_discrete(labels = c("1\n ","2\n  "))+
  scale_fill_manual(values = c("gray90","black"))
Fig1d

#create layout for panel
lay <- rbind(c(1,1,1,1),
             c(1,1,1,1),
             c(1,1,1,1),
             c(1,1,1,1),
             c(1,1,1,1),
             c(1,1,1,1),
             c(NA,4,4,NA),
             c(NA,NA,NA,NA),
             c(2,2,3,3),
             c(2,2,3,3),
             c(2,2,3,3),
             c(2,2,3,3),
             c(2,2,3,3))
        
myplots<-list(Fig1a,Fig1b,Fig1d)
grid.arrange(Fig1a,Fig1b,Fig1d,mylegend, layout_matrix=lay) #put panel together


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
xc_isa = multipatt(cover.relrow2, clustmem$clust2, control=how(nperm=999))
summary(xc_isa)

#indicator species by treatment
trt_isa = multipatt(cover.relrow2, data2$treatment, control=how(nperm=999))
summary(trt_isa)

#test for cluster differences in community composition
permanova2<-adonis(cover.Bio2~clustmem$clust2*data2$year, perm=1000, method="bray")
permanova2

##within cluster 1, do compensatory dynamics occur?
ggplot(data2, aes(y=propForb, x=treatment, fill=as.factor(clustmem$clust2)))+
  geom_boxplot()

ggplot(data2, aes(y=RaoQ, x=as.factor(clust)))+
  geom_boxplot()

mR<-lme(RaoQ ~treatment*as.factor(clust), random=~1|year/shelterBlock, data2, na.action=na.exclude)
summary(mR)
anova(mR) #treatment is significant
r.squaredGLMM(mR) 
qqnorm(residuals(mR))
qqline(residuals(mR))
shapiro.test(residuals(mR))
LSmR<-lsmeans(mR, ~treatment)
contrast(LSmR, "pairwise")

Fig2a<- ggplot(data=data2, aes(x=as.factor(clust), y=RaoQ, fill=as.factor(clust), alpha=0.8))+
  geom_boxplot()+
  ggtitle("a)")+
  #theme(legend.position="none")+
  labs(x="Cluster", y="Rao's Q")+
  annotate("text", x= c("1", "2"), y = c(10, 10), label = c("a", "b"), color = "black")+
  #facet_wrap(~clust)+
  theme(legend.position="none")+
  scale_fill_manual(values = c("grey90","black"))
Fig2a

mA<-lme(Avena ~treatment*as.factor(clust), random=~1|year/shelterBlock, data2, na.action=na.exclude)
summary(mA)
anova(mA) #treatment is significant
r.squaredGLMM(mA) 
qqnorm(residuals(mA))
qqline(residuals(mA))
shapiro.test(residuals(mA))
LSmA<-lsmeans(mA, ~treatment)
contrast(LSmA, "pairwise")

Fig2b<- ggplot(data=data2, aes(x=as.factor(clust), y=Avena, fill=as.factor(clust), alpha=0.8))+
  geom_boxplot()+
  ggtitle("b)")+
  #theme(legend.position="none")+
  labs(x="Cluster", y="Avena barbata (% cover)")+
  annotate("text", x= c("1", "2"), y = c(110, 110), label = c("a", "b"), color = "black")+
  #facet_wrap(~clust)+
  theme(legend.position="none")+
  scale_fill_manual(values = c("grey90","black"))
Fig2b

#make figure 2 for manuscript
grid.arrange(Fig2a,Fig2b, ncol=2)

ggplot(data2, aes(y=FDis, x=as.factor(clust), fill=as.factor(clust)))+
  geom_boxplot()

ggplot(data2, aes(y=H, x=as.factor(clust), fill=as.factor(clust)))+
  geom_boxplot()

ggplot(data2, aes(y=ANPPgm, x=as.factor(clust),  color=as.factor(clust)))+
  geom_boxplot()+
  geom_smooth(method="lm", se=F)

data2 <- data2 %>% mutate(treatment=ordered(treatment, levels = c(controlRain="controlRain",  fallDry="fallDry", springDry="springDry", consistentDry="consistentDry")))
ggplot(data2, aes(y=forbCover, x=as.factor(clust2),  fill=treatment))+
  geom_boxplot()+
  scale_fill_manual(values = c("royalblue2","peachpuff","lightsteelblue1",  "sienna"), guide = guide_legend(title = "Treatment")) +
  labs(x="Cluster Membership", y="Forb Cover (%)")+
  theme_classic()

data2$clust15 <- clust_xc$clust[match(data2$TB, clust_xc$TB)]
data2_cv<-data2 %>% group_by(treatment, clust, shelterBlock) %>% summarize(CVanpp=sd(ANPPgm)/mean(ANPPgm)*100) #%>% 
data2_cv<- na.omit(data2_cv) %>% group_by(clust15, treatment)%>%summarize(CV=mean(CVanpp), se=sd(CVanpp)/sqrt(length(CVanpp)))

#supplementary figure for CV ANPP
ggplot(data2_cv, aes(y=CVanpp, x=as.factor(clust), fill=as.factor(clust)))+
  #geom_point(cex=6, pch=21)+
  #geom_errorbar(aes(ymin=CV-se, ymax=CV+se), width=1)+
  geom_boxplot()+
  scale_fill_manual(values = c("white","black"), guide = guide_legend(title = "Cluster")) +
  labs(x="Cluster Membership", y="CV ANPP")+
  theme_classic()

c<-lme(ANPPgm ~treatment2*year, random=~1|shelterBlock, subset(data2, clust==1), na.action=na.exclude)
summary(c)
anova(c) #treatment is significant
r.squaredGLMM(c) #12% of variation explained by fixed effects, 57% by whole model (interannual variation?)
qqnorm(residuals(c))
qqline(residuals(c))
shapiro.test(residuals(c))
#normally distributed, continue
LS.all2<-lsmeans(c, ~treatment2)
contrast(LS.all2, "pairwise")

d<-lme(ANPPgm ~treatment2*year, random=~1|shelterBlock, subset(data2, clust=="2"), na.action=na.exclude)
summary(d)
anova(d) #treatment is significant
r.squaredGLMM(d) #12% of variation explained by fixed effects, 57% by whole model (interannual variation?)
qqnorm(residuals(d))
qqline(residuals(d))
shapiro.test(residuals(d))
#normally distributed, continue
LS.all2d<-lsmeans(d, ~treatment)
contrast(LS.all2d, "pairwise")

May_all_XC<-May_all_XC%>%arrange(treatment, shelterBlock)
May_ANPP_XC<-May_ANPP_XC%>%mutate(treatment=ordered(treatment, levels = c(controlRain="controlRain", springDry="springDry", fallDry="fallDry", consistentDry="consistentDry")))%>%
  arrange(treatment, shelterBlock) #from Cover_dataanalyses.R

ggplot(data2, aes(y=nbsp, x=as.factor(clust)))+
  geom_boxplot()+
  scale_color_manual(values = c("royalblue2","lightsteelblue1", "peachpuff", "sienna"), guide = guide_legend(title = "Treatment")) +
  labs(x="Cluster Membership", y="Species Richness")+
  theme_classic()

ggplot(data2, aes(y=RaoQ, x=as.factor(clust)))+
  geom_boxplot()+
  scale_color_manual(values = c("royalblue2","lightsteelblue1", "peachpuff", "sienna"), guide = guide_legend(title = "Treatment")) +
  labs(x="Cluster Membership", y="Functional Diversity (RaoQ)")+
  theme_classic()


##Figure 3
library(tidyverse)
library(nlme)
library(ggplot2)
library(multcomp)
library(lsmeans)
library(MuMIn)
library(RColorBrewer)
setwd("~/Dropbox/ClimVar/DATA/")
May_ANPP<-read.csv("Plant_composition_data/ANPP/ANPP_CleanedData/ClimVar_ANPP-peak.csv")

#change plots, years, shelter to factors
May_ANPP[,'plot'] <- as.factor(as.character(May_ANPP[,'plot']))
May_ANPP[,'year'] <- as.factor(as.character(May_ANPP[,'year']))
May_ANPP[,'shelter'] <- as.factor(as.character(May_ANPP[,'shelter']))

#create subset with no species manipulations (control community) only
May_XC<-filter(May_ANPP, subplot=='XC')

m1<-lme(weight_g_m ~treatment, random=~1|shelterBlock/year, May_XC, na.action=na.exclude)
summary(m1)
anova(m1) #treatment is significant
r.squaredGLMM(m1) #12% of variation explained by fixed effects, 57% by whole model (interannual variation?)
qqnorm(residuals(m1))
qqline(residuals(m1))
shapiro.test(residuals(m1))
#normally distributed, continue
LS1<-lsmeans(m1, ~treatment)
contrast(LS1, "pairwise")
#control rain ANPP is significantly greater than all the drought except spring dry, no surprise
#control rain is most similar to spring dry

m1b<-lme(ANPPgm ~treatment*as.factor(clust), random=~1|year/shelterBlock, data2, na.action=na.exclude)
summary(m1b)
anova(m1b) #treatment is significant
r.squaredGLMM(m1b) 
qqnorm(residuals(m1b))
qqline(residuals(m1b))
shapiro.test(residuals(m1b))
LS1b<-lsmeans(m1b, ~treatment)
contrast(LS1b, "pairwise")


May_XC$treatment2 <- as.character(May_XC$treatment)
#Then turn it back into a factor with the levels in the correct order
May_XC$treatment2 <- factor(May_XC$treatment2, levels = c("controlRain", "fallDry","springDry","consistentDry"))

label1 <- paste("Treatment: F['3,33']==4.59")
label2<-paste("Year: F['2,33']==12.36")
label3<-paste("Treatment*Year: F['6,3'3]== 2.10")
label4<- "p==0.001"
Fig3a<-ggplot(May_XC, aes(x = treatment2, y = weight_g_m)) + 
  geom_boxplot(aes(fill = treatment2)) + 
  ggtitle("a) Overall treatment effects on ANPP")+
  scale_fill_manual(values = c("Black", "gray40","lightgrey", "white"), guide = guide_legend(title = "Treatment")) +
  #geom_jitter(position=position_jitter(0.2), aes(color=year)) +
  scale_x_discrete(labels = c("Control\n ", "Early\nDrought\n", "Late\n Drought\n","Consistent\n Drought\n"))+
  #geom_jitter(position=position_jitter(0.2), aes(color=May_XC$shelterBlock, shape=as.factor(year))) +
  scale_color_manual(values = c("gray80","gray50", "black"), guide = guide_legend(title = "Year")) +
  #facet_wrap(~year)+
  annotate("text", x= c("controlRain","consistentDry", "fallDry","springDry"), y = c( 975, 900, 900,900), label = c("a", "b", "b", "ab"), color = "black", cex=6) +
  #geom_label(aes(x=3, y = max(weight_g_m - 2)), label =  c(label1,label2,label3), parse = T)+ #x=c("springDry","consistentDry"), y=c(1025,1025), parse=TRUE)+
  annotate("text", label= paste(c("Treatment: F ['3,33']==4.59", "p==0.001","Year: F['2,33']==12.36", "p < 0.001","TreatmentxYear: F['6,33'] ==  2.10", "p == 0.079")),
           x=c(3.5,4,3.5,4,3.5,4), y=c(1025,1025,1000,1000,975,975), cex=4, hjust = 1, parse=T)+
  theme(legend.position="none")+
  labs(x="Rainfall Treatment", y = expression(paste("ANPP ", g/m^{2})))
Fig3a

Fig3b<-ggplot(d=May_XC, aes(x=as.factor(year), y=weight_g_m)) +
  ggtitle("b) Overall year effects")+
  geom_boxplot(aes(y=weight_g_m), shape=16)+
  #scale_fill_manual(values = c("gray99","gray80", "gray50"), guide = guide_legend(title = "Year")) +
  #geom_jitter(position=position_jitter(0.2), aes(color=treatment2)) +
  #scale_color_manual(values = c("royalblue2", "peachpuff","lightsteelblue1", "sienna"), guide = guide_legend(title = "Treatment")) +
  labs(x="Year", y = expression(paste("ANPP ", g/m^{2})))+
  theme(legend.position="none")+
  annotate("text", x= c("2015", "2016","2017"), y = c(900, 975, 975), label = c("a", "b", "b"), color = "black", cex=6) 
Fig3b

data2$clust <- factor(data2$clust, levels = c("1", "2"), 
                                labels = c("a) Cluster 1: Functionally Diverse", "b) Cluster 2: Avena dominated"))
ann_text <- data.frame(y=c( 950, 950, 950,950),treatment2=c("controlRain", "fallDry","springDry","consistentDry"),label = c("a", "b", "ab", "b"), clust = factor("a) Cluster 1: Functionally Diverse",levels = c("a) Cluster 1: Functionally Diverse", "b) Cluster 2: Avena dominated")))
data2$treatment2 <- as.character(data2$treatment)
#Then turn it back into a factor with the levels in the correct order
data2$treatment2 <- factor(data2$treatment2, levels = c("controlRain", "fallDry","springDry","consistentDry"))

Fig3c<-ggplot(data2, aes(y=ANPPgm, x=treatment2,  fill=treatment2))+
  #ggtitle("c) Treatment response by community cluster")+
  geom_boxplot()+
  facet_wrap(~as.factor(clust))+
  scale_x_discrete(labels = c("Control\n ", "Early\nDrought\n", "Late\n Drought\n","Consistent\n Drought\n"))+
  scale_fill_manual(values = c("Black", "gray40","lightgrey", "white"), guide = guide_legend(title = "Treatment"), labels=c("Control", "Early Drought", "Late Drought", "Consistent Drought")) +
  #annotate("text", x= c(0.7, 0.9, 1.1,1.3), y = c( 950, 950, 950,950), label = c("a", "b", "ab", "b"), color = "black", cex=6) +
  labs(x="", y = expression(paste("ANPP ", g/m^{2})))+
  theme(legend.position="none")
Fig3c

Fig3c<-Fig3c+geom_text(data= ann_text,mapping = aes(x = treatment2, y = y, label = label))
Fig3c #add figure 3c to phenology for clusters to make new figure

#create layout for panel
lay2 <- rbind(c(1,1,2,2),
             c(1,1,2,2),
             c(3,3,3,3),
             c(3,3,3,3))
             
grid.arrange(Fig3a, Fig3b,Fig3c, layout_matrix=lay2) #put panel together

#Figure 4
##functional group ANPP
FG_ANPP<-read.csv("ANPP/ANPP_CleanedData/ClimVar_ANPP-func-groups.csv")

#change plots, years, shelter to factors
FG_ANPP[,'plot'] <- as.factor(as.character(FG_ANPP[,'plot']))
FG_ANPP[,'year'] <- as.factor(as.character(FG_ANPP[,'year']))
FG_ANPP[,'shelter'] <- as.factor(as.character(FG_ANPP[,'shelter']))

#rename func groups
FG_ANPP$func[FG_ANPP$func=="N"] <- "N-fixer"
FG_ANPP$func[FG_ANPP$func=="G"] <- "Grass"
FG_ANPP$func[FG_ANPP$func=="F"] <- "Forb"

#create subset with no species manipulations (control community) only
FG_XC<-filter(FG_ANPP, harvest=="Second", subplot=="XC") %>% dplyr::select(-harvest) %>% arrange(treatment,shelterBlock)

#create subset for grass
FG_XC_grass<-filter(FG_ANPP, subplot=='XC', harvest=="Second", func=="Grass") %>% dplyr::select(-harvest) %>% arrange(treatment,shelterBlock)
FG_XC_grass$clust<-clustmem$clust2

m2<-lme(weight_g_m ~treatment, random=~1|shelterBlock/year, FG_XC_grass, na.action=na.exclude)
summary(m2)
anova(m2) #treatment is significant
r.squaredGLMM(m2) #10% of variation explained by fixed effects, 59% by whole model (interannual variation?)
qqnorm(residuals(m2))
qqline(residuals(m2))
shapiro.test(residuals(m2))
#normally distributed, continue
LS2<-lsmeans(m2, ~treatment)
contrast(LS2, "pairwise")

m2a<-lme(weight_g_m ~treatment*as.factor(clust), random=~1|year/shelterBlock, FG_XC_grass, na.action=na.exclude)
summary(m2a)
anova(m2a) #treatment is significant
r.squaredGLMM(m2a) #10% of variation explained by fixed effects, 59% by whole model (interannual variation?)
qqnorm(residuals(m2a))
qqline(residuals(m2a))
shapiro.test(residuals(m2a))
#normally distributed, continue
LS2a<-lsmeans(m2a, ~treatment*clust)
contrast(LS2a, "pairwise")

FG_XC_grass$treatment2 <- as.character(FG_XC_grass$treatment)
#Then turn it back into a factor with the levels in the correct order
FG_XC_grass$treatment2 <- factor(FG_XC_grass$treatment2, levels = c("controlRain",  "fallDry", "springDry","consistentDry"))

Fig4a<-ggplot(FG_XC_grass, aes(x = treatment2, y = weight_g_m)) +
  ylim(0,1000)+
  ggtitle("a)")+
  geom_boxplot(aes(fill = treatment2)) + 
  scale_fill_manual(values = c("Black", "gray40","lightgrey", "white"), guide = guide_legend(title = "Treatment")) +
  #geom_jitter(position=position_jitter(0.2), aes(color=year)) +
  #geom_jitter(position=position_jitter(0.2), aes(color=May_XC$shelterBlock, shape=as.factor(year))) +
  #scale_color_manual(values = c("gray80","gray50", "black"), guide = guide_legend(title = "Year")) +
  theme(legend.position = "none")+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
  scale_x_discrete(labels = c("Control", "Early\nDrought", "Late\nDrought","Consistent\nDrought"))+
  #facet_wrap(~year)+
  annotate("text", x= c("controlRain","consistentDry", "fallDry","springDry"), y = c( 990, 900, 900,900), label = c("a", "b", "b", "ab"), color = "black", cex=6) +
  xlab("Rainfall Treatment") +
  ylab("Grass ANPP (g/m2)")
Fig4a

Fig4b<-ggplot(FG_XC_grass, aes(x = treatment2, y = weight_g_m)) +
  ylim(0,1000)+
  ggtitle("b)")+
  geom_boxplot(aes(fill = as.factor(clust))) + 
  scale_fill_manual(values = c("white", "black"), guide = guide_legend(title = "Cluster"), label=c("Diverse (1)","Avena (2)")) +
  #geom_jitter(position=position_jitter(0.2), aes(color=year)) +
  #geom_jitter(position=position_jitter(0.2), aes(color=May_XC$shelterBlock, shape=as.factor(year))) +
  #scale_color_manual(values = c("gray80","gray50", "black"), guide = guide_legend(title = "Year")) +
  scale_x_discrete(labels = c("Control", "Early\nDrought", "Late\nDrought","Consistent\nDrought"))+
  theme(legend.position = c(0.85,0.85), legend.background = element_blank(),legend.box.background = element_rect(colour = "black"))+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  #facet_wrap(~year)+
  #annotate("text", x= c("controlRain","consistentDry", "fallDry","springDry"), y = c( 990, 900, 900,900), label = c("a", "b", "b", "ab"), color = "black") +
  xlab("Rainfall Treatment") +
  ylab("Grass ANPP (g/m2)")
Fig4b

#add to Fig 2
data2$year<-as.factor(data2$year)
FG_XC_2<-left_join(data2,FG_XC, by=c("treatment", "shelterBlock","year"))
Fig2c<-ggplot(data=subset(FG_XC_2, func!="Grass"), aes(x = as.factor(clust), y = weight_g_m.y)) +
  #ylim(0,100)+
  ggtitle("c)")+
  geom_boxplot(aes(fill = as.factor(func))) + 
  scale_fill_manual(values = c("white", "black"), guide = guide_legend(title = "Functional\nGroup")) +
  theme(legend.position = c(0.8, 0.7),legend.background = element_blank(),legend.box.background = element_rect(colour = "black"))+
  #annotate("text", x= c("controlRain","consistentDry", "fallDry","springDry"), y = c( 990, 900, 900,900), label = c("a", "b", "b", "ab"), color = "black") +
  xlab("Cluster") +
  ylab("ANPP (g/m2)")+
  annotate("text", x= c("1", "2"), y = c(155, 155), label = c("a         a", "b         b"), color = "black")
Fig2c

#make figure 2 for manuscript
grid.arrange(Fig2a,Fig2b,Fig2c, ncol=3)

FG_XC_forb<-filter(FG_ANPP, subplot=='XC', harvest=="Second", func=="Forb") %>% dplyr::select(-harvest)%>% arrange(treatment,shelterBlock)
FG_XC_forb$clust<-clustmem$clust2

m3<-lme(log(weight_g_m+1) ~treatment, random=~1|year/shelterBlock, FG_XC_forb, na.action=na.exclude)
summary(m3)
anova(m3)#no treatment effect
r.squaredGLMM(m3) #<1% of variation explained by fixed effects, 35% by whole model (interannual variation?)
qqnorm(residuals(m3))
qqline(residuals(m3))
shapiro.test(residuals(m3))
#normally distributed, continue
LS3<-lsmeans(m3, ~treatment)
contrast(LS3, "pairwise") #no differences

m3<-lme(weight_g_m ~treatment*as.factor(clust), random=~1|year/shelterBlock, FG_XC_forb, na.action=na.exclude)
summary(m3)
anova(m3)#no treatment effect
r.squaredGLMM(m3) #<1% of variation explained by fixed effects, 35% by whole model (interannual variation?)
qqnorm(residuals(m3))
qqline(residuals(m3))
shapiro.test(residuals(m3))
LS3<-lsmeans(m3, ~treatment*clust)
contrast(LS3, "pairwise") #no differences

FG_XC_forb$treatment2 <- as.character(FG_XC_forb$treatment)
#Then turn it back into a factor with the levels in the correct order
FG_XC_forb$treatment2 <- factor(FG_XC_forb$treatment2, levels = c("controlRain",  "fallDry","springDry", "consistentDry"))

Fig4c<-ggplot(FG_XC_forb, aes(x = treatment2, y = weight_g_m)) + 
  ylim(0,160)+
  ggtitle("c)")+
  geom_boxplot(aes(fill = treatment2)) + 
  scale_fill_manual(values = c("Black", "gray40","lightgrey", "white"), guide = guide_legend(title = "Treatment")) +
  #geom_jitter(position=position_jitter(0.2), aes(color=year)) +
  #geom_jitter(position=position_jitter(0.2), aes(color=May_XC$shelterBlock, shape=as.factor(year))) +
  theme(legend.position = "none")+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
  scale_x_discrete(labels = c("Control", "Early\nDrought", "Late\nDrought","Consistent\nDrought"))+
  #facet_wrap(~year)+
  annotate("text", x= c("controlRain","consistentDry", "fallDry","springDry"), y = c( 155, 155, 155, 155), label = c("a", "a", "a", "a"), color = "black", cex=6) +
  xlab("Rainfall Treatment") +
  ylab("Forb ANPP (g/m2)")
Fig4c

Fig4d<-ggplot(FG_XC_forb, aes(x = treatment2, y = weight_g_m)) + 
  ylim(0,160)+
  ggtitle("d)")+
  geom_boxplot(aes(fill = as.factor(clust))) + 
  scale_fill_manual(values = c("white", "black"), guide = guide_legend(title = "Cluster"), label=c("Diverse (1)","Avena (2)")) +
  #geom_jitter(position=position_jitter(0.2), aes(color=year)) +
  #geom_jitter(position=position_jitter(0.2), aes(color=May_XC$shelterBlock, shape=as.factor(year))) +
  theme(legend.position = "none")+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
  scale_x_discrete(labels = c("Control", "Early\nDrought", "Late\nDrought","Consistent\nDrought"))+
  #facet_wrap(~year)+
  annotate("text", x= c("controlRain","consistentDry", "fallDry","springDry"), y = c( 160, 150, 150, 150), label = c("", "", "", "*"), color = "black", cex=14) +
  xlab("Rainfall Treatment") +
  ylab("Forb ANPP (g/m2)")
Fig4d

FG_XC_nfix<-filter(FG_ANPP, subplot=='XC', harvest=="Second", func=="N-fixer") %>% dplyr::select(-harvest)%>% arrange(treatment,shelterBlock)
FG_XC_nfix$clust<-clustmem$clust2

m4<-lme(log(weight_g_m+1) ~treatment, random=~1|year/shelterBlock, FG_XC_nfix, na.action=na.exclude)
summary(m4)
anova(m4) #treatment is significant
r.squaredGLMM(m4) #28% of variation explained by fixed effects, 37% by whole model (interannual variation?)
qqnorm(residuals(m4))
qqline(residuals(m4))
shapiro.test(residuals(m4))
#normally distributed, continue
LS4<-lsmeans(m4, ~treatment)
contrast(LS4, "pairwise") #consistent drought is different from fall and spring, but not control rain

m4a<-lme(log(weight_g_m+1) ~treatment*as.factor(clust), random=~1|year/shelterBlock, FG_XC_nfix, na.action=na.exclude)
summary(m4a)
anova(m4a) #treatment is significant
r.squaredGLMM(m4a) #36% of variation explained by fixed effects, 44% by whole model (interannual variation?)
qqnorm(residuals(m4a))
qqline(residuals(m4a))
shapiro.test(residuals(m4a))
#normally distributed, continue
LS4a<-lsmeans(m4a, ~treatment*clust)
contrast(LS4a, "pairwise") #consistent drought is different from fall and spring, but not control rain

FG_XC_nfix$treatment2 <- as.character(FG_XC_nfix$treatment)
#Then turn it back into a factor with the levels in the correct order
FG_XC_nfix$treatment2 <- factor(FG_XC_nfix$treatment2, levels = c("controlRain", "fallDry", "springDry","consistentDry"))

Fig4e<-ggplot(FG_XC_nfix, aes(x = treatment2, y = weight_g_m)) +
  ylim(0,150)+
  ggtitle("e)")+
  geom_boxplot(aes(fill = treatment2)) + 
  scale_fill_manual(values = c("Black", "gray40","lightgrey", "white"), guide = guide_legend(title = "Treatment")) +
  #geom_jitter(position=position_jitter(0.2), aes(color=year)) +
  #geom_jitter(position=position_jitter(0.2), aes(color=May_XC$shelterBlock, shape=as.factor(year))) +
  #scale_color_manual(values = c("gray80","gray50", "black"), guide = guide_legend(title = "Year")) +
  scale_x_discrete(labels = c("Control", "Early\nDrought", "Late\nDrought","Consistent\nDrought"))+
  theme(legend.position = "none")+
  #facet_wrap(~year)+
  annotate("text", x= c("controlRain","consistentDry", "fallDry","springDry"), y = c( 149, 149, 149, 149), label = c("ab", "b", "a", "a"), color = "black", cex=6) +
  xlab("Rainfall Treatment") +
  ylab("N-Fixer ANPP g/m2")
Fig4e

Fig4f<-ggplot(FG_XC_nfix, aes(x = treatment2, y = weight_g_m)) + 
  ggtitle("f)")+
  ylim(0,150)+
  geom_boxplot(aes(fill = as.factor(clust))) + 
  scale_fill_manual(values = c("white","black"), guide = guide_legend(title = "Cluster"), label=c("Diverse (1)","Avena (2)")) +
  #geom_jitter(position=position_jitter(0.2), aes(color=year)) +
  #geom_jitter(position=position_jitter(0.2), aes(color=May_XC$shelterBlock, shape=as.factor(year))) +
  #scale_color_manual(values = c("gray80","gray50", "black"), guide = guide_legend(title = "Year")) +
  scale_x_discrete(labels = c("Control", "Early\nDrought", "Late\nDrought","Consistent\nDrought"))+
  theme(legend.position = "none")+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  #facet_wrap(~year)+
  #annotate("text", x= c("controlRain","consistentDry", "fallDry","springDry"), y = c( 175, 175, 175, 175), label = c("ab", "b", "a", "a"), color = "black") +
  xlab("Rainfall Treatment") +
  ylab("N-Fixer ANPP (g/m2)")
Fig4f

#create layout for panel
lay3 <- rbind(c(1,1,2,2,2),
              c(3,3,4,4,4),
              c(5,5,6,6,6))

grid.arrange(Fig4a, Fig4b,Fig4c, Fig4d, Fig4e, Fig4f, layout_matrix=lay3, heights = c(1, 1, 1.1)) #put panel together

#run Phenology_cleaning.R script first
pheno_2015_xc$treatment2 <- factor(pheno_2015_xc$treatment, levels = c("controlRain", "fallDry","springDry","consistentDry"))
Fig5a<-ggplot(data=pheno_2015_xc, aes(x=time, y=meanPG))+
  ggtitle("a) Overall treatment effects on phenology")+
  geom_errorbar(aes(ymax = meanPG+sePG, ymin = meanPG-sePG), width=.25)+
  geom_line(aes(x=time, y=meanPG, group=treatment2))+
  geom_point(aes(fill=treatment2), pch=21, cex=5)+
  scale_fill_manual(values = c("Black", "gray40","lightgrey", "white"), guide = guide_legend(title = "Treatment:"), labels=c("Control", "Early\nDrought", "Late\nDrought","Consistent\nDrought"))+
  labs(x="Time", y="Greenness (%)")+
  theme(legend.position="bottom")
Fig5a

#combine phenology with overall ANPP graph for new figure 3 
Fig3b_v2<-ggplot(data=pheno_2015_xc, aes(x=time, y=meanPG))+
  ggtitle("b) Overall treatment effects on phenology")+
  geom_errorbar(aes(ymax = meanPG+sePG, ymin = meanPG-sePG), width=.25)+
  geom_line(aes(x=time, y=meanPG, group=treatment2))+
  geom_point(aes(fill=treatment2, shape = treatment2), cex=5)+
  scale_fill_manual(values = c("Black", "gray40","lightgrey", "white"), 
                    guide = guide_legend(title = "Treatment:"), 
                    labels=c("Control", "Early\nDrought", "Late\nDrought","Consistent\nDrought"))+
  scale_shape_manual(values=c(21,24,22,23),guide = guide_legend(title = "Treatment:"), #change legend title
                     labels=c("Control", "Early\nDrought", "Late\nDrought", "Consistent\nDrought"))+
  labs(x="Time", y="Greenness (%)")+
  annotate("text", label= paste(c("Treatment: F['3,357']==81.02", "p < 0.001","Time: F['3,357']==145.64", "p < 0.001","Treatment*x*Time: F['9,357']==3.60", "p < 0.001")),
           x=c(3.5,4,3.5,4,3.5,4), y=c(110,110,106,106,102,102), cex=4, hjust = 1, parse=T)+
  theme(legend.position="bottom")
Fig3b_v2

#create layout for panel
lay2 <- rbind(c(1,1,2,2),
              c(1,1,2,2))

grid.arrange(Fig3a, Fig3b_v2, layout_matrix=lay2) #put panel together

pheno_clust_sum<-pheno_clust %>% filter(subplot=="XC") %>%
  group_by(time, clust, treatment) %>%
  summarize(meanPG=mean(Percent.Green), sePG=sd(Percent.Green)/sqrt(length(Percent.Green)))

pheno_clust_sum$clust <- factor(pheno_clust_sum$clust, levels = c("a) Cluster 1: Functionally Diverse", "b) Cluster 2: Avena dominated"), 
                  labels = c("c) Cluster 1: Functionally Diverse", "d) Cluster 2: Avena dominated"))
pheno_clust_sum$treatment2 <- factor(pheno_clust_sum$treatment, levels = c("controlRain", "fallDry","springDry","consistentDry"))

Fig5b<-ggplot(data=pheno_clust_sum, aes(x=time, y=meanPG))+
  facet_wrap(~as.factor(clust), ncol=2)+
  geom_errorbar(aes(ymax = meanPG+sePG, ymin = meanPG-sePG), width=.25)+
  geom_line(aes(x=time, y=meanPG, group=treatment2))+
  geom_point(aes(fill=treatment2, shape=treatment2), cex=5)+
  scale_fill_manual(values = c("Black", "gray40","lightgrey", "white"), 
                    guide = guide_legend(title = "Treatment:"), 
                    labels=c("Control", "Early\nDrought", "Late\nDrought","Consistent\nDrought"))+
  scale_shape_manual(values=c(21,24,22,23),guide = guide_legend(title = "Treatment:"), #change legend title
                     labels=c("Control", "Early\nDrought", "Late\nDrought", "Consistent\nDrought"))+
  labs(x="Time", y="Greenness (%)") +
  theme(legend.position="bottom")
Fig5b

grid.arrange(Fig5a,Fig5b, ncol=2)

##new figure 4
grid.arrange(Fig3c, Fig5b, ncol=1)
treatment_names <- c(
  `controlRain` = "Control",
  `consistentDry` = "Consistent Drought",
  `fallDry` = "Early Drought",
  `springDry` = "Late Drought"
)

#supplementary figure, run Cover_dataanalyses.R first
gf_graphic2 <- gf %>%
  group_by(subplot, func, treatment, year) %>%
  summarize(meancover=mean(cover), secover=sd(cover)/sqrt(length(cover)))
gf$TB<-paste(gf$treatment,gf$shelterBlock,sep=".")
gf$clust15 <- clust_xc$clust[match(gf$TB, clust_xc$TB)]
gf_graphic2$treatment2 <- as.character(gf_graphic2$treatment)
#Then turn it back into a factor with the levels in the correct order
gf_graphic2$treatment2 <- factor(gf_graphic2$treatment2, levels = c("controlRain",  "fallDry", "springDry","consistentDry"))
ggplot(data = subset(gf_graphic2, subplot %in% c("XC")), aes(x=year, y=meancover, shape=func)) + 
  geom_rect(data = subset(gf_graphic2,treatment2 == 'consistentDry'),aes(fill = treatment2),xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf,alpha = 0.3) +
  geom_rect(data = subset(gf_graphic2,treatment2 == 'controlRain'),aes(fill = treatment2),xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf,alpha = 0.3) +
  geom_rect(data = subset(gf_graphic2,treatment2 == 'fallDry'),aes(fill = treatment2),xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf,alpha = 0.3) +
  geom_rect(data = subset(gf_graphic2,treatment2 == 'springDry'),aes(fill = treatment2),xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf,alpha = 0.3) +
  scale_fill_manual(values=c("sienna", "royalblue2",  "peachpuff2", "lightsteelblue1"),
                    guide = guide_legend(title = "Treatment:"), 
                    labels=c("Consistent Drought","Control",  "Early Drought", "Late Drought"))+
  geom_point(cex=3.5) +
  geom_errorbar(aes(ymax = meancover+secover, ymin = meancover-secover), width=.25) + 
  labs(x="Year", y="Percent cover", shape="Functional Group") +
  facet_wrap(~treatment2, ncol=4, labeller=as_labeller(treatment_names))+ 
  theme_bw()+ 
  theme(text = element_text(size=14))+
  #theme(strip.background = element_blank(),strip.text.x = element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

treatment_names2 <- c(
  `controlRain` = "Control",
  `consistentDry` = "Consistent Drought",
  `fallDry` = "Early Drought",
  `springDry` = "Late Drought",
  `1` = "1: Functionally Diverse",
  `2` = "2: Avena dominated"
)

gf_graphic3 <- gf %>% filter(subplot=="XC")%>%
  group_by(clust15, func, treatment, year) %>%
  summarize(meancover=mean(cover), secover=sd(cover)/sqrt(length(cover)))
gf_graphic3$treatment2 <- as.character(gf_graphic3$treatment)
#Then turn it back into a factor with the levels in the correct order
gf_graphic3$treatment2 <- factor(gf_graphic3$treatment2, levels = c("controlRain",  "fallDry", "springDry","consistentDry"))
ggplot(data = gf_graphic3, aes(x=year, y=meancover, shape=func)) + 
  geom_rect(data = subset(gf_graphic3,treatment2 == 'consistentDry'),aes(fill = treatment2),xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf,alpha = 0.3) +
  geom_rect(data = subset(gf_graphic3,treatment2 == 'controlRain'),aes(fill = treatment2),xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf,alpha = 0.3) +
  geom_rect(data = subset(gf_graphic3,treatment2 == 'fallDry'),aes(fill = treatment2),xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf,alpha = 0.3) +
  geom_rect(data = subset(gf_graphic3,treatment2 == 'springDry'),aes(fill = treatment2),xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf,alpha = 0.3) +
  scale_fill_manual(values=c("sienna", "royalblue2",  "peachpuff2", "lightsteelblue1"))+
  geom_point(cex=3.5) +
  geom_errorbar(aes(ymax = meancover+secover, ymin = meancover-secover), width=.25) + 
  labs(x="Year", y="Percent cover", shape="Functional group") +
  scale_fill_manual(values=c("sienna", "royalblue2",  "peachpuff2", "lightsteelblue1"),
                    guide = guide_legend(title = "Treatment:"), 
                    labels=c("Consistent Drought","Control",  "Early Drought", "Late Drought"))+
  facet_wrap(~clust15*treatment2, ncol=4, labeller=as_labeller(treatment_names2))+ 
  theme_bw()+ 
  theme(text = element_text(size=14))+
  #theme(strip.background = element_blank(), strip.text.x = element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
