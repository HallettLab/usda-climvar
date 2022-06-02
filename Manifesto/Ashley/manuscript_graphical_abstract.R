#the purpose of this script is to modify Figure 1 to fit as a graphical abstract
### Data Import
source("./manuscript_figures.R")

library(cowplot)
Fig1a<-ggplot(spscoresall_2, aes(x=NMDS1, y=NMDS2))+
  xlim(-1,1.3)+
  geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,fill=as.factor(clust2),group=as.factor(clust2)),alpha=0.6) + # add the convex hulls
  #scale_fill_manual(values=c("lightgrey", "black"), guide = guide_legend(title = "Cluster:"), #change legend title
  #labels=c("1", "2")) #change labels in the legend)
  geom_path(data=spscoresall_2, arrow=arrow(), mapping=aes(x=NMDS1, y=NMDS2, group=TB))+
  geom_point(cex=4, aes(shape=as.factor(treatment), fill=treatment))+
  ggtitle("Community composition over time from 2015 to 2017")+
  scale_fill_manual(values=c("gray90", "black","white","Black", "gray40","lightgrey" ), guide = guide_legend(title = "Treatment:"), #change legend title
                    labels=c("Cluster1","Cluster2","Consistent\nDrought","Control", "Early\nDrought", "Late\nDrought"))+ #change labels in the legend)+
  scale_shape_manual(values=c(23,21,24,22),guide = guide_legend(title = "Treatment:"), #change legend title
                     labels=c("Consistent\nDrought", "Control", "Early\nDrought", "Late\nDrought"))+
  theme(legend.position="none")
#theme(legend.position=c(0.7, 0.1), legend.direction="horizontal",legend.key.size = unit(0.1,"cm"),legend.title=element_text(size=11), legend.text=element_text(size=9))
#theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))
Fig1a
Fig1a<-Fig1a+geom_label(data=spc, mapping=aes(x=MDS1, y=MDS2, label=name), position=position_nudge(x=0.2), cex=3)
Fig1a<-Fig1a+annotate("text", y=c(1.2, 1.2), x=c(0.5,-0.7), label=c("Cluster 2","Cluster 1"), cex=6)
Fig1a

Fig1d<- ggplot(data=veclength2, aes(x=as.factor(clust), y=veclength, fill=as.factor(clust), alpha=0.8))+
  ylim(0,1.5)+
  geom_boxplot()+
  #ggtitle("c) Community change by cluster")+
  #theme(legend.position="none")+
  labs(x="", y="Community Change\n(Vector Length)")+
  annotate("text", x= c("1", "2"), y = c(1.5, 1.5), label = c("a", "b"), color = "black")+
  #facet_wrap(~clust)+
  theme(legend.position="none",  text = element_text(size = 12), panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  scale_x_discrete(labels = c("Functionally diverse\nCluster 1 ","Avena dominant\nCluster 2  "))+
  scale_fill_manual(values = c("gray90","black"))
Fig1d

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(Fig1a)
mylegend

graphical_abstract<-ggdraw(Fig1a) +
  draw_plot(Fig1d, .6, 0.06, .38, .3)
graphical_abstract

ggdraw(plot_grid(graphical_abstract, mylegend, ncol=1, rel_heights = c(1, .1)))
