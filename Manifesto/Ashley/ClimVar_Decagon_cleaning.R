require(gdata)
library(tidyr)
library(dplyr)
library(lubridate)
library(ggplot2)

## Set working directory
setwd("~/Dropbox/ClimVar/DATA/Decagon data")

##A FUNCTION TO IMPORT AND CLEAN DECAGON FILES##
cleanDecagon<-function(X){
  mydat<-read.xls(X, sheet=1, header=T, na.strings="#N/A!")
  mydf<-tbl_df(mydat[-c(1,2),])
  mydf$plot<-as.character(names(mydf[1]))
  names(mydf)=c("time", "B", "C",  "F", "G", "XC", "plot")
  mydf<-mydf%>%
    mutate(plot=extract_numeric(plot))%>%
    mutate(date=parse_date_time(time, "m%d%Y I%M p")) %>%
    mutate(year=year(date), month=month(date), day=day(date), hour=hour(date),julian=yday(date))%>%
    mutate(time=julian+(hour/24)) %>%
    mutate(date=as.character(date))
}

##LIST ALL DECAGON FILES##
allfilenames<-list.files()
excelfilenames<-as.matrix(subset(allfilenames, grepl(".xls", allfilenames)==T))
#double check if there is an excel version of the master; if so remove from list
#no master (in separate folder)

#remove Mar 25 2017 files from the list until after we talk to Caitlin about the dates & data
excelfilenames<- excelfilenames %>% subset(excelfilenames[,1]!="PLOT1 25Mar17-0135.xls") 
excelfilenames<- excelfilenames %>%subset(excelfilenames[,1]!="PLOT10 25Mar17-0124.xls")
excelfilenames<- excelfilenames %>% subset(excelfilenames[,1]!="PLOT11 25Mar17-0122.xls")
  excelfilenames<- excelfilenames %>%subset(excelfilenames[,1]!="PLOT12 25Mar17-0120.xls") 
  excelfilenames<- excelfilenames %>%subset(excelfilenames[,1]!="PLOT13 25Mar17-0120.xls") 
    excelfilenames<- excelfilenames %>%subset(excelfilenames[,1]!="PLOT14 25Mar17-0118.xls") 
    excelfilenames<- excelfilenames %>%subset(excelfilenames[,1]!="PLOT15 25Mar17-0117.xls") 
    excelfilenames<- excelfilenames %>%subset(excelfilenames[,1]!="PLOT16 25Mar17-0115.xls") 
    excelfilenames<- excelfilenames %>%subset(excelfilenames[,1]!="PLOT9 25Mar17-0127.xls") 
    excelfilenames<- excelfilenames %>%subset(excelfilenames[,1]!="PLOT8 25Mar17-0128.xls") 
    excelfilenames<- excelfilenames %>%subset(excelfilenames[,1]!="PLOT7 25Mar17-0129.xls") 
    excelfilenames<- excelfilenames %>%subset(excelfilenames[,1]!="PLOT6 25Mar17-0130.xls") 
    excelfilenames<- excelfilenames %>%subset(excelfilenames[,1]!="PLOT5 25Mar17-0131.xls") 
    excelfilenames<- excelfilenames %>%subset(excelfilenames[,1]!="PLOT4 25Mar17-0132.xls") 
    excelfilenames<- excelfilenames %>%subset(excelfilenames[,1]!="PLOT3 25Mar17-0133.xls")
    excelfilenames<- excelfilenames %>%subset(excelfilenames[,1]!="PLOT2 25Mar17-0134.xls")

##CLEAN ALL DECAGON FILES##
myupdates_list<-apply(excelfilenames, 1, cleanDecagon)
#create a single data frame of new files
myupdates_df <- tbl_df(do.call("rbind", myupdates_list))


##FORMAT FOR R##
#tidy myupdates_df
dat<-myupdates_df
dat<-dat%>%
  gather(subplot, sm, B:XC)
key<-read.csv("Shelter_key.csv")
dat2<-tbl_df(merge(key, dat))
dat2$plot<-as.factor(as.character(dat2$plot))
dat2$sm<-as.numeric(dat2$sm)

#check for duplicate records
duplicated(dat2)
duplicates<-dat2[duplicated(dat2),]
duplicates$year<-as.factor(as.character(duplicates$year))
levels(duplicates$year)
#there seems to be a problem with plot 10 from Mar 24 2017 data. For some reason, data coded as year 2000 and 2044. I will remove these & remove duplicates
dat2<- unique(dat2)
#remove year 2000 and 2044 from record
dat2<-dat2 %>% filter(year!="2000")%>%filter(year!="2044")
#check again for duplicates from same plot/time but with dif sm values
duplicates3<-dat2[duplicated(dat2[,c(1:12)]),]
duplicates3$year<-as.factor(as.character(duplicates3$year))
#no dups, let's move on

#pdf("smXsppXtrt.pdf")
biggraph<-ggplot(data=dat2, aes(x=time,
                      y=sm, group=subplot, color=subplot)) + geom_line() + facet_grid(shelterBlock~treatment)
biggraph
#print(biggraph)
#dev.off()

#pdf("smXtrt_controlspp.pdf")
controldat<-subset(dat2, subplot=="XC")

controlgraph<-ggplot(data=controldat, aes(x=time,
                     y=sm, color=treatment, group=treatment)) + geom_line() + facet_grid(~shelterBlock)  + 
  scale_y_continuous(breaks=c(seq(-.4,.8,.1)))
controlgraph
controlAC<-subset(controldat, shelterBlock=="A" | shelterBlock=="C")
controlAC$shelterBlock<-as.character(controlAC$shelterBlock)
controlgraph2<-ggplot(data=controlAC, aes(x=time,
                                          y=sm, color=treatment, group=treatment)) + geom_line() + facet_grid(~shelterBlock)  + 
  scale_y_continuous(breaks=c(seq(-.4,.8,.1)))
#print(controlgraph)
#print(controlgraph2)
#dev.off()

# pull out focal plots, an
smdat2 <- dat2 %>%
  # filter(subplot == "B" | subplot == "G" | subplot == "F" ) %>%
  # filter(subplot == "XC" | subplot == "C") %>%
  #  filter(subplot == "XC" | subplot == "B") %>%
  #filter( subplot == "G" | subplot == "F" ) %>%
  filter(subplot !="C") %>%
  tbl_df() %>%
  mutate(doy3= julian/365,
         doy4 = year + doy3) 

smdat2 %>%
  group_by(treatment, shelterBlock) %>%
  summarize(meansm = mean(sm, na.rm=T))

getOption("device")

ggplot(subset(smdat2), aes(x=doy4, y=sm,  color = subplot, group = interaction(subplot, plot))) +
  geom_line(size = .5) + theme_bw() + facet_wrap(~treatment) # facet_grid(shelterBlock~treatment) 

ggplot(subset(smdat2), aes(x=doy4, y=sm,  color = subplot)) + geom_line() + facet_wrap(~treatment)
ggplot(subset(smdat2, subplot != "B"), aes(x=doy4, y=sm,  color = treatment)) + geom_line() + facet_wrap(~subplot)

# aggregate across treatments
smdat3 <- smdat2 %>%
  group_by(subplot, treatment, year, doy4, doy3) %>%
  summarize(sm = mean(sm, na.rm=T))

dat2%>%
  group_by(treatment, shelterBlock)%>%
  summarise(meansm=mean(sm, na.rm=T))


ggplot(subset(smdat3), aes(x=doy4, y=sm,  color = treatment, group = (treatment))) +
  geom_line(size = .5) + theme_bw() + facet_wrap(~subplot) # facet_grid(shelterBlock~treatment) 

ggplot(subset(smdat2), aes(x=doy4, y=sm,  color = treatment, group = (treatment))) +
  geom_line(size = .5) + theme_bw() + facet_wrap(~plot) # facet_grid(shelterBlock~treatment) 

ggplot(subset(smdat2), aes(x=doy4, y=sm,  color = subplot, group = (subplot))) +
  geom_line(size = .5) + theme_bw() + facet_wrap(subplot~plot) # facet_grid(shelterBlock~treatment) 

#really weird data in plot 10, subplot B - did a sensor malfunction? remove plot 10, subplot B, year 2016 from the dataset
smdat4<-smdat2 %>% filter(plot!="10")
  #mutate_all(.funs = function(x) replace(x, which(x < 0 ), NA))

smdat5 <- smdat4 %>%
  group_by(subplot, treatment, year, doy4, doy3) %>%
  summarize(sm = mean(sm, na.rm=T))

#create a plot showing sm data by treatment
ggplot(subset(smdat5), aes(x=doy4, y=sm,  color = treatment, group = (treatment))) +
  geom_line(size = .5) + theme_bw() # facet_wrap(~subplot) # facet_grid(shelterBlock~treatment) 

smdat4 %>%
  group_by(treatment) %>%
  summarize(meansm = mean(sm, na.rm=T))

#create a new variable for growing season?
smdat4<-smdat4 %>% mutate(season=ifelse(doy4 < 2015.5, "2015", ifelse(doy4 >2015.8 & doy4 < 2016.5, "2016", ifelse(doy4 > 2016.8 & doy4 < 2017.5, "2017", "summer"))))

#summarize by season
sm_season<-smdat4 %>%
  group_by(treatment, season) %>% filter(season!="summer") %>%
  summarize(meansm = mean(sm, na.rm=T))

CV <- function(x){(sd(x)/mean(x))*100}
moistCV<-aggregate(sm ~ treatment*shelterBlock*subplot*season, data= smdat4, FUN = CV) %>% filter(season!="summer")
colnames(moistCV)[colnames(moistCV)=="sm"] <- "sm_cv"
sm<-smdat4%>%group_by(treatment, season,shelterBlock,subplot) %>% filter(season!="summer") %>% summarize(meansm=mean(sm), sesm=(sd(sm)/sqrt(length(sm)))) 
sm<-merge(sm,moistCV)
sm_XC<-filter(sm, subplot=="XC")
May_ANPP_XC<-merge(sm_XC, May_ANPP_XC)

ggplot(May_ANPP_XC, x=treatment, y=meansm, fill=treatment)+
  annotate("text", x= c("consistentDry", "controlRain","fallDry","springDry"), y = c(0.4, 0.47, 0.45,0.45), label = c("a", "b", "c", "c"), color = "black") +
  geom_boxplot(aes(x=treatment, y=meansm, fill = treatment), shape=16)+
  theme_bw()+
  scale_fill_manual(values = c("sienna", "royalblue2","lightsteelblue1", "peachpuff" ), guide = guide_legend(title = "Treatment"))+
  xlab("Rainfall Treatment") +
  ylab("Volumetric soil moisture")

ggplot(May_ANPP_XC, aes(x=reorder(treatment, sm_cv, FUN=mean), y=sm_cv, fill=treatment))+
  annotate("text", x= c("consistentDry", "controlRain","fallDry","springDry"), y = c(0.4, 0.47, 0.45,0.45), label = c("a", "b", "c", "c"), color = "black") +
  geom_boxplot(aes(fill = treatment))+
  theme_bw()+
  scale_fill_manual(values = c("sienna", "royalblue2","lightsteelblue1", "peachpuff" ), guide = guide_legend(title = "Treatment"))+
  xlab("Rainfall Treatment") +
  ylab("Soil moisture CV")

ggplot(May_ANPP_XC, x=shelterBlock, y=meansm )+
  #annotate("text", x= c("consistentDry", "controlRain","fallDry","springDry"), y = c(0.4, 0.47, 0.45,0.45), label = c("a", "b", "c", "c"), color = "black") +
  geom_boxplot(aes(x=shelterBlock, y=meansm), shape=16)+
  theme_bw()+
  scale_fill_manual(values = c("sienna", "royalblue2","lightsteelblue1", "peachpuff" ), guide = guide_legend(title = "Treatment"))+
  xlab("Rainfall Treatment") +
  ylab("Volumetric soil moisture")

ggplot(May_ANPP_XC, aes(x=meansm, y=percentGrass, color=treatment))+
  geom_point()+
  geom_smooth(method="lm", se=F)

sm1<-lme(meansm ~ treatment, random=~1|season/shelterBlock, data=May_ANPP_XC, na.action=na.exclude)
summary(sm1)
anova(sm1)
r.squaredGLMM(sm1) #26% of variation explained by fixed effects, 37% by whole model 
qqnorm(residuals(sm1))
qqline(residuals(sm1))
shapiro.test(residuals(sm1))
#normally distributed, continue
LSsm1<-lsmeans(sm1, ~treatment)
contrast(LSsm1, "pairwise")

sm2<-lme(log(sm_cv+1) ~ treatment, random=~1|season/shelterBlock, data=May_ANPP_XC, na.action=na.exclude)
summary(sm2)
anova(sm2)
r.squaredGLMM(sm2) #26% of variation explained by fixed effects, 37% by whole model 
qqnorm(residuals(sm2))
qqline(residuals(sm2))
shapiro.test(residuals(sm2))
#normally distributed, continue
LSsm2<-lsmeans(sm2, ~treatment)
contrast(LSsm2, "pairwise")


#create a new variable for growing season?
smdat4<-smdat4 %>% mutate(season=ifelse(doy4 < 2015.5, "2015", ifelse(doy4 >2015.8 & doy4 < 2016.5, "2016", ifelse(doy4 > 2016.8 & doy4 < 2017.5, "2017", "summer"))))

#summarize by season
sm_2015<-smdat4 %>% filter(year=="2015")%>%
  group_by(treatment, subplot, shelterBlock) %>% filter(subplot!="C", subplot !="XC") %>%
  summarize(meansm = mean(sm, na.rm=T))

sm_anpp_2015<- merge(sm_2015, May_ANPP_2015)

m4<-lme(weight_g_m ~subplot*meansm, random=~1|shelterBlock, sm_anpp_2015, na.action=na.exclude)
summary(m4)
anova(m4)
r.squaredGLMM(m4)
qqnorm(residuals(m4))
qqline(residuals(m4))
shapiro.test(residuals(m4)) #normal
LS4<-lsmeans(m4, ~subplot*meansm)
contrast(LS4, "pairwise")

ggplot(sm_anpp_2015, aes(x=meansm, y=weight_g_m, color=subplot))+ #color=subplot, shape=AvDom))+
  geom_point()+
  facet_wrap(~AvDom)+
  geom_smooth(method="lm", se=F)
