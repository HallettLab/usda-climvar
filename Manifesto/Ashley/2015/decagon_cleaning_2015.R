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

sm2015<-dat2 %>% filter(year=="2014"|year=="2015")
write.csv(sm2015, file ="~/Dropbox/ClimVar/DATA/Decagon data/ClimVar_sm_2015.csv")

pdf("smXsppXtrt.pdf")
biggraph<-ggplot(data=dat2, aes(x=time,
                                y=sm, group=subplot, color=subplot)) + geom_line() + facet_grid(shelterBlock~treatment)
biggraph
print(biggraph)
dev.off()

pdf("smXtrt_controlspp.pdf")
controldat<-subset(dat2, subplot=="XC")

controlgraph<-ggplot(data=controldat, aes(x=time,
                                          y=sm, color=treatment, group=treatment)) + geom_line() + facet_grid(~shelterBlock)  + 
  scale_y_continuous(breaks=c(seq(-.4,.8,.1)))

controlAC<-subset(controldat, shelterBlock=="A" | shelterBlock=="C")
controlAC$shelterBlock<-as.character(controlAC$shelterBlock)
controlgraph2<-ggplot(data=controlAC, aes(x=time,
                                          y=sm, color=treatment, group=treatment)) + geom_line() + facet_grid(~shelterBlock)  + 
  scale_y_continuous(breaks=c(seq(-.4,.8,.1)))
print(controlgraph)
print(controlgraph2)
dev.off()

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
smdat4<-smdat2 %>% 
  mutate_all(.funs = function(x) replace(x, which(x < 0 ), NA))

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
smdat4<-smdat4 %>% mutate( season=ifelse(doy4 %in% 2014:2015.5, "one", ifelse(doy4 %in% 2015.8:2016.5, "two", ifelse(doy4 %in% 2016.8:2017.5, "three", "summer"))))

CV <- function(x){(sd(x)/mean(x))*100}
moistCV<-aggregate(sm ~ treatment*shelterBlock*subplot*year, data= smdat4, FUN = CV)
colnames(moistCV)[colnames(moistCV)=="sm"] <- "sm_cv"


