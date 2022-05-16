library(tidyverse)


# Accumulation curves for BNPP analysis -----------------------------------

#read in data
bnppdat <- read.csv("./Dropbox/ClimVar/DATA/Plant_composition_data/BNPP/BNPP_EnteredData/ClimVar_BNPP_20150520.csv") %>%
  as_tibble() %>%
  mutate(uniquePlot=paste(Plot, Subplot, sep="_")) %>%
  mutate(uniqueSample=paste(uniquePlot, Depth, sep="_"))
names(bnppdat)[14]="rmass"

#subset data for accumulation
dat <- bnppdat %>%
  #select(uniqueSample, Interval, rmass) %>%
  arrange(uniqueSample, Interval)

#create dat2 with accumulated mass
dat2 <- as.data.frame(cbind(uniqueSample=as.character(), Interval=as.numeric(), rmass=as.numeric(),
                          allrmass=as.numeric()))

uniqueSamples<-unique(dat$uniqueSample)
for (i in 1: length(uniqueSamples)) {
  subber <- subset(dat, uniqueSample==uniqueSamples[i])
  subber$rmass2<-cumsum(subber$rmass)
 dat2 <- rbind(dat2, subber)
}

#add rmass2 (the cummulative mass) back to the full dataset
dat3 <- merge(bnppdat, dat2)



#A function to return fitted values for each interval based on a logistic growth model             
Logmodel<-function(mydat){
  wilson<-nls(rmass2~a*(log(Interval*10)) + b, start=list(a=.5, b=.5), data=mydat,trace=TRUE, control = list(maxiter = 500))
  a<-coef(wilson)[1]
  b<-coef(wilson)[2]
  phi3<-coef(wilson)[3]
  #Interval<-c(min(mydat$Interval):max(mydat$Interval)) #construct a range of x values bounded by the data
  Interval<-c(1:9)
  predictmass<-a*log(Interval*10) + b #predicted count
  predict<-data.frame(Interval,predictmass) #create the prediction data frame#And add a nice plot (I cheated and added the awesome inset jpg in another program)
  predict$uniqueSample<-unique(mydat$uniqueSample)
  predict$a<-a
  predict$b<-b
  return(predict)
}

# Create a dataframe with fit values using all Interval data
logdata <- data.frame(Interval=numeric(), predictmass=numeric(), uniqueSample=character(), a=numeric(), b=numeric())
uniqueSamples<-unique(dat3$uniqueSample)
for (i in 1:length(uniqueSamples)){
  subber<-subset(dat3, uniqueSample==uniqueSamples[i])
  subout<-Logmodel(subber)
  logdata<-rbind(logdata, subout)
}

# Create a dataframe with fit values using a subset of Interval data
logdatamin <- data.frame(Interval=numeric(), predictmass=numeric(), uniqueSample=character(), a=numeric(), b=numeric())
uniqueSamples<-unique(dat3$uniqueSample)
for (i in 1:length(uniqueSamples)){
  subber<-subset(dat3, uniqueSample==uniqueSamples[i] & Interval<5)
  subout<-Logmodel(subber)
  logdatamin<-rbind(logdatamin, subout)
}

#rename logdatamin columns
str(logdatamin)
names(logdatamin)[2:5]=c("minpredictmass", "uniqueSample", "mina", "minb")

#merge everything together
dat4 <- merge(dat3, logdata)
dat5 <- merge(dat4, logdatamin) %>%
  group_by(uniqueSample) %>%
  mutate(totbmass=max(rmass2)) %>%
  mutate(pertot=rmass/totbmass*100)

proptot_all <- ggplot(subset(dat5, Interval>0), aes(x=Interval*10, y=pertot)) + geom_point() + theme_bw() +
  labs(x="Time interval (min)", y="Percentage of total root biomass") +   facet_grid(Depth~uniquePlot, scale="fixed") 

proptot_no1 <- ggplot(subset(dat5, Interval>1), aes(x=Interval*10, y=pertot)) + geom_point() + theme_bw() +
  labs(x="Time interval (min)", y="Percentage of total root biomass") +   facet_grid(Depth~uniquePlot, scale="fixed") 


# pdf("BNPP_Percentoftotal_MarchHarvest.pdf")
# proptot_all
# proptot_no1
# dev.off()

#Plot it
myplot2 <- ggplot((dat5), aes(x=Interval*10, y=rmass2, color=uniquePlot)) + geom_point() + 
  facet_wrap(~Depth, scale="fixed") + 
  theme_bw() +
  labs(x="Time interval (min)", y="Root biomass (g)", color="Plot") +
  geom_line(aes(x=Interval*10, y=predictmass), size=.5)  + 
  geom_line(aes(x=Interval*10, y=minpredictmass), size=.5, linetype="dotted") 
  
#Plot it
myplot3 <- ggplot((dat5), aes(x=Interval*10, y=rmass2, color=uniquePlot)) + geom_point() + 
  facet_wrap(uniquePlot~Depth, scale="free") + 
  theme_bw() +
  labs(x="Time interval (min)", y="Root biomass (g)", color="Plot") +
  geom_line(aes(x=Interval*10, y=predictmass), size=.5)  + 
  geom_line(aes(x=Interval*10, y=minpredictmass), size=.5, linetype="dotted") 

# pdf("BNPP_calibration_fixedscales_MayHarvest.pdf")
# myplot2
# dev.off()

#pdf("BNPP_calibration_MarchHarvest.pdf")
#myplot+ggtitle("30 minutes effort")
#myplot2 + ggtitle("40 minutes effort")
#myplot3 + ggtitle("50 minutes effort")
#dev.off()












