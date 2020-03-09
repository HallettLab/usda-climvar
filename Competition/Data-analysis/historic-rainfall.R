### Historical projects for Browns Valley ###
library(tidyverse)
library(lubridate)


# set pathway to climvar dropbox competition data folder

# specify dropbox pathway (varies by user -- tweak when share with Caitlin)
if(file.exists("~/Dropbox/Shared/ClimVar/Competition/Data/")){
  # LGS
  datpath <- "~/Dropbox/Shared/ClimVar/Competition/Data/"
  # LMH
}else{
  datpath <- "~/Dropbox/ClimVar/Competition/Data/"
}


## Pull in the prism data and clean
rain <- read_csv(paste0(datpath, "Competition_CleanedData/PRISM_brownsvalley_long.csv"), skip = 10) %>%
  mutate(ppt = `ppt (inches)`*2.54*10) %>%
  separate(Date, c("year", "month")) %>%
  mutate(year = as.numeric(year),
         month = as.numeric(month)) %>%
  mutate(year = ifelse(month == 12 | month == 11 | month == 10 | month == 9, year + 1, year)) %>%
  mutate(season = "Early",
         season = ifelse(month == 2 | month == 3 | month == 4, "Late", season)) %>%
  filter(month != 5, month != 6, month!= 7, month != 8)

## Summarize by year 
## Using 50% as the cutoff 
rainsummary <-  rain %>%
  group_by(year, season) %>%
  summarize(ppt = sum(ppt)) %>%
  spread(season, ppt) %>%
  mutate(Total = Early + Late) 

## check rainfall patterns for site description
last50 <- rainsummary %>%
  filter(year%in%c(1967:2016))

rainsummary <- rainsummary %>%
  mutate(raintype = "fallWet",
         raintype = ifelse(Early < quantile(rainsummary$Early, .5), "fallDry", raintype),
         raintype = ifelse(Late < quantile(rainsummary$Late, .5), "fallWet", raintype),
         raintype = ifelse(Total < quantile(rainsummary$Total, .5), "fallDry", raintype)) 


## Visualize the rainfall scenarios
#pdf("Rainfal history.pdf", width = 10, height = 8)
ggplot(rainsummary, aes(x=year, y=Total))  + geom_line()+
  geom_point(aes(color = raintype), size = 4) + theme_bw() + labs(x="Year", y="Total rainfall (mm)")
#dev.off()