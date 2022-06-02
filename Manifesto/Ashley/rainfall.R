### Historical projects for Browns Valley ###
library(tidyverse)
library(lubridate)


# set pathway to climvar dropbox data folder

# specify dropbox pathway (varies by user -- tweak when share with Caitlin)
datpath <- "~/Dropbox/ClimVar/DATA/Met Station CIMIS/"


## Pull in the prism data and clean
rain <- read_csv(paste0(datpath, "cimis_brownsvalley_daily_20190320.csv")) %>%
  separate(Date, c("month", "day", "year")) %>%
  mutate(year = as.numeric(year),
         month = as.numeric(month)) %>%
  mutate(year = ifelse(month == 12 | month == 11 | month == 10 | month == 9, year + 1, year)) %>%
  mutate(season = "Early",
         season = ifelse(month == 2 | month == 3 | month == 4 | month == 5, "Late", season)) %>%
  filter( month != 6, month!= 7, month != 8)

## Summarize by year 
## Using 50% as the cutoff 
rainsummary <-  rain %>%
  group_by(year, season) %>%
  summarize(ppt = sum(`Precip (mm)`, na.rm=TRUE)) #%>%
  
rainsummary2<-rainsummary %>% spread(season, ppt) %>%
  mutate(Total = Early + Late) 

## check rainfall patterns for site description
last50 <- rainsummary2 %>%
  filter(year%in%c(1967:2017))

rainsummary2 <- rainsummary2 %>%
  mutate(raintype = "fallWet",
         raintype = ifelse(Early < quantile(rainsummary2$Early, .5), "fallDry", raintype),
         raintype = ifelse(Late < quantile(rainsummary2$Late, .5), "fallWet", raintype),
         raintype = ifelse(Total < quantile(rainsummary2$Total, .5), "fallDry", raintype)) 


## Visualize the rainfall scenarios
#pdf("Rainfal history.pdf", width = 10, height = 8)
ggplot(rainsummary2, aes(x=year, y=Total))  + geom_line()+
  geom_point(aes(color = raintype), size = 4) + theme_bw() + labs(x="Year", y="Total rainfall (mm)")
#dev.off()

rainsummary$season <- as.character(rainsummary$season)
#Then turn it back into a factor with the levels in the correct order
rainsummary$season <- factor(rainsummary$season, levels = c("Late",  "Early"))


ggplot(subset(rainsummary, year=="2015"|year=="2016"|year=="2017"), aes(x=year, y=ppt, fill=season))+
  geom_bar( stat="identity")+
  ggtitle("a) Total water-year precipitation")+
  theme_bw()+
  theme(legend.position = c(0.2, 0.8), strip.background = element_blank(), 
        text = element_text(size = 20), 
        strip.text.x = element_text(size = 16, face = "italic"), strip.text.y = element_text(size = 16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  labs(y = "Total Precipitation (mm)", x="Year") +
  scale_fill_manual(values = c( "gray40","lightgrey"), guide = guide_legend(title = "Season: "))

ggplot(subset(rainsummary, year=="2015"|year=="2016"), aes(x=as.factor(year), y=ppt, fill=season))+
  geom_bar( stat="identity")+
  #ggtitle("a) Total water-year precipitation")+
  theme_bw()+
  theme(legend.position = c(0.2, 0.8), strip.background = element_blank(), 
        text = element_text(size = 20), 
        strip.text.x = element_text(size = 16, face = "italic"), strip.text.y = element_text(size = 16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  labs(y = "Total Precipitation (mm)", x="Year") +
  scale_fill_manual(values = c( "gray40","lightgrey"), guide = guide_legend(title = "Season: "))
