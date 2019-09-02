# explore weather with climvar phenology and soil moisture
# created: mar 2019
# contact: caitlin.t.white@colorado.edu

# script purpose:
# compare phenology with SFREC weather in prep for compost project
# ctw also curious to compare CIMIS weather against soil moisture data and run storm compilation script on precip dat (but TBD later..)

# note: daily weather data from CIMIS weather station at browns valley
# > In the future, stream weather data from CIMIS API. ctw too time crunched for now so using static downloaded dataset in script..
# > CIMIS qa flags described here: https://cimis.water.ca.gov/Content/PDF/CurrentFlags2.pdf


# -- SETUP -----
rm(list = ls())
library(tidyverse)
library(zoo)
library(lubridate)
library(cowplot)
options(stringsAsFactors = F)
theme_set(theme_bw())
na_vals <- c("" , " ", "NA", NA)

# set path to data in Climvar Dropbox
datpath <- "~/DropBox/ClimVar/DATA/"

#read in data
#daily weather
weather <- read.csv(paste0(datpath, "Met\ Station\ CIMIS/cimis_brownsvalley_daily_20190320.csv"),
                    na.strings = na_vals)

# phenology
pheno <- read.csv(paste0(datpath, "Plant_composition_data/Phenology/Phenology_CleanedData/ClimVar_Phenology_clean.csv"),
                  na.strings = na_vals)
#shelterkey
treatments <- read.csv(paste0(datpath, "Plant_composition_data/Shelter_key.csv"))

# -- DATA PREP ----
# 1) daily weather
str(weather)
summary(weather)
#what are the unique flag codes?
lapply(weather[grepl("qc", colnames(weather))], unique) #M = missing; R = outside 98% CI; Y = outside 96% CI; H = hourly value was flagged; A or E = infill with historical average; P = pending
#how many of each QA flag?
weather[grepl("qc", colnames(weather))] %>%
  mutate_all(as.factor) %>%
  summary
#which precip dates are missing?
weather$Date[is.na(weather$Precip..mm.)] #"3/5/2016", "10/9/2017", "11/14/2017", NA
#which date is NA?
weather[is.na(weather$Date),] #.. all NAs, can remove (also is the last row)
weather <- weather[!is.na(weather$Date),]

# i'm going to be lazy and ignore QA flags bc this isn't for analysis..
weather$Date <- as.Date(weather$Date, format = "%m/%d/%Y")
weather$yr <- lubridate::year(weather$Date)

#quick plots..
#precip
qplot(x = Date, y = `Precip..mm.`, data = weather)
#plot all vars
weather[!grepl("qc", colnames(weather))] %>%
  gather(var, val, `ETo..mm.`:`Avg.Soil.Temp..C.`) %>%
  ggplot(aes(Date, val)) +
  geom_point(alpha = 0.5) +
  #stat_summary(fun.y = "mean", geom = "line", col = "skyblue")+
  facet_wrap(~var, scales = "free_y") 
# there's some weird min temps at 0 (flat line).. broken sensor?


# 2) phenology
# check data structure
str(pheno)
# change date from character to date
pheno$date <- as.Date(pheno$date, format = "%Y-%m-%d")


# -- VISUALIZE WEATHER AND PHENO -----
# plot air temp, soil temp, min air temp (check frost), precip, and solar radiation with pct green

tidy_weather <- weather[!grepl("qc", colnames(weather))] %>%
  # create water year (Oct 1 - Sep 30) 
  mutate(WY = ifelse(month(Date) %in% 10:12, yr+1, yr)) %>%
  # create water year day of year sequence (won't be technically right for 2013 since dataset only starts Aug 1 2013, but still in order)
  group_by(WY) %>%
  mutate(doWY = seq(1:length(Date))) %>%
  # correct doWY for 2013 dates (dataset begins at Aug 1 2013)
  mutate(doWY = ifelse(WY == 2013, (doWY-max(doWY))+365, doWY)) %>%
  ungroup() %>%
  # tidy weather variables (i.e. make long form)
  gather(var, val, `ETo..mm.`:`Avg.Soil.Temp..C.`) %>%
  filter(grepl("Preci|Avg[.]Air|Min[.]Air|Soil|Sol", var)) %>%
  mutate(var = gsub("[.][.]",".",var)) %>%
  group_by(WY,var) %>%
  # add rolling mean for panel plot
  mutate(rollmean_7day = zoo::rollmean(val, 7, fill = NA))

weather_plot <- ggplot() +
  geom_point(data = subset(tidy_weather, WY >2013 & !grepl("Precip", var)), aes(doWY, val), alpha = 0.4, col = "dodgerblue3") +
  geom_line(data = subset(tidy_weather, WY>2013 & !grepl("Precip", var)), 
               aes(x = doWY, y = rollmean_7day), col = "black", lwd = 1) +
  geom_line(data = subset(tidy_weather, WY >2013 & grepl("Precip", var)), aes(doWY, val)) +
  labs(x = "Day of water year (e.g. Oct 1 = 1, Sep 30 = 365)") +
  scale_x_continuous(limits = c(0,366), expand = c(0,0)) +
  facet_grid(var~WY, scales = "free_y")


pheno_plot <- subset(pheno, plot != "All" & subplot != "All") %>%
  merge(treatments) %>%
  #add in water year and day of WY to plot with weather data
  full_join(unique(tidy_weather[c("Date", "WY", "doWY")]), by = c("date" = "Date")) %>%
  gather(var, val, percent_green:percent_bare) %>%
  filter(var == "percent_green"& WY !=2013) %>%
  # assign pct green for 2017 photo dates (general estimates based on greenness in photos..)
  # not splitting up estimates by drought treatment type
  mutate(val = ifelse(date == "2017-04-13", 100,
                      ifelse(date == "2017-05-10", 65,
                             ifelse(date == "2017-05-26", 10, val)))) %>%
  ggplot(aes(doWY, val)) +
  geom_point(aes(col = treatment), alpha = 0.5) +
  stat_summary(aes(label = gsub("[0-9]{4}-", "",date)), fun.y = "mean", geom = "text", size = 3) +
  scale_x_continuous(limits = c(0,366), expand = c(0,0))+
  labs(title = "ClimVar phenology with SFREC CIMIS weather data, date labels at mean % green for given date, spring 2017 estimated from photos",
       subtitle = paste("Black line = 7-day rolling mean for each CIMIS variable except daily precip; date range: 2013-10-01 to", max(weather$Date))) +
  theme(legend.position = c(0.08,0.5),
        legend.title = element_text(size = 10),
        legend.background = element_rect(color="black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(size = 12),
        plot.subtitle = element_text(size = 10, face = "italic")) +
  facet_grid(var~WY, scales = "free_y")

# stack weather and phenology figs
weather_pheno <- plot_grid(pheno_plot, weather_plot,
          nrow = 2,
          align = "v",
          rel_heights = c(0.5,1))    

# write out to desktop for now.. (not sure where best to save)
ggsave("~/Desktop/weather_pheno.pdf",
       weather_pheno,
       width = 8, height = 6,
       scale = 1.4)  
