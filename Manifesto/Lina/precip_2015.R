library(tidyr)
library(dplyr)
library(lubridate)
library(ggplot2)

#CIMIS Precip data of Browns Valley in 2014-2015
ppt_data <- read.csv("./BNPP/BNPP_CleanedData/CIMISdata.csv")
ppt_data$Date <- mdy(ppt_data$Date)

mydf_ppt <- ppt_data %>%
  mutate(year=year(Date), month=month(Date), day=day(Date), julian=yday(Date))
class(ppt_data$Date)

#Plot Precip of 2014-2015
ggplot(mydf_ppt, aes(x = Date, y = Precip_mm))+
  geom_line()+
  theme_classic()+
  labs(x = "Day of 2014-2015 growing season", y = "Precipitation (mm)")


#soil moisture data in 2014-2015
sm_data <- read.csv("./Dropbox/ClimVar/DATA/Decagon data/ClimVar_sm_2015.csv")

mydf_sm <- sm_data %>%
  mutate_all(.funs = function(x) replace(x, which(x < 0 ), NA)) %>%
  #mutate(date=parse_date_time(time, "m%d%Y I%M p")) %>%
  mutate(year=year(date), month=month(date), day=day(date), hour=hour(date),julian=yday(date))%>%
  group_by(treatment, date) %>%
  summarise(mean = mean(sm, na.rm = TRUE))

ggplot(mydf_sm, aes(x = date, y = mean, col = treatment, group = (treatment)))+
  geom_line()+
  theme_classic()+
  labs(x = "Day of 2014-2015 growing season", y = "Soil moisture (m3 m-3)")
