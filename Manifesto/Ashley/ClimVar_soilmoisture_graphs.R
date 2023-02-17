source("ClimVar_Decagon_cleaning.R")
library(lubridate)

dat3 <- dat2 %>%
  filter(year == 2014 | year == 2015) %>%
  mutate(date = ymd_hms(date)) %>%
  mutate(season = ifelse(month%in%c(9:12, 1), "early", "late")) %>%
  mutate(complete_date = (paste(year, month, day, sep="-")),
         complete_date = ymd(complete_date)) %>%
  filter(month%in%c(9:12, 1:5)) %>%
  mutate(toremove = ifelse(plot == 6 & subplot == "G", 1, 0)) %>%
  filter(date < '2015-05-10 18:00:00',
         toremove == 0,
         date > '2014-11-12 22:00:00') %>%
  group_by(treatment, plot, subplot, shelterBlock, complete_date, season) %>%
  summarize(sm = mean(sm))
  
ggplot(dat3, aes(x=complete_date, y=sm, 
                 group = interaction(plot, subplot),
                 color = shelterBlock)) + facet_wrap(~treatment) + geom_line()

dat4 <- dat3 %>%
  tbl_df() %>%
  mutate(treatment=ordered(treatment, levels = c( consistentDry="consistentDry", fallDry="fallDry",springDry="springDry", controlRain="controlRain"))) %>%
  group_by(treatment, complete_date) %>%
  summarize(sm = mean(sm, na.rm=T)) 
  

ggplot(dat4, aes(x=complete_date, y=sm, color = treatment))  + geom_line(lwd = 1.5) +
  theme_bw() + 
  theme(legend.position = "none", strip.background = element_blank(), 
        text = element_text(size = 20), 
        strip.text.x = element_text(size = 16, face = "italic"), strip.text.y = element_text(size = 16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  labs(y = expression(paste("Soil moisture (m"^3, " m"^-3,")")), x = "Day of 2014-2015 growing season") +
  scale_color_manual(values = c( "indianred4", "indianred",  "dodgerblue1", "dodgerblue4"))
ggsave("FigS0_sm.pdf", width = 8, height = 5)



## JUST THE LAST YEAR

dat3 <- dat2 %>%
  filter(year == 2016 | year == 2017) %>%
  mutate(date = ymd_hms(date)) %>%
  mutate(season = ifelse(month%in%c(9:12, 1), "early", "late")) %>%
  mutate(complete_date = (paste(year, month, day, sep="-")),
         complete_date = ymd(complete_date)) %>%
  filter(month%in%c(9:12, 1:5)) %>%
  mutate(toremove = ifelse(plot == 6 & subplot == "G", 1, 0)) %>%
  filter(date < '2017-05-10 18:00:00',
         toremove == 0,
         date > '2016-11-12 22:00:00') %>%
  group_by(treatment, plot, subplot, shelterBlock, complete_date, season) %>%
  summarize(sm = mean(sm))

ggplot(dat3, aes(x=complete_date, y=sm, 
                 group = interaction(plot, subplot),
                 color = shelterBlock)) + facet_wrap(~treatment) + geom_line()

dat4 <- dat3 %>%
  tbl_df() %>%
  mutate(treatment2 = ifelse(treatment == "consistentDry", "fallDry", treatment),
         treatment2 = ifelse(treatment == "springDry", "controlRain", treatment2)) %>%
 # mutate(treatment=ordered(treatment, levels = c( consistentDry="consistentDry", fallDry="fallDry",springDry="springDry", controlRain="controlRain"))) %>%
  group_by(treatment2, complete_date) %>%
  summarize(sm = mean(sm, na.rm=T)) 


ggplot(dat4, aes(x=complete_date, y=sm, color = treatment2))  + geom_line(lwd = 1.5) +
  theme_bw() + 
  theme(legend.position = "none", strip.background = element_blank(), 
        text = element_text(size = 20), 
        strip.text.x = element_text(size = 16, face = "italic"), strip.text.y = element_text(size = 16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  labs(y = expression(paste("Soil moisture (m"^3, " m"^-3,")")), x = "Day of 2016-2017 growing season") +
  scale_color_manual(values = c( "dodgerblue4",  "indianred"))
ggsave("Soilmoisture_201617.pdf", width = 8, height = 5)
