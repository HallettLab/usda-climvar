source("./Manifesto/Ashley/ClimVar_Decagon_cleaning.R")
library(lubridate)

dat3 <- dat2 %>%
  filter(year == 2016 | year == 2017) %>%
  mutate(date = ymd_hms(date)) %>%
  mutate(season = ifelse(month%in%c(9:12, 1), "early", "late")) %>%
  mutate(complete_date = (paste(year, month, day, sep="-")),
         complete_date = ymd(complete_date)) %>%
  filter(month%in%c(9:12, 1:5)) %>%
  #mutate(toremove = ifelse(plot == 6 & subplot == "G", 1, 0)) %>%
  filter(date < '2017-05-10 18:00:00',
         subplot != 'C',
         date > '2016-10-10 22:00:00') %>%
  group_by(treatment, plot, subplot, shelterBlock, complete_date, season) %>%
  summarize(sm = mean(sm))
  
ggplot(dat3, aes(x=complete_date, y=sm, 
                 group = interaction(plot, subplot),
                 color = shelterBlock)) + facet_wrap(~treatment) + geom_line()

dat4 <- dat3 %>%
  tbl_df() %>%
  #mutate(falltreatment = ifelse(treatment =="consistentDry" | treatment == "fallDry", "dry", "wet")) %>%

  mutate(treatment=ordered(treatment, levels = c(controlRain="controlRain",  fallDry="fallDry", 
                                                 springDry="springDry", consistentDry="consistentDry")))%>%
  
  group_by(treatment, complete_date) %>%
  summarize(sm = mean(sm, na.rm=T)) 
  

p3<-ggplot(dat4, aes(x=complete_date, y=sm, color = treatment))  + geom_line(lwd = 1) +
  theme_bw() + 
  ylim(0,0.57)+
  scale_x_date(limits = as.Date(c("2016-10-11","2017-05-10")))+
  scale_color_manual(values = c("royalblue2", "peachpuff", "lightsteelblue1","sienna"), 
                     guide = guide_legend(title = "Treatment"),
                     labels=c("Control", "Early Drought", "Late Drought","Consistent Drought" )) +
  theme(strip.background = element_blank(), 
        text = element_text(size = 12), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  geom_segment(aes(x = as.Date("2016-10-11"), y = 0.55, xend = as.Date("2017-05-10"), yend = 0.55), color="sienna", size=5.5)+
  geom_segment(aes(x = as.Date("2016-10-11"), y = 0.53, xend = as.Date("2017-02-01"), yend = 0.53), color="peachpuff", size=5)+
  geom_segment(aes(x = as.Date("2017-02-01"), y = 0.53, xend = as.Date("2017-05-10"), yend = 0.53), color="lightsteelblue1", size=5)+
    labs(y = expression(paste("Soil moisture (m"^3, " m"^-3,")")), x = "Day of 2016-2017 growing season") 
  #scale_color_manual(values = c( "indianred3", "dodgerblue2"))
#ggsave("Competition/Figures/soil-moisture_201617.pdf", width = 8, height = 5)
p3

dat5 <- dat2 %>%
  filter(year == 2015 | year == 2016) %>%
  mutate(date = ymd_hms(date)) %>%
  mutate(season = ifelse(month%in%c(9:12, 1), "early", "late")) %>%
  mutate(complete_date = (paste(year, month, day, sep="-")),
         complete_date = ymd(complete_date)) %>%
  filter(month%in%c(9:12, 1:5)) %>%
  #mutate(toremove = ifelse(plot == 6 & subplot == "G", 1, 0)) %>%
  filter(date < '2016-05-10 18:00:00',
         subplot != 'C',
         date > '2015-10-10 22:00:00') %>%
  group_by(treatment, plot, subplot, shelterBlock, complete_date, season) %>%
  summarize(sm = mean(sm))

ggplot(dat5, aes(x=complete_date, y=sm, 
                 group = interaction(plot, subplot),
                 color = shelterBlock)) + facet_wrap(~treatment) + geom_line()

dat6 <- dat5 %>%
  tbl_df() %>%
  #mutate(falltreatment = ifelse(treatment =="consistentDry" | treatment == "fallDry", "dry", "wet")) %>%
  mutate(treatment=ordered(treatment, levels = c( consistentDry="consistentDry", fallDry="fallDry",
                                                   springDry="springDry", controlRain="controlRain"))) %>%
  group_by(treatment, complete_date) %>%
  summarize(sm = mean(sm, na.rm=T)) 


p2<-ggplot(dat6, aes(x=complete_date, y=sm, color = treatment))  + geom_line(lwd = 1) +
  theme_bw() + 
  scale_x_date(limits = as.Date(c("2015-10-11","2016-05-10")))+
  ylim(0,0.57)+
  scale_color_manual(values = c("sienna", "peachpuff", "lightsteelblue1","royalblue2"), 
                     guide = guide_legend(title = "Treatment"),
                     labels=c("Consistent Drought", "Early Drought", "Late Drought","Control" )) +  
  theme(strip.background = element_blank(), 
        text = element_text(size = 12), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  geom_segment(aes(x = as.Date("2015-10-11"), y = 0.55, xend = as.Date("2016-05-10"), yend = 0.55), color="sienna", size=5.5)+
  geom_segment(aes(x = as.Date("2015-10-11"), y = 0.53, xend = as.Date("2016-02-01"), yend = 0.53), color="peachpuff", size=5)+
  geom_segment(aes(x = as.Date("2016-02-01"), y = 0.53, xend = as.Date("2016-05-10"), yend = 0.53), color="lightsteelblue1", size=5)+
  labs(y = expression(paste("Soil moisture (m"^3, " m"^-3,")")), x = "Day of 2015-2016 growing season") 
  #scale_color_manual(values = c( "indianred3", "dodgerblue2"))
#ggsave("Competition/Figures/soil-moisture_201617.pdf", width = 8, height = 5)
p2

dat7 <- dat2 %>%
  filter(year == 2014 | year == 2015) %>%
  mutate(date = ymd_hms(date)) %>%
  mutate(season = ifelse(month%in%c(9:12, 1), "early", "late")) %>%
  mutate(complete_date = (paste(year, month, day, sep="-")),
         complete_date = ymd(complete_date)) %>%
  filter(month%in%c(9:12, 1:5)) %>%
  #mutate(toremove = ifelse(plot == 6 & subplot == "G", 1, 0)) %>%
  filter(date < '2015-05-10 18:00:00',
         subplot != 'C',
         date > '2014-10-10 22:00:00') %>%
  group_by(treatment, subplot, shelterBlock, season, complete_date) %>%
  summarize(sm = mean(sm))

ggplot(dat7, aes(x=complete_date, y=sm, 
                 #group = interaction(plot, subplot),
                 color = shelterBlock)) + facet_wrap(~treatment) + geom_line()


dat8 <- dat7 %>%
  tbl_df() %>%
  #mutate(falltreatment = ifelse(treatment =="consistentDry" | treatment == "fallDry", "dry", "wet")) %>%
  mutate(treatment=ordered(treatment, levels = c( consistentDry="consistentDry", fallDry="fallDry",
                                                  springDry="springDry", controlRain="controlRain"))) %>%
  group_by(treatment, complete_date) %>%
  summarize(sm = mean(sm, na.rm=T)) 


p1<-ggplot(dat8, aes(x=complete_date, y=sm, color = treatment))  + geom_line(lwd = 1) +
  theme_bw() + 
  ylim(0,0.57)+
  scale_x_date(limits = as.Date(c("2014-10-11","2015-05-10")))+
  scale_color_manual(values = c("sienna", "peachpuff", "lightsteelblue1","royalblue2"), 
                     guide = guide_legend(title = "Treatment"),
                     labels=c("Consistent Drought", "Early Drought", "Late Drought","Control" )) +
  theme(strip.background = element_blank(), 
      text = element_text(size = 12), 
       panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  geom_segment(aes(x = as.Date("2014-10-11"), y = 0.55, xend = as.Date("2015-05-10"), yend = 0.55), color="sienna", size=5.5)+
  geom_segment(aes(x = as.Date("2014-10-11"), y = 0.53, xend = as.Date("2015-02-01"), yend = 0.53), color="peachpuff", size=5)+
  geom_segment(aes(x = as.Date("2015-02-01"), y = 0.53, xend = as.Date("2015-05-10"), yend = 0.53), color="lightsteelblue1", size=5)+
  labs(y = expression(paste("Soil moisture (m"^3, " m"^-3,")")), x = "Day of 2014-2015 growing season") 
  #scale_color_manual(values = c( "indianred3", "dodgerblue2"))
#ggsave("Competition/Figures/soil-moisture_201617.pdf", width = 8, height = 5)
p1

library(gridExtra)
#create layout for panel
lay <- rbind(c(1,1),
            c(1,1),
            c(2,2),
            c(2,2),
            c(3,3),
            c(3,3))
grid.arrange(p1, p2, p3, layout_matrix=lay) #put panel together

library(ggpubr)
ggarrange (p1, p2, p3, labels="auto", common.legend=T, legend="bottom",ncol=1)
