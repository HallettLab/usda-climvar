dat <- dat.cover_with0 %>%
  mutate(treatment=ordered(treatment, levels = c(consistentDry="consistentDry", springDry="springDry", fallDry="fallDry", controlRain="controlRain"))) %>%
  mutate(status = as.character(status))

ggplot(subset(dat, status!="" & status=="native"), aes(x=treatment, y=cover, fill=func)) + geom_boxplot() + facet_wrap(~year)

dat %>%
  filter(status=="")

unique(dat$status)
