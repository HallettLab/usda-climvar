library(dplyr)
library(tidyr)
library(ggplot2)

dat.nat <- read.csv("ClimVar_MasterCover_Native_1516.csv") %>%
  tbl_df() %>%
  mutate(treatment=ordered(treatment, levels = c(consistentDry="consistentDry", springDry="springDry", fallDry="fallDry", controlRain="controlRain"))) %>%
  mutate(status = as.character(status))

# Looks like Clarkia, Eschsholzia and Achillea are the main drivers of "natives do better in wet falls" pattern
# Might want to reconsider Mimulus given it's low recruitment
ggplot(subset(dat.nat, status!="" & status=="native"), aes(x=treatment, y=cover, fill=as.factor(year))) + geom_boxplot() + facet_wrap(~species_name, scales = "free")
#ggplot(subset(dat.nat, status!="" & status=="non-native"), aes(x=treatment, y=cover, fill=as.factor(year))) + geom_boxplot() + facet_wrap(~species_name, scales = "free")

dat.all <- read.csv("ClimVar_MasterCover_1516.csv") %>%
  tbl_df() %>%
  mutate(treatment=ordered(treatment, levels = c(consistentDry="consistentDry", springDry="springDry", fallDry="fallDry", controlRain="controlRain"))) %>%
  mutate(status = as.character(status))

ggplot(subset(dat.all, status!="" & subplot!="C"), aes(x=treatment, y=cover, fill=as.factor(year))) + geom_boxplot() + facet_wrap(~species_name, scales = "free")
