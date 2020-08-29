## Some preliminary scripts analyzing pcode4meter growth under different competition scenarios

library(tidyverse)
library(readr)

background <- read_csv("~/Dropbox/ClimVar/Competition/Data/Competition_CleanedData/old/ClimVar_Comp_background-biomass-2.csv")
comp.dat <- read_csv("~/Dropbox/ClimVar/Competition/Data/Competition_CleanedData/ClimVar_Comp_combined-biomass.csv")



# take the average across blocks
comp.dat2 <- select(comp.dat, plot:insitu_pdisturbed) %>% # remove duplicated rows due to different prediction methods for pcode4 fit
  distinct() %>%
  group_by(background, bdensity, pcode4, falltreatment) %>%
  # filter(!is.na(disturbance)) %>%
  summarize(mean_weight = mean(p.ind.wgt.g, na.rm = T))


# some graphs
ggplot(background, aes(x=falltreatment, y=density)) + geom_boxplot() + facet_grid(backgrounddensity~backgroundspp)

ggplot(comp.dat2, aes(x=bdensity, y = mean_weight, color = falltreatment, group = falltreatment)) + geom_point() +
  geom_line() + 
  facet_grid(pcode4~background, scales = "free")

ggplot(subset(comp.dat2, !is.na(bdensity)), aes(x=falltreatment, y = mean_weight, color = background, group = background)) + geom_point() +
  geom_line() + 
  facet_grid(pcode4~bdensity, scales = "free")
ggsave("~/Dropbox/ClimVar/Competition/Figures/Exploratory/all-competition.pdf", width = 8, height = 10)


grass.dat <- comp.dat2 %>%
  filter(grepl("Av|Bro|Vulpia", background),
         pcode4%in%c("AVFA", "BRHO", "VUMY"))

ggplot(subset(grass.dat, !is.na(bdensity)), aes(x=falltreatment, y = mean_weight, color = background, group = background)) + geom_point() +
  geom_line() + 
  facet_grid(pcode4~bdensity, scales = "free")
ggsave("~/Dropbox/ClimVar/Competition/Figures/Exploratory/grass-competition.pdf", width = 8, height = 10)


