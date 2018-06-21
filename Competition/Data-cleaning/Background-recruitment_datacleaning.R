require(gdata)
library(tidyverse)
library(nlme)

# Read in data
dat.background <-read.xls("~/Dropbox/ClimVar/Competition/Data/Competition_EnteredData/Competition_background_fallrecruitment_20161229.xlsx", sheet=2, header=T, na.strings="#N/A!")
names(dat.background) = c("page", "recorder", "date", "plot", "subplot", "subplotname", "spptreatment",
                        "position", "todelete", "density5x5cm", "conversion", "notes", "QA", "todelete2")

# Convert for area
dat.background2 <- dat.background %>%
  tbl_df() %>%
  filter(spptreatment != "x Control x") %>%
  dplyr::select(plot, spptreatment, density5x5cm, conversion) %>%
  mutate(density = as.numeric(as.character(density5x5cm)),
         ## density should be multiplied by 100 for the full plot; but because the 
         ## seeds were concentrated in the center i'm estimating that it was in fact a smaller area 
         density = ifelse(conversion != "WP (0.5x0.5m)", density*25, density)) %>%
  separate(spptreatment, c("backgroundspp", "backgrounddensity"), sep = "_") 

# Read in keys  
shelter.key <- read.csv("~/Dropbox/ClimVar/Competition/Data/Competition_EnteredData/Shelter_key.csv") %>%
  tbl_df()

seed.key <- read.csv("~/Dropbox/ClimVar/Competition/Data/Competition_EnteredData/CompExpt_seedingKey.csv") %>%
  tbl_df()

# Join it all
dat.backclean <-left_join(shelter.key, dat.background2) %>%
  mutate(falltreatment = "wet",
         falltreatment = ifelse(treatment == "fallDry" | treatment == "consistentDry", "dry", falltreatment))

dat.backclean <-left_join(dat.backclean, seed.key) %>%
  mutate(perRecruit = density/seedsAdded)

write.csv(dat.backclean, "~/Dropbox/ClimVar/Competition/Data/Competition_CleanedData/ClimVar_Comp_background-spp-recruitment.csv")


# Quick visual
ggplot((dat.backclean), aes(x=backgrounddensity, y = perRecruit, fill=treatment)) + geom_boxplot() + 
  facet_wrap( ~backgroundspp, scales = "free") #+ scale_y_log10()


ggplot((dat.backclean), aes(x=backgrounddensity, y = perRecruit, fill=falltreatment)) + geom_boxplot() + 
  facet_wrap( ~backgroundspp, scales = "free") #+ scale_y_log10()

## Quick analysis of recruitment patterns
# Lasthenia
l <- lme(perRecruit ~ backgrounddensity + falltreatment, random = ~1|shelterBlock, data = subset(dat.backclean, backgroundspp == "Lasthenia"))
summary(l)

l <- lme(perRecruit ~ backgrounddensity*falltreatment, random = ~1|shelterBlock, data = subset(dat.backclean, backgroundspp == "Lasthenia"))
summary(l)

# m1<-lmer(perRecruit ~ backgrounddensity +  falltreatment + (1|shelterBlock), data=subset(dat.backclean, backgroundspp == "Lasthenia"),
#          contrasts=list(backgrounddensity = contr.treatment, falltreatment=contr.treatment)) #contrasts must be treatment for this
# summary(m1)

# Avena
l <- lme(perRecruit ~ backgrounddensity + falltreatment, random = ~1|shelterBlock, data = subset(dat.backclean, backgroundspp == "Avena"))
summary(l)

l <- lme(perRecruit ~ backgrounddensity *falltreatment, random = ~1|shelterBlock, data = subset(dat.backclean, backgroundspp == "Avena"))
summary(l)


# Bromus
l <- lme(perRecruit ~ backgrounddensity + falltreatment, random = ~1|shelterBlock, data = subset(dat.backclean, backgroundspp == "Bromus"))
summary(l)

l <- lme(perRecruit ~ backgrounddensity *falltreatment, random = ~1|shelterBlock, data = subset(dat.backclean, backgroundspp == "Bromus"))
summary(l)


# Escholtzia
l <- lme(perRecruit ~ backgrounddensity + falltreatment, random = ~1|shelterBlock, data = subset(dat.backclean, backgroundspp == "Eschscholzia"))
summary(l)

l <- lme(perRecruit ~ backgrounddensity *falltreatment, random = ~1|shelterBlock, data = subset(dat.backclean, backgroundspp == "Eschscholzia"))
summary(l)

# Trifolium
l <- lme(perRecruit ~ backgrounddensity + falltreatment, random = ~1|shelterBlock, data = subset(dat.backclean, backgroundspp == "Trifolium"))
summary(l)

l <- lme(perRecruit ~ backgrounddensity *falltreatment, random = ~1|shelterBlock, data = subset(dat.backclean, backgroundspp == "Trifolium"))
summary(l)


# Vulpia
l <- lme(perRecruit ~ backgrounddensity + falltreatment, random = ~1|shelterBlock, data = subset(dat.backclean, backgroundspp == "Vulpia"))
summary(l)

l <- lme(perRecruit ~ backgrounddensity *falltreatment, random = ~1|shelterBlock, data = subset(dat.backclean, backgroundspp == "Vulpia"))
summary(l)
