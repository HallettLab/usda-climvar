# -- SETUP -----
rm(list = ls()) # clean environment
library(tidyverse)
options(stringsAsFactors = F)
theme_set(theme_bw())
na_vals = c(" ","", NA, "NA")

# set pathway to climvar dropbox competition data folder
datpath <- "~/Dropbox/ClimVar/Competition/Data/"

# se function
se <-function(x, na.rm = T){
  x <- x[is.na(x)==F]
  sd(x)/sqrt(length(x))
}

# read in data
# cleaned, combined competition dataset (cleaned background, phytometers, predicted vals, plot treatments)
comp.dat <- read.csv(paste0(datpath, "Competition_CleanedData/Competition_combined_clean.csv"),
                     na.strings = na_vals, strip.white = T) %>%
  tbl_df()

# comp.datagg <- comp.dat %>%
#   group_by(falltreatment, background, bdensity, pcode4) %>%
#   summarize(allweight = mean(p_totwgt), indweight = mean(p.ind.wgt.g)) %>%
#   filter(pcode4%in%c("AVFA", "BRHO", "VUMY") & 
#   background%in%c("Avena fatua", "Bromus hordeaceus", "Vulpia myuros"))
# 
# ggplot(comp.datagg, aes(x=falltreatment, y=allweight, color = background, group = background)) + geom_point() + geom_line() +
#   facet_grid(pcode4~bdensity, scales = "free")
# 
# 
# ggplot(comp.datagg, aes(x=falltreatment, y=indweight, color = background, group = background)) + geom_point() + geom_line() +
#   facet_grid(pcode4~bdensity, scales = "free")

# phyto.dat <- comp.dat %>%
#   select(plot:bdensity, insitu_bdisturbed, phytometer, pcode4, insitu_pstems, p_seedfit, pfit_source) %>%
#   filter(!is.na(p_seedfit)) %>%
#   mutate(p_seedfit = ifelse(insitu_pstems == 0, 0, p_seedfit)) %>%
#   select(-phytometer, -bcode4) %>%
#   rename(seedsIn = insitu_pstems, seedsOut = p_seedfit) %>%
#   mutate(R = seedsOut/seedsIn,
#          R = ifelse(is.na(R), 0, R))


# isolate the phytometer data
phyto.dat0 <- comp.dat %>%
  select(plot:bdensity, insitu_bdisturbed, phytometer, pcode4, insitu_pstems, p_seedfit, p.ind.wgt.g, p_totwgt, pfit_source) %>%
  filter(!is.na(p_seedfit)) %>%
 # mutate(p_seedfit = ifelse(insitu_pstems == 0, 0, p_seedfit)) %>%
  select(-phytometer, -bcode4, -shelter) %>%
  rename(stemsIn = insitu_pstems,  disturbed = insitu_bdisturbed,
         block = shelterBlock) %>%
  mutate(seedsIn = ifelse(pcode4 == "AVFA" | pcode4 == "TRHI", 10, 12),
         seedsIn = ifelse(pcode4 == "LACA" | pcode4 == "ESCA", 15, seedsIn),
         seedsIn = ifelse(stemsIn > seedsIn, stemsIn, seedsIn)) %>%
  # code to subset seedfit based on individual pwgt or total pwgt
  subset(grepl("ind", pfit_source)) %>% # uncomment this line if want based on indidivual phyto wgt
  #subset(grepl("tot", pfit_source)) %>% # uncomment this line if want based on total phyto wgt
  select(-pfit_source) 

# and bring back allodat
allodat <- read.csv(paste0(datpath,"Competition_CleanedData/Competition_allometric_clean.csv"),
                    na.strings = na_vals, strip.white = T) %>%
  rename(bcode4 = species)

allodat2 <- allodat %>%
  rename(pcode4 = bcode4)

phyto.dat <- left_join(phyto.dat0, allodat2) %>%
  mutate(seedsOut = (intercept + p.ind.wgt.g*slope)*stemsIn) %>%
  select(-p_seedfit, -p.ind.wgt.g, -p_totwgt, -intercept:-slope_se) %>%
  mutate(R = seedsOut/seedsIn)

# ggplot(phyto.dat, aes(x= seedsOut, y = p_seedfit, color = pcode4)) + geom_point()+ facet_wrap(~pcode4)



# islate the background data
back.dat0 <- comp.dat %>%
  filter(!is.na(background)) %>%
  select(plot:bdensity,insitu_bdisturbed, b.ind.wgt.g, seedsAdded, insitu_plot_bdensity) %>%
  unique() 

back.dat <- left_join(back.dat0, allodat) %>%
  mutate(p_seedfit = (intercept + b.ind.wgt.g*slope)*insitu_plot_bdensity) %>%
  mutate(p_seedfit = ifelse(insitu_plot_bdensity == 0, 0, p_seedfit)) %>%
  select(plot:bdensity,insitu_bdisturbed,seedsAdded, p_seedfit, insitu_plot_bdensity) %>%
  rename(seedsIn = seedsAdded, seedsOut = p_seedfit, pcode4 = bcode4,
         stemsIn = insitu_plot_bdensity, disturbed = insitu_bdisturbed,
         block = shelterBlock) %>%
  filter(!is.na(seedsOut)) %>%
  select(plot, falltreatment, treatment, block, background, bdensity, disturbed, pcode4, stemsIn, seedsOut, seedsIn) %>%
  mutate(R = seedsOut/seedsIn)

phyto.datall <- rbind(phyto.dat, back.dat) %>%
  filter(pcode4 != "TRHI")

phytoagg_focal <- phyto.datall %>%
  group_by(falltreatment, pcode4, bdensity) %>%
  summarize(meanR = mean(R), seR = se(R)) %>%
  mutate(bdensity = ifelse(is.na(bdensity), "none", bdensity)) %>%
  mutate(bdensity = as.factor(bdensity))

ggplot(phytoagg_focal, aes(x=falltreatment, y=meanR, group = pcode4, color = pcode4)) + geom_point() + geom_line() +
  facet_wrap(~bdensity) + 
  geom_errorbar(aes(ymin = meanR - seR, ymax = meanR + seR), width = .2) + 
  scale_y_log10() + geom_hline(yintercept = 1, lty = "dashed")


phytoagg_background <- phyto.datall %>%
  group_by(falltreatment, background, bdensity) %>%
  summarize(meanR = mean(R), seR = se(R)) %>%
  tbl_df() %>%
  mutate(bdensity = ifelse(is.na(bdensity), "none", bdensity)) %>%
  mutate(bdensity = as.factor(bdensity))

ggplot(phytoagg_background, aes(x=falltreatment, y=meanR, group = background, color = background)) + geom_point() + geom_line() +
  facet_wrap(~bdensity) + 
  geom_errorbar(aes(ymin = meanR - seR, ymax = meanR + seR), width = .2) + 
  scale_y_log10() + geom_hline(yintercept = 1, lty = "dashed")


phytoagg_both <- phyto.datall %>%
  group_by(falltreatment, background, bdensity, pcode4) %>%
  summarize(meanR = mean(R), seR = se(R),
            meanseeds = mean(seedsOut)) %>%
  tbl_df() %>%
  mutate(bdensity = ifelse(is.na(bdensity), "none", bdensity)) %>%
  mutate(bdensity = as.factor(bdensity),
         meanR = ifelse(meanR==0, .0001, meanR)) %>%
  mutate(grouping = paste(background, pcode4))

ggplot(phytoagg_both, aes(x=falltreatment, y=meanR, group = grouping, color = grouping)) + geom_point() + geom_line() +
  facet_wrap(~bdensity) + 
 # geom_errorbar(aes(ymin = meanR - seR, ymax = meanR + seR), width = .2) + 
  scale_y_log10() + geom_hline(yintercept = 1, lty = "dashed")

ggplot(phytoagg_both, aes(x=falltreatment, y=meanR, group = grouping, color = pcode4)) + geom_point() + geom_line() +
  facet_wrap(~bdensity) + 
  # geom_errorbar(aes(ymin = meanR - seR, ymax = meanR + seR), width = .2) + 
  scale_y_log10() + geom_hline(yintercept = 1, lty = "dashed")

ggplot(phytoagg_both, aes(x=falltreatment, y = meanR, color = background, group = background)) + geom_point() + geom_line() + 
  facet_grid(pcode4~bdensity, scales = "free") 


ggplot(subset(phytoagg_both, pcode4%in%c("AVFA", "BRHO", "VUMY") & 
              background%in%c("Avena fatua", "Bromus hordeaceus", "Vulpia myuros")), 
       aes(x=falltreatment, y = meanR, color = background, group = background)) + geom_point() + geom_line() + 
  facet_grid(pcode4~bdensity, scales = "free") 



ggplot(subset(phytoagg_both, pcode4%in%c("AVFA", "BRHO", "VUMY") & 
                background%in%c("Avena fatua", "Bromus hordeaceus", "Vulpia myuros")), 
       aes(x=falltreatment, y = meanR, color = background, group = background)) + geom_point() + geom_line() + 
  facet_grid(pcode4~bdensity, scales = "free") 
