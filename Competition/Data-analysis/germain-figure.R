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


phyto.dat <- comp.dat %>%
  select(plot:bdensity, insitu_bdisturbed, phytometer, pcode4, insitu_pstems, p_totwgt_seedfit ) %>%
  filter(!is.na(p_totwgt_seedfit)) %>%
  mutate(p_totwgt_seedfit = ifelse(insitu_pstems == 0, 0, p_totwgt_seedfit)) %>%
  select(-phytometer, -bcode4) %>%
  rename(seedsIn = insitu_pstems, seedsOut = p_totwgt_seedfit) %>%
  mutate(R = seedsOut/seedsIn,
         R = ifelse(is.na(R), 0, R))



phytoagg_focal <- phyto.dat %>%
  group_by(falltreatment, pcode4, bdensity) %>%
  summarize(meanR = mean(R), seR = se(R)) %>%
  mutate(bdensity = ifelse(is.na(bdensity), "none", bdensity)) %>%
  mutate(bdensity = as.factor(bdensity))

ggplot(phytoagg_focal, aes(x=falltreatment, y=meanR, group = pcode4, color = pcode4)) + geom_point() + geom_line() +
  facet_wrap(~bdensity) + 
  geom_errorbar(aes(ymin = meanR - seR, ymax = meanR + seR), width = .2) + 
  scale_y_log10() + geom_hline(yintercept = 1, lty = "dashed")


phytoagg_background <- phyto.dat %>%
  group_by(falltreatment, background, bdensity) %>%
  summarize(meanR = mean(R), seR = se(R)) %>%
  tbl_df() %>%
  mutate(bdensity = ifelse(is.na(bdensity), "none", bdensity)) %>%
  mutate(bdensity = as.factor(bdensity))

ggplot(phytoagg_background, aes(x=falltreatment, y=meanR, group = background, color = background)) + geom_point() + geom_line() +
  facet_wrap(~bdensity) + 
  geom_errorbar(aes(ymin = meanR - seR, ymax = meanR + seR), width = .2) + 
  scale_y_log10() + geom_hline(yintercept = 1, lty = "dashed")


phytoagg_both <- phyto.dat %>%
  group_by(falltreatment, background, bdensity, pcode4) %>%
  summarize(meanR = mean(R), seR = se(R)) %>%
  tbl_df() %>%
  mutate(bdensity = ifelse(is.na(bdensity), "none", bdensity)) %>%
  mutate(bdensity = as.factor(bdensity),
         meanR = ifelse(meanR==0, .0001, meanR)) %>%
  mutate(grouping = paste(background, pcode4))

ggplot(phytoagg_both, aes(x=falltreatment, y=meanR, group = grouping, color = grouping)) + geom_point() + geom_line() +
  facet_wrap(~bdensity) + 
 # geom_errorbar(aes(ymin = meanR - seR, ymax = meanR + seR), width = .2) + 
  scale_y_log10() + geom_hline(yintercept = 1, lty = "dashed")


