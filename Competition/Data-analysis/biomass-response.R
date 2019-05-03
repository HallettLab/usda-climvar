# preliminary analysis of phytometer growth and fecundity under different competition scenarios
# authors: LMH, CTW
# initiated: Oct 2018 (modified over time)

# script purpose:
# ...

# notes:
# ...



# -- SETUP -----
rm(list = ls()) # clean environment
library(tidyverse)
options(stringsAsFactors = F)
theme_set(theme_bw())
na_vals = c(" ","", NA, "NA")

# set pathway to climvar dropbox competition data folder
datpath <- "~/Dropbox/ClimVar/Competition/Data/"

# read in data
# cleaned, combined competition dataset (cleaned background, phytometers, predicted vals, plot treatments)
comp.dat <- read.csv(paste0(datpath, "Competition_CleanedData/Competition_combined_clean.csv"),
                     na.strings = na_vals, strip.white = T)
# spp list
spplist <- read.csv(paste0(datpath, "Competition_SpeciesKey.csv"),
                    na.strings = na_vals, strip.white = T)


# -- VISUALIZE TRENDS -----
# take the average across blocks
comp.dat2 <- comp.dat %>%
  group_by(background, bcode4, bdensity, phytometer, pcode4, falltreatment) %>%
  # filter(!is.na(disturbance)) %>%
  summarize(phyto_mean_weight = mean(p.ind.wgt.g, na.rm = T)) %>%
  ungroup() %>%
  mutate(background = ifelse(is.na(background), "Control", background),
         bdensity = ifelse(is.na(bdensity), "none", bdensity))


#### phyto mean weight by fall drought treatment, background competitors, and background density ####
ggplot(comp.dat2, aes(x=falltreatment, y=phyto_mean_weight)) + 
  geom_boxplot() +  
  ggtitle("Mean phytometer weight by background competitor and background density treatment") +
  facet_grid(bdensity~background)

ggplot(comp.dat2, aes(x=bdensity, y = phyto_mean_weight, color = falltreatment, group = falltreatment)) + geom_point() +
  geom_line() + 
  ggtitle("Mean phytometer weight by fall drought treatment, background competitor (x strip), and phytometer (y strip)") +
  facet_grid(pcode4~bcode4, scales = "free")

ggplot(comp.dat2, aes(x=falltreatment, y = phyto_mean_weight, color = background, group = background)) + 
  geom_point() +
  geom_line() + 
  labs(title = "Mean phytometer weight by fall drought treatment, background density and competitor",
       subtitle = "y strip = phytometers, colored in by background competitor plot") +
  facet_grid(phytometer~bdensity, scales = "free_y", labeller = label_wrap_gen(width = 10))
ggsave("~/Dropbox/ClimVar/Competition/Figures/Exploratory/all-competition_May2019.pdf", width = 8, height = 10)

comp.dat2 %>%
  filter(bcode4 %in% c(spplist$code4[spplist$fxnl_grp =="Grass"],"Control"),
         pcode4 %in% spplist$code4[spplist$fxnl_grp == "Grass"]) %>%
  ggplot(aes(x=falltreatment, y = phyto_mean_weight, color = background, group = background)) + geom_point() +
  geom_line() + 
  ggtitle("Mean phytometer weight by fall treatment and background density, grasses only") +
  facet_grid(phytometer~bdensity, scales = "free")
ggsave("~/Dropbox/ClimVar/Competition/Figures/Exploratory/grass-competition_May2019.pdf", width = 8, height = 10)


#### fecundity trends #####
# take the average across blocks
comp.dat3 <- ungroup(comp.dat) %>%
  group_by(background, bdensity, phytometer, falltreatment, seedsAdded) %>%
  # filter(!is.na(disturbance)) %>%
  summarize(b_mean_density = mean(insitu_plot_bdensity, na.rm = T),
            b_mean_weight = mean(b.ind.wgt.g, na.rm = T),
            phyto_mean_weight = mean(p.ind.wgt.g, na.rm = T),
            phyto_mean_seed = mean(p_totwgt_seedfit, na.rm = T)) %>%
  ungroup() %>%
  mutate(background = ifelse(is.na(background), "Control", background),
         bdensity = ifelse(is.na(bdensity), "none", bdensity)) %>%
  filter(!is.nan(phyto_mean_seed))


# some graphs of density and biomass
comp.dat3[c("falltreatment", "b_mean_density", "bdensity", "background")] %>%
  filter(background!= "Control") %>%
  distinct() %>%
  ggplot(aes(x=falltreatment, y=b_mean_density)) + 
  ggtitle("Mean background density by competitor and density treatment") +
  geom_point() + facet_grid(background~bdensity, scales = "free_y")

ggplot(subset(comp.dat3, bdensity != "none"), 
       aes(x=bdensity, y = phyto_mean_weight, color = background, group = background)) + geom_point() +
  geom_line() + 
  facet_grid(phytometer~falltreatment, scales = "free")

ggplot(comp.dat3, aes(x=bdensity, y = phyto_mean_weight, color = falltreatment, group = falltreatment)) + geom_point() +
  geom_line() + 
  facet_grid(phytometer~background, scales = "free")


# some graphs of seeds
ggplot(comp.dat3, aes(x=bdensity, y = phyto_mean_seed, color = background, group = background)) + geom_point() +
  geom_line() + 
  ggtitle("Mean fecundity, by fall treatment, background density and competitor") +
  facet_grid(phytometer~falltreatment, scales = "free")

ggplot(comp.dat3, aes(x=bdensity, y = phyto_mean_seed, color = falltreatment, group = falltreatment)) + geom_point() +
  geom_line() +
  ggtitle("Phytomer mean fecundity") +
  facet_grid(phytometer~background, scales = "free")

ggplot(subset(comp.dat, !is.na(insitu_plot_bdensity) & !pcode4 %in% c("ESCA", "TRHI")), aes(insitu_plot_bdensity, p_totwgt_seedfit)) +
  geom_point(aes(col = background)) +
  labs(x = "Competitor plot density (# plants)",
       y = "Phytometer fecundity (# seeds)",
       title = "Projected phytometer fecundity and competitor density, by density treatment and phytometer") +
  scale_color_manual(values = c("Avena fatua" = "springgreen4", "Bromus hordeaceus" = "springgreen1", "Vulpia myuros" = "green3",
                                "Eschscholzia californica" = "chocolate2", "Lasthenia californica" = "darkgoldenrod1",  "Trifolium hirtum" = "orchid")) +
  facet_grid(phytometer~bdensity, scales = "free")
ggsave("~/Dropbox/ClimVar/Competition/Figures/Exploratory/all-fecundity.pdf", width = 8, height = 10)


# seems like:
# Avena is always the superior competitor
# Brome vs Vulpia is weather dependent: Brome can increase when rare under wet, Vuplia under dry
# Lasthenia is always the weakest competitor