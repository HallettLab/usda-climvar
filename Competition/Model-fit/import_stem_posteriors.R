# extract posteriors
library(rstan)
library(bayesplot)
library(ggplot2)
library(tidyverse)

# set working directory to project directory, if not already there!

load("./Competition/Model-fit/stems_avfa_dry_posteriors.rdata")
avfa_dry <- rstan::extract(stems_avfa_dry)
remove(stems_avfa_dry)

load("./Competition/Model-fit/stems_avfa_wet_posteriors.rdata")
avfa_wet <- rstan::extract(stems_avfa_wet)
remove(stems_avfa_wet)

load("./Competition/Model-fit/stems_brho_wet_posteriors.rdata")
brho_wet <- rstan::extract(stems_brho_wet)
remove(stems_brho_wet)

load("./Competition/Model-fit/stems_brho_dry_posteriors.rdata")
brho_dry <- rstan::extract(stems_brho_dry)
remove(stems_brho_dry)

load("./Competition/Model-fit/stems_vumy_dry_posteriors.rdata")
vumy_dry <- rstan::extract(stems_vumy_dry)
remove(stems_vumy_dry)

load("./Competition/Model-fit/stems_vumy_wet_posteriors.rdata")
vumy_wet <- rstan::extract(stems_vumy_wet)
remove(stems_vumy_wet)

#vumy_dry$lambda <- exp(vumy_dry$lambda)
#vumy_dry$alpha_avfa <- exp(vumy_dry$alpha_avfa)
#vumy_dry$alpha_brho <- exp(vumy_dry$alpha_brho)
#vumy_dry$alpha_laca <- exp(vumy_dry$alpha_laca)
#vumy_dry$alpha_vumy <- exp(vumy_dry$alpha_vumy)

#vumy_wet$lambda <- exp(vumy_wet$lambda)
#vumy_wet$alpha_avfa <- exp(vumy_wet$alpha_avfa)
#vumy_wet$alpha_brho <- exp(vumy_wet$alpha_brho)
#vumy_wet$alpha_laca <- exp(vumy_wet$alpha_laca)
#vumy_wet$alpha_vumy <- exp(vumy_wet$alpha_vumy)

# unlist each of the species/rainfall combinations
vw <- as_tibble(do.call("cbind", vumy_wet)) %>%
  mutate(species = "vumy", treatment = "wet") %>%
  select(-lp__)

vd <- as_tibble(do.call("cbind", vumy_dry)) %>%
  mutate(species = "vumy", treatment = "dry") %>%
  select(-lp__)

aw <- as_tibble(do.call("cbind", avfa_wet)) %>%
  mutate(species = "avfa", treatment = "wet") %>%
  select(-lp__)

ad <- as_tibble(do.call("cbind", avfa_dry)) %>%
  mutate(species = "avfa", treatment = "dry") %>%
  select(-lp__)


bw <- as_tibble(do.call("cbind", brho_wet)) %>%
  mutate(species = "brho", treatment = "wet") %>%
  select(-lp__)

bd <- as_tibble(do.call("cbind", brho_dry)) %>%
  mutate(species = "brho", treatment = "dry") %>%
  select(-lp__)

# put them together
posteriors <- rbind(vw, vd, aw, ad, bw, bd) %>%
  mutate(speciestreatment = paste(species, treatment, sep = " "))

#rm(list=setdiff(ls(), "posteriors")) 
rm(vw, vd, aw, ad, bw, bd)

ggplot(posteriors, aes(x=lambda, fill = treatment, line = treatment)) + geom_density() + 
  facet_wrap(~species, ncol = 1, scales = "free_y")


ggplot(posteriors, aes(x=alpha_avfa, fill = treatment, line = treatment)) + geom_density() + 
  facet_wrap(~species, ncol = 1, scales = "free_y")


ggplot(posteriors, aes(x=alpha_brho, fill = treatment, line = treatment)) + geom_density() + 
  facet_wrap(~species, ncol = 1, scales = "free_y")


ggplot(posteriors, aes(x=alpha_vumy, fill = treatment, line = treatment)) + geom_density() + 
  facet_wrap(~species, ncol = 1, scales = "free_y")

