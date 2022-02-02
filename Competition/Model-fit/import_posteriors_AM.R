# extract posteriors
library(rstan)
library(bayesplot)
library(ggplot2)
library(tidyverse)

## AVFA
load("./Competition/Model-fit/posteriors/avfa_dry_posteriors_updated_germ.rdata")
avfa_dry <- rstan::extract(no_dist_seeds_avfa_dry)
#stan_dens(no_dist_seeds_avfa_dry, pars = c("lambda", "alpha_avfa", "alpha_brho", "alpha_esca", "alpha_laca", "alpha_vumy", "alpha_trhi"))
remove(no_dist_seeds_avfa_dry)

load("./Competition/Model-fit/posteriors/avfa_wet_posteriors_updated_germ.rdata")
avfa_wet <- rstan::extract(no_dist_seeds_avfa_wet)
#stan_dens(no_dist_seeds_avfa_wet, pars = c("lambda", "alpha_avfa", "alpha_brho", "alpha_esca", "alpha_laca", "alpha_vumy", "alpha_trhi"))
remove(no_dist_seeds_avfa_wet)

# Convert to data frame
params_avfa <- rbind(as.data.frame(lapply(avfa_dry, mean))[-8],
                     as.data.frame(lapply(avfa_wet, mean))[-8])
params_avfa$species <- "avfa"
params_avfa$treatment <- c("fallDry","fallWet")

## BRHO
load("./Competition/Model-fit/posteriors/brho_wet_posteriors_updated_germ.rdata")
brho_wet <- rstan::extract(no_dist_seeds_brho_wet)
#stan_dens(no_dist_seeds_brho_wet, pars = c("lambda", "alpha_avfa", "alpha_brho", "alpha_esca", "alpha_laca", "alpha_vumy", "alpha_trhi"))
remove(no_dist_seeds_brho_wet)

load("./Competition/Model-fit/posteriors/brho_dry_posteriors_updated_germ.rdata")
brho_dry <- rstan::extract(no_dist_seeds_brho_dry)
#stan_dens(no_dist_seeds_brho_dry, pars = c("lambda", "alpha_avfa", "alpha_brho", "alpha_esca", "alpha_laca", "alpha_vumy", "alpha_trhi"))
remove(no_dist_seeds_brho_dry)

# Convert to data frame
params_brho <- rbind(as.data.frame(lapply(brho_dry, mean))[-8],
                     as.data.frame(lapply(brho_wet, mean))[-8])
params_brho$species <- "brho"
params_brho$treatment <- c("fallDry","fallWet")

## VUMY
load("./Competition/Model-fit/posteriors/vumy_dry_posteriors_updated_germ.rdata")
vumy_dry <- rstan::extract(no_dist_seeds_vumy_dry)
#stan_dens(no_dist_seeds_vumy_dry, pars = c("lambda", "alpha_avfa", "alpha_brho", "alpha_esca", "alpha_laca", "alpha_vumy", "alpha_trhi"))
remove(no_dist_seeds_vumy_dry)

load("./Competition/Model-fit/posteriors/vumy_wet_posteriors_updated_germ.rdata")
vumy_wet <- rstan::extract(no_dist_seeds_vumy_wet)
#stan_dens(no_dist_seeds_vumy_wet, pars = c("lambda", "alpha_avfa", "alpha_brho", "alpha_esca", "alpha_laca", "alpha_vumy", "alpha_trhi"))
remove(no_dist_seeds_vumy_wet)

# Convert to data frame
params_vumy <- rbind(as.data.frame(lapply(vumy_dry, mean))[-8],
                     as.data.frame(lapply(vumy_wet, mean))[-8])
params_vumy$species <- "vumy"
params_vumy$treatment <- c("fallDry","fallWet")

## LACA
load("./Competition/Model-fit/posteriors/laca_dry_posteriors_updated_germ.rdata")
laca_dry <- rstan::extract(no_dist_seeds_laca_dry)
#stan_dens(no_dist_seeds_laca_dry, pars = c("lambda", "alpha_avfa", "alpha_brho", "alpha_esca", "alpha_laca", "alpha_vumy", "alpha_trhi"))
remove(no_dist_seeds_laca_dry)

load("./Competition/Model-fit/posteriors/laca_wet_posteriors_updated_germ.rdata")
laca_wet <- rstan::extract(no_dist_seeds_laca_wet)
#stan_dens(no_dist_seeds_laca_wet, pars = c("lambda", "alpha_avfa", "alpha_brho", "alpha_esca", "alpha_laca", "alpha_vumy", "alpha_trhi"))
remove(no_dist_seeds_laca_wet)

# Convert to data frame
params_laca <- rbind(as.data.frame(lapply(laca_dry, mean))[-8],
                     as.data.frame(lapply(laca_wet, mean))[-8])
params_laca$species <- "laca"
params_laca$treatment <- c("fallDry","fallWet")

## ESCA
load("./Competition/Model-fit/posteriors/esca_dry_posteriors_updated_germ.rdata")
esca_dry <- rstan::extract(no_dist_seeds_esca_dry)
#stan_dens(no_dist_seeds_esca_dry, pars = c("lambda", "alpha_avfa", "alpha_brho", "alpha_esca", "alpha_laca", "alpha_vumy", "alpha_trhi"))
remove(no_dist_seeds_esca_dry)

load("./Competition/Model-fit/posteriors/esca_wet_posteriors_updated_germ.rdata")
esca_wet <- rstan::extract(no_dist_seeds_esca_wet)
#stan_dens(no_dist_seeds_esca_wet, pars = c("lambda", "alpha_avfa", "alpha_brho", "alpha_esca", "alpha_laca", "alpha_vumy", "alpha_trhi"))
remove(no_dist_seeds_esca_wet)

# Convert to data frame
params_esca <- rbind(as.data.frame(lapply(esca_dry, mean))[-8],
                     as.data.frame(lapply(esca_wet, mean))[-8])
params_esca$species <- "esca"
params_esca$treatment <- c("fallDry","fallWet")

## TRHI
load("./Competition/Model-fit/posteriors/trhi_dry_posteriors_updated_germ.rdata")
trhi_dry <- rstan::extract(no_dist_seeds_trhi_dry)
stan_dens(no_dist_seeds_trhi_dry, pars = c("lambda", "alpha_avfa", "alpha_brho", "alpha_esca", "alpha_laca", "alpha_vumy", "alpha_trhi"))
remove(no_dist_seeds_trhi_dry)

load("./Competition/Model-fit/posteriors/trhi_wet_posteriors_updated_germ.rdata")
trhi_wet <- rstan::extract(no_dist_seeds_trhi_wet)
stan_dens(no_dist_seeds_trhi_wet, pars = c("lambda", "alpha_avfa", "alpha_brho", "alpha_esca", "alpha_laca", "alpha_vumy", "alpha_trhi"))
remove(no_dist_seeds_trhi_wet)

# Convert to data frame
params_trhi <- rbind(as.data.frame(lapply(trhi_dry, mean))[-8],
                     as.data.frame(lapply(trhi_wet, mean))[-8])
params_trhi$species <- "trhi"
params_trhi$treatment <- c("fallDry","fallWet")

# Create total dataframe
params <- rbind(params_avfa, params_vumy, params_brho, params_esca, params_laca, params_trhi)
rm(list = ls()[grep("params_",ls())])

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

