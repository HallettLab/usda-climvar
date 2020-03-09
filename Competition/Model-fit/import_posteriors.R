# run coexistence models
library(rstan)
library(bayesplot)
library(ggplot2)

load("./Competition/Model-fit/avfa_dry_posteriors.rdata")
avfa_dry <- extract(seeds_avfa_dry)
stan_dens(seeds_avfa_dry, pars = c("lambda", "alpha_avfa", "alpha_brho", "alpha_laca", "alpha_vumy"))
remove(seeds_avfa_dry)

load("./Competition/Model-fit/avfa_wet_posteriors.rdata")
avfa_wet <- extract(seeds_avfa_wet)
stan_dens(seeds_avfa_wet, pars = c("lambda", "alpha_avfa", "alpha_brho", "alpha_laca", "alpha_vumy"))
remove(seeds_avfa_wet)

load("./Competition/Model-fit/brho_wet_posteriors.rdata")
brho_wet <- extract(seeds_brho_wet)
stan_dens(seeds_brho_wet, pars = c("lambda", "alpha_avfa", "alpha_brho", "alpha_laca", "alpha_vumy"))
remove(seeds_brho_wet)

load("./Competition/Model-fit/brho_dry_posteriors.rdata")
brho_dry <- extract(seeds_brho_dry)
stan_dens(seeds_brho_dry, pars = c("lambda", "alpha_avfa", "alpha_brho", "alpha_laca", "alpha_vumy"))
remove(seeds_brho_dry)

load("./Competition/Model-fit/vumy_dry_posteriors.rdata")
vumy_dry <- extract(seeds_vumy_dry)
stan_dens(seeds_vumy_dry, pars = c("lambda", "alpha_avfa", "alpha_brho", "alpha_laca", "alpha_vumy"))
remove(seeds_vumy_dry)

load("./Competition/Model-fit/vumy_wet_posteriors.rdata")
vumy_wet <- extract(seeds_vumy_wet)
stan_dens(seeds_vumy_wet, pars = c("lambda", "alpha_avfa", "alpha_brho", "alpha_laca", "alpha_vumy"))
remove(seeds_vumy_wet)

vumy_dry$lambda <- exp(vumy_dry$lambda)
vumy_dry$alpha_avfa <- exp(vumy_dry$alpha_avfa)
vumy_dry$alpha_brho <- exp(vumy_dry$alpha_brho)
vumy_dry$alpha_laca <- exp(vumy_dry$alpha_laca)
vumy_dry$alpha_vumy <- exp(vumy_dry$alpha_vumy)

vumy_wet$lambda <- exp(vumy_wet$lambda)
vumy_wet$alpha_avfa <- exp(vumy_wet$alpha_avfa)
vumy_wet$alpha_brho <- exp(vumy_wet$alpha_brho)
vumy_wet$alpha_laca <- exp(vumy_wet$alpha_laca)
vumy_wet$alpha_vumy <- exp(vumy_wet$alpha_vumy)

