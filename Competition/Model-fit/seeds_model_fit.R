# Run experimental data with fecundity model

# source data
source("./Competition/Data-analysis/coexistence-model_formatting.R")

library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

library(here)

setwd(here("Competition", "Model-fit"))
disturbed <- which(seedsin.seedsout$disturbed == 1)
data <- seedsin.seedsout[-disturbed,]

initials <- list(lambda=exp(10), alpha_avfa=exp(0.03), alpha_brho=exp(0.03), alpha_laca=exp(0.03), 
                 alpha_vumy=exp(0.03))
initials1<- list(initials, initials, initials)
# ---------------------------------------------------------------------------------------------
# avfa model fitting for wet

dat <- subset(data, species == "AVFA")
dat <- subset(dat, falltreatment == "wet")

Fecundity <- as.integer(round(dat$seedsOut))

avfa <- as.integer(dat$AVFA_seedsIn)
brho <- as.integer(dat$BRHO_seedsIn)
laca <- as.integer(dat$LACA_seedsIn)
vumy <- as.integer(dat$VUMY_seedsIn)

N <- as.integer(length(Fecundity))

intra <- avfa

no_dist_seeds_avfa_wet <- stan(file = "Generic_four_species_BH_model.stan", data = c("N", "Fecundity", "intra", "avfa", "brho",
                                                                             "laca", "vumy"),
                       iter = 9000, chains = 3, thin = 3, control = list(adapt_delta = 0.85, max_treedepth = 10),
                       init = initials1)

save(no_dist_seeds_avfa_wet, file = "no_dist_avfa_wet_posteriors.rdata")
# ---------------------------------------------------------------------------------------------
# avfa model fitting for dry

dat <- subset(data, species == "AVFA")
dat <- subset(dat, falltreatment == "dry")

Fecundity <- as.integer(round(dat$seedsOut))

avfa <- as.integer(dat$AVFA_seedsIn)
brho <- as.integer(dat$BRHO_seedsIn)
laca <- as.integer(dat$LACA_seedsIn)
vumy <- as.integer(dat$VUMY_seedsIn)

N <- as.integer(length(Fecundity))

intra <- avfa

no_dist_seeds_avfa_dry <- stan(file = "Generic_four_species_BH_model.stan", data = c("N", "Fecundity", "intra", "avfa", "brho",
                                                                             "laca", "vumy"),
                       iter = 9000, chains = 3, thin = 3, control = list(adapt_delta = 0.85, max_treedepth = 20),
                       init = initials1)

save(no_dist_seeds_avfa_dry, file = "no_dist_avfa_dry_posteriors.rdata")
# ---------------------------------------------------------------------------------------------
# brho model fitting for wet

dat <- subset(data, species == "BRHO")
dat <- subset(dat, falltreatment == "wet")

Fecundity <- as.integer(round(dat$seedsOut))

avfa <- as.integer(dat$AVFA_seedsIn)
brho <- as.integer(dat$BRHO_seedsIn)
laca <- as.integer(dat$LACA_seedsIn)
vumy <- as.integer(dat$VUMY_seedsIn)

N <- as.integer(length(Fecundity))

intra <- brho

no_dist_seeds_brho_wet <- stan(file = "Generic_four_species_BH_model.stan", data = c("N", "Fecundity", "intra", "avfa", "brho",
                                                                             "laca", "vumy"),
                       iter = 9000, chains = 3, thin = 3, control = list(adapt_delta = 0.85, max_treedepth = 10),
                       init = initials1)

save(no_dist_seeds_brho_wet, file = "no_dist_brho_wet_posteriors.rdata")
# ---------------------------------------------------------------------------------------------
# brho model fitting for dry

dat <- subset(data, species == "BRHO")
dat <- subset(dat, falltreatment == "dry")

Fecundity <- as.integer(round(dat$seedsOut))

avfa <- as.integer(dat$AVFA_seedsIn)
brho <- as.integer(dat$BRHO_seedsIn)
laca <- as.integer(dat$LACA_seedsIn)
vumy <- as.integer(dat$VUMY_seedsIn)

N <- as.integer(length(Fecundity))

intra <- brho

no_dist_seeds_brho_dry <- stan(file = "Generic_four_species_BH_model.stan", data = c("N", "Fecundity", "intra", "avfa", "brho",
                                                                             "laca", "vumy"),
                       iter = 9000, chains = 3, thin = 3, control = list(adapt_delta = 0.85, max_treedepth = 20),
                       init = initials1)

save(no_dist_seeds_brho_dry, file = "no_dist_brho_dry_posteriors.rdata")
# ---------------------------------------------------------------------------------------------
# laca model fitting for dry

dat <- subset(data, species == "LACA")
dat <- subset(dat, falltreatment == "dry")

Fecundity <- as.integer(round(dat$seedsOut))

avfa <- as.integer(dat$AVFA_seedsIn)
brho <- as.integer(dat$BRHO_seedsIn)
laca <- as.integer(dat$LACA_seedsIn)
vumy <- as.integer(dat$VUMY_seedsIn)

N <- as.integer(length(Fecundity))

intra <- laca

seeds_laca_dry <- stan(file = "Generic_four_species_BH_model.stan", data = c("N", "Fecundity", "intra", "avfa", "brho",
                                                                       "laca", "vumy"),
                       iter = 6000, chains = 3, thin = 2, control = list(adapt_delta = 0.95, max_treedepth = 20), 
                       init = initials1)

save(seeds_laca_dry, file = "laca_dry_posteriors.rdata")
# ---------------------------------------------------------------------------------------------
# laca model fitting for wet

dat <- subset(data, species == "LACA")
dat <- subset(dat, falltreatment == "wet")

Fecundity <- as.integer(round(dat$seedsOut))

avfa <- as.integer(dat$AVFA_seedsIn)
brho <- as.integer(dat$BRHO_seedsIn)
laca <- as.integer(dat$LACA_seedsIn)
vumy <- as.integer(dat$VUMY_seedsIn)

N <- as.integer(length(Fecundity))

intra <- laca

seeds_laca_wet <- stan(file = "Generic_four_species_BH_model.stan", data = c("N", "Fecundity", "intra", "avfa", "brho",
                                                                       "laca", "vumy"),
                       iter = 6000, chains = 3, thin = 2, control = list(adapt_delta = 0.95, max_treedepth = 20), 
                       init = initials1)

save(seeds_laca_wet, file = "laca_wet_posteriors.rdata")
# ---------------------------------------------------------------------------------------------
# vumy model fitting for wet

dat <- subset(data, species == "VUMY")
dat <- subset(dat, falltreatment == "wet")

Fecundity <- as.integer(round(dat$seedsOut))

avfa <- as.integer(dat$AVFA_seedsIn)
brho <- as.integer(dat$BRHO_seedsIn)
laca <- as.integer(dat$LACA_seedsIn)
vumy <- as.integer(dat$VUMY_seedsIn)

N <- as.integer(length(Fecundity))

intra <- vumy

no_dist_seeds_vumy_wet <- stan(file = "Generic_four_species_BH_model.stan", data = c("N", "Fecundity", "intra", "avfa", "brho",
                                                                             "laca", "vumy"),
                       iter = 9000, chains = 3, thin = 3, control = list(adapt_delta = 0.85, max_treedepth = 20),
                       init = initials1)

save(no_dist_seeds_vumy_wet, file = "no_dist_vumy_wet_posteriors.rdata")
# ---------------------------------------------------------------------------------------------
# vumy model fitting for dry

dat <- subset(data, species == "VUMY")
dat <- subset(dat, falltreatment == "dry")

Fecundity <- as.integer(round(dat$seedsOut))

avfa <- as.integer(dat$AVFA_seedsIn)
brho <- as.integer(dat$BRHO_seedsIn)
laca <- as.integer(dat$LACA_seedsIn)
vumy <- as.integer(dat$VUMY_seedsIn)

N <- as.integer(length(Fecundity))

intra <- vumy

no_dist_seeds_vumy_dry <- stan(file = "Generic_four_species_BH_model.stan", data = c("N", "Fecundity", "intra", "avfa", "brho",
                                                                             "laca", "vumy"),
                       iter = 9000, chains = 3, thin = 3, control = list(adapt_delta = 0.85, max_treedepth = 20),
                       init = initials1)

save(no_dist_seeds_vumy_dry, file = "no_dist_vumy_dry_posteriors.rdata")
