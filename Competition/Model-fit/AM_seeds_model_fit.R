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

initials <- list(lambda=1000, alpha_avfa=1, alpha_brho=1, alpha_esca=1,
                 alpha_laca=1, alpha_vumy=1)
initials1<- list(initials, initials, initials, initials)
# ---------------------------------------------------------------------------------------------
# avfa model fitting for wet

dat <- subset(data, species == "AVFA")
dat <- subset(dat, falltreatment == "wet")

Fecundity <- as.integer(round(dat$seedsOut))

avfa <- as.integer(dat$AVFA_seedsIn)
brho <- as.integer(dat$BRHO_seedsIn)
laca <- as.integer(dat$LACA_seedsIn)
vumy <- as.integer(dat$VUMY_seedsIn)
esca <- as.integer(dat$ESCA_seedsIn)

N <- as.integer(length(Fecundity))

intra <- avfa

no_dist_seeds_avfa_wet <- stan(file = "Five_species_BH_model.stan", data = c("N", "Fecundity", "intra", "avfa", "brho",
                                                                             "esca", "laca", "vumy"),
                       iter = 40000, chains = 4, thin = 3, control = list(adapt_delta = 0.99, max_treedepth = 20),
                       init = initials1)

save(no_dist_seeds_avfa_wet, file = "./posteriors/AM_avfa_wet_posteriors_tot.rdata")
# ---------------------------------------------------------------------------------------------
# avfa model fitting for dry

dat <- subset(data, species == "AVFA")
dat <- subset(dat, falltreatment == "dry")

Fecundity <- as.integer(round(dat$seedsOut))

avfa <- as.integer(dat$AVFA_seedsIn)
brho <- as.integer(dat$BRHO_seedsIn)
laca <- as.integer(dat$LACA_seedsIn)
vumy <- as.integer(dat$VUMY_seedsIn)
esca <- as.integer(dat$ESCA_seedsIn)

N <- as.integer(length(Fecundity))

intra <- avfa

no_dist_seeds_avfa_dry <- stan(file = "Five_species_BH_model.stan", data = c("N", "Fecundity", "intra", "avfa", "brho", "esca",
                                                                             "laca", "vumy"),
                       iter = 50000, chains = 4, thin = 3, control = list(adapt_delta = 0.99, max_treedepth = 20),
                       init = initials1)

save(no_dist_seeds_avfa_dry, file = "./posteriors/AM_avfa_dry_posteriors_tot.rdata")
# ---------------------------------------------------------------------------------------------
# brho model fitting for wet

dat <- subset(data, species == "BRHO")
dat <- subset(dat, falltreatment == "wet")

Fecundity <- as.integer(round(dat$seedsOut))

avfa <- as.integer(dat$AVFA_seedsIn)
brho <- as.integer(dat$BRHO_seedsIn)
laca <- as.integer(dat$LACA_seedsIn)
vumy <- as.integer(dat$VUMY_seedsIn)
esca <- as.integer(dat$ESCA_seedsIn)

N <- as.integer(length(Fecundity))

intra <- brho

no_dist_seeds_brho_wet <- stan(file = "Five_species_BH_model.stan", data = c("N", "Fecundity", "intra", "avfa", "brho", "esca",
                                                                                     "laca", "vumy"),
                               iter = 40000, chains = 4, thin = 3, control = list(adapt_delta = 0.99, max_treedepth = 20),
                               init = initials1)

save(no_dist_seeds_brho_wet, file = "./posteriors/AM_brho_wet_posteriors_tot.rdata")
# ---------------------------------------------------------------------------------------------
# brho model fitting for dry

dat <- subset(data, species == "BRHO")
dat <- subset(dat, falltreatment == "dry")

Fecundity <- as.integer(round(dat$seedsOut))

avfa <- as.integer(dat$AVFA_seedsIn)
brho <- as.integer(dat$BRHO_seedsIn)
laca <- as.integer(dat$LACA_seedsIn)
vumy <- as.integer(dat$VUMY_seedsIn)
esca <- as.integer(dat$ESCA_seedsIn)

N <- as.integer(length(Fecundity))

intra <- brho

no_dist_seeds_brho_dry <- stan(file = "Five_species_BH_model.stan", data = c("N", "Fecundity", "intra", "avfa", "brho", "esca",
                                                                             "laca", "vumy"),
                       iter = 40000, chains = 4, thin = 3, control = list(adapt_delta = 0.99, max_treedepth = 20),
                       init = initials1)

save(no_dist_seeds_brho_dry, file = "./posteriors/AM_brho_dry_posteriors_tot.rdata")
# ---------------------------------------------------------------------------------------------
# laca model fitting for dry

dat <- subset(data, species == "LACA")
dat <- subset(dat, falltreatment == "dry")

Fecundity <- as.integer(round(dat$seedsOut))

avfa <- as.integer(dat$AVFA_seedsIn)
brho <- as.integer(dat$BRHO_seedsIn)
laca <- as.integer(dat$LACA_seedsIn)
vumy <- as.integer(dat$VUMY_seedsIn)
esca <- as.integer(dat$ESCA_seedsIn)

N <- as.integer(length(Fecundity))

intra <- laca

no_dist_seeds_laca_dry <- stan(file = "Five_species_BH_model.stan", data = c("N", "Fecundity", "intra", "avfa", "brho","esca",
                                                                       "laca", "vumy"),
                       iter = 40000, chains = 4, thin = 3, control = list(adapt_delta = 0.99, max_treedepth = 20), 
                       init = initials1)

save(no_dist_seeds_laca_dry, file = "./posteriors/AM_laca_dry_posteriors_tot.rdata")
# ---------------------------------------------------------------------------------------------
# laca model fitting for wet

dat <- subset(data, species == "LACA")
dat <- subset(dat, falltreatment == "wet")

Fecundity <- as.integer(round(dat$seedsOut))

avfa <- as.integer(dat$AVFA_seedsIn)
brho <- as.integer(dat$BRHO_seedsIn)
laca <- as.integer(dat$LACA_seedsIn)
vumy <- as.integer(dat$VUMY_seedsIn)
esca <- as.integer(dat$ESCA_seedsIn)

N <- as.integer(length(Fecundity))

intra <- laca

no_dist_seeds_laca_wet <- stan(file = "Five_species_BH_model.stan", data = c("N", "Fecundity", "intra", "avfa", "brho","esca",
                                                                       "laca", "vumy"),
                       iter = 40000, chains = 4, thin = 3, control = list(adapt_delta = 0.99, max_treedepth = 20), 
                       init = initials1)

save(no_dist_seeds_laca_wet, file = "./posteriors/AM_laca_wet_posteriors_tot.rdata")
# ---------------------------------------------------------------------------------------------
# vumy model fitting for wet

dat <- subset(data, species == "VUMY")
dat <- subset(dat, falltreatment == "wet")

Fecundity <- as.integer(round(dat$seedsOut))

avfa <- as.integer(dat$AVFA_seedsIn)
brho <- as.integer(dat$BRHO_seedsIn)
laca <- as.integer(dat$LACA_seedsIn)
vumy <- as.integer(dat$VUMY_seedsIn)
esca <- as.integer(dat$ESCA_seedsIn)

N <- as.integer(length(Fecundity))

intra <- vumy

no_dist_seeds_vumy_wet <- stan(file = "Five_species_BH_model.stan", data = c("N", "Fecundity", "intra", "avfa", "brho","esca",
                                                                             "laca", "vumy"),
                       iter = 40000, chains = 4, thin = 3, control = list(adapt_delta = 0.99, max_treedepth = 20),
                       init = initials1)

save(no_dist_seeds_vumy_wet, file = "./posteriors/AM_vumy_wet_posteriors_tot.rdata")
# ---------------------------------------------------------------------------------------------
# vumy model fitting for dry

dat <- subset(data, species == "VUMY")
dat <- subset(dat, falltreatment == "dry")

Fecundity <- as.integer(round(dat$seedsOut))

avfa <- as.integer(dat$AVFA_seedsIn)
brho <- as.integer(dat$BRHO_seedsIn)
laca <- as.integer(dat$LACA_seedsIn)
vumy <- as.integer(dat$VUMY_seedsIn)
esca <- as.integer(dat$ESCA_seedsIn)

N <- as.integer(length(Fecundity))

intra <- vumy

no_dist_seeds_vumy_dry <- stan(file = "Five_species_BH_model.stan", data = c("N", "Fecundity", "intra", "avfa", "brho","esca",
                                                                             "laca", "vumy"),
                       iter = 40000, chains = 4, thin = 3, control = list(adapt_delta = 0.99, max_treedepth = 20),
                       init = initials1)

save(no_dist_seeds_vumy_dry, file = "./posteriors/AM_vumy_dry_posteriors_tot.rdata")

# ---------------------------------------------------------------------------------------------
# esca model fitting for wet

dat <- subset(data, species == "ESCA")
dat <- subset(dat, falltreatment == "wet")

Fecundity <- as.integer(round(dat$seedsOut))

avfa <- as.integer(dat$AVFA_seedsIn)
brho <- as.integer(dat$BRHO_seedsIn)
laca <- as.integer(dat$LACA_seedsIn)
vumy <- as.integer(dat$VUMY_seedsIn)
esca <- as.integer(dat$ESCA_seedsIn)

N <- as.integer(length(Fecundity))

intra <- esca

no_dist_seeds_esca_wet <- stan(file = "Five_species_BH_model.stan", data = c("N", "Fecundity", "intra", "avfa", "brho","esca",
                                                                             "laca", "vumy"),
                               iter = 40000, chains = 4, thin = 3, control = list(adapt_delta = 0.99, max_treedepth = 20),
                               init = initials1)

save(no_dist_seeds_esca_wet, file = "./posteriors/AM_esca_wet_posteriors_tot.rdata")

# ---------------------------------------------------------------------------------------------
# esca model fitting for dry

dat <- subset(data, species == "ESCA")
dat <- subset(dat, falltreatment == "dry")

Fecundity <- as.integer(round(dat$seedsOut))

avfa <- as.integer(dat$AVFA_seedsIn)
brho <- as.integer(dat$BRHO_seedsIn)
laca <- as.integer(dat$LACA_seedsIn)
vumy <- as.integer(dat$VUMY_seedsIn)
esca <- as.integer(dat$ESCA_seedsIn)

N <- as.integer(length(Fecundity))

intra <- esca

no_dist_seeds_esca_dry <- stan(file = "Five_species_BH_model.stan", data = c("N", "Fecundity", "intra", "avfa", "brho","esca",
                                                                             "laca", "vumy"),
                               iter = 40000, chains = 4, thin = 3, control = list(adapt_delta = 0.99, max_treedepth = 20),
                               init = initials1)

save(no_dist_seeds_esca_dry, file = "./posteriors/AM_esca_dry_posteriors_tot.rdata")
