# Run experimental data with fecundity model
setwd('./UO Hallett/Projects/usda-climvar')
# source data
source("./Competition/Data-analysis/coexistence-model_formatting.R")

library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

library(here)

setwd(here("Competition", "Model-fit"))
disturbed <- which(seedsin.seedsout$disturbed == 1)
data <- seedsin.seedsout[-disturbed,]

initials <- list(lambda=250, alpha_avfa=1, alpha_brho=1, alpha_esca=1,
                 alpha_laca=1, alpha_vumy=1, alpha_trhi=1)
initials1<- list(initials, initials, initials, initials)

### AVFA ----
# wet

dat <- subset(data, species == "AVFA")
dat <- subset(dat, falltreatment == "wet")

Fecundity <- as.integer(round(dat$seedsOut))

avfa <- as.integer(dat$AVFA_seedsIn)
brho <- as.integer(dat$BRHO_seedsIn)
laca <- as.integer(dat$LACA_seedsIn)
vumy <- as.integer(dat$VUMY_seedsIn)
esca <- as.integer(dat$ESCA_seedsIn)
trhi <- as.integer(dat$TRHI_seedsIn)

N <- as.integer(length(Fecundity))

intra <- avfa

intra_g <- 0.78

no_dist_seeds_avfa_wet <- stan(file = "Six_species_BH_model.stan", data = c("N", "Fecundity", "intra", "intra_g",
                                                                            "avfa", "brho","esca", "laca", "vumy", "trhi"),
                       iter = 50000, chains = 4, thin = 3, control = list(adapt_delta = 0.95, max_treedepth = 20),
                       init = initials1)

save(no_dist_seeds_avfa_wet, file = "./posteriors/avfa_wet_posteriors_updated_germ.rdata")
# dry

dat <- subset(data, species == "AVFA")
dat <- subset(dat, falltreatment == "dry")

Fecundity <- as.integer(round(dat$seedsOut))

avfa <- as.integer(dat$AVFA_seedsIn)
brho <- as.integer(dat$BRHO_seedsIn)
laca <- as.integer(dat$LACA_seedsIn)
vumy <- as.integer(dat$VUMY_seedsIn)
esca <- as.integer(dat$ESCA_seedsIn)
trhi <- as.integer(dat$TRHI_seedsIn)

N <- as.integer(length(Fecundity))

intra <- avfa

no_dist_seeds_avfa_dry <- stan(file = "Six_species_BH_model.stan", data = c("N", "Fecundity", "intra", "intra_g",
                                                                            "avfa", "brho","esca", "laca", "vumy", "trhi"),
                       iter = 50000, chains = 4, thin = 3, control = list(adapt_delta = 0.95, max_treedepth = 20),
                       init = initials1)

save(no_dist_seeds_avfa_dry, file = "./posteriors/avfa_dry_posteriors_updated_germ.rdata")
### BRHO ----
# wet

dat <- subset(data, species == "BRHO")
dat <- subset(dat, falltreatment == "wet")

Fecundity <- as.integer(round(dat$seedsOut))

avfa <- as.integer(dat$AVFA_seedsIn)
brho <- as.integer(dat$BRHO_seedsIn)
laca <- as.integer(dat$LACA_seedsIn)
vumy <- as.integer(dat$VUMY_seedsIn)
esca <- as.integer(dat$ESCA_seedsIn)
trhi <- as.integer(dat$TRHI_seedsIn)

N <- as.integer(length(Fecundity))

intra <- brho

intra_g <- 0.8

no_dist_seeds_brho_wet <- stan(file = "Six_species_BH_model.stan", data = c("N", "Fecundity", "intra", "intra_g",
                                                                            "avfa", "brho","esca", "laca", "vumy", "trhi"),
                               iter = 50000, chains = 4, thin = 3, control = list(adapt_delta = 0.97, max_treedepth = 20),
                               init = initials1)

save(no_dist_seeds_brho_wet, file = "./posteriors/brho_wet_posteriors_updated_germ.rdata")
# dry

dat <- subset(data, species == "BRHO")
dat <- subset(dat, falltreatment == "dry")

Fecundity <- as.integer(round(dat$seedsOut))

avfa <- as.integer(dat$AVFA_seedsIn)
brho <- as.integer(dat$BRHO_seedsIn)
laca <- as.integer(dat$LACA_seedsIn)
vumy <- as.integer(dat$VUMY_seedsIn)
esca <- as.integer(dat$ESCA_seedsIn)
trhi <- as.integer(dat$TRHI_seedsIn)

N <- as.integer(length(Fecundity))

intra <- brho

intra_g <- 0.8

no_dist_seeds_brho_dry <- stan(file = "Six_species_BH_model.stan", data = c("N", "Fecundity", "intra", "intra_g",
                                                                            "avfa", "brho","esca", "laca", "vumy", "trhi"),
                       iter = 40000, chains = 4, thin = 3, control = list(adapt_delta = 0.99, max_treedepth = 20),
                       init = initials1)

save(no_dist_seeds_brho_dry, file = "./posteriors/brho_dry_posteriors_updated_germ.rdata")
### LACA ----
# dry

dat <- subset(data, species == "LACA")
dat <- subset(dat, falltreatment == "dry")

Fecundity <- as.integer(round(dat$seedsOut))

avfa <- as.integer(dat$AVFA_seedsIn)
brho <- as.integer(dat$BRHO_seedsIn)
laca <- as.integer(dat$LACA_seedsIn)
vumy <- as.integer(dat$VUMY_seedsIn)
esca <- as.integer(dat$ESCA_seedsIn)
trhi <- as.integer(dat$TRHI_seedsIn)

N <- as.integer(length(Fecundity))

intra <- laca

intra_g <- 0.8

no_dist_seeds_laca_dry <- stan(file = "Six_species_BH_model.stan", data = c("N", "Fecundity", "intra", "intra_g",
                                                                            "avfa", "brho","esca", "laca", "vumy", "trhi"),
                       iter = 40000, chains = 4, thin = 3, control = list(adapt_delta = 0.95, max_treedepth = 20), 
                       init = initials1)

save(no_dist_seeds_laca_dry, file = "./posteriors/laca_dry_posteriors_updated_germ.rdata")

# wet

dat <- subset(data, species == "LACA")
dat <- subset(dat, falltreatment == "wet")

Fecundity <- as.integer(round(dat$seedsOut))

avfa <- as.integer(dat$AVFA_seedsIn)
brho <- as.integer(dat$BRHO_seedsIn)
laca <- as.integer(dat$LACA_seedsIn)
vumy <- as.integer(dat$VUMY_seedsIn)
esca <- as.integer(dat$ESCA_seedsIn)
trhi <- as.integer(dat$TRHI_seedsIn)

N <- as.integer(length(Fecundity))

intra <- laca

intra_g <- 0.8

no_dist_seeds_laca_wet <- stan(file = "Six_species_BH_model.stan", data = c("N", "Fecundity", "intra", "intra_g",
                                                                            "avfa", "brho","esca", "laca", "vumy", "trhi"),
                       iter = 40000, chains = 4, thin = 3, control = list(adapt_delta = 0.95, max_treedepth = 20), 
                       init = initials1)

save(no_dist_seeds_laca_wet, file = "./posteriors/laca_wet_posteriors_updated_germ.rdata")
### VUMY ----
# wet

dat <- subset(data, species == "VUMY")
dat <- subset(dat, falltreatment == "wet")

Fecundity <- as.integer(round(dat$seedsOut))

avfa <- as.integer(dat$AVFA_seedsIn)
brho <- as.integer(dat$BRHO_seedsIn)
laca <- as.integer(dat$LACA_seedsIn)
vumy <- as.integer(dat$VUMY_seedsIn)
esca <- as.integer(dat$ESCA_seedsIn)
trhi <- as.integer(dat$TRHI_seedsIn)

N <- as.integer(length(Fecundity))

intra <- vumy

intra_g <- 0.6

no_dist_seeds_vumy_wet <- stan(file = "Six_species_BH_model.stan", data = c("N", "Fecundity", "intra", "intra_g",
                                                                            "avfa", "brho","esca", "laca", "vumy", "trhi"),
                       iter = 40000, chains = 4, thin = 3, control = list(adapt_delta = 0.97, max_treedepth = 20),
                       init = initials1)

save(no_dist_seeds_vumy_wet, file = "./posteriors/vumy_wet_posteriors_updated_germ.rdata")
# dry

dat <- subset(data, species == "VUMY")
dat <- subset(dat, falltreatment == "dry")

Fecundity <- as.integer(round(dat$seedsOut))

avfa <- as.integer(dat$AVFA_seedsIn)
brho <- as.integer(dat$BRHO_seedsIn)
laca <- as.integer(dat$LACA_seedsIn)
vumy <- as.integer(dat$VUMY_seedsIn)
esca <- as.integer(dat$ESCA_seedsIn)
trhi <- as.integer(dat$TRHI_seedsIn)

N <- as.integer(length(Fecundity))

intra <- vumy

intra_g <- 0.6

no_dist_seeds_vumy_dry <- stan(file = "Six_species_BH_model.stan", data = c("N", "Fecundity", "intra", "intra_g",
                                                                            "avfa", "brho","esca", "laca", "vumy", "trhi"),
                       iter = 40000, chains = 4, thin = 3, control = list(adapt_delta = 0.97, max_treedepth = 20),
                       init = initials1)

save(no_dist_seeds_vumy_dry, file = "./posteriors/vumy_dry_posteriors_updated_germ.rdata")

### ESCA ----
# wet

dat <- subset(data, species == "ESCA")
dat <- subset(dat, falltreatment == "wet")

Fecundity <- as.integer(round(dat$seedsOut))

avfa <- as.integer(dat$AVFA_seedsIn)
brho <- as.integer(dat$BRHO_seedsIn)
laca <- as.integer(dat$LACA_seedsIn)
vumy <- as.integer(dat$VUMY_seedsIn)
esca <- as.integer(dat$ESCA_seedsIn)
trhi <- as.integer(dat$TRHI_seedsIn)

N <- as.integer(length(Fecundity))

intra <- esca

intra_g <- 0.95

no_dist_seeds_esca_wet <- stan(file = "Six_species_BH_model.stan", data = c("N", "Fecundity", "intra", "intra_g",
                                                                            "avfa", "brho","esca", "laca", "vumy", "trhi"),
                               iter = 40000, chains = 4, thin = 3, control = list(adapt_delta = 0.99, max_treedepth = 20),
                               init = initials1)

save(no_dist_seeds_esca_wet, file = "./posteriors/esca_wet_posteriors_updated_germ.rdata")

# dry

dat <- subset(data, species == "ESCA")
dat <- subset(dat, falltreatment == "dry")

Fecundity <- as.integer(round(dat$seedsOut))

avfa <- as.integer(dat$AVFA_seedsIn)
brho <- as.integer(dat$BRHO_seedsIn)
laca <- as.integer(dat$LACA_seedsIn)
vumy <- as.integer(dat$VUMY_seedsIn)
esca <- as.integer(dat$ESCA_seedsIn)
trhi <- as.integer(dat$TRHI_seedsIn)

N <- as.integer(length(Fecundity))

intra <- esca

intra_g <- 0.95

no_dist_seeds_esca_dry <- stan(file = "Six_species_BH_model.stan", data = c("N", "Fecundity", "intra", "intra_g",
                                                                            "avfa", "brho","esca", "laca", "vumy", "trhi"),
                               iter = 40000, chains = 4, thin = 3, control = list(adapt_delta = 0.97, max_treedepth = 20),
                               init = initials1)

save(no_dist_seeds_esca_dry, file = "./posteriors/esca_dry_posteriors_updated_germ.rdata")

### TRHI ----
# wet

dat <- subset(data, species == "TRHI")
dat <- subset(dat, falltreatment == "wet")

Fecundity <- as.integer(round(dat$seedsOut))

avfa <- as.integer(dat$AVFA_seedsIn)
brho <- as.integer(dat$BRHO_seedsIn)
laca <- as.integer(dat$LACA_seedsIn)
vumy <- as.integer(dat$VUMY_seedsIn)
esca <- as.integer(dat$ESCA_seedsIn)
trhi <- as.integer(dat$TRHI_seedsIn)

N <- as.integer(length(Fecundity))

intra <- trhi

intra_g <- 0.2

no_dist_seeds_trhi_wet <- stan(file = "Six_species_BH_model.stan", data = c("N", "Fecundity", "intra", "intra_g",
                                                                            "avfa", "brho","esca", "laca", "vumy", "trhi"),
                               iter = 40000, chains = 4, thin = 3, control = list(adapt_delta = 0.98, max_treedepth = 20),
                               init = initials1)

save(no_dist_seeds_trhi_wet, file = "./posteriors/trhi_wet_posteriors_updated_germ.rdata")

# dry

dat <- subset(data, species == "TRHI")
dat <- subset(dat, falltreatment == "dry")

Fecundity <- as.integer(round(dat$seedsOut))

avfa <- as.integer(dat$AVFA_seedsIn)
brho <- as.integer(dat$BRHO_seedsIn)
laca <- as.integer(dat$LACA_seedsIn)
vumy <- as.integer(dat$VUMY_seedsIn)
esca <- as.integer(dat$ESCA_seedsIn)
trhi <- as.integer(dat$TRHI_seedsIn)

N <- as.integer(length(Fecundity))

intra <- trhi

intra_g <- 0.2

initials <- list(lambda=250, alpha_avfa=2, alpha_brho=0.5, alpha_esca=0.5,
                 alpha_laca=0.5, alpha_vumy=0.5, alpha_trhi=0.5)
initials1<- list(initials, initials, initials, initials)

no_dist_seeds_trhi_dry <- stan(file = "TRHI_Six_species_BH_model.stan", data = c("N", "Fecundity", "intra", "intra_g",
                                                                                 "avfa", "brho","esca", "laca", "vumy", "trhi"),
                               iter = 200000, chains = 4, thin = 3, control = list(adapt_delta = 0.995, max_treedepth = 20),
                               init = initials1)

save(no_dist_seeds_trhi_dry, file = "./posteriors/trhi_dry_posteriors_updated_germ.rdata")

