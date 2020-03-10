# Run experimental data with fecundity model

# source data
source("./Competition/Data-analysis/coexistence-model_formatting.R")

library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

library(here)

setwd(here("Competition", "Model-fit"))
data <- stemsin.seedsout

initials <- list(lambda=exp(10), alpha_avfa=exp(0.03), alpha_brho=exp(0.03), alpha_laca=exp(0.03), 
                 alpha_vumy=exp(0.03))
initials1<- list(initials, initials, initials)
# ---------------------------------------------------------------------------------------------
# avfa model fitting for wet

dat <- subset(data, species == "AVFA")
dat <- subset(dat, falltreatment == "wet")

Fecundity <- as.integer(round(dat$seedsOut))

avfa <- as.integer(dat$AVFA_stemsIn)
brho <- as.integer(dat$BRHO_stemsIn)
laca <- as.integer(dat$LACA_stemsIn)
vumy <- as.integer(dat$VUMY_stemsIn)

intra <- avfa
remove <- which(intra==0)

Fecundity <- Fecundity[-remove]
N <- as.integer(length(Fecundity))

avfa <- avfa[-remove]
brho <- brho[-remove]
laca <- laca[-remove]
vumy <- vumy[-remove]

intra <- avfa

stems_avfa_wet <- stan(file = "Generic_four_species_BH_model.stan", data = c("N", "Fecundity", "intra", "avfa", "brho",
                                                                             "laca", "vumy"),
                       iter = 9000, chains = 3, thin = 3, control = list(adapt_delta = 0.85, max_treedepth = 10),
                       init = initials1)

save(stems_avfa_wet, file = "stems_avfa_wet_posteriors.rdata")
# ---------------------------------------------------------------------------------------------
# avfa model fitting for dry

dat <- subset(data, species == "AVFA")
dat <- subset(dat, falltreatment == "dry")

Fecundity <- as.integer(round(dat$seedsOut))

avfa <- as.integer(dat$AVFA_stemsIn)
brho <- as.integer(dat$BRHO_stemsIn)
laca <- as.integer(dat$LACA_stemsIn)
vumy <- as.integer(dat$VUMY_stemsIn)

intra <- avfa
remove <- which(intra==0)

Fecundity <- Fecundity[-remove]
N <- as.integer(length(Fecundity))

avfa <- avfa[-remove]
brho <- brho[-remove]
laca <- laca[-remove]
vumy <- vumy[-remove]

intra <- avfa

stems_avfa_dry <- stan(file = "Generic_four_species_BH_model.stan", data = c("N", "Fecundity", "intra", "avfa", "brho",
                                                                       "laca", "vumy"),
                 iter = 9000, chains = 3, thin = 3, control = list(adapt_delta = 0.85, max_treedepth = 10),
                 init = initials1)

save(stems_avfa_dry, file = "stems_avfa_dry_posteriors.rdata")
# ---------------------------------------------------------------------------------------------
# brho model fitting for wet

dat <- subset(data, species == "BRHO")
dat <- subset(dat, falltreatment == "wet")

Fecundity <- as.integer(round(dat$seedsOut))

avfa <- as.integer(dat$AVFA_stemsIn)
brho <- as.integer(dat$BRHO_stemsIn)
laca <- as.integer(dat$LACA_stemsIn)
vumy <- as.integer(dat$VUMY_stemsIn)

intra <- brho
remove <- which(intra==0)

Fecundity <- Fecundity[-remove]
N <- as.integer(length(Fecundity))

avfa <- avfa[-remove]
brho <- brho[-remove]
laca <- laca[-remove]
vumy <- vumy[-remove]

intra <- brho

stems_brho_wet <- stan(file = "Generic_four_species_BH_model.stan", data = c("N", "Fecundity", "intra", "avfa", "brho",
                                                                              "laca", "vumy"),
                        iter = 9000, chains = 3, thin = 3, control = list(adapt_delta = 0.85, max_treedepth = 10),
                        init = initials1)

save(stems_brho_wet, file = "stems_brho_wet_posteriors.rdata")
# ---------------------------------------------------------------------------------------------
# brho model fitting for dry

dat <- subset(data, species == "BRHO")
dat <- subset(dat, falltreatment == "dry")

Fecundity <- as.integer(round(dat$seedsOut))

avfa <- as.integer(dat$AVFA_stemsIn)
brho <- as.integer(dat$BRHO_stemsIn)
laca <- as.integer(dat$LACA_stemsIn)
vumy <- as.integer(dat$VUMY_stemsIn)

intra <- brho
remove <- which(intra==0)

Fecundity <- Fecundity[-remove]
N <- as.integer(length(Fecundity))

avfa <- avfa[-remove]
brho <- brho[-remove]
laca <- laca[-remove]
vumy <- vumy[-remove]

intra <- brho

stems_brho_dry <- stan(file = "Generic_four_species_BH_model.stan", data = c("N", "Fecundity", "intra", "avfa", "brho",
                                                                             "laca", "vumy"),
                       iter = 9000, chains = 3, thin = 3, control = list(adapt_delta = 0.85, max_treedepth = 10),
                       init = initials1)

save(stems_brho_dry, file = "stems_brho_dry_posteriors.rdata")
# ---------------------------------------------------------------------------------------------
# laca model fitting for dry

dat <- subset(data, species == "LACA")
dat <- subset(dat, falltreatment == "dry")

Fecundity <- as.integer(round(dat$seedsOut))

avfa <- as.integer(dat$AVFA_stemsIn)
brho <- as.integer(dat$BRHO_stemsIn)
laca <- as.integer(dat$LACA_stemsIn)
vumy <- as.integer(dat$VUMY_stemsIn)

intra <- laca
remove <- which(intra==0)

Fecundity <- Fecundity[-remove]
N <- as.integer(length(Fecundity))

avfa <- avfa[-remove]
brho <- brho[-remove]
laca <- laca[-remove]
vumy <- vumy[-remove]

intra <- laca

stems_laca_dry <- stan(file = "Generic_four_species_BH_model.stan", data = c("N", "Fecundity", "intra", "avfa", "brho",
                                                                             "laca", "vumy"),
                       iter = 9000, chains = 3, thin = 3, control = list(adapt_delta = 0.85, max_treedepth = 10),
                       init = initials1)

save(stems_laca_dry, file = "stems_laca_dry_posteriors.rdata")

# ---------------------------------------------------------------------------------------------
# laca model fitting for wet

dat <- subset(data, species == "LACA")
dat <- subset(dat, falltreatment == "wet")

Fecundity <- as.integer(round(dat$seedsOut))

avfa <- as.integer(dat$AVFA_stemsIn)
brho <- as.integer(dat$BRHO_stemsIn)
laca <- as.integer(dat$LACA_stemsIn)
vumy <- as.integer(dat$VUMY_stemsIn)

intra <- laca
remove <- which(intra==0)

Fecundity <- Fecundity[-remove]
N <- as.integer(length(Fecundity))

avfa <- avfa[-remove]
brho <- brho[-remove]
laca <- laca[-remove]
vumy <- vumy[-remove]

intra <- laca

stems_laca_wet <- stan(file = "Generic_four_species_BH_model.stan", data = c("N", "Fecundity", "intra", "avfa", "brho",
                                                                             "laca", "vumy"),
                       iter = 9000, chains = 3, thin = 3, control = list(adapt_delta = 0.85, max_treedepth = 15),
                       init = initials1)

save(stems_laca_wet, file = "stems_laca_wet_posteriors.rdata")

# ---------------------------------------------------------------------------------------------
# vumy model fitting for wet

dat <- subset(data, species == "VUMY")
dat <- subset(dat, falltreatment == "wet")

Fecundity <- as.integer(round(dat$seedsOut))

avfa <- as.integer(dat$AVFA_stemsIn)
brho <- as.integer(dat$BRHO_stemsIn)
laca <- as.integer(dat$LACA_stemsIn)
vumy <- as.integer(dat$VUMY_stemsIn)

intra <- vumy
remove <- which(intra==0)

Fecundity <- Fecundity[-remove]
N <- as.integer(length(Fecundity))

avfa <- avfa[-remove]
brho <- brho[-remove]
laca <- laca[-remove]
vumy <- vumy[-remove]

intra <- vumy

stems_vumy_wet <- stan(file = "Generic_four_species_BH_model.stan", data = c("N", "Fecundity", "intra", "avfa", "brho",
                                                                             "laca", "vumy"),
                       iter = 9000, chains = 3, thin = 3, control = list(adapt_delta = 0.85, max_treedepth = 10),
                       init = initials1)

save(stems_vumy_wet, file = "stems_vumy_wet_posteriors.rdata")
# ---------------------------------------------------------------------------------------------
# vumy model fitting for dry

dat <- subset(data, species == "VUMY")
dat <- subset(dat, falltreatment == "dry")

Fecundity <- as.integer(round(dat$seedsOut))

avfa <- as.integer(dat$AVFA_stemsIn)
brho <- as.integer(dat$BRHO_stemsIn)
laca <- as.integer(dat$LACA_stemsIn)
vumy <- as.integer(dat$VUMY_stemsIn)

intra <- vumy
remove <- which(intra==0)

Fecundity <- Fecundity[-remove]
N <- as.integer(length(Fecundity))

avfa <- avfa[-remove]
brho <- brho[-remove]
laca <- laca[-remove]
vumy <- vumy[-remove]

intra <- vumy

stems_vumy_dry <- stan(file = "Generic_four_species_BH_model.stan", data = c("N", "Fecundity", "intra", "avfa", "brho",
                                                                             "laca", "vumy"),
                       iter = 9000, chains = 3, thin = 3, control = list(adapt_delta = 0.85, max_treedepth = 10),
                       init = initials1)

save(stems_vumy_dry, file = "stems_vumy_dry_posteriors.rdata")
