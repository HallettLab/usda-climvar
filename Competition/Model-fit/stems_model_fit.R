# Run experimental data with fecundity model

# source data
setwd("/Users/hallett/Repositories/usda-climvar/Competition/Data-analysis")
source("coexistence-model_formatting.R")


library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

setwd("/Users/hallett/Repositories/usda-climvar/Competition/Model-fit")

data <- stemsin.seedsout
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

initials <- list(lambda=10, alpha_avfa=0.03, alpha_brho=0.03, alpha_laca=0.03, alpha_vumy=0.03)
initials1<- list(initials, initials, initials)

stems_avfa_wet <- stan(file = "Generic_four_species_BH_model.stan", data = c("N", "Fecundity", "intra", "avfa", "brho",
                                                   "laca", "vumy"),
                  iter = 3000, chains = 3, thin = 2, control = list(adapt_delta = 0.8, max_treedepth = 10))

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
                 iter = 3000, chains = 3, thin = 2, control = list(adapt_delta = 0.8, max_treedepth = 10))

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
                 iter = 9000, chains = 3, thin = 3, control = list(adapt_delta = 0.95, max_treedepth = 20))

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
                 iter = 9000, chains = 3, thin = 3, control = list(adapt_delta = 0.95, max_treedepth = 20))


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
                 iter = 9000, chains = 3, thin = 3, control = list(adapt_delta = 0.95, max_treedepth = 20))


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
                 iter = 3000, chains = 3, thin = 1, control = list(adapt_delta = 0.95, max_treedepth = 20))


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
                 iter = 9000, chains = 3, thin = 3, control = list(adapt_delta = 0.95, max_treedepth = 20))


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
                 iter = 9000, chains = 3, thin = 3, control = list(adapt_delta = 0.98, max_treedepth = 20))

