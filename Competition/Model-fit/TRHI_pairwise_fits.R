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

initials <- list(lambda=250, alpha_trhi=1, alpha_inter=1)
initials1<- list(initials, initials, initials, initials)

dat <- data[data$species == "TRHI" & data$falltreatment == "dry",]
dat2 <- data[data$species == "TRHI" & data$falltreatment == "wet",]

# AVFA 

dat_avfa <- dat[dat$background == "Trifolium hirtum" | dat$background == "Avena fatua",]
dat_avfa <- dat_avfa[!is.na(dat_avfa$seedsOut),]

Fecundity <- as.integer(round(dat_avfa$seedsOut))
trhi <- as.integer(dat_avfa$TRHI_seedsIn)
inter <- as.integer(dat_avfa$AVFA_seedsIn)

N <- as.integer(length(Fecundity))

trhi_avfa <- stan(file = "TRHI_pairwise.stan", data = c("N", "Fecundity", "trhi", "inter"),
                               iter = 40000, chains = 4, thin = 3, control = list(adapt_delta = 0.99, max_treedepth = 20),
                               init = initials1)

# BRHO 

dat_brho <- dat[dat$background == "Trifolium hirtum" | dat$background == "Bromus hordeaceus",]
dat_brho <- dat_brho[!is.na(dat_brho$seedsOut),]

Fecundity <- as.integer(round(dat_brho$seedsOut))
trhi <- as.integer(dat_brho$TRHI_seedsIn)
inter <- as.integer(dat_brho$BRHO_seedsIn)

N <- as.integer(length(Fecundity))

trhi_brho <- stan(file = "TRHI_pairwise.stan", data = c("N", "Fecundity", "trhi", "inter"),
                  iter = 500000, chains = 4, thin = 3, control = list(adapt_delta = 0.99, max_treedepth = 20),
                  init = initials1)

# VUMY 

dat_vumy <- dat[dat$background == "Trifolium hirtum" | dat$background == "Vulpia myuros",]
dat_vumy <- dat_vumy[!is.na(dat_vumy$seedsOut),]

Fecundity <- as.integer(round(dat_vumy$seedsOut))
trhi <- as.integer(dat_vumy$TRHI_seedsIn)
inter <- as.integer(dat_vumy$VUMY_seedsIn)

N <- as.integer(length(Fecundity))

trhi_vumy <- stan(file = "TRHI_pairwise.stan", data = c("N", "Fecundity", "trhi", "inter"),
                  iter = 1000000, chains = 4, thin = 3, control = list(adapt_delta = 0.99, max_treedepth = 20),
                  init = initials1)

# LACA 

dat_laca <- dat[dat$background == "Trifolium hirtum" | dat$background == "Lasthenia californica",]
dat_laca <- dat_laca[!is.na(dat_laca$seedsOut),]

Fecundity <- as.integer(round(dat_laca$seedsOut))
trhi <- as.integer(dat_laca$TRHI_seedsIn)
inter <- as.integer(dat_laca$LACA_seedsIn)

N <- as.integer(length(Fecundity))

trhi_laca <- stan(file = "TRHI_pairwise.stan", data = c("N", "Fecundity", "trhi", "inter"),
                  iter = 40000, chains = 4, thin = 3, control = list(adapt_delta = 0.99, max_treedepth = 20),
                  init = initials1)

# ESCA 

dat_esca <- dat[dat$background == "Trifolium hirtum" | dat$background == "Eschscholzia californica",]
dat_esca <- dat_esca[!is.na(dat_esca$seedsOut),]

Fecundity <- as.integer(round(dat_esca$seedsOut))
trhi <- as.integer(dat_esca$TRHI_seedsIn)
inter <- as.integer(dat_esca$ESCA_seedsIn)

N <- as.integer(length(Fecundity))

trhi_esca <- stan(file = "TRHI_pairwise.stan", data = c("N", "Fecundity", "trhi", "inter"),
                  iter = 40000, chains = 4, thin = 3, control = list(adapt_delta = 0.99, max_treedepth = 20),
                  init = initials1)