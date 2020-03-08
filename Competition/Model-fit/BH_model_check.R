# generate data

do.BH <- function(lambda, alpha_avfa, alpha_brho, alpha_esca, alpha_laca, alpha_trhi, alpha_vumy,
                  avfa, brho, esca, laca, trhi, vumy, focals) {
  Ntp1_avfa <- rpois(length(avfa), lambda*focals[1:50]/
                       (1+ alpha_avfa * (avfa + focals[1:50])))
  Ntp1_brho <- rpois(length(brho), lambda*focals[51:100]/
                       (1+ alpha_avfa * focals[51:100] + alpha_brho * brho))
  Ntp1_esca <- rpois(length(esca), lambda*focals[101:150]/
                       (1+ alpha_avfa * focals[101:150] + alpha_esca * esca))
  Ntp1_laca <- rpois(length(esca), lambda*focals[151:200]/
                       (1+ alpha_avfa * focals[151:200] + alpha_laca * laca))
  Ntp1_trhi <- rpois(length(trhi), lambda*focals[201:250]/
                       (1+ alpha_avfa * focals[201:250] + alpha_trhi * trhi))
  Ntp1_vumy <- rpois(length(vumy), lambda*focals[251:300]/
                       (1+ alpha_avfa * focals[251:300] + alpha_vumy * vumy))
  
  return(c(Ntp1_avfa, Ntp1_brho, Ntp1_esca, Ntp1_laca, Ntp1_trhi, Ntp1_vumy))
}

lambda <- 15
alpha_avfa <- .1
alpha_brho <- .02
alpha_esca <- .1
alpha_laca <- .01
alpha_trhi <- .01
alpha_vumy <- .09

focals <- sample.int(3, size=50*6, replace=TRUE)
avfa_int <- sample.int(200, size=50, replace=TRUE)
brho_int <- sample.int(200, size=50, replace=TRUE)
esca_int <- sample.int(200, size=50, replace=TRUE)
laca_int <- sample.int(200, size=50, replace=TRUE)
trhi_int <- sample.int(200, size=50, replace=TRUE)
vumy_int <- sample.int(200, size=50, replace=TRUE)

#avfa is focal here
intra <- focals

Fecundity <- do.BH(lambda = lambda, alpha_avfa = alpha_avfa, alpha_brho = alpha_brho, 
                   alpha_esca = alpha_esca, alpha_laca = alpha_laca, alpha_trhi = alpha_trhi, 
                   alpha_vumy = alpha_vumy, avfa = avfa_int, brho = brho_int, esca = esca_int, 
                   laca = laca_int, trhi = trhi_int, vumy = vumy_int, focals = focals)

Fecundity <- as.integer(Fecundity)

avfa <- intra
avfa[1:50] <- avfa[1:50] + avfa_int

brho <- as.integer(rep(0, length(focals)))
brho[51:100] <- brho_int

esca <- as.integer(rep(0, length(focals)))
esca[101:150] <- esca_int

laca <- as.integer(rep(0, length(focals)))
laca[151:200] <- laca_int

trhi <- as.integer(rep(0, length(focals)))
trhi[201:250] <- trhi_int

vumy <- as.integer(rep(0, length(focals)))
vumy[251:300] <- vumy_int

# -----------------------------------------
# Check stan model
N <- as.integer(300)

library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

setwd("/Users/lash1937/Documents/GitHub/Usda-climvar/Competition/Model-fit")
PrelimFit <- stan(file = "BH_model.stan", data = c("N", "Fecundity", "intra", "avfa", "brho",
                                                   "esca", "laca", "trhi", "vumy"),
                  iter = 3000, chains = 1, thin = 1)
  # control = list(adapt_delta = 0.8, max_treedepth = 10)
