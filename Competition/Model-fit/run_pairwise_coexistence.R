# run coexistence models for dry versus wet

source("./Competition/Model-fit/import_posteriors_AM.R")

run.to.equilibrium <- function(surv, germ, lambda, alpha_intra, Nt) {
  Ntp1 <- (1-germ)*surv*Nt + germ*lambda*Nt/(1+ alpha_intra * Nt)
  return(Ntp1)

}

run.invader <- function(surv, germ, lambda, alpha_inter, resid_abund, invader_abund) {
  Ntp1 <- (1-germ)*surv*invader_abund + germ*lambda*invader_abund/(1+ alpha_inter * resid_abund)
  LDGR <- log(Ntp1/invader_abund)
  return(LDGR)
  
}

avfa_dry$germ <- .78
avfa_wet$germ <- .78
avfa_dry$surv <- .01
avfa_wet$surv <- .01

brho_dry$germ <- .8
brho_wet$germ <- .8
brho_dry$surv <- .013
brho_wet$surv <- .013

vumy_dry$germ <- .6
vumy_wet$germ <- .6
vumy_dry$surv <- .045
vumy_wet$surv <- .045

laca_dry$germ <- .8
laca_wet$germ <- .8
laca_dry$surv <- .01
laca_wet$surv <- .01

esca_dry$germ <- .95
esca_wet$germ <- .95
esca_dry$surv <- .01
esca_wet$surv <- .01

trhi_dry$germ <- .2
trhi_wet$germ <- .2
trhi_dry$surv <- .01
trhi_wet$surv <- .01

all_datset <- list(avfa_dry, avfa_wet, brho_dry, brho_wet, vumy_dry, vumy_wet, 
                   laca_dry, laca_wet, esca_dry, esca_wet, trhi_dry, trhi_wet) # set which dataset we want to do
all_intra <- c("alpha_avfa", "alpha_avfa", "alpha_brho", "alpha_brho", "alpha_vumy", 
               "alpha_vumy", "alpha_laca", "alpha_laca", "alpha_esca", "alpha_esca", 
               "alpha_trhi", "alpha_trhi")
options <- length(all_intra)

time <- 200
runs <- 200

N <- array(NA, c(options, runs, time))
N[,,1] <- 100

for (x in 1:options) {
  datset <- all_datset[[x]]
  intra <- all_intra[[x]]
  
  post_length <- length(datset$lambda)
  
  all_intras <- datset[[intra]]
  
  for (t in 1:(time-1)) {
    posts <- sample(post_length, runs, replace=TRUE)
    lambda <- datset$lambda[posts]
    alpha_intra <- all_intras[posts]
    N[x, ,t+1] <- run.to.equilibrium(surv=datset$surv, germ=datset$germ, 
                                     lambda=lambda, alpha_intra=alpha_intra, Nt=N[x, ,t])
  }
}

species <- c("avfa", "brho", "vumy", "laca", "esca", "trhi")

# make into dry versus wet conditions dataframes
residents_dry <-  data.frame(N[1,,200], N[3,,200], N[5,,200], N[7,,200], N[9,,200], N[11,,200])
names(residents_dry) <- species

residents_wet <-  data.frame(N[2,,200], N[4,,200], N[6,,200], N[8,,200], N[10,,200], N[12,,200])
names(residents_wet) <- species

# invade into residents
# -------------------------------------------------------------------------------------------------------
# invader into dry  
reps <- 200
avfa_into_brho_dry <- matrix(NA, reps, runs)
avfa_into_vumy_dry <- matrix(NA, reps, runs)
avfa_into_laca_dry <- matrix(NA, reps, runs)
avfa_into_esca_dry <- matrix(NA, reps, runs)
avfa_into_trhi_dry <- matrix(NA, reps, runs)

brho_into_avfa_dry <- matrix(NA, reps, runs)
brho_into_vumy_dry <- matrix(NA, reps, runs)
brho_into_laca_dry <- matrix(NA, reps, runs)
brho_into_esca_dry <- matrix(NA, reps, runs)
brho_into_trhi_dry <- matrix(NA, reps, runs)

vumy_into_avfa_dry <- matrix(NA, reps, runs)
vumy_into_brho_dry <- matrix(NA, reps, runs)
vumy_into_laca_dry <- matrix(NA, reps, runs)
vumy_into_esca_dry <- matrix(NA, reps, runs)
vumy_into_trhi_dry <- matrix(NA, reps, runs)

laca_into_avfa_dry <- matrix(NA, reps, runs)
laca_into_brho_dry <- matrix(NA, reps, runs)
laca_into_vumy_dry <- matrix(NA, reps, runs)
laca_into_esca_dry <- matrix(NA, reps, runs)
laca_into_trhi_dry <- matrix(NA, reps, runs)

esca_into_avfa_dry <- matrix(NA, reps, runs)
esca_into_brho_dry <- matrix(NA, reps, runs)
esca_into_laca_dry <- matrix(NA, reps, runs)
esca_into_vumy_dry <- matrix(NA, reps, runs)
esca_into_trhi_dry <- matrix(NA, reps, runs)

trhi_into_avfa_dry <- matrix(NA, reps, runs)
trhi_into_brho_dry <- matrix(NA, reps, runs)
trhi_into_laca_dry <- matrix(NA, reps, runs)
trhi_into_vumy_dry <- matrix(NA, reps, runs)
trhi_into_esca_dry <- matrix(NA, reps, runs)




invader_abund <- 1

# avena invades
post_length <- length(avfa_dry$lambda)
for (r in 1:reps) {
  
  posts <- sample(post_length, runs, replace=TRUE)
  
  avfa_into_brho_dry[r,] <- run.invader(surv=avfa_dry$surv, germ= avfa_dry$germ, 
                                        lambda=avfa_dry$lambda[posts], alpha_inter=avfa_dry$alpha_brho[posts],
                                        resid_abund=residents_dry$brho, invader_abund=invader_abund)
  
  avfa_into_vumy_dry[r,] <- run.invader(surv=avfa_dry$surv, germ= avfa_dry$germ, 
                                        lambda=avfa_dry$lambda[posts], alpha_inter=avfa_dry$alpha_vumy[posts],
                                        resid_abund=residents_dry$vumy, invader_abund=invader_abund)
  
  avfa_into_laca_dry[r,] <- run.invader(surv=avfa_dry$surv, germ= avfa_dry$germ, 
                                        lambda=avfa_dry$lambda[posts], alpha_inter=avfa_dry$alpha_laca[posts],
                                        resid_abund=residents_dry$laca, invader_abund=invader_abund)
  
  avfa_into_esca_dry[r,] <- run.invader(surv=avfa_dry$surv, germ= avfa_dry$germ, 
                                        lambda=avfa_dry$lambda[posts], alpha_inter=avfa_dry$alpha_esca[posts],
                                        resid_abund=residents_dry$esca, invader_abund=invader_abund)
  
  avfa_into_trhi_dry[r,] <- run.invader(surv=avfa_dry$surv, germ= avfa_dry$germ, 
                                        lambda=avfa_dry$lambda[posts], alpha_inter=avfa_dry$alpha_trhi[posts],
                                        resid_abund=residents_dry$trhi, invader_abund=invader_abund)
}

# brho invades
post_length <- length(brho_dry$lambda)
for (r in 1:reps) {
  posts <- sample(post_length, runs, replace=TRUE)
  
  brho_into_avfa_dry[r,] <- run.invader(surv=brho_dry$surv, germ= brho_dry$germ, 
                                        lambda=brho_dry$lambda[posts], alpha_inter=brho_dry$alpha_avfa[posts],
                                        resid_abund=residents_dry$avfa, invader_abund=invader_abund)
  
  brho_into_vumy_dry[r,] <- run.invader(surv=brho_dry$surv, germ= brho_dry$germ, 
                                        lambda=brho_dry$lambda[posts], alpha_inter=brho_dry$alpha_vumy[posts],
                                        resid_abund=residents_dry$vumy, invader_abund=invader_abund)
  
  brho_into_laca_dry[r,] <- run.invader(surv=brho_dry$surv, germ= brho_dry$germ, 
                                        lambda=brho_dry$lambda[posts], alpha_inter=brho_dry$alpha_laca[posts],
                                        resid_abund=residents_dry$laca, invader_abund=invader_abund)
  
  brho_into_esca_dry[r,] <- run.invader(surv=brho_dry$surv, germ= brho_dry$germ, 
                                        lambda=brho_dry$lambda[posts], alpha_inter=brho_dry$alpha_esca[posts],
                                        resid_abund=residents_dry$esca, invader_abund=invader_abund)
  
  brho_into_trhi_dry[r,] <- run.invader(surv=brho_dry$surv, germ= brho_dry$germ, 
                                        lambda=brho_dry$lambda[posts], alpha_inter=brho_dry$alpha_trhi[posts],
                                        resid_abund=residents_dry$trhi, invader_abund=invader_abund)
  
}

# Vumy invades
post_length <- length(vumy_dry$lambda)
for (r in 1:reps) {
  posts <- sample(post_length, runs, replace=TRUE)
  
  vumy_into_avfa_dry[r,] <- run.invader(surv=vumy_dry$surv, germ= vumy_dry$germ, 
                                        lambda=vumy_dry$lambda[posts], alpha_inter=vumy_dry$alpha_avfa[posts],
                                        resid_abund=residents_dry$avfa, invader_abund=invader_abund)
  
  vumy_into_brho_dry[r,] <- run.invader(surv=vumy_dry$surv, germ= vumy_dry$germ, 
                                        lambda=vumy_dry$lambda[posts], alpha_inter=vumy_dry$alpha_brho[posts],
                                        resid_abund=residents_dry$brho, invader_abund=invader_abund)
  
  vumy_into_laca_dry[r,] <- run.invader(surv=vumy_dry$surv, germ= vumy_dry$germ, 
                                        lambda=vumy_dry$lambda[posts], alpha_inter=vumy_dry$alpha_laca[posts],
                                        resid_abund=residents_dry$laca, invader_abund=invader_abund)
  
  vumy_into_esca_dry[r,] <- run.invader(surv=vumy_dry$surv, germ= vumy_dry$germ, 
                                        lambda=vumy_dry$lambda[posts], alpha_inter=vumy_dry$alpha_esca[posts],
                                        resid_abund=residents_dry$esca, invader_abund=invader_abund)
  
  vumy_into_trhi_dry[r,] <- run.invader(surv=vumy_dry$surv, germ= vumy_dry$germ, 
                                        lambda=vumy_dry$lambda[posts], alpha_inter=vumy_dry$alpha_trhi[posts],
                                        resid_abund=residents_dry$trhi, invader_abund=invader_abund)
  
}

# Laca invades
post_length <- length(laca_dry$lambda)
for (r in 1:reps) {
  posts <- sample(post_length, runs, replace=TRUE)
  
  laca_into_avfa_dry[r,] <- run.invader(surv=laca_dry$surv, germ= laca_dry$germ, 
                                        lambda=laca_dry$lambda[posts], alpha_inter=laca_dry$alpha_avfa[posts],
                                        resid_abund=residents_dry$avfa, invader_abund=invader_abund)
  
  laca_into_brho_dry[r,] <- run.invader(surv=laca_dry$surv, germ= laca_dry$germ, 
                                        lambda=laca_dry$lambda[posts], alpha_inter=laca_dry$alpha_brho[posts],
                                        resid_abund=residents_dry$brho, invader_abund=invader_abund)
  
  laca_into_vumy_dry[r,] <- run.invader(surv=laca_dry$surv, germ= laca_dry$germ, 
                                        lambda=laca_dry$lambda[posts], alpha_inter=laca_dry$alpha_vumy[posts],
                                        resid_abund=residents_dry$vumy, invader_abund=invader_abund)
  
  laca_into_esca_dry[r,] <- run.invader(surv=laca_dry$surv, germ= laca_dry$germ, 
                                        lambda=laca_dry$lambda[posts], alpha_inter=laca_dry$alpha_esca[posts],
                                        resid_abund=residents_dry$esca, invader_abund=invader_abund)
  
  laca_into_trhi_dry[r,] <- run.invader(surv=laca_dry$surv, germ= laca_dry$germ, 
                                        lambda=laca_dry$lambda[posts], alpha_inter=laca_dry$alpha_trhi[posts],
                                        resid_abund=residents_dry$trhi, invader_abund=invader_abund)
  
}

# Esca invades
post_length <- length(esca_dry$lambda)
for (r in 1:reps) {
  posts <- sample(post_length, runs, replace=TRUE)
  
  esca_into_avfa_dry[r,] <- run.invader(surv=esca_dry$surv, germ= esca_dry$germ, 
                                        lambda=esca_dry$lambda[posts], alpha_inter=esca_dry$alpha_avfa[posts],
                                        resid_abund=residents_dry$avfa, invader_abund=invader_abund)
  
  esca_into_brho_dry[r,] <- run.invader(surv=esca_dry$surv, germ= esca_dry$germ, 
                                        lambda=esca_dry$lambda[posts], alpha_inter=esca_dry$alpha_brho[posts],
                                        resid_abund=residents_dry$brho, invader_abund=invader_abund)
  
  esca_into_vumy_dry[r,] <- run.invader(surv=esca_dry$surv, germ= esca_dry$germ, 
                                        lambda=esca_dry$lambda[posts], alpha_inter=esca_dry$alpha_vumy[posts],
                                        resid_abund=residents_dry$vumy, invader_abund=invader_abund)
  
  esca_into_laca_dry[r,] <- run.invader(surv=esca_dry$surv, germ= esca_dry$germ, 
                                        lambda=esca_dry$lambda[posts], alpha_inter=esca_dry$alpha_laca[posts],
                                        resid_abund=residents_dry$laca, invader_abund=invader_abund)
  
  esca_into_trhi_dry[r,] <- run.invader(surv=esca_dry$surv, germ= esca_dry$germ, 
                                        lambda=esca_dry$lambda[posts], alpha_inter=esca_dry$alpha_trhi[posts],
                                        resid_abund=residents_dry$trhi, invader_abund=invader_abund)
  
}

# trhi invades
post_length <- length(trhi_dry$lambda)
for (r in 1:reps) {
  posts <- sample(post_length, runs, replace=TRUE)
  
  trhi_into_avfa_dry[r,] <- run.invader(surv=trhi_dry$surv, germ= trhi_dry$germ, 
                                        lambda=trhi_dry$lambda[posts], alpha_inter=trhi_dry$alpha_avfa[posts],
                                        resid_abund=residents_dry$avfa, invader_abund=invader_abund)
  
  trhi_into_brho_dry[r,] <- run.invader(surv=trhi_dry$surv, germ= trhi_dry$germ, 
                                        lambda=trhi_dry$lambda[posts], alpha_inter=trhi_dry$alpha_brho[posts],
                                        resid_abund=residents_dry$brho, invader_abund=invader_abund)
  
  trhi_into_vumy_dry[r,] <- run.invader(surv=trhi_dry$surv, germ= trhi_dry$germ, 
                                        lambda=trhi_dry$lambda[posts], alpha_inter=trhi_dry$alpha_vumy[posts],
                                        resid_abund=residents_dry$vumy, invader_abund=invader_abund)
  
  trhi_into_laca_dry[r,] <- run.invader(surv=trhi_dry$surv, germ= trhi_dry$germ, 
                                        lambda=trhi_dry$lambda[posts], alpha_inter=trhi_dry$alpha_laca[posts],
                                        resid_abund=residents_dry$laca, invader_abund=invader_abund)
  
  trhi_into_esca_dry[r,] <- run.invader(surv=trhi_dry$surv, germ= trhi_dry$germ, 
                                        lambda=trhi_dry$lambda[posts], alpha_inter=trhi_dry$alpha_esca[posts],
                                        resid_abund=residents_dry$esca, invader_abund=invader_abund)
  
}

# -------------------------------------------------------------------------------------------------------
# invader into wet  
reps <- 200
avfa_into_brho_wet <- matrix(NA, reps, runs)
avfa_into_vumy_wet <- matrix(NA, reps, runs)
avfa_into_laca_wet <- matrix(NA, reps, runs)
avfa_into_esca_wet <- matrix(NA, reps, runs)
avfa_into_trhi_wet <- matrix(NA, reps, runs)

brho_into_avfa_wet <- matrix(NA, reps, runs)
brho_into_vumy_wet <- matrix(NA, reps, runs)
brho_into_laca_wet <- matrix(NA, reps, runs)
brho_into_esca_wet <- matrix(NA, reps, runs)
brho_into_trhi_wet <- matrix(NA, reps, runs)

vumy_into_avfa_wet <- matrix(NA, reps, runs)
vumy_into_brho_wet <- matrix(NA, reps, runs)
vumy_into_laca_wet <- matrix(NA, reps, runs)
vumy_into_esca_wet <- matrix(NA, reps, runs)
vumy_into_trhi_wet <- matrix(NA, reps, runs)

laca_into_avfa_wet <- matrix(NA, reps, runs)
laca_into_brho_wet <- matrix(NA, reps, runs)
laca_into_vumy_wet <- matrix(NA, reps, runs)
laca_into_esca_wet <- matrix(NA, reps, runs)
laca_into_trhi_wet <- matrix(NA, reps, runs)

esca_into_avfa_wet <- matrix(NA, reps, runs)
esca_into_brho_wet <- matrix(NA, reps, runs)
esca_into_laca_wet <- matrix(NA, reps, runs)
esca_into_vumy_wet <- matrix(NA, reps, runs)
esca_into_trhi_wet <- matrix(NA, reps, runs)

trhi_into_avfa_wet <- matrix(NA, reps, runs)
trhi_into_brho_wet <- matrix(NA, reps, runs)
trhi_into_laca_wet <- matrix(NA, reps, runs)
trhi_into_vumy_wet <- matrix(NA, reps, runs)
trhi_into_esca_wet <- matrix(NA, reps, runs)


invader_abund <- 1

# avena invades
post_length <- length(avfa_wet$lambda)
for (r in 1:reps) {
  
  posts <- sample(post_length, runs, replace=TRUE)
  
  avfa_into_brho_wet[r,] <- run.invader(surv=avfa_wet$surv, germ= avfa_wet$germ, 
                                        lambda=avfa_wet$lambda[posts], alpha_inter=avfa_wet$alpha_brho[posts],
                                        resid_abund=residents_wet$brho, invader_abund=invader_abund)
  
  avfa_into_vumy_wet[r,] <- run.invader(surv=avfa_wet$surv, germ= avfa_wet$germ, 
                                        lambda=avfa_wet$lambda[posts], alpha_inter=avfa_wet$alpha_vumy[posts],
                                        resid_abund=residents_wet$vumy, invader_abund=invader_abund)
  
  avfa_into_laca_wet[r,] <- run.invader(surv=avfa_wet$surv, germ= avfa_wet$germ, 
                                        lambda=avfa_wet$lambda[posts], alpha_inter=avfa_wet$alpha_laca[posts],
                                        resid_abund=residents_wet$laca, invader_abund=invader_abund)
  
  avfa_into_esca_wet[r,] <- run.invader(surv=avfa_wet$surv, germ= avfa_wet$germ, 
                                        lambda=avfa_wet$lambda[posts], alpha_inter=avfa_wet$alpha_esca[posts],
                                        resid_abund=residents_wet$esca, invader_abund=invader_abund)
  
  avfa_into_trhi_wet[r,] <- run.invader(surv=avfa_wet$surv, germ= avfa_wet$germ, 
                                        lambda=avfa_wet$lambda[posts], alpha_inter=avfa_wet$alpha_trhi[posts],
                                        resid_abund=residents_wet$trhi, invader_abund=invader_abund)
}

# brho invades
post_length <- length(brho_wet$lambda)
for (r in 1:reps) {
  posts <- sample(post_length, runs, replace=TRUE)
  
  brho_into_avfa_wet[r,] <- run.invader(surv=brho_wet$surv, germ= brho_wet$germ, 
                                        lambda=brho_wet$lambda[posts], alpha_inter=brho_wet$alpha_avfa[posts],
                                        resid_abund=residents_wet$avfa, invader_abund=invader_abund)
  
  brho_into_vumy_wet[r,] <- run.invader(surv=brho_wet$surv, germ= brho_wet$germ, 
                                        lambda=brho_wet$lambda[posts], alpha_inter=brho_wet$alpha_vumy[posts],
                                        resid_abund=residents_wet$vumy, invader_abund=invader_abund)
  
  brho_into_laca_wet[r,] <- run.invader(surv=brho_wet$surv, germ= brho_wet$germ, 
                                        lambda=brho_wet$lambda[posts], alpha_inter=brho_wet$alpha_laca[posts],
                                        resid_abund=residents_wet$laca, invader_abund=invader_abund)
  
  brho_into_esca_wet[r,] <- run.invader(surv=brho_wet$surv, germ= brho_wet$germ, 
                                        lambda=brho_wet$lambda[posts], alpha_inter=brho_wet$alpha_esca[posts],
                                        resid_abund=residents_wet$esca, invader_abund=invader_abund)
  
  brho_into_trhi_wet[r,] <- run.invader(surv=brho_wet$surv, germ= brho_wet$germ, 
                                        lambda=brho_wet$lambda[posts], alpha_inter=brho_wet$alpha_trhi[posts],
                                        resid_abund=residents_wet$trhi, invader_abund=invader_abund)
}

# Vumy invades
post_length <- length(vumy_wet$lambda)
for (r in 1:reps) {
  posts <- sample(post_length, runs, replace=TRUE)
  
  vumy_into_avfa_wet[r,] <- run.invader(surv=vumy_wet$surv, germ= vumy_wet$germ, 
                                        lambda=vumy_wet$lambda[posts], alpha_inter=vumy_wet$alpha_avfa[posts],
                                        resid_abund=residents_wet$avfa, invader_abund=invader_abund)
  
  vumy_into_brho_wet[r,] <- run.invader(surv=vumy_wet$surv, germ= vumy_wet$germ, 
                                        lambda=vumy_wet$lambda[posts], alpha_inter=vumy_wet$alpha_brho[posts],
                                        resid_abund=residents_wet$brho, invader_abund=invader_abund)
  
  vumy_into_laca_wet[r,] <- run.invader(surv=vumy_wet$surv, germ= vumy_wet$germ, 
                                        lambda=vumy_wet$lambda[posts], alpha_inter=vumy_wet$alpha_laca[posts],
                                        resid_abund=residents_wet$laca, invader_abund=invader_abund)
  
  vumy_into_esca_wet[r,] <- run.invader(surv=vumy_wet$surv, germ= vumy_wet$germ, 
                                        lambda=vumy_wet$lambda[posts], alpha_inter=vumy_wet$alpha_esca[posts],
                                        resid_abund=residents_wet$esca, invader_abund=invader_abund)
  
  vumy_into_trhi_wet[r,] <- run.invader(surv=vumy_wet$surv, germ= vumy_wet$germ, 
                                        lambda=vumy_wet$lambda[posts], alpha_inter=vumy_wet$alpha_trhi[posts],
                                        resid_abund=residents_wet$trhi, invader_abund=invader_abund)
}

# Laca invades
post_length <- length(laca_wet$lambda)
for (r in 1:reps) {
  posts <- sample(post_length, runs, replace=TRUE)
  
  laca_into_avfa_wet[r,] <- run.invader(surv=laca_wet$surv, germ= laca_wet$germ, 
                                        lambda=laca_wet$lambda[posts], alpha_inter=laca_wet$alpha_avfa[posts],
                                        resid_abund=residents_wet$avfa, invader_abund=invader_abund)
  
  laca_into_brho_wet[r,] <- run.invader(surv=laca_wet$surv, germ= laca_wet$germ, 
                                        lambda=laca_wet$lambda[posts], alpha_inter=laca_wet$alpha_brho[posts],
                                        resid_abund=residents_wet$brho, invader_abund=invader_abund)
  
  laca_into_vumy_wet[r,] <- run.invader(surv=laca_wet$surv, germ= laca_wet$germ, 
                                        lambda=laca_wet$lambda[posts], alpha_inter=laca_wet$alpha_vumy[posts],
                                        resid_abund=residents_wet$vumy, invader_abund=invader_abund)
  
  laca_into_esca_wet[r,] <- run.invader(surv=laca_wet$surv, germ= laca_wet$germ, 
                                        lambda=laca_wet$lambda[posts], alpha_inter=laca_wet$alpha_esca[posts],
                                        resid_abund=residents_wet$esca, invader_abund=invader_abund)
  
  laca_into_trhi_wet[r,] <- run.invader(surv=laca_wet$surv, germ= laca_wet$germ, 
                                        lambda=laca_wet$lambda[posts], alpha_inter=laca_wet$alpha_trhi[posts],
                                        resid_abund=residents_wet$trhi, invader_abund=invader_abund)
}

# Esca invades
post_length <- length(esca_wet$lambda)
for (r in 1:reps) {
  posts <- sample(post_length, runs, replace=TRUE)
  
  esca_into_avfa_wet[r,] <- run.invader(surv=esca_wet$surv, germ= esca_wet$germ, 
                                        lambda=esca_wet$lambda[posts], alpha_inter=esca_wet$alpha_avfa[posts],
                                        resid_abund=residents_wet$avfa, invader_abund=invader_abund)
  
  esca_into_brho_wet[r,] <- run.invader(surv=esca_wet$surv, germ= esca_wet$germ, 
                                        lambda=esca_wet$lambda[posts], alpha_inter=esca_wet$alpha_brho[posts],
                                        resid_abund=residents_wet$brho, invader_abund=invader_abund)
  
  esca_into_vumy_wet[r,] <- run.invader(surv=esca_wet$surv, germ= esca_wet$germ, 
                                        lambda=esca_wet$lambda[posts], alpha_inter=esca_wet$alpha_vumy[posts],
                                        resid_abund=residents_wet$vumy, invader_abund=invader_abund)
  
  esca_into_laca_wet[r,] <- run.invader(surv=esca_wet$surv, germ= esca_wet$germ, 
                                        lambda=esca_wet$lambda[posts], alpha_inter=esca_wet$alpha_laca[posts],
                                        resid_abund=residents_wet$laca, invader_abund=invader_abund)
  
  esca_into_trhi_wet[r,] <- run.invader(surv=esca_wet$surv, germ= esca_wet$germ, 
                                        lambda=esca_wet$lambda[posts], alpha_inter=esca_wet$alpha_trhi[posts],
                                        resid_abund=residents_wet$trhi, invader_abund=invader_abund)
}

# trhi invades
post_length <- length(trhi_wet$lambda)
for (r in 1:reps) {
  posts <- sample(post_length, runs, replace=TRUE)
  
  trhi_into_avfa_wet[r,] <- run.invader(surv=trhi_wet$surv, germ= trhi_wet$germ, 
                                        lambda=trhi_wet$lambda[posts], alpha_inter=trhi_wet$alpha_avfa[posts],
                                        resid_abund=residents_wet$avfa, invader_abund=invader_abund)
  
  trhi_into_brho_wet[r,] <- run.invader(surv=trhi_wet$surv, germ= trhi_wet$germ, 
                                        lambda=trhi_wet$lambda[posts], alpha_inter=trhi_wet$alpha_brho[posts],
                                        resid_abund=residents_wet$brho, invader_abund=invader_abund)
  
  trhi_into_vumy_wet[r,] <- run.invader(surv=trhi_wet$surv, germ= trhi_wet$germ, 
                                        lambda=trhi_wet$lambda[posts], alpha_inter=trhi_wet$alpha_vumy[posts],
                                        resid_abund=residents_wet$vumy, invader_abund=invader_abund)
  
  trhi_into_laca_wet[r,] <- run.invader(surv=trhi_wet$surv, germ= trhi_wet$germ, 
                                        lambda=trhi_wet$lambda[posts], alpha_inter=trhi_wet$alpha_laca[posts],
                                        resid_abund=residents_wet$laca, invader_abund=invader_abund)
  
  trhi_into_esca_wet[r,] <- run.invader(surv=trhi_wet$surv, germ= trhi_wet$germ, 
                                        lambda=trhi_wet$lambda[posts], alpha_inter=trhi_wet$alpha_esca[posts],
                                        resid_abund=residents_wet$esca, invader_abund=invader_abund)
  
}

invasion_dry <- data.frame(as.vector(avfa_into_brho_dry), as.vector(avfa_into_vumy_dry), as.vector(avfa_into_laca_dry), as.vector(avfa_into_esca_dry), as.vector(avfa_into_trhi_dry),
                           as.vector(brho_into_avfa_dry), as.vector(brho_into_vumy_dry), as.vector(brho_into_laca_dry), as.vector(brho_into_esca_dry), as.vector(brho_into_trhi_dry),
                           as.vector(vumy_into_avfa_dry), as.vector(vumy_into_brho_dry), as.vector(vumy_into_laca_dry), as.vector(vumy_into_esca_dry), as.vector(vumy_into_trhi_dry),
                           as.vector(laca_into_avfa_dry), as.vector(laca_into_brho_dry), as.vector(laca_into_vumy_dry), as.vector(laca_into_esca_dry), as.vector(laca_into_trhi_dry),
                           as.vector(esca_into_avfa_dry), as.vector(esca_into_brho_dry), as.vector(esca_into_vumy_dry), as.vector(esca_into_laca_dry), as.vector(esca_into_trhi_dry),
                           as.vector(trhi_into_avfa_dry), as.vector(trhi_into_brho_dry), as.vector(trhi_into_vumy_dry), as.vector(trhi_into_laca_dry), as.vector(trhi_into_esca_dry))
names(invasion_dry) <- c("avfa_into_brho", "avfa_into_vumy", "avfa_into_laca", "avfa_into_esca", "avfa_into_trhi",
                         "brho_into_avfa", "brho_into_vumy", "brho_into_laca", "brho_into_esca", "brho_into_trhi",
                         "vumy_into_avfa", "vumy_into_brho", "vumy_into_laca", "vumy_into_esca", "vumy_into_trhi",
                         "laca_into_avfa", "laca_into_brho", "laca_into_vumy", "laca_into_esca", "laca_into_trhi",
                         "esca_into_avfa", "esca_into_brho", "esca_into_vumy", "esca_into_laca", "esca_into_trhi",
                         "trhi_into_avfa", "trhi_into_brho", "trhi_into_vumy", "trhi_into_laca", "trhi_into_esca")

invasion_wet <- data.frame(as.vector(avfa_into_brho_wet), as.vector(avfa_into_vumy_wet), as.vector(avfa_into_laca_wet), as.vector(avfa_into_esca_wet), as.vector(avfa_into_trhi_wet),
                           as.vector(brho_into_avfa_wet), as.vector(brho_into_vumy_wet), as.vector(brho_into_laca_wet), as.vector(brho_into_esca_wet), as.vector(brho_into_trhi_wet),
                           as.vector(vumy_into_avfa_wet), as.vector(vumy_into_brho_wet), as.vector(vumy_into_laca_wet), as.vector(vumy_into_esca_wet), as.vector(vumy_into_trhi_wet),
                           as.vector(laca_into_avfa_wet), as.vector(laca_into_brho_wet), as.vector(laca_into_vumy_wet), as.vector(laca_into_esca_wet), as.vector(laca_into_trhi_wet),
                           as.vector(esca_into_avfa_wet), as.vector(esca_into_brho_wet), as.vector(esca_into_vumy_wet), as.vector(esca_into_laca_wet), as.vector(esca_into_trhi_wet),
                           as.vector(trhi_into_avfa_wet), as.vector(trhi_into_brho_wet), as.vector(trhi_into_vumy_wet), as.vector(trhi_into_laca_wet), as.vector(trhi_into_esca_wet))
names(invasion_wet) <- c("avfa_into_brho", "avfa_into_vumy", "avfa_into_laca", "avfa_into_esca", "avfa_into_trhi",
                         "brho_into_avfa", "brho_into_vumy", "brho_into_laca", "brho_into_esca", "brho_into_trhi",
                         "vumy_into_avfa", "vumy_into_brho", "vumy_into_laca", "vumy_into_esca", "vumy_into_trhi",
                         "laca_into_avfa", "laca_into_brho", "laca_into_vumy", "laca_into_esca", "laca_into_trhi",
                         "esca_into_avfa", "esca_into_brho", "esca_into_vumy", "esca_into_laca", "esca_into_trhi",
                         "trhi_into_avfa", "trhi_into_brho", "trhi_into_vumy", "trhi_into_laca", "trhi_into_esca")

## Grab mean equilibrium abundances for both rainfall conditions
equil_abund <- as.data.frame(rbind(apply(residents_dry, 2, mean),
                   apply(residents_wet, 2, mean)))
equil_abund$trt <- c("dry","wet")

equil_abund <- equil_abund %>%
  pivot_longer(cols = avfa:trhi, names_to = "species", values_to = "abundance")
  

rm(list=setdiff(ls(), c("invasion_dry", "invasion_wet","equil_abund","params")))

### Plot ----

dry_means <- invasion_dry %>% 
  summarise_all(list(mean))

wet_means <- invasion_wet %>% 
  summarise_all(list(mean))

invasion_means <- rbind(dry_means, wet_means)
invasion_means$trt <- c("dry","wet")

invasion_means <- invasion_means %>%
 pivot_longer(cols = avfa_into_brho:trhi_into_esca, 
              names_to = "invasion", values_to = "growth")

invasion_means$pair <- NA
invasion_means$pair[grepl("avfa", invasion_means$invasion) & grepl("brho", invasion_means$invasion)] <- "avfa_brho"
invasion_means$pair[grepl("avfa", invasion_means$invasion) & grepl("vumy", invasion_means$invasion)] <- "avfa_vumy"
invasion_means$pair[grepl("avfa", invasion_means$invasion) & grepl("laca", invasion_means$invasion)] <- "avfa_laca"
invasion_means$pair[grepl("avfa", invasion_means$invasion) & grepl("esca", invasion_means$invasion)] <- "avfa_esca"
invasion_means$pair[grepl("avfa", invasion_means$invasion) & grepl("trhi", invasion_means$invasion)] <- "avfa_trhi"

invasion_means$pair[grepl("brho", invasion_means$invasion) & grepl("vumy", invasion_means$invasion)] <- "brho_vumy"
invasion_means$pair[grepl("brho", invasion_means$invasion) & grepl("laca", invasion_means$invasion)] <- "brho_laca"
invasion_means$pair[grepl("brho", invasion_means$invasion) & grepl("esca", invasion_means$invasion)] <- "brho_esca"
invasion_means$pair[grepl("brho", invasion_means$invasion) & grepl("trhi", invasion_means$invasion)] <- "brho_trhi"

invasion_means$pair[grepl("vumy", invasion_means$invasion) & grepl("laca", invasion_means$invasion)] <- "vumy_laca"
invasion_means$pair[grepl("vumy", invasion_means$invasion) & grepl("esca", invasion_means$invasion)] <- "vumy_esca"
invasion_means$pair[grepl("vumy", invasion_means$invasion) & grepl("trhi", invasion_means$invasion)] <- "vumy_trhi"

invasion_means$pair[grepl("laca", invasion_means$invasion) & grepl("esca", invasion_means$invasion)] <- "laca_esca"
invasion_means$pair[grepl("laca", invasion_means$invasion) & grepl("trhi", invasion_means$invasion)] <- "laca_trhi"

invasion_means$pair[grepl("esca", invasion_means$invasion) & grepl("trhi", invasion_means$invasion)] <- "esca_trhi"

invasion_means$plot <- c(0.5, 0.5, 0.5, 0.5, 0.5, 
                         1.3, 0.5, 0.5, 0.5, 0.5,
                         1.3, 1.3, 0.5, 0.5, 0.5,
                         1.3, 1.3, 1.3, 0.5, 0.5,
                         1.3, 1.3, 1.3, 1.3, 0.5,
                         1.3, 1.3, 1.3, 1.3, 1.3,
                         0.7, 0.7, 0.7, 0.7, 0.7,
                         1.5, 0.7, 0.7, 0.7, 0.7,
                         1.5, 1.5, 0.7, 0.7, 0.7,
                         1.5, 1.5, 1.5, 0.7, 0.7,
                         1.5, 1.5, 1.5, 1.5, 0.7,
                         1.5, 1.5, 1.5, 1.5, 1.5)

invasion_means$pair <- factor(invasion_means$pair, 
                              levels = c("avfa_brho", "avfa_esca", "avfa_laca", "avfa_vumy", "avfa_trhi",
                                         "brho_esca", "brho_laca", "brho_vumy", "brho_trhi",
                                         "vumy_esca", "vumy_laca","vumy_trhi",
                                         "laca_esca", "laca_trhi",
                                         "esca_trhi"))

pdf("./Competition/Figures/Invader_LDGR_plots.pdf", height = 6, width = 8)
ggplot(invasion_means, aes(x = plot, y = growth, color = factor(trt))) + 
  geom_point(size = 3) + 
  geom_hline(yintercept = 0, linetype = "dashed") + ylim(-6.5,6) + ylab("Invader Log LDGR") + xlab("Invader-Resident") + 
  facet_wrap(~pair) + 
  scale_x_continuous(breaks = c(0.6, 1.4), labels = c("1-2", "2-1"), limits = c(.2,1.8)) +
  labs(color = "Fall Treatment") + 
  scale_color_manual(values = c("#f4a261", "#2a9d8f")) + 
  theme_bw()
dev.off()
