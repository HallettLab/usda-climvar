# run coexistence models for dry versus wet

source("./Competition/Model-fit/import_posteriors.R")

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

vumy_dry$germ <- .8
vumy_wet$germ <- .8
vumy_dry$surv <- .045
vumy_wet$surv <- .045

all_datset <- list(avfa_dry, avfa_wet, brho_dry, brho_wet, vumy_dry, vumy_wet) # set which dataset we want to do
all_intra <- c("alpha_avfa", "alpha_avfa", "alpha_brho", "alpha_brho", "alpha_vumy", "alpha_vumy")
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

# make into dry versus wet conditions dataframes
residents_dry <-  data.frame(N[1,,200], N[3,,200], N[5,,200])
names(residents_dry) <- c("AVFA", "BRHO", "VUMY")

residents_wet <-  data.frame(N[2,,200], N[4,,200], N[6,,200])
names(residents_wet) <- c("AVFA", "BRHO", "VUMY")

# invade into residents
# -------------------------------------------------------------------------------------------------------
# invader into dry  
reps <- 200
avfa_into_brho_dry <- matrix(NA, reps, runs)
avfa_into_vumy_dry <- matrix(NA, reps, runs)
brho_into_avfa_dry <- matrix(NA, reps, runs)
brho_into_vumy_dry <- matrix(NA, reps, runs)
vumy_into_avfa_dry <- matrix(NA, reps, runs)
vumy_into_brho_dry <- matrix(NA, reps, runs)

invader_abund <- 1

# avena invades
post_length <- length(avfa_dry$lambda)
for (r in 1:reps) {

  posts <- sample(post_length, runs, replace=TRUE)
  
  avfa_into_brho_dry[r,] <- run.invader(surv=avfa_dry$surv, germ= avfa_dry$germ, 
                                    lambda=avfa_dry$lambda[posts], alpha_inter=avfa_dry$alpha_brho[posts],
                                    resid_abund=residents_dry$BRHO, invader_abund=invader_abund)
  
  avfa_into_vumy_dry[r,] <- run.invader(surv=avfa_dry$surv, germ= avfa_dry$germ, 
                                        lambda=avfa_dry$lambda[posts], alpha_inter=avfa_dry$alpha_vumy[posts],
                                        resid_abund=residents_dry$VUMY, invader_abund=invader_abund)
}

# brho invades
post_length <- length(brho_dry$lambda)
for (r in 1:reps) {
  posts <- sample(post_length, runs, replace=TRUE)
  
  brho_into_avfa_dry[r,] <- run.invader(surv=brho_dry$surv, germ= brho_dry$germ, 
                                        lambda=brho_dry$lambda[posts], alpha_inter=brho_dry$alpha_avfa[posts],
                                        resid_abund=residents_dry$AVFA, invader_abund=invader_abund)
  
  brho_into_vumy_dry[r,] <- run.invader(surv=brho_dry$surv, germ= brho_dry$germ, 
                                        lambda=brho_dry$lambda[posts], alpha_inter=brho_dry$alpha_vumy[posts],
                                        resid_abund=residents_dry$VUMY, invader_abund=invader_abund)
  
}

# Vumpy invades
post_length <- length(vumy_dry$lambda)
for (r in 1:reps) {
  posts <- sample(post_length, runs, replace=TRUE)
  
  vumy_into_avfa_dry[r,] <- run.invader(surv=vumy_dry$surv, germ= vumy_dry$germ, 
                                        lambda=vumy_dry$lambda[posts], alpha_inter=vumy_dry$alpha_avfa[posts],
                                        resid_abund=residents_dry$AVFA, invader_abund=invader_abund)
  
  vumy_into_brho_dry[r,] <- run.invader(surv=vumy_dry$surv, germ= vumy_dry$germ, 
                                        lambda=vumy_dry$lambda[posts], alpha_inter=vumy_dry$alpha_brho[posts],
                                        resid_abund=residents_dry$BRHO, invader_abund=invader_abund)
  
}

# -------------------------------------------------------------------------------------------------------
# invader into wet  
reps <- 20
avfa_into_brho_wet <- matrix(NA, reps, runs)
avfa_into_vumy_wet <- matrix(NA, reps, runs)
brho_into_avfa_wet <- matrix(NA, reps, runs)
brho_into_vumy_wet <- matrix(NA, reps, runs)
vumy_into_avfa_wet <- matrix(NA, reps, runs)
vumy_into_brho_wet <- matrix(NA, reps, runs)

invader_abund <- 1

# avena invades
post_length <- length(avfa_wet$lambda)
for (r in 1:reps) {
  
  posts <- sample(post_length, runs, replace=TRUE)
  
  avfa_into_brho_wet[r,] <- run.invader(surv=avfa_wet$surv, germ= avfa_wet$germ, 
                                        lambda=avfa_wet$lambda[posts], alpha_inter=avfa_wet$alpha_brho[posts],
                                        resid_abund=residents_wet$BRHO, invader_abund=invader_abund)
  
  avfa_into_vumy_wet[r,] <- run.invader(surv=avfa_wet$surv, germ= avfa_wet$germ, 
                                        lambda=avfa_wet$lambda[posts], alpha_inter=avfa_wet$alpha_vumy[posts],
                                        resid_abund=residents_wet$VUMY, invader_abund=invader_abund)
}

# brho invades
post_length <- length(brho_wet$lambda)
for (r in 1:reps) {
  posts <- sample(post_length, runs, replace=TRUE)
  
  brho_into_avfa_wet[r,] <- run.invader(surv=brho_wet$surv, germ= brho_wet$germ, 
                                        lambda=brho_wet$lambda[posts], alpha_inter=brho_wet$alpha_avfa[posts],
                                        resid_abund=residents_wet$AVFA, invader_abund=invader_abund)
  
  brho_into_vumy_wet[r,] <- run.invader(surv=brho_wet$surv, germ= brho_wet$germ, 
                                        lambda=brho_wet$lambda[posts], alpha_inter=brho_wet$alpha_vumy[posts],
                                        resid_abund=residents_wet$VUMY, invader_abund=invader_abund)
  
}

# Vumpy invades
post_length <- length(vumy_wet$lambda)
for (r in 1:reps) {
  posts <- sample(post_length, runs, replace=TRUE)
  
  vumy_into_avfa_wet[r,] <- run.invader(surv=vumy_wet$surv, germ= vumy_wet$germ, 
                                        lambda=vumy_wet$lambda[posts], alpha_inter=vumy_wet$alpha_avfa[posts],
                                        resid_abund=residents_wet$AVFA, invader_abund=invader_abund)
  
  vumy_into_brho_wet[r,] <- run.invader(surv=vumy_wet$surv, germ= vumy_wet$germ, 
                                        lambda=vumy_wet$lambda[posts], alpha_inter=vumy_wet$alpha_brho[posts],
                                        resid_abund=residents_wet$BRHO, invader_abund=invader_abund)
  
}

invasion_dry <- data.frame(as.vector(avfa_into_brho_dry), as.vector(avfa_into_vumy_dry),
                           as.vector(brho_into_avfa_dry), as.vector(brho_into_vumy_dry),
                           as.vector(vumy_into_avfa_dry), as.vector(vumy_into_brho_dry))
names(invasion_dry) <- c("avfa_into_brho_dry", "avfa_into_vumy_dry", "brho_into_avfa_dry",
                         "brho_into_vumy_dry", "vumy_into_avfa_dry", "vumy_into_brho_dry")

invasion_wet <- data.frame(as.vector(avfa_into_brho_wet), as.vector(avfa_into_vumy_wet),
                           as.vector(brho_into_avfa_wet), as.vector(brho_into_vumy_wet),
                           as.vector(vumy_into_avfa_wet), as.vector(vumy_into_brho_wet))
names(invasion_wet) <- c("avfa_into_brho", "avfa_into_vumy", "brho_into_avfa",
                         "brho_into_vumy", "vumy_into_avfa", "vumy_into_brho")

rm(list=setdiff(ls(), c("invasion_dry", "invasion_wet")))

