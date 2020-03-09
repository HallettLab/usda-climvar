# run coexistence models for dry versus wet

source("./Competition/Model-fit/import_posteriors.R")

run.to.equilibrium <- function(surv, germ, lambda, alpha_intra, Nt) {
  Ntp1 <- (1-germ)*surv*Nt + germ*lambda*Nt/(1+ alpha_intra * Nt)
  return(Ntp1)

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


