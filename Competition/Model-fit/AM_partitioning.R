### Setup
#setwd('./UO Hallett/Projects/usda-climvar')
library(tidyverse)
library(rstan)

### Create historical rainfall conditions ----
# First determine how common each environmental type is
# what about for what we actually see in terms of the number of years in each env. condition
# first read in the data
rain <- read_csv("../avena-erodium/Data/PRISM_brownsvalley_long.csv", skip = 10) %>%
  mutate(ppt = `ppt (inches)`*2.54*10) %>%
  separate(Date, c("year", "month")) %>%
  mutate(year = as.numeric(year),
         month = as.numeric(month)) %>%
  mutate(year = ifelse(month == 12 | month == 11 | month == 10 | month == 9, year + 1, year)) %>%
  mutate(season = "Early",
         season = ifelse(month == 2 | month == 3 | month == 4, "Late", season)) %>%
  filter(month != 5, month != 6, month!= 7, month != 8)

## Summarize by year 
## Using 50% as the cutoff 
rainsummary <-  rain %>%
  group_by(year, season) %>%
  summarize(ppt = sum(ppt)) %>%
  spread(season, ppt) %>%
  mutate(Total = Early + Late) 

rainsummary <- rainsummary %>%
  mutate(raintype = "fallWet",
         raintype = ifelse(Early < quantile(rainsummary$Early, .5), "fallDry", raintype))

fall.dry <- length(which(rainsummary$raintype == "fallDry")) / nrow(rainsummary)
fall.wet <- length(which(rainsummary$raintype == "fallWet")) / nrow(rainsummary)

### Load or set model parameters ----
source('./Competition/Model-fit/import_posteriors_AM.R')

# Set germination and survival rates
params$germ[params$species == "avfa"] <- 0.78
params$germ[params$species == "brho"] <- 0.8
params$germ[params$species == "vumy"] <- 0.8
params$germ[params$species == "laca"] <- 0.8
params$germ[params$species == "esca"] <- 0.92

params$surv[params$species == "avfa"] <- 0.01
params$surv[params$species == "brho"] <- 0.013
params$surv[params$species == "vumy"] <- 0.045
params$surv[params$species == "laca"] <- 0.013
params$surv[params$species == "esca"] <- 0.013

# Create weighted parameters for mechanism partitioning
params_weighted <- params[params$treatment == "fallDry",1:6]*fall.dry + params[params$treatment == "fallWet",1:6]*fall.wet
params_weighted$species <- c("avfa","vumy","brho","esca","laca")
params_weighted <- left_join(params_weighted, unique(params[,c("species","germ","surv")]))

### Coexistence functions ----
# Determine equilibrium conditions for each species in isolation 
pop.equilibrium <- function (N0, s, g, a_intra, lambda) {
  # to run for only a single timestep
  N <- s*(1-g)*N0 + N0*(lambda*g)/(1+a_intra*N0)
  return(N)
}

# invader population growth rate one time step forward
pop.invade <- function (N0, resident, s, g, a_inter, lambda) {
  # to run for only a single timestep
  N <- s*(1-g)*N0 + N0*(lambda*g)/(1+a_inter*resident)
  return(N)
}

# resident population growth rate one time step forward
pop.resident <- function (N0, resident, s, g, a_intra, a_inter, lambda) {
  # to run for only a single timestep
  N <- s*(1-g)*resident + resident*(lambda*g)/(1+a_intra*resident+a_inter*N0)
  return(N)
}
### Coexistence - No partitioning ----

## Run all species to equilibrium in isolation
# Species to loop across and time in years
species <- c("avfa","brho","vumy","laca","esca")
time <- length(rainsummary$raintype)

# Pre-allocate empty object to store results, and set initial conditions
N <- matrix(NA, nrow = time+1, ncol = length(species))
colnames(N) <- species
N[1,] <- 300

## Run to equilibrium given historical rainfall conditions
for(i in 1:length(species)){

  resident <- params[params$species == species[i],]
  resident_a <- paste0("alpha_", species[i])
  
  for (t in 1:time) {
    
    sim <- resident[resident$treatment == rainsummary$raintype[t],]
    
    N[t+1,i] <- pop.equilibrium(N0 = N[t,i], s = sim$surv, g = sim$germ, lambda = sim$lambda,
                                a_intra = as.numeric(sim[resident_a]))
  
    }
}

## Run pairwise invasions into equilibrium resident 

# Pre-allocate empty object to store results, a separate matrix for each invader
coexist_out <- list(avfa = matrix(NA, nrow = 72, ncol = length(species)-1),
                    brho = matrix(NA, nrow = 72, ncol = length(species)-1),
                    vumy = matrix(NA, nrow = 72, ncol = length(species)-1),
                    laca = matrix(NA, nrow = 72, ncol = length(species)-1),
                    esca = matrix(NA, nrow = 72, ncol = length(species)-1))
colnames(coexist_out[['avfa']]) <- species[species != 'avfa']
colnames(coexist_out[['brho']]) <- species[species != 'brho']
colnames(coexist_out[['vumy']]) <- species[species != 'vumy']
colnames(coexist_out[['laca']]) <- species[species != 'laca']
colnames(coexist_out[['esca']]) <- species[species != 'esca']

# Loop once for each invader
for (i in 1:length(species)){
  
  invader <- params[params$species == species[i],]
  invader_out <- coexist_out[[species[i]]]
  
  inter <- species[species != species[i]]
  
  # Loop once for each potential resident that isn't the invader
  for (j in 1:length(inter)){
  
    resident <- inter[j]
    resident_a <- paste0("alpha_",resident)
    
    temp <- 1
    
    # Calculate GRWR for each year after burn-in
    for (t in 50:time) {
      
      sim <- invader[invader$treatment == rainsummary$raintype[t],] 
      
      invader_out[temp,j] <- pop.invade(N0 = 1, s = sim$surv, g = sim$germ, lambda = sim$lambda,
                                        resident = as.numeric(N[t,resident]),
                                        a_inter = as.numeric(sim[resident_a]))
      
      temp  <- temp + 1 
    }
  }
  
  # Save in output object
  coexist_out[[species[i]]] <- invader_out
  
}

# Calculate LGDR
coexist_out <- lapply(coexist_out, colMeans)
coexist_out <- lapply(coexist_out, log)

### Coexistence - Variation-independent mechanisms ----

## Run all species to equilibrium with average rainfall conditions
# Species to loop across and time in years
species <- c("avfa","brho","vumy","laca","esca")
time <- length(rainsummary$raintype)

# Pre-allocate empty object to store results, and set initial conditions
N_no_var <- matrix(NA, nrow = time+1, ncol = length(species))
colnames(N_no_var) <- species
N_no_var[1,] <- 300

## Run to equilibrium given historical rainfall conditions
for(i in 1:length(species)){
  
  resident <- params_weighted[params_weighted$species == species[i],]
  resident_a <- paste0("alpha_", species[i])
  
  for (t in 1:time) {
    
    N_no_var[t+1,i] <- pop.equilibrium(N0 = N_no_var[t,i], s = resident$surv, g = resident$germ, lambda = resident$lambda,
                                a_intra = as.numeric(resident[resident_a]))
  }
}

## Run pairwise invasions into equilibrium resident 

# Pre-allocate empty object to store results, a separate matrix for each invader
coexist_no_var <- list(avfa = matrix(NA, nrow = 72, ncol = length(species)-1),
                    brho = matrix(NA, nrow = 72, ncol = length(species)-1),
                    vumy = matrix(NA, nrow = 72, ncol = length(species)-1),
                    laca = matrix(NA, nrow = 72, ncol = length(species)-1),
                    esca = matrix(NA, nrow = 72, ncol = length(species)-1))
colnames(coexist_no_var[['avfa']]) <- species[species != 'avfa']
colnames(coexist_no_var[['brho']]) <- species[species != 'brho']
colnames(coexist_no_var[['vumy']]) <- species[species != 'vumy']
colnames(coexist_no_var[['laca']]) <- species[species != 'laca']
colnames(coexist_no_var[['esca']]) <- species[species != 'esca']

# Loop once for each invader
for (i in 1:length(species)){
  
  invader <- params_weighted[params_weighted$species == species[i],]
  invader_out <- coexist_no_var[[species[i]]]
  
  inter <- species[species != species[i]]
  
  # Loop once for each potential resident that isn't the invader
  for (j in 1:length(inter)){
    
    resident <- inter[j]
    resident_a <- paste0("alpha_",resident)
    
    temp <- 1
    
    # Calculate GRWR for each year after burn-in
    for (t in 50:time) {
      
      invader_out[temp,j] <- pop.invade(N0 = 1, s = invader$surv, g = invader$germ, lambda = invader$lambda,
                                        resident = as.numeric(N_no_var[t,resident]),
                                        a_inter = as.numeric(invader[resident_a]))
      
      temp  <- temp + 1 
    }
  }
  
  # Save in output object
  coexist_no_var[[species[i]]] <- invader_out
  
}

# Calculate LGDR
coexist_no_var <- lapply(coexist_no_var, colMeans)
coexist_eps_0 <- lapply(coexist_no_var, log)

### Coexistence - Variable lambda ----

## Run all species to equilibrium with variable lambda
# Species to loop across and time in years
species <- c("avfa","brho","vumy","laca","esca")
time <- length(rainsummary$raintype)

# Pre-allocate empty object to store results, and set initial conditions
N_var_lamb <- matrix(NA, nrow = time+1, ncol = length(species))
colnames(N_var_lamb) <- species
N_var_lamb[1,] <- 300

## Run to equilibrium given historical rainfall conditions
for(i in 1:length(species)){
  
  resident <- params_weighted[params_weighted$species == species[i],]
  resident_var <- params[params$species == species[i],]
  resident_a <- paste0("alpha_", species[i])
  
  for (t in 1:time) {
    
    sim <- resident_var[resident_var$treatment == rainsummary$raintype[t],]
    
    N_var_lamb[t+1,i] <- pop.equilibrium(N0 = N_var_lamb[t,i], s = resident$surv, g = resident$germ, lambda = sim$lambda,
                                       a_intra = as.numeric(resident[resident_a]))
  }
}

## Run pairwise invasions into equilibrium resident 

# Pre-allocate empty object to store results, a separate matrix for each invader
coexist_var_lamb <- list(avfa = matrix(NA, nrow = 72, ncol = length(species)-1),
                       brho = matrix(NA, nrow = 72, ncol = length(species)-1),
                       vumy = matrix(NA, nrow = 72, ncol = length(species)-1),
                       laca = matrix(NA, nrow = 72, ncol = length(species)-1),
                       esca = matrix(NA, nrow = 72, ncol = length(species)-1))
colnames(coexist_var_lamb[['avfa']]) <- species[species != 'avfa']
colnames(coexist_var_lamb[['brho']]) <- species[species != 'brho']
colnames(coexist_var_lamb[['vumy']]) <- species[species != 'vumy']
colnames(coexist_var_lamb[['laca']]) <- species[species != 'laca']
colnames(coexist_var_lamb[['esca']]) <- species[species != 'esca']

# Loop once for each invader
for (i in 1:length(species)){
  
  invader <- params_weighted[params_weighted$species == species[i],]
  invader_var <- params[params$species == species[i],]
  
  invader_out <- coexist_var_lamb[[species[i]]]
  
  inter <- species[species != species[i]]
  
  # Loop once for each potential resident that isn't the invader
  for (j in 1:length(inter)){
    
    resident <- inter[j]
    resident_a <- paste0("alpha_",resident)
    
    temp <- 1
    
    # Calculate GRWR for each year after burn-in
    for (t in 50:time) {
      
      sim <- invader_var[invader_var$treatment == rainsummary$raintype[t],]
      
      invader_out[temp,j] <- pop.invade(N0 = 1, s = invader$surv, g = invader$germ, lambda = sim$lambda,
                                        resident = as.numeric(N_var_lamb[t,resident]),
                                        a_inter = as.numeric(invader[resident_a]))
      
      temp  <- temp + 1 
    }
  }
  
  # Save in output object
  coexist_var_lamb[[species[i]]] <- invader_out
  
}

# Calculate LGDR
coexist_var_lamb <- lapply(coexist_var_lamb, colMeans)
coexist_var_lamb <- lapply(coexist_var_lamb, log)
coexist_eps_lamb <- Map('-', coexist_var_lamb, coexist_eps_0)

### Coexistence - variable alpha ----

## Run all species to equilibrium with variable lambda
# Species to loop across and time in years
species <- c("avfa","brho","vumy","laca","esca")
time <- length(rainsummary$raintype)

# Pre-allocate empty object to store results, and set initial conditions
N_var_alph <- matrix(NA, nrow = time+1, ncol = length(species))
colnames(N_var_alph) <- species
N_var_alph[1,] <- 300

## Run to equilibrium given historical rainfall conditions
for(i in 1:length(species)){
  
  resident <- params_weighted[params_weighted$species == species[i],]
  resident_var <- params[params$species == species[i],]
  resident_a <- paste0("alpha_", species[i])
  
  for (t in 1:time) {
    
    sim <- resident_var[resident_var$treatment == rainsummary$raintype[t],]
    
    N_var_alph[t+1,i] <- pop.equilibrium(N0 = N_var_alph[t,i], s = resident$surv, g = resident$germ, lambda = resident$lambda,
                                         a_intra = as.numeric(sim[resident_a]))
  }
}

## Run pairwise invasions into equilibrium resident 

# Pre-allocate empty object to store results, a separate matrix for each invader
coexist_var_alph <- list(avfa = matrix(NA, nrow = 72, ncol = length(species)-1),
                         brho = matrix(NA, nrow = 72, ncol = length(species)-1),
                         vumy = matrix(NA, nrow = 72, ncol = length(species)-1),
                         laca = matrix(NA, nrow = 72, ncol = length(species)-1),
                         esca = matrix(NA, nrow = 72, ncol = length(species)-1))
colnames(coexist_var_alph[['avfa']]) <- species[species != 'avfa']
colnames(coexist_var_alph[['brho']]) <- species[species != 'brho']
colnames(coexist_var_alph[['vumy']]) <- species[species != 'vumy']
colnames(coexist_var_alph[['laca']]) <- species[species != 'laca']
colnames(coexist_var_alph[['esca']]) <- species[species != 'esca']

# Loop once for each invader
for (i in 1:length(species)){
  
  invader <- params_weighted[params_weighted$species == species[i],]
  invader_var <- params[params$species == species[i],]
  
  invader_out <- coexist_var_alph[[species[i]]]
  
  inter <- species[species != species[i]]
  
  # Loop once for each potential resident that isn't the invader
  for (j in 1:length(inter)){
    
    resident <- inter[j]
    resident_a <- paste0("alpha_",resident)
    
    temp <- 1
    
    # Calculate GRWR for each year after burn-in
    for (t in 50:time) {
      
      sim <- invader_var[invader_var$treatment == rainsummary$raintype[t],]
      
      invader_out[temp,j] <- pop.invade(N0 = 1, s = invader$surv, g = invader$germ, lambda = invader$lambda,
                                        resident = as.numeric(N_var_alph[t,resident]),
                                        a_inter = as.numeric(sim[resident_a]))
      
      temp  <- temp + 1 
    }
  }
  
  # Save in output object
  coexist_var_alph[[species[i]]] <- invader_out
  
}

# Calculate LGDR
coexist_var_alph <- lapply(coexist_var_alph, colMeans)
coexist_var_alph <- lapply(coexist_var_alph, log)
coexist_eps_alpha <- Map('-', coexist_var_alph, coexist_eps_0)

### Coexistence - interaction term ----
## Add all partitioned terms together
coexist_eps_int <- Map('+', coexist_eps_0, coexist_eps_alpha)
coexist_eps_int <- Map('+', coexist_eps_int, coexist_eps_lamb)

## Subtract partitioned terms from overall invader growth rate
coexist_eps_int <- Map('-', coexist_out, coexist_eps_int)
