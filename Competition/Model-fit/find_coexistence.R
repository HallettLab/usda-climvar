source('./Competition/Model-fit/import_posteriors_AM.R')
source('./Competition/Model-fit/run_pairwise_coexistence.R')

### Invasion setup ----

## Define invasion function
pop.invade <- function (N0, resident, s, g, a_inter, lambda) {
  # to run for only a single timestep
  N <- s*(1-g)*N0 + N0*(lambda*g)/(1+a_inter*resident)
  return(N)
}

## Set seedbank parameters
params$germ[params$species == "avfa"] <- 0.78
params$germ[params$species == "brho"] <- 0.8
params$germ[params$species == "vumy"] <- 0.6
params$germ[params$species == "laca"] <- 0.8
params$germ[params$species == "esca"] <- 0.92
params$germ[params$species == "trhi"] <- 0.20

params$surv[params$species == "avfa"] <- 0.01
params$surv[params$species == "brho"] <- 0.013
params$surv[params$species == "vumy"] <- 0.045
params$surv[params$species == "laca"] <- 0.013
params$surv[params$species == "esca"] <- 0.013
params$surv[params$species == "trhi"] <- 0.01

species <- c("avfa","brho","vumy","laca","esca","trhi")

## Create output object
coexist_out <- list(avfa = matrix(NA, nrow = 2, ncol = length(species)),
                    brho = matrix(NA, nrow = 2, ncol = length(species)),
                    vumy = matrix(NA, nrow = 2, ncol = length(species)),
                    laca = matrix(NA, nrow = 2, ncol = length(species)),
                    esca = matrix(NA, nrow = 2, ncol = length(species)),
                    trhi = matrix(NA, nrow = 2, ncol = length(species)))
colnames(coexist_out[['avfa']]) <- species
colnames(coexist_out[['brho']]) <- species
colnames(coexist_out[['vumy']]) <- species
colnames(coexist_out[['laca']]) <- species
colnames(coexist_out[['esca']]) <- species
colnames(coexist_out[['trhi']]) <- species

### Run invasions for every species pair until coexistence is found ----

## For every species as an invader
for(i in 1:length(species)){
  
  # Grab relevant output object, set invader parameters in wet and dry treatments
  invader_out <- coexist_out[[species[i]]]
  invader <- species[i]
  inv_param_dry <- params[params$species == species[i] & params$treatment == "fallDry",]
  inv_param_wet <- params[params$species == species[i] & params$treatment == "fallWet",]
  
  # Create vector of potential resident species 
  inter <- species[species != species[i]]
  
  ## For every potential resident
  for(k in 1:length(inter)) {
    
    # Grab equilibrium abundances in each treatment
    resident <- inter[k]
    res_abund_dry <- equil_abund$abundance[equil_abund$species == inter[k] & equil_abund$trt == "dry"]
    res_abund_wet <- equil_abund$abundance[equil_abund$species == inter[k] & equil_abund$trt == "wet"]
    resident_a <- paste0("alpha_",resident)
    
    # Initialize growth variables to start while loop
    growth_dry <- 0
    growth_wet <- 0
    
    ## While dry treatment population growth rate isn't increasing 
    while(growth_dry <= 1) {
      
      # Calculate invader GRWR at current resident population, beginning with equilibrium population
      growth_dry <- pop.invade(N0 = 1, s = inv_param_dry$surv, g = inv_param_dry$germ, lambda = inv_param_dry$lambda,
                               resident = res_abund_dry,
                               a_inter = as.numeric(inv_param_dry[resident_a]))
      
      # Iterate the resident abundance downward by a small step
      res_abund_dry <- res_abund_dry - 0.1
      
    }
    
    ## While wet treatment population growth rate isn't increasing 
    while(growth_wet <= 1) {
      
      # Calculate invader GRWR at current resident population, beginning with equilibrium population
      growth_wet <- pop.invade(N0 = 1, s = inv_param_wet$surv, g = inv_param_wet$germ, lambda = inv_param_wet$lambda,
                               resident = res_abund_wet,
                               a_inter = as.numeric(inv_param_wet[resident_a]))
      
      # Iterate the resident abundance downward by a small step
      res_abund_wet <- res_abund_wet - 0.1
      
    }
    
    # Save abundances that allowed positive invader GRWR
    invader_out[1,resident] <- res_abund_dry
    invader_out[2,resident] <- res_abund_wet
  
  }
  
  # Save output obejct
  coexist_out[[species[i]]] <- invader_out
  
}

## Widen equilibrium abundance data for easier comparison with successful invasion abundances
abund_compare <- equil_abund %>%
  pivot_wider(names_from = species, values_from = abundance) %>%
  select(-trt)

coexist_comp <- coexist_out

## For each species
for(i in 1:length(species)) {
  
  # Divide successful invasion abundance by equilibrium abundance to find as fraction of equilibrium
  compare <- as.data.frame(coexist_comp[[species[i]]]/abund_compare)
  compare$trt <- c("dry","wet")
  
  # Make longer for plotting
  compare <- compare %>%
    pivot_longer(cols = !trt, names_to = "resident", values_to = "prop_abundance")
  
  compare$invader <- species[i]
  
  coexist_comp[[species[i]]] <- compare
  
}

## Reformat list output as a data frame
coexist_comp <- bind_rows(coexist_comp)

# Reorder x axis residents by competitive ability
coexist_comp <- coexist_comp %>%
  mutate(resident = factor(resident, levels = c("avfa","vumy","brho","laca","trhi","esca")))

## Plot

pdf('./Competition/Figures/figures/all_coexist_abundance.pdf', height = 7, width = 6)

ggplot(coexist_comp, aes(x = rep(1:6, 12), y = prop_abundance, fill = as.factor(trt))) + 
  geom_bar(position = "dodge", stat = "identity") + 
  facet_grid(rows = vars(invader)) + 
  scale_x_continuous(breaks = 1:6, labels = species) + 
  scale_fill_manual(values = c("#f4a261","#2a9d8f"), name = "Fall Rainfall") + 
  xlab("Resident") + ylab("Relative Resident Equilibrium Carrying Capacity") + 
  labs(color = "Treatment") +
  theme_bw()

dev.off()