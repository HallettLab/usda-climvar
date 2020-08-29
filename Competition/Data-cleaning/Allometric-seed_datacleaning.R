# calculate biomass-fecundity allometric relationship by competiton species
# authors: LMH, CTW
# initate: Jan 2019

# script purpose:
# 1) derive allometric relationship of biomass to seed production (fecundity) by species by fall drought treatment
# 2) compile predicted values, confidence intervals and prediction intervals (using `predict`)
# write out (1) and (2) datasets to ClimVar Dropbox Competition_CleanedData subfolder 

# notes: 
# AVBA, BRHO, LACA, and VUMY specimens collected by fall dry and fall wet. check to see if treatment makes a difference on fecundity.
# more than 20 specimens collected for some species (e.g. april and late-season samples for forbs)
# ESCA specimens collected from Mariposa, CA, spring 2020, in ambient conditions (from CTW family property, by CTW's mom, sent to Boulder, CO, where CTW processing)
# TRHI specimens collected from USDA Compost project, spring 2020 (collectred by Nikolai, AS et al. processing in Eugene, OR)
# no ESCA or TRHI specimens with seeds available for collection like the other species at ClimVar project Spring 2017 (sampling trip #2 = too late for TRHI, too early for ESCA)

# dependencies:
## scripts dependent on file generated in this script (i.e. changes made to csv written out here may affect code in dependent scripts):
# 1. in Competition/Data-cleaning
# > Combine_datacleaning.R; model_prep.R
# >> other scripts depdendent on Combine_datacleaning.R: germain_figure.R
# 2. in Competition/Data-analysis
# > coexistence_model_formatting.R * <-- many scripts in Competition/Model-fit dependent on this script




# -- SETUP -----
rm(list = ls()) # clean environment
library(readxl)
library(tidyverse)
library(broom)
options(stringsAsFactors = F)
theme_set(theme_bw())
na_vals = c("", " ", NA, "NA")
datpath <- "~/Dropbox/ClimVar/Competition/Data/"

# read in cleaned phytometer dataset
phytodat <- read.csv(paste0(datpath, "Competition_CleanedData/Competition_phytometers_clean.csv"),
                     na.strings = na_vals)

# iterate through each tab in competition specimens 2017 workbook, read in and compile data.
# all tabs have same headers, each species x dry/wet on its own tab
# get spreadsheets
specsheets <- excel_sheets(paste0(datpath, "Competition_EnteredData/Competition_specimens_spring2017.xlsx"))
specsheets # first tab is metadata, tabs 2 through end are data tabs

# get colname order from metadata spreadsheet for ordering master compiled data frame (below)
specmeta <- read_excel(paste0(datpath, "Competition_EnteredData/Competition_specimens_spring2017.xlsx"), 
                       sheet = specsheets[1],
                       na = na_vals)
# clean up read in metadata table
# preserve variable descriptions and order
specmeta <- specmeta[grep("^Variable", specmeta[[1]]):nrow(specmeta),]
# set row one as colnames then remove
colnames(specmeta) <- specmeta[1,]; specmeta <- specmeta[-1,]

# update sep 2020: read in and append ESCA + TRHI dat after CTW and LMH process 2020 samples
esca <- read.csv(paste0(datpath, "Competition_EnteredData/Competition_ESCAspecimens_spring2020.csv"), na.strings = na_vals) %>%
  # convert esca dates to date
  mutate_at(c("clip_date", "wgh_date"), function(x) as.Date(x, format = "%m/%d/%y"))



# -- READ IN AND COMPILE ALL SPECIMEN TABS -----
# initiate df for storing specimen data
specdat <- data.frame()
for(i in specsheets[-1]){ # don't read in metadata tab
  print(paste("Reading in",i))
  tempdat <- read_excel(paste0(datpath,"Competition_EnteredData/Competition_specimens_spring2017.xlsx"), 
                        sheet = i, na = na_vals, trim_ws = T)
  #print(paste("Colnames for", i, "are:"))
  #print(colnames(tempdat))
  if(ncol(specdat) == 0){
    specdat <- rbind(specdat, tempdat)
  }else{
    #compare colnames in master vs read in dataset
    master_names <- colnames(specdat); tempdat_names <- colnames(tempdat)
    missing_master <- master_names[!master_names %in% tempdat_names]
    missing_tempdat <- tempdat_names[!tempdat_names %in% master_names]
    if(length(missing_master)>0){
      print(paste("Master colnames missing in", i, "being added to bind:"))
      print(missing_master)
      addtoTemp <- as.data.frame(matrix(ncol=length(missing_master), nrow = nrow(tempdat)))
      colnames(addtoTemp) <- missing_master
      tempdat <- cbind(tempdat, addtoTemp)
    }
    if(length(missing_tempdat)>0){
      print(paste(i, "colnames missing in master specimen data frame being added to bind:"))
      print(missing_tempdat)
      addtoMaster <- as.data.frame(matrix(ncol=length(missing_tempdat), nrow = nrow(specdat)))
      colnames(addtoMaster) <- missing_tempdat
      specdat <- cbind(specdat, addtoMaster)
    }
    # alphabetize temporarily to rbind; reorder cols when at last dataset read in
    tempdat <- tempdat[sort(colnames(tempdat))]
    specdat <- specdat[sort(colnames(specdat))]
    specdat <- rbind(specdat, tempdat)
  }
  if(i == specsheets[length(specsheets)]){
    # reorder columns
    specdat <- specdat[specmeta$Variable]
    print("All specimens imported and compiled in master data frame!")
    # clean up environment
    rm(addtoMaster, addtoTemp, tempdat, master_names, tempdat_names, missing_master, missing_tempdat, i)
  }
}


# -- SENSITIVITY CHECKS -----
# dry vs wet (does it matter? or can pool all specimens)
ggplot(specdat, aes(wgt_g, seeds, col = trt)) +
  geom_point(alpha = 0.5) +
  geom_smooth() +
  facet_grid(trt~species, scales = "free")

# april vs may samples for forbs (LACA)
subset(specdat, species == "LACA") %>%
  mutate(clip_month = ifelse(lubridate::month(clip_date)==4, "April", "May"),
         clip_month = ifelse(is.na(clip_month), "May", clip_month)) %>%
  ggplot(aes(wgt_g, seeds, col = clip_month)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se=F) +
  facet_wrap(~species, scales = "free") #hm.. april looks like better time to clip? maybe plants pooping out by may
# potentially should exclude may samples?
# > note: 8 missing rows are 8 individuals bagged together in april (CTW misunderstood initial instructions), so no individual weights. okay to ignore

# check how slope changes if exclude may
# everything
summary(lm(seeds~wgt_g, data = subset(specdat, species == "LACA")))
# april only
summary(lm(seeds~wgt_g, data = subset(specdat, species == "LACA" & clip_date == as.Date("2017-04-19"))))


# quick plot of esca
## compare flower, pod, wgt relationships to seed
select(esca, specimen, flowers, pods, seeds, wgt_g) %>%
  gather(met, val, flowers, pods, wgt_g) %>%
  ggplot(aes(val, seeds)) +
  geom_point(aes(col = as.factor(specimen))) +
  geom_smooth(method = "lm", col = "grey20") +
  facet_wrap(~met)
# what is the relationship of flowers to pods?
plot(esca$flowers, esca$pods)



# -- SIMPLIFY (PER LMH) AND JOIN WITH CLEANED COMPETITION BIOMASS -----
# > note: ctw cleaned up code here significantly since all datasets compiled in for loop
# > LMH previously simplified each into separate datasets bc read in as separate list items
# put it together!
allo.tog <- specdat %>%
  select(species, trt, specimen, clip_date, seeds, wgt_g) %>%
  filter(!is.na(wgt_g)) %>% # removes LACA plants that were bagged together before weighing
  rbind(data.frame(cbind(esca, trt = "ambient2020"))[names(.)])

# graph it all together!
ggplot(allo.tog, aes(x=wgt_g, y=seeds, color = trt)) + geom_point() + geom_smooth(method = "lm", se =F) + 
  facet_wrap(~species, scales = "free") + geom_smooth(data = allo.tog, aes(x=wgt_g, y=seeds), method = "lm", se = F, color = "black")
# treatment doesn't seem to have any effect, can group all together in models


# -- RUN LM LOOP -----
# build loop that:
# stores simple linear regression
# captures lm results
# captures predicted values based on field data with CIs and prediction intervals

species <- sort(unique(allo.tog$species))
# initiate data frames for storing:
allo.out <- data.frame()
predict_out <- data.frame()

for(i in species){
  print(paste("Iterating through", i, "regression"))
  
  # create lm predicting seeds from mass
  # > mar 2020: LS and LMH decide makes sense to force intercept through 0 since no biomass, no seeds
  # l <- lm(seeds~wgt_g, data = subset(allo.tog, species == i))  # let lm determine intercept
  l <- lm(seeds~0 + wgt_g, data = subset(allo.tog, species == i)) # force intercept at 0
  
  # store lm results
  temp_lm <- cbind(phytometer = i, tidy(l))
  temp_lm <- temp_lm[!names(temp_lm) == "statistic"]
  temp_lm <- rename(temp_lm, "est" = estimate, "se" = std.error, "pval" = p.value)
  # if intercept forced through 0, add in intercept term as 0 so it's clear that's what it is
  if(!any(grepl("intercept", temp_lm$term))){
    temp_lm <- rbind(temp_lm, data.frame(phytometer = i, term = "Intercept", est = 0, se = NA, pval = NA))
  }
  # spread out model results by model term (e.g. intercept, beta term), so each beta estimate, se, and pval for each predictor-side component in its own col
  # > this works whether the intercept goes through 0 or not
  temp_lm <- temp_lm %>%
    mutate(term = casefold(gsub("[(]|[)]", "", term))) %>%
    gather(met,val, est:pval) %>%
    unite(gosharks, term, met) %>%
    spread(gosharks, val)
  
  # add to allo.out
  allo.out <- rbind(allo.out, temp_lm)
  
  # compile predict results
  # > update 2020 aug 28 (after group mtg): want to predict BOTH based on inidividual weight, and total phytometer weight
  # > for forb phytometers, we typically only clipped up to 3 individuals but more may have grown in (clipped in apr, so left 4th and greater individuals in ground to see how they looked in may)
  # > for grasses, we typically clipped all phyto individuals (bc harvesting in may)
  # for-loop through both types of weight to append to master predicted dataset
  for(w in c("p.ind.wgt.g", "p_totwgt")){
    temp_predict <- phytodat[phytodat$phytometer == i, c("plot", "backgroundspp", "backgrounddensity", "phytometer", w)] # <- laurens changed this to dry wgt of the specimens.. not of the total phytometers in the plot..
    # rename wgt column as "wgt_g" so can run predict function on same lm created above
    colnames(temp_predict)[ncol(temp_predict)] <- "wgt_g"
    # crunch 95% confidence interval
    CIs <- predict.lm(l, newdata = temp_predict, interval = "confidence")
    colnames(CIs)[2:3] <- paste0(colnames(CIs)[2:3],"CI.95") # prefix to "upper" and "lower"
    # crunch 95% prediction interval
    PIs <- predict.lm(l, newdata = temp_predict, interval = "predict")
    colnames(PIs)[2:3] <- paste0(colnames(PIs)[2:3],"PI.95") # prefix to "upper" and "lower"
    # bind dataset used for prediction, CIs, PIs, and type of phyto dry mass (individual or total) used for prediction
    temp_predict <- cbind(temp_predict, CIs, PIs[,2:3], pfit_source = w)
    # append to predict_out
    predict_out <- rbind(predict_out, temp_predict)
  }
  
  # if last species, clean up and print done
  if(i == species[length(species)]){
    # lower case colnames, remove parentheses and "_est" from beta colnames, renumber rownames
    colnames(allo.out) <- casefold(gsub("[(]|[)]|_est", "", colnames(allo.out)))
    colnames(allo.out) <- gsub("wgt_g", "slope", colnames(allo.out))
    rownames(allo.out) <- seq(1,nrow(allo.out),1)
    allo.out <- data.frame(allo.out)
    # clarify wgt in prediction dataset is for phytometers, renumber rownames
    rownames(predict_out) <- seq(1,nrow(predict_out),1)
    predict_out <- rename(predict_out, pwgt_g = wgt_g) %>%
      data.frame()
    
    print("lm results and predictions compiled!")
  }
}

# preliminary look (NA = no background density bc control plot)
ggplot(predict_out, aes(x=backgroundspp, y = fit)) + geom_boxplot() +  
  geom_point(aes(col = pfit_source), alpha = 0.5, position = position_dodge(width = 0.3)) +
  facet_grid(phytometer~backgrounddensity, scales = "free_y")
# > note: NAs in prediction dataset are plots where data are missing (phytos not clipped/collected .. oopsies)


# -- FINISHING -----
# merge phyodat and predicted values
# > LMH requests wide form for prediction dataset (thinks that makes the most sense so phytodat not repeated.. can just run two different models based on column)
phytodat2 <- left_join(phytodat, select(predict_out, -pwgt_g)) %>%
  rename(p_seedfit = fit) %>%
  tbl_df()

# graph it all together!
ggplot(allo.tog) + 
  geom_point(aes(x=wgt_g, y=seeds, color = trt)) + 
  geom_abline(data = allo.out, aes(intercept = 0, slope = slope, lty = species)) +
  facet_wrap(~species, scales = "free")

# change spp colname in allometric table to more generic
allo.out <- rename(allo.out, species = phytometer)

# write out allometric table
write.csv(allo.out, paste0(datpath,"Competition_CleanedData/Competition_allometric_clean.csv"), row.names = F)
# write out phytodat with predicted values
write.csv(phytodat2, paste0(datpath,"Competition_CleanedData/Competition_phytometers_predicted.csv"), row.names = F)
