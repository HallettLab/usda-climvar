# calculate biomass-fecundity allometric relationship by competiton species
# authors: LMH, CTW
# initate: Jan 2019

# script purpose:
# derive allometric relationship of biomass to seed production (fecundity) by species by fall drought treatment
# compile predicted values, confidence intervals and prediction intervals (using `predict`)
# write out to competition cleaned data 

# notes: 
# specimens collected by fall dry and fall wet. check to see if treatment makes a difference on fecundity.
# more than 20 specimens collected for some species (e.g. april and late-season samples for forbs)

# > as of 2019-05-02, still need to address +20specimen collections for LACA, but writing out data table so we can move on for now

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
specmeta <- specmeta[23:nrow(specmeta),]
# set row one as colnames then remove
colnames(specmeta) <- specmeta[1,]; specmeta <- specmeta[-1,]


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

# check how slope changes if exclude may
# everything
summary(lm(seeds~wgt_g, data = subset(specdat, species == "LACA")))
# april only
summary(lm(seeds~wgt_g, data = subset(specdat, species == "LACA" & clip_date == as.Date("2017-04-19"))))


# -- SIMPLIFY (PER LMH) AND JOIN WITH CLEANED COMPETITION BIOMASS -----
# > note: ctw cleaned up code here significantly since all datasets compiled in for loop
# > LMH previously simplified each into separate datasets bc read in as separate list items
# put it together!
allo.tog <- specdat %>%
  select(species, trt, specimen, clip_date, seeds, wgt_g) %>%
  filter(!is.na(wgt_g)) # removes LACA plants that were bagged together before weighing

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
  l <- lm(seeds~wgt_g, data = subset(allo.tog, species == i))
  # store lm results
  temp_lm <- cbind(phytometer = c(rep(i,2)), tidy(l))
  temp_lm <- temp_lm[,-5]
  colnames(temp_lm)[3:5] <- c("est", "se","pval")
  temp_lm <- temp_lm %>%
    #mutate(term = casefold(gsub("[(]|[)]", "", term))) %>%
    gather(met,val, est:pval) %>%
    unite(gosharks, term, met) %>%
    spread(gosharks, val)
  
  # add to allo.out
  allo.out <- rbind(allo.out, temp_lm)
  
  # compile predict results
  temp_predict <- phytodat[phytodat$phytometer == i, c("plot", "backgroundspp", "backgrounddensity", "phytometer", "p_totwgt")]
  colnames(temp_predict)[ncol(temp_predict)] <- "wgt_g"
  CIs <- predict.lm(l, newdata = temp_predict ,interval = "confidence")
  colnames(CIs)[2:3] <- paste0(colnames(CIs)[2:3],"CI.95") 
  PIs <- predict.lm(l, newdata = temp_predict ,interval = "predict")
  colnames(PIs)[2:3] <- paste0(colnames(PIs)[2:3],"PI.95") 
  temp_predict <- cbind(temp_predict, CIs, PIs[,2:3]) %>%
    rename(p_totwgt = wgt_g)
  
  # add to predict_out
  predict_out <- rbind(predict_out, temp_predict)
  
  # if last species, clean up and print done
  if(i == species[length(species)]){
    allo.out <- allo.out[,c(1:2,4,3,5,7,6)]
    # lower case colnames, remove parentheses and "_est" from beta colnames
    colnames(allo.out) <- casefold(gsub("[(]|[)]|_est", "", colnames(allo.out)))
    colnames(allo.out) <- gsub("wgt_g", "slope", colnames(allo.out))
    rownames(allo.out) <- seq(1,nrow(allo.out),1)
    rownames(predict_out) <- seq(1,nrow(predict_out),1)
    print("lm results and predictions compiled!")
  }
}


# -- FINISHING -----
# merge phyodat and predicted values
phytodat2 <- left_join(phytodat, predict_out) %>%
  rename(p_totwgt_seedfit = fit)

# write out allometric table
write.csv(allo.out, paste0(datpath,"Competition_CleanedData/Competition_allometric_clean.csv"), row.names = F)
# write out phytodat with predicted values
write.csv(phytodat2, paste0(datpath,"Competition_CleanedData/Competition_phytometers_predicted.csv"), row.names = F)
