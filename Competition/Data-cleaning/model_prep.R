

# -- SETUP -----
rm(list = ls()) # clean environment
library(tidyverse)
options(stringsAsFactors = F)
theme_set(theme_bw())
na_vals = c(" ","", NA, "NA")

# set pathway to climvar dropbox competition data folder
datpath <- "~/Dropbox/ClimVar/Competition/Data/"

# read in data
# cleaned, combined competition dataset (cleaned background, phytometers, predicted vals, plot treatments)
comp.dat <- read.csv(paste0(datpath, "Competition_CleanedData/Competition_combined_clean.csv"),
                     na.strings = na_vals, strip.white = T) %>%
  tbl_df()

# and bring back allodat
allodat <- read.csv(paste0(datpath,"Competition_CleanedData/Competition_allometric_clean.csv"),
                    na.strings = na_vals, strip.white = T) %>%
  rename(bcode4 = phytometer)

phyto.dat <- comp.dat %>%
  select(plot:bdensity, insitu_bdisturbed, phytometer, pcode4, insitu_pstems, p_totwgt_seedfit ) %>%
  filter(!is.na(p_totwgt_seedfit)) %>%
  mutate(p_totwgt_seedfit = ifelse(insitu_pstems == 0, 0, p_totwgt_seedfit)) %>%
  select(-phytometer, -bcode4) %>%
  rename(seedsIn = insitu_pstems, seedsOut = p_totwgt_seedfit) 

back.dat0 <- comp.dat %>%
  filter(!is.na(background)) %>%
  select(plot:bdensity,insitu_bdisturbed, b.ind.wgt.g, seedsAdded, insitu_plot_bdensity) %>%
  unique() 

back.dat <- left_join(back.dat0, allodat) %>%
  mutate(p_totwgt_seedfit = intercept + b.ind.wgt.g*insitu_plot_bdensity*slope) %>%
  mutate(p_totwgt_seedfit = ifelse(insitu_plot_bdensity == 0, 0, p_totwgt_seedfit)) %>%
  select(plot:bdensity,insitu_bdisturbed,seedsAdded, p_totwgt_seedfit) %>%
  rename(seedsIn = seedsAdded, seedsOut = p_totwgt_seedfit, pcode4 = bcode4) %>%
  filter(!is.na(seedsOut))

tog <- rbind(phyto.dat, back.dat) %>%
  filter(pcode4!= "TRHI") %>%
  gather(var, val, seedsIn:seedsOut) %>%
  mutate(varnew = paste(pcode4, var, sep = "_")) %>%
  select(-pcode4, -var) %>%
  spread(varnew, val)
  
brvu <- tog %>%
  filter(background == "Bromus hordeaceus" | background == "Vulpia myuros" | is.na(background))
  
m1E <- as.formula(log(BRHO_seedsOut +1) ~  log((BRHO_seedsIn+1)*exp(log(lambda)-
                                                                      log((1+aiE*(BRHO_seedsIn+1)+aiA*(VUMY_seedsIn+1))))))


treatments <- unique(brvu$falltreatment)

ERoutput <- as.data.frame(matrix(nrow = 0, ncol = 7))
names(ERoutput) = c("estimate", "se", "t", "p", "params", "treatment", "species")
for (i in 1:length(treatments)){
  m1out <- nlsLM(m1E, start=list(lambda=1, aiE = .01, aiA=.01),
                 lower = c(0, 0, 0), upper = c(200, 1, 1),
                 control=nls.lm.control(maxiter=500), trace=T,
                 data = subset(brvu, falltreatment == treatments[i]))
  
  outreport <- as.data.frame(summary(m1out)$coef[1:3, 1:4])
  names(outreport) = c("estimate", "se", "t", "p")
  outreport$params <- row.names(outreport)
  outreport$treatment <- treatments[i]
  outreport$species <- "Bromus"
  ERoutput <- rbind(ERoutput, outreport)
}


brvu <- tog %>%
  filter(background == "Bromus hordeaceus" | background == "Vulpia myuros" | is.na(background))

m1E <- as.formula(log(VUMY_seedsOut +1) ~  log((VUMY_seedsIn+1)*exp(log(lambda)-
                                                                      log((1+aiE*(BRHO_seedsIn+1)+aiA*(VUMY_seedsIn+1))))))


treatments <- unique(brvu$falltreatment)

ERoutput <- as.data.frame(matrix(nrow = 0, ncol = 7))
names(ERoutput) = c("estimate", "se", "t", "p", "params", "treatment", "species")
for (i in 1:length(treatments)){
  m1out <- nlsLM(m1E, start=list(lambda=1, aiE = .01, aiA=.01),
                 lower = c(0, 0, 0), upper = c(200, 1, 1),
                 control=nls.lm.control(maxiter=500), trace=T,
                 data = subset(brvu, falltreatment == treatments[i]))
  
  outreport <- as.data.frame(summary(m1out)$coef[1:3, 1:4])
  names(outreport) = c("estimate", "se", "t", "p")
  outreport$params <- row.names(outreport)
  outreport$treatment <- treatments[i]
  outreport$species <- "Vulpia"
  ERoutput <- rbind(ERoutput, outreport)
}

# plot, falltreat, treatment, shelterBlock, uniqueplot (background), species seedsin seeds out
  
  ## graph just phyto - seeds
phyto.dat <- comp.dat %>%
  select(plot:bdensity, insitu_bdisturbed, phytometer, pcode4, insitu_pstems, p_totwgt_seedfit ) %>%
  filter(!is.na(p_totwgt_seedfit)) %>%
  mutate(p_totwgt_seedfit = ifelse(insitu_pstems == 0, 0, p_totwgt_seedfit)) %>%
  select(-phytometer, -bcode4) %>%
  rename(seedsIn = insitu_pstems, seedsOut = p_totwgt_seedfit) 

back.dat0 <- comp.dat %>%
  filter(!is.na(background)) %>%
  select(plot:bdensity,insitu_bdisturbed, b.ind.wgt.g, insitu_plot_bdensity) %>%
  unique()

back.dat2 <- left_join(back.dat0, allodat) %>%
  mutate(p_totwgt_seedfit = intercept + b.ind.wgt.g*slope) %>%
  mutate(p_totwgt_seedfit = ifelse(insitu_plot_bdensity == 0, 0, p_totwgt_seedfit)) %>%
  select(plot:bdensity,insitu_bdisturbed, p_totwgt_seedfit) %>%
  rename(seedsOut = p_totwgt_seedfit, pcode4 = bcode4) %>%
  filter(!is.na(seedsOut)) %>%
  mutate(seedsIn = 1)

phytotog <- rbind(phyto.dat, back.dat2) %>%
  mutate(R = seedsOut/seedsIn,
         R = ifelse(is.na(R), 0, R)) %>%
  filter(pcode4 != "TRHI" & pcode4 != "ESCA")

phytotog2 <- phytotog %>%
  group_by(falltreatment, background, bdensity, pcode4) %>%
  summarize(meanR = mean(R), seR = sd(R)/sqrt(length(R))) %>%
  filter(!is.na(background))

ggplot(phytotog2, aes(x= falltreatment, y = meanR, color = background, group = background)) +
  facet_grid(pcode4~bdensity, scales = "free") + geom_point() + geom_line()

ggplot(filter(phytotog2, background%in%c("Avena fatua", "Bromus hordeaceus", "Vulpia myuros") & pcode4 != "LACA"), 
       aes(x= falltreatment, y = meanR, color = background, group = background)) +
  facet_grid(pcode4~bdensity, scales = "free") + geom_point() + geom_line() + 
  geom_errorbar(aes(ymin = meanR - seR, ymax = meanR + seR), width = .2)


# graph just phyto - bmass

## graph just phyto - seeds
phyto.dat2 <- comp.dat %>%
  select(plot:bdensity, insitu_bdisturbed, phytometer, pcode4, p.ind.wgt.g ) %>%
   select(-phytometer, -bcode4) 

back.dat3 <- comp.dat %>%
  filter(!is.na(background)) %>%
  select(plot:bdensity,insitu_bdisturbed, b.ind.wgt.g) %>%
  unique() %>%
  rename(pcode4 = bcode4, p.ind.wgt.g = b.ind.wgt.g)


phytotog.bmass <- rbind(phyto.dat2, back.dat3) %>%
  filter(!is.na(background))

phytotog.bmass.mean <- phytotog.bmass %>%
  group_by(falltreatment, background, bdensity, pcode4) %>%
  summarize(meanbmass = mean(p.ind.wgt.g, na.rm=T), sebmass = sd(p.ind.wgt.g)/sqrt(length(p.ind.wgt.g))) %>%
  filter(!is.na(background))

ggplot(phytotog.bmass.mean, aes(x= falltreatment, y = meanbmass, color = background, group = background)) + 
  facet_grid(pcode4~bdensity, scales = "free") + geom_point() + geom_line() 
