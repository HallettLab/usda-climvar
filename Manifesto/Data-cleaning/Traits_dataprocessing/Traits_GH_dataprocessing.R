library(tidyverse)
library(FD)
library(readr)


## Read in trait data
trts <- read.csv("~/Dropbox/ClimVar/DATA/Plant_composition_data/Traits/Traits_GHscreening_Brad/Bradtraits_cleaned.csv") %>%
  tbl_df()

## Read in cover data
veg <- read.csv("~/Dropbox/ClimVar/DATA/Plant_composition_data/Cover/Cover_CleanedData/ClimVar_species-cover.csv") %>%
  tbl_df() %>%
  mutate(species_name = as.character(species_name))


## Clean up species names and group similar species for trait analyses
veg$species_name <- ifelse(veg$species_name == "Bromus madritensis madritensis", 
                           "Bromus madritensis", 
                           veg$species_name)

veg$species_name <- ifelse(veg$species_name == "Avena barbata", 
                           "Avena sp", 
                           veg$species_name)

veg$species_name <- ifelse(veg$species_name == "Avena fatua", 
                           "Avena sp", 
                           veg$species_name)


veg$species_name <- ifelse(veg$species_name == "Bromus sterilis", 
                           "Bromus diandrus", 
                           veg$species_name)

veg$species_name <- ifelse(veg$species_name == "Cynosaurus echinatus", 
                           "Cynosurus echinatus", 
                           veg$species_name)

veg$species_name <- ifelse(veg$species_name == "Hordeum marinum", 
                           "Hordeum murinum", 
                           veg$species_name)

veg$species_name <- ifelse(veg$species_name == "Erodium moschatum", 
                           "Erodium botrys", 
                           veg$species_name)


veg$species_name <- ifelse(veg$species_name == "Vulpia myorus" | veg$species_name == "Vulpia bromoides", 
                           "Vulpia myuros", 
                           veg$species_name)

veg <- veg %>%
  group_by(plot, species_name, subplot, year, genus, func, status, treatment,
           shelterBlock, func2) %>%
  summarize(cover = sum(cover)) %>%
  tbl_df() %>%
  group_by(plot, subplot) %>%
  mutate(totcover = sum(cover)) %>%
  mutate(relcover = cover/totcover) %>%
  filter(subplot != "L")

##  Parallel these changes in the trait data
trts$Taxon <- as.character(trts$Taxon)
trts$Taxon <- ifelse(trts$Taxon == "Avena fatua", 
                     "Avena sp", 
                     trts$Taxon)

trts$species_name  <- trts$Taxon


## Compare species overlap between trts and veg
sort(unique(trts$Taxon))
sort(unique(veg$species_name))

sppoverlap <- intersect(trts$Taxon, veg$species_name)


### clean up veg data
vegdat <- veg[which(veg$species_name%in%sppoverlap),] %>%
  arrange(plot, subplot, year, treatment, shelterBlock, species_name) %>%
  filter(species_name != "Lupinus bicolor") %>%
  select(plot, subplot, year, treatment, shelterBlock, species_name, cover) %>%
  spread(species_name, cover, fill = 0) %>%
  filter(subplot != "L")

vegdat2 <- veg[which(veg$species_name%in%sppoverlap),] %>%
  filter(subplot != "L")


traitdata <- trts[which(trts$species_name%in%sppoverlap),] %>%
  filter(species_name != "Lupinus bicolor") %>%
  arrange(species_name) %>%
  select(species_name, Ht, LDMC, SLA, RMF, Total, RGR, Dens, DiamC, SRLC, SRLF, PropF) %>%
  unique() 

## format for FD
abundance <- as.matrix(vegdat[1:nrow(vegdat), 6:ncol(vegdat)])
dim(abundance)

dim(traitdata)

traits <- as.matrix(traitdata[1:nrow(traitdata), 2:12])
row.names(traits) <- traitdata$species_name

#T=number of traits
T<-dim(traits)[2]
t2<-dim(traits)[1]
t2


#c=number of communities
C<-dim(abundance)[1]
c2<-dim(abundance)[2]
c2


#check coherence of number of species in 'traits' and 'abundances'
if(dim(abundance)[2]!=dim(traits)[1])stop("error:differentnumberofspeciesin'traits'and'abundances'matrices")

##Species names must be in same order in both files

# calculate trait diversity
results <- dbFD(traits, abundance, corr="cailliez")
results <- as.data.frame(results) %>%
  tbl_df()

mykey <- vegdat[,1:5] %>%
  tbl_df()

# put together with the key
trtout <- cbind(mykey, results) %>%
  tbl_df() %>%
  mutate(treatment=ordered(treatment, levels = c(consistentDry="consistentDry", springDry="springDry", fallDry="fallDry", controlRain="controlRain"))) 

## add in species columns
trtout2 <- merge(trtout, vegdat) %>%
  tbl_df()
names(trtout2)
names(trtout2)[26] = "Avena"

### and grass to forb ratios ###

gf <- veg %>%
  group_by(plot, subplot, treatment, year, func) %>%
  summarize(cover = sum(cover)) %>%
  spread(func, cover) %>%
  mutate(propForb = forb/(grass + forb))
names(gf)[5:6] = c("forbCover", "grassCover")

trtout2.5 <- left_join(trtout2, gf)

## add in ANPP
anpp <- read.csv(("~/Dropbox/ClimVar/DATA/Plant_composition_data/ANPP/ANPP_CleanedData/ClimVar_ANPP-Peak.csv")) %>%
  mutate(ANPPgm = weight_g_m) %>%
  select(-weight_g_m)
trtout3 <- left_join(trtout2.5, anpp) %>%
  mutate(treatment2=ordered(treatment, levels = c(consistentDry="consistentDry", springDry="springDry", fallDry="fallDry", controlRain="controlRain"))) %>%
  mutate(spring = as.factor(ifelse(treatment == "consistentDry" | treatment == "springDry", 0, 1)),
         fall = as.factor(ifelse(treatment == "consistentDry" | treatment == "fallDry", 0, 1)))

write.csv(trtout3, "~/Dropbox/ClimVar/DATA/Plant_composition_data/Traits/Traits_ProcessedData-GH/ClimVar_trait-diversity-GH.csv", row.names = F)
