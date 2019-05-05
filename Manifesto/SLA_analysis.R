library(tidyr)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)

se <- function(x) {
  sd(x)/sqrt(length(x))
}

datpath <- "~/Dropbox/ClimVar/DATA/Plant_composition_data/"
key <- read.csv(paste0(datpath, "Shelter_key.csv"))
sla <- read.csv(paste0(datpath, "SLA/SLA_CleanedData/SLA_leafwgts_2015-2016_cleaned.csv"))

SLA <- merge(key, sla) %>%
  tbl_df() %>%
  mutate(treatment=ordered(treatment, levels = c(consistentDry="consistentDry", springDry="springDry", fallDry="fallDry", controlRain="controlRain"))) 


SLA2 <- SLA %>%
  group_by(treatment, species) %>%
  summarize(meanSLA = mean(SLA), seSLA = se(SLA), 
            meanArea = mean(leaf_area_cm2), seArea = se(leaf_area_cm2))


A <- ggplot(SLA2, aes(x=treatment, y=meanSLA, fill=species)) + 
  geom_bar(stat="identity", position="dodge") + 
  geom_errorbar(aes(ymin = meanSLA - seSLA, ymax = meanSLA + seSLA)) + 
  facet_wrap(~species) + theme_bw() + theme(legend.position = "none") +
  xlab("Treatment") + ylab ("SLA (cm2/ g)")


B <- ggplot(SLA2, aes(x=treatment, y=meanArea, fill=species)) + 
  geom_bar(stat="identity", position="dodge") + 
  geom_errorbar(aes(ymin = meanArea - seArea, ymax = meanArea + seArea)) + 
  facet_wrap(~species) + theme_bw() + theme(legend.position = "none") +
  xlab("Treatment") + ylab ("Leaf area (cm2)")

pdf(paste0(datpath, "SLA/SLA_Figures/SLA_leafArea_bytrt.pdf"))
grid.arrange(A, B)
dev.off()
