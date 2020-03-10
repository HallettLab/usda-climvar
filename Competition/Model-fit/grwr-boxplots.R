source("./Competition/Model-fit/run_pairwise_coexistence.R")

library(tidyverse)

dry <- invasion_dry %>%
  tbl_df() %>%
  gather("trt", "grwr", 1:6) %>%
  separate(trt, c("invader", "into", "resident", "treatment")) %>%
  select(-into)

wet <- invasion_wet %>%
  tbl_df() %>%
  gather("trt", "grwr", 1:6) %>%
  separate(trt, c("invader", "into", "resident")) %>%
  select(-into) %>%
  mutate(treatment = "wet")

grwr <- rbind(dry, wet)

ggplot(grwr, aes(x=treatment, y = grwr, color = resident)) + geom_boxplot() + facet_wrap(~invader) + 
  geom_hline(yintercept = 0, lty = "dashed")
