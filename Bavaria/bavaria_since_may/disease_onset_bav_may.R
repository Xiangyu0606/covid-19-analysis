

library(tidyverse)
library(segmented)
library(surveillance)
library(future)
library(future.apply)


now=lubridate::ymd("2020-09-24")
#' Pick date and estimated number of disease onsets that day as the only two columns
covid = read_tsv("D:/Bavaria/before_may (quasi-poisson)/data/nowcasting_results_2020-09-24.csv") %>%
  rename(Date=date, Est_Onsets=nowcast_med) %>%
  dplyr::select(Date, Est_Onsets)


# Stop 21 days before now, in order to avoid any nowcasting effects.
covid_ts <- covid %>% filter(Date <= now - 21) %>% arrange(Date)


source("D:/Bavaria/before_may (quasi-poisson)/breakpoint_fun_onset.R")
# Get necessary data (date and  case numbers)

dat_bp_may =  covid_ts %>%
  select(date=Date, backpro=Est_Onsets) %>%
  filter(date>=lubridate::ymd("2020-05-01"))



# Estimate segmented regression with 3 breakpoints based on discrete optimization,
# can take several minutes, we load saved results
disease_3bp_may = estimate_bp_disc_optim(data = dat_bp_may, bp = 3)
save(disease_3bp_may, file = "D:/Bavaria/bavaria_since_may/disease_3bp_may.RData")
load("D:/Bavaria/bavaria_since_may/disease_3bp_may.RData")

# Estimate segmented regression with 4 Breakpoints based on discrete optimization,
# can take several minutes, we load saved results
disease_4bp_may = estimate_bp_disc_optim(data = dat_bp_may, 4)
save(disease_4bp_may, file = "D:/Bavaria/bavaria_since_may/disease_4bp_may.RData")
load("D:/Bavaria/bavaria_since_may/disease_4bp_may.RData")


# Estimate segmented regression with 3 or 4 Breakpoints based on 'segmented' package with startpoints from
# discrete optimization
disease_3_seg_may = estimate_bp_segmented(
  data           = dat_bp_may,
  bp             = 3,
  start_bp       = disease_3bp_may$bps,
  segmented_seed = 0520)
disease_4_seg_may = estimate_bp_segmented(
  data           = dat_bp_may,
  bp             = 4,
  start_bp       = disease_4bp_may$bps-1,
  segmented_seed = 0520)

disease_5_seg_may = estimate_bp_segmented(
  data           = dat_bp_may,
  bp             = 5,
  start_bp       = c(11,36,75,99,103),
  segmented_seed = 0520)
disease_6_seg_may = estimate_bp_segmented(
  data           = dat_bp_may,
  bp             = 6,
  start_bp       = c(11,36,75,92,108,122),
  segmented_seed = 0520) 

# Compare deviance and dispersion of models
c (dis_3bp = summary(disease_3_seg_may$segmented_model)$dispersion,
   dev_3bp = summary(disease_3_seg_may$segmented_model)$deviance,
   dis_4bp = summary(disease_4_seg_may$segmented_model)$dispersion,
   dev_4bp = summary(disease_4_seg_may$segmented_model)$deviance,
   dis_5bp = summary(disease_5_seg_may$segmented_model)$dispersion,
   dev_5bp = summary(disease_5_seg_may$segmented_model)$deviance,
   dis_6bp = summary(disease_6_seg_may$segmented_model)$dispersion,
   dev_6bp = summary(disease_6_seg_may$segmented_model)$deviance)
############################## Results #########################################
# Plot of estimated segmented regression model

png("D:/Bavaria/bavaria_since_may/bavaria-disease-onset-may-6bp.png", width = 1200, height = 480)
disease_6_seg_may$plot
dev.off()



# More info on estimated breakpoints
# Breakpoints + 95%-CIs
knitr::kable(disease_6_seg_may$breakpoints, "html") %>%
  readr::write_file("D:/covid-19Analysis/bavaria_since_may/disease-onset-3change-points-may.html")

