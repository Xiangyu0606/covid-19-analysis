library(tidyverse)
library(segmented)
library(surveillance)
library(future)
library(future.apply)

now=lubridate::ymd("2020-09-24")
#' Pick date and estimated number of disease onsets that day as the only two columns
covid = read_tsv("D:/bavaria/before_may (quasi-poisson)/data/nowcasting_results_2020-09-24.csv") %>%
  rename(Date=date, Est_Onsets=nowcast_med) %>%
  dplyr::select(Date, Est_Onsets)


# Stop 21 days before now, in order to avoid any nowcasting effects.
covid_ts <- covid %>% filter(Date <= now - 21) %>% arrange(Date)


#'
source("D:/Bavaria/before_may_nb/breakpoint_fun_onset_bav_nb.R")
# Get necessary data (date and  case numbers)

dat_bp_may_nb =  covid_ts %>%
  select(date=Date, backpro=Est_Onsets) %>%
  filter(date>=lubridate::ymd("2020-05-01"))

# Estimate segmented regression with 3 breakpoints based on discrete optimization,
# can take several minutes, we load saved results
disease_bav_may_nb_3bp = estimate_bp_disc_optim(data = dat_bp_may_nb, bp = 3)
save(disease_bav_may_nb_3bp, file = "D:/Bavaria/bavaria_since_may/NegBinomial/disease_bav_may_nb_3bp.RData")
load("D:/Bavaria/bavaria_since_may/NegBinomial/disease_bav_may_nb_3bp.RData")  

# Estimate segmented regression with 4 Breakpoints based on discrete optimization,
# can take several minutes, we load saved results
disease_bav_may_nb_4bp = estimate_bp_disc_optim(data = dat_bp_may_nb, bp = 4)
save(disease_bav_may_nb_4bp, file = "D:/Bavaria/bavaria_since_may/NegBinomial/disease_bav_may_nb_4bp.RData")
load("D:/Bavaria/bavaria_since_may/NegBinomial/disease_bav_may_nb_4bp.RData")   

disease_3_seg_may_nb = estimate_bp_segmented(
  data           = dat_bp_may_nb,
  bp             = 3,
  start_bp       = disease_bav_may_nb_3bp$bps,
  segmented_seed = 0520)  
disease_4_seg_may_nb = estimate_bp_segmented(
  data           = dat_bp_may_nb,
  bp             = 4,
  start_bp       = disease_bav_may_nb_4bp$bps-1,
  segmented_seed = 0520) 

disease_5_seg_may_nb = estimate_bp_segmented(
  data           = dat_bp_may_nb,
  bp             = 5,
  start_bp       = c(33,44,61,66,81),
  segmented_seed = 0520)
disease_6_seg_may_nb = estimate_bp_segmented(
  data           = dat_bp_may_nb,
  bp             = 6,
  start_bp       = c(33,44,51,73,93,105),
  segmented_seed = 0520)
# Compare deviance and dispersion of models
c(
  dis_3bp = summary(disease_3_seg_may_nb$segmented_model)$dispersion,
  dev_3bp = summary(disease_3_seg_may_nb$segmented_model)$deviance,
  aic_3bp = summary(disease_3_seg_may_nb$segmented_model)$aic,
  dis_4bp = summary(disease_4_seg_may_nb$segmented_model)$dispersion,
  dev_4bp = summary(disease_4_seg_may_nb$segmented_model)$deviance,
  aic_4bp = summary(disease_4_seg_may_nb$segmented_model)$aic,
  dis_5bp = summary(disease_5_seg_may_nb$segmented_model)$dispersion,
  dev_5bp = summary(disease_5_seg_may_nb$segmented_model)$deviance,
  aic_5bp = summary(disease_5_seg_may_nb$segmented_model)$aic,
  dis_6bp = summary(disease_6_seg_may_nb$segmented_model)$dispersion,
  dev_6bp = summary(disease_6_seg_may_nb$segmented_model)$deviance,
  aic_6bp = summary(disease_6_seg_may_nb$segmented_model)$aic)


############################## Results #########################################
# Plot of estimated segmented regression model
png("D:/bavaria/bavaria_since_may/NegBinomial/disease-onset_bav_may_nb_4BP.png", width = 1200,height = 480)
disease_4_seg_may_nb$plot
dev.off()

# More info on estimated breakpoints
# Breakpoints + 95%-CIs
knitr::kable(disease_4_seg_may_nb$breakpoints, "html") %>%
  readr::write_file("D:/covid-19Analysis/bavaria_since_may/NegBinomial/disease-4change-points-may-nb.html")

