
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

dat_bp =  covid_ts %>%
  select(date=Date, backpro=Est_Onsets) %>%
  filter(date<=lubridate::ymd("2020-05-01"))

# Estimate segmented regression with 3 breakpoints based on discrete optimization,
# can take several minutes, we load saved results
res_3_disc = estimate_bp_disc_optim(data = dat_bp, bp = 3)
save(res_3_disc, file = "D:/Bavaria/before_may (quasi-poisson)/results/res_disc_opt_3bp_onset.RData")
load("D:/Bavaria/before_may (quasi-poisson)/results/res_disc_opt_3bp_onset.RData")

# Estimate segmented regression with 4 Breakpoints based on discrete optimization,
# can take several minutes, we load saved results
res_4_disc = estimate_bp_disc_optim(data = dat_bp, bp = 4)
save(res_4_disc, file = "D:/Bavaria/before_may (quasi-poisson)/results/res_disc_opt_4bp_onset.RData")
load("D:/Bavaria/before_may (quasi-poisson)/results/res_disc_opt_4bp_onset.RData")

# Estimate segmented regression with 3 or 4 Breakpoints based on 'segmented' package with startpoints from
# discrete optimization
res_3_seg = estimate_bp_segmented(
  data           = dat_bp,
  bp             = 3,
  start_bp       = res_3_disc$bps,
  segmented_seed = 0520)
res_4_seg = estimate_bp_segmented(
  data           = dat_bp,
  bp             = 4,
  start_bp       = res_4_disc$bps-1,
  segmented_seed = 0520)

# Compare deviance and dispersion of models
round(
  c(
    dev_3_bp = summary(res_3_seg$segmented_model)$deviance,
    overdisp_3_bp = summary(res_3_seg$segmented_model)$dispersion,
    dev_4_bp = summary(res_4_seg$segmented_model)$deviance,
    overdisp_4_bp = summary(res_4_seg$segmented_model)$dispersion),
  2)

############################## Results #########################################
# Plot of estimated segmented regression model
png("D:/Bavaria/before_may (quasi-poisson)/results/disease-onset-3bp.png", width = 1000, height = 460)
res_3_seg$plot
dev.off()
png("D:/Bavaria/before_may (quasi-poisson)/results/disease-onset-4bp.png", width = 1000, height = 460)
res_4_seg$plot
dev.off()

# More info on estimated breakpoints
# Breakpoints + 95%-CIs
knitr::kable(res_4_seg$breakpoints, "html") %>%
  readr::write_file("D:/covid-19 analysis/bavarian/results/disease-onset-change-points.html")
# Daily multiplicative change in infection numbers
knitr::kable(res_4_seg$coef, "html", digits = 3) %>%
  readr::write_file("D:/covid-19 analysis/bavarian/results/disease-onset-factors.html")
