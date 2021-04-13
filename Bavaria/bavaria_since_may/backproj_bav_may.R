# Load packages
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


quantiles_incu <- data.frame(q=c(0.5, 0.975), value=c(5.1, 11.5))
# Fit distribution
f_target2 <- function(theta) {
  qlnorm(quantiles_incu$q[c(1,2)], theta[1], theta[2]) - quantiles_incu$value[c(1,2)]
}
incu_lnorm <- nleqslv::nleqslv(c(1,1), f_target2)$x
# Compare observed and fitted value
data.frame(
  q        = quantiles_incu$q,
  observed = quantiles_incu$value,
  fitted   = qlnorm( quantiles_incu$q, incu_lnorm[1], incu_lnorm[2]))

inc_pdf = data.frame(t = seq(0,15, length=1000)) %>%
  mutate(pdf = dlnorm(t, incu_lnorm[1], incu_lnorm[2]))

# Discretize incubation period distribution
cdf <- plnorm(seq(0,14,by=1), incu_lnorm[1], incu_lnorm[2])
pmf <- structure(c(0,diff(cdf)), names=seq(0,14,by=1))
# Normalize the discrete incubation period distribution
pmf <- pmf/sum(pmf)
df <- data.frame(days=as.numeric(names(pmf)), pmf=pmf)

#' The backprojection can be done using the function `surveillance::backprojNP`
#' The backprojected curve shows the number of infections per day and can be
#' compared to interventions similar to Werber et al. (2013) (https://doi.org/10.1093/aje/kwt069)

# Extract data from nowcast
sts_symp <- sts(
  epoch       = covid_ts$Date,
  observed    = matrix(covid_ts$Est_Onsets, ncol = 1, nrow = nrow(covid_ts)),
  epochAsDate = TRUE)

# Perform back projection with smoothing to adjust for weekday effects (k=6)
bp <- backprojNP(sts_symp, incu.pmf=pmf, control=list(k=6, eq3a.method="C"))
# Reduce to relevant subset of the data making it comparable with Dehning et al. (2020)
bp <- bp[epoch(bp) >= as.Date("2020-05-01")]

bpdf <- bp %>%
  as.data.frame() %>%
  mutate(epoch_numeric=as.numeric(epoch), t = epoch - max(epoch) + 1)



#'
source("D:/Bavaria/before_may (quasi-poisson)/breakpoint_fun.R")
# Get necessary data (date and backprojected case numbers)
dat_bp_back_may = bpdf %>%
  select(date=epoch, backpro = upperbound) %>%
  filter(date>=lubridate::ymd("2020-05-01"))
# Estimate segmented regression with 3 breakpoints based on discrete optimization,
# can take several minutes, we load saved results
backproj_bav_may_3bp = estimate_bp_disc_optim(data = dat_bp_back_may, bp = 3)
save(backproj_bav_may_3bp, file = "D:/Bavaria/bavaria_since_may/backproj_bav_may_3bp.RData")
load("D:/Bavaria/bavaria_since_may/backproj_bav_may_3bp.RData")

# Estimate segmented regression with 4 Breakpoints based on discrete optimization,
# can take several minutes, we load saved results
backproj_bav_may_4bp= estimate_bp_disc_optim(data = dat_bp_back_may, bp = 4)
save(backproj_bav_may_4bp, file = "D:/bavaria/bavaria_since_may/backproj_bav_may_4bp.RData")
load("D:/Bavaria/bavaria_since_may/backproj_bav_may_4bp.RData")

# Estimate segmented regression with 3 or 4 Breakpoints based on 'segmented' package with startpoints from
# discrete optimization
backproj_3_seg_bav_may = estimate_bp_segmented(
  data           = dat_bp_back_may,
  bp             = 3,
  start_bp       = backproj_bav_may_3bp$bps,
  segmented_seed = 0520)
backproj_4_seg_bav_may = estimate_bp_segmented(
  data           = dat_bp_back_may,
  bp             = 4,
  start_bp       = backproj_bav_may_4bp$bps,
  segmented_seed = 0520)

backproj_5_seg_bav_may = estimate_bp_segmented(
  data           = dat_bp_back_may,
  bp             = 5,
  start_bp       = c(30,60,90,100,123),
  segmented_seed = 0520)#70.01
backproj_6_seg_bav_may = estimate_bp_segmented(
  data           = dat_bp_back_may,
  bp             = 6,
  start_bp       = c(16,28,62,90,99,123),
  segmented_seed = 0520)#50.08
backproj_7_seg_bav_may = estimate_bp_segmented(
  data           = dat_bp_back_may,
  bp             = 7,
  start_bp       = c(16,28,56,62,90,99,123),
  segmented_seed = 0520)

# Compare deviance and dispersion of models
c(
    dev_3_bp      = summary(backproj_3_seg_bav_may)$deviance,
    dis_3_bp = summary(backproj_3_seg_bav_may$segmented_model)$dispersion,
    dev_4_bp      = summary(backproj_4_seg_bav_may$segmented_model)$deviance,
    dis_4_bp = summary(backproj_4_seg_bav_may$segmented_model)$dispersion)
  
c(
  dev_5_bp      = summary(backproj_5_seg_bav_may$segmented_model)$deviance,
  dis_5_bp = summary(backproj_5_seg_bav_may$segmented_model)$dispersion,
  dev_6_bp      = summary(backproj_6_seg_bav_may$segmented_model)$deviance,
  dis_6_bp = summary(backproj_6_seg_bav_may$segmented_model)$dispersion,
  dev_7_bp      = summary(backproj_7_seg_bav_may$segmented_model)$deviance,
  dis_7_bp = summary(backproj_7_seg_bav_may$segmented_model)$dispersion)

############################## Results #########################################
# Plot of estimated segmented regression model
### Use model with 3 change points
png("D:/Bavaria/bavaria_since_may/backproj_4_seg_bav_may.png", width = 1200, height = 480)
backproj_4_seg_bav_may$plot
dev.off()


# More info on estimated breakpoints
# Breakpoints + 95%-CIs
knitr::kable(backproj_7_seg_bav_may$breakpoints, "html") %>%
  readr::write_file("D:/covid-19Analysis/bavaria_since_may/backproj-7change-points-may.html")


