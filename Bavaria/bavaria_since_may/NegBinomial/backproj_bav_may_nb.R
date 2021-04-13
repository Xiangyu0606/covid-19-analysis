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


#' ## Backprojection of the epidemic curve
#'
#' Non-parametric back-projection as in Becker et al. (1991). The exposure time is
#' the relevant time scale to assess interventions.
#'
#' We take a literature based approach to deduce an incubation time distribution.
#'
#' Lauer et al. (2020) - log normal distribution - same as in Dehning et al. (2020)
#' Source: [Lauer et al. (2020)](https://www.ncbi.nlm.nih.gov/pubmed/32150748).
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




# Get necessary data (date and backprojected case numbers)
dat_bp_bav_back_may_nb = bpdf %>%
  select(date=epoch, backpro = upperbound) %>%
  filter(date>=lubridate::ymd("2020-05-01")) %>%mutate(backpro = round(backpro))

source("D:/Bavaria/before_may_nb/breakpoint_fun_nb.R")
# Estimate segmented regression with 3 breakpoints based on discrete optimization,
# can take several minutes, we load saved results

backproj_bav_may_3bp_nb = estimate_bp_disc_optim(data = dat_bp_bav_back_may_nb, bp = 3)
save(backproj_bav_may_3bp_nb, file = "D:/covid-19Analysis/bavaria_since_may/NegBinomial/backproj_bav_may_3bp_nb.RData")
load("D:/Bavaria/bavaria_since_may/NegBinomial/backproj_bav_may_3bp_nb.RData")

# Estimate segmented regression with 4 Breakpoints based on discrete optimization,
# can take several minutes, we load saved results
backproj_bav_may_4bp_nb = estimate_bp_disc_optim(data = dat_bp_bav_back_may_nb, bp = 4)
save(backproj_bav_may_4bp_nb, file = "D:/covid-19Analysis/bavaria_since_may/NegBinomial/backproj_bav_may_4bp_nb.RData")
load("D:/Bavaria/bavaria_since_may/NegBinomial/backproj_bav_may_4bp_nb.RData")

# Estimate segmented regression with 3 or 4 Breakpoints based on 'segmented' package with startpoints from
# discrete optimization

backproj_3_seg_may_nb = estimate_bp_segmented(
  data           = dat_bp_bav_back_may_nb,
  bp             = 3,
  start_bp       = backproj_bav_may_3bp_nb$bps,
  segmented_seed = 0520) 
backproj_4_seg_may_nb = estimate_bp_segmented(
  data           = dat_bp_bav_back_may_nb,
  bp             = 4,
  start_bp       = backproj_bav_may_4bp_nb$bps,
  segmented_seed = 0520)   

backproj_5_seg_may_nb = estimate_bp_segmented(
  data           = dat_bp_bav_back_may_nb,
  bp             = 5,
  start_bp       = c(30,60,90,100,123),
  segmented_seed = 0520)  
backproj_6_seg_may_nb = estimate_bp_segmented(
  data           = dat_bp_bav_back_may_nb,
  bp             = 6,
  start_bp       = c(21,26,65,89,99,123),
  segmented_seed = 0520) 

# Compare deviance and dispersion of models
c(
  dis_3bp = summary(backproj_3_seg_may_nb$segmented_model)$dispersion,
  dev_3bp = summary(backproj_3_seg_may_nb$segmented_model)$deviance,
  aic_3bp = summary(backproj_3_seg_may_nb$segmented_model)$aic,
  dis_4bp = summary(backproj_4_seg_may_nb$segmented_model)$dispersion,
  dev_4bp = summary(backproj_4_seg_may_nb$segmented_model)$deviance,
  aic_4bp = summary(backproj_4_seg_may_nb$segmented_model)$aic,
  dis_5bp = summary(backproj_5_seg_may_nb$segmented_model)$dispersion,
  dev_5bp = summary(backproj_5_seg_may_nb$segmented_model)$deviance,
  aic_5bp = summary(backproj_5_seg_may_nb$segmented_model)$aic,
  dis_6bp = summary(backproj_6_seg_may_nb$segmented_model)$dispersion,
  dev_6bp = summary(backproj_6_seg_may_nb$segmented_model)$deviance,
  aic_6bp = summary(backproj_6_seg_may_nb$segmented_model)$aic)

png("D:/Bavaria/bavaria_since_may/NegBinomial/backproj_bav_4BP_may_nb.png", width = 1200, height = 480)
backproj_4_seg_may_nb$plot
dev.off()

list(summary(backproj_6_seg_may_nb$segmented_model)$deviance,summary(backproj_6_seg_may_nb$segmented_model)$dispersion,
     summary(backproj_6_seg_may_nb$segmented_model)$aic)

knitr::kable(backproj_4_seg_may_nb$breakpoints, "html") %>%
  readr::write_file("D:/bavaria/bavaria_since_may/NegBinomial/disease-4change-points-may-nb.html")
