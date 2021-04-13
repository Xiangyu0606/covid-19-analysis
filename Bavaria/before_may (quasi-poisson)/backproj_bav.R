
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



# Extract data from nowcast
sts_symp <- sts(
  epoch       = covid_ts$Date,
  observed    = matrix(covid_ts$Est_Onsets, ncol = 1, nrow = nrow(covid_ts)),
  epochAsDate = TRUE)

# Perform back projection with smoothing to adjust for weekday effects (k=6)
bp <- backprojNP(sts_symp, incu.pmf=pmf, control=list(k=6, eq3a.method="C"))
# Reduce to relevant subset of the data making it comparable with Dehning et al. (2020)
bp <- bp[epoch(bp) <= as.Date("2020-05-01") & epoch(bp) >= as.Date("2020-02-27"),]

bpdf <- bp %>%
  as.data.frame() %>%
  mutate(epoch_numeric=as.numeric(epoch), t = epoch - max(epoch) + 1)

source("D:/Bavaria/before_may (quasi-poisson)/breakpoint_fun.R")
# Get necessary data (date and backprojected case numbers)
dat_bp = bpdf %>%
  select(date=epoch, backpro = upperbound) %>%
  filter(date>=lubridate::ymd("2020-02-27"))
# Estimate segmented regression with 3 breakpoints based on discrete optimization,
# can take several minutes, we load saved results
res_3_disc = estimate_bp_disc_optim(data = dat_bp, bp = 3)
save(res_3_disc, file = "D:/Bavaria/before_may (quasi-poisson)/results/res_disc_opt_3bp.RData")
load("D:/Bavaria/before_may (quasi-poisson)/results/res_disc_opt_3bp.RData")

# Estimate segmented regression with 4 Breakpoints based on discrete optimization,
# can take several minutes, we load saved results
res_4_disc = estimate_bp_disc_optim(data = dat_bp, bp = 4)

save(res_4_disc, file = "D:/covid-19Analysis/bavarian/results/res_disc_opt_4bp.RData")
load("D:/Bavaria/before_may (quasi-poisson)/results/res_disc_opt_4bp.RData")

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
  start_bp       = res_4_disc$bps,
  segmented_seed = 0520)

# Compare deviance and dispersion of models

round(
  c(
    dev_3_bp      = summary(res_3_seg$segmented_model)$deviance,
    overdisp_3_bp = summary(res_3_seg$segmented_model)$dispersion,
    dev_4_bp      = summary(res_4_seg$segmented_model)$deviance,
    overdisp_4_bp = summary(res_4_seg$segmented_model)$dispersion),
  2)


############################## Results #########################################
# Plot of estimated segmented regression model
### Use model with 3 change points
png("D:/Bavaria/before_may (quasi-poisson)/results/backprojections-3bp.png", width = 1000, height = 460)
res_3_seg$plot
dev.off()
png("D:/Bavaria/before_may (quasi-poisson)/results/backprojections-4bp.png", width = 1000, height = 460)
res_4_seg$plot
dev.off()

# More info on estimated breakpoints
# Breakpoints + 95%-CIs
knitr::kable(res_3_seg$breakpoints, "html") %>%
  readr::write_file("D:/Bavaria/before_may (quasi-poisson)/results/back-projection-change-points.html")
# Daily multiplicative change in infection numbers
knitr::kable(res_3_seg$coef, "html", digits = 3) %>%
  readr::write_file("D:/Bavaria/before_may (quasi-poisson)/results/back-projection-factors.html")
