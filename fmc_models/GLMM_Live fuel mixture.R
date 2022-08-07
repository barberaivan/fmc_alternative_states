# GLMM structure

# location parameter model:
# mu ~ vpd +
#   site +   // random effect
#   ds +     // random effect, with sigma_ds varyng by site as fixed effects.
#   point_id // only for live fuel mixture

# dispersion parameter model: 
# phi ~ 
#   site +   // random effect
#   ds       // random effect, with sigma_ds constant across sites.


# Packages  ---------------------------------------------------------------

library(tidyverse)
library(rstan)
library(arm)
library(bayestestR) # hdi
library(mgcv) 
theme_set(theme_bw(base_family = "Ubuntu", base_size = 12))

# Functions ---------------------------------------------------------------

hdint <- function(x, ci = 0.95, name = "mu") {
  ci <- hdi(x, ci = ci)
  result <- c(ci$CI_low, ci$CI_high)
  names(result) <- paste(c(name, name), c("lower", "upper"), sep = "_")
  result
}

hdmean <- function(x, ci = 0.95, name = "mu") {
  ci <- hdi(x, ci = ci)
  result <- c(ci$CI_low, mean(x), ci$CI_high)
  names(result) <- paste(c(name, name, name), c("lower", "mean", "upper"), sep = "_")
  result
}

etimean <- function(x, ci = 0.95, name = "mu") {
  out <- 1 - ci
  q <- quantile(x, probs = c(out / 2, 1 - out / 2))
  result <- c(q[1], mean(x), q[2])
  names(result) <- paste(rep(name, 3), c("lower", "mean", "upper"), sep = "_")
  result
}

normalize <- function(x) x / sum(x)

# Subset fuel type --------------------------------------------------------

# fuel_type names (prior means)

### Dead
# Litter(50)
# 1 h fuel sticks(50)     
# 10 h fuel sticks(50)

### Live 
# Live fuel mixture(150) 
# Caña(100)
# Coihue (150)
# Laura (150)

f <- "Live fuel mixture"
prior_mean <- 150

# Data -------------------------------------------------------------------

setwd("/home/ivan/Insync/Humedad y EEA/fmc_alternative_states") # to run in raw R.

fdata <- readRDS("data_fmc_dataframe.R")
fdata <- fdata[fdata$fuel_type == f, ]

# Bring vapour pressure deficit data (vpd)
air_data <- readRDS("microclimate/predictions_vpd.R")
air_temp <- air_data$together
air_temp$site <- paste(air_temp$Community, air_temp$transect, sep = " / ")

# Order data
sites_present <- levels(fdata$site)[which(levels(fdata$site) %in% unique(fdata$site))]

fdata$site <- factor(fdata$site, levels = sites_present)
fdata <- fdata[order(fdata$site, fdata$time, fdata$point_id), ]
fdata$date_num <- factor(as.numeric(fdata$time), 
                         levels = as.character(unique(fdata$time)),
                         labels = as.character(1:length(unique(fdata$time))))

# datesite (ds)
fdata$ds <- paste(fdata$site, fdata$date_num, sep = " / ")
fdata$ds <- factor(fdata$ds, levels = unique(fdata$ds))

# point
fdata$point_id <- factor(fdata$point_id, levels = unique(fdata$point_id %>% as.numeric))

# table(fdata$site, fdata$site_num)

site_matrix <- model.matrix(moisture ~ site - 1, fdata)
date_matrix <- model.matrix(moisture ~ date_num - 1, fdata)
ds_matrix <- model.matrix(moisture ~ ds - 1, fdata)
point_matrix <- model.matrix(moisture ~ point_id - 1, fdata)

unique_id <- which(!duplicated(fdata[, c("site", "date_num", "ds")]))
unique_id_point <- which(!duplicated(fdata[, c("site", "date_num", "point_id")]))

# Merge fuel moisture and vpd data
fdata_site <- aggregate(moisture ~ transect + Community + site, fdata, mean)
fdata_site <- left_join(fdata_site, air_temp, by = c("transect", "Community", "site"))

# vpd sequence for prediction
N_pred <- 200
vpd_pred <- seq(min(fdata_site$vpd_lower), max(fdata_site$vpd_upper), 
                length.out = N_pred)

# scale vpd
vpd_scaled <- scale(fdata_site$vpd_max)  # "vpd_max" is observed mean by site
# vpd for prediction
vpd_pred_scaled <- (vpd_pred - attr(vpd_scaled, "scaled:center")) / attr(vpd_scaled, "scaled:scale")


# Prior predictive check --------------------------------------------------
N_sites <- 10
global_mean <- 150

prior_phi_int_mean = log(0.05)
prior_phi_int_sd = 1
prior_phi_sigma_sites_sd = 1
prior_phi_sigma_ds_sd = 1

# sample parameters
sigma_sites <- abs(rnorm(1, 0, prior_phi_sigma_sites_sd))
error_s <- rnorm(1, 0, sigma_sites)
sigma_ds <- abs(rnorm(1, 0, prior_phi_sigma_ds_sd))
error_d <- rnorm(1, 0, sigma_ds)

phi_intercept <- rnorm(1, prior_phi_int_mean, prior_phi_int_sd)

phi_sample <- exp(phi_intercept + error_s + error_d)

# curve(dgamma(x, shape = 1 / phi_sample, rate = 1 / (global_mean * phi_sample)),
#       from = 0, to = 500, #add = TRUE,
#       col = rgb(0, 0, 0, 1))
#       #ylim = c(0, NA))
# # Priors are OK. (very uninformative)

# Stan data ---------------------------------------------------------------

# data
sdata <- list(
  # Variables
  N = nrow(fdata),
  N_dates = length(unique(fdata$time)),
  N_sites = length(unique(fdata$site)),
  N_points = length(unique(fdata$point_id)),
  N_unique = length(unique_id),
  N_pred = N_pred,
  unique_id = unique_id,

  y = fdata$moisture,
  site_matrix = site_matrix,
  ds_matrix = ds_matrix,
  point_matrix = point_matrix,
  
  vpd = vpd_scaled %>% as.numeric,
  vpd_pred = vpd_pred_scaled,
  
  # Priors 
  # (log scale)
  
  # mu
  
  # fixed effects
  prior_mu_eta_mean = log(prior_mean),
  prior_mu_eta_sd = 1,
  prior_mu_lambda_sd = 1, # 2.5 if a spline is used, 1 for the linear term.
  
  # random effects
  prior_mu_sigma_sites_sd = 0.7,
  prior_mu_sigma_ds_sd = 0.7,
  prior_mu_sigma_points_sd = 0.7,    
  
  # phi
  
  prior_phi_eta_mean = log(0.05), # notice it's lower than for dead fuels, 
                                  # for this is a Gamma dispersion (1 / shape), 
                                  # while dead fuels have a truncated normal sigma.
  prior_phi_eta_sd = 1,
  
  # random effects
  prior_phi_sigma_sites_sd = 1,
  prior_phi_sigma_ds_sd = 1
)
str(sdata)

# Sampling settings
nc <- 10    # chains and cores
ns <- 1500  # samples 
nw <- 1500  # warm up iterations

N_samples <- nc * ns # to use later

# Model fit ----------------------------------------------------------------

# Compile, sample and save
# stan_code <- stan_model("fmc_models/GLMM_Live fuel mixture.stan", verbose = TRUE)
# glm1 <- sampling(
#   stan_code, data = sdata, seed = 564, refresh = 10,
#   chains = nc, cores = nc, iter = ns + nw, warmup = nw,
#   #chains = 1, cores = 1, iter = 10,
#   control = list(adapt_delta = 0.999, max_treedepth = 15),
#   pars = c("mu_sigma_sites", "mu_sigma_ds", "mu_sigma_points",
#            "mu_lambda",
#            "phi_sigma_sites", "phi_sigma_ds", # the last one is a scalar here, a vector in dead model.
#            "sites_phi", "pred_phi",
#            "phi",
#            "sites_mu", "pred_mu",
#            "mu",
#            "ds_means")
# )
# saveRDS(glm1, paste("fmc_models/model samples_", f, ".R", sep = ""))
# sglm1 <- summary(glm1)[[1]]
# saveRDS(sglm1, paste("fmc_models/model summary_", f, ".R", sep = ""))
# head(sglm1)
# min(sglm1[, "n_eff"]); max(sglm1[, "Rhat"])
# Done.

# 1: There were 5 divergent transitions after warmup. See
# https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
# to find out why this is a problem and how to eliminate them. 
# 2: Examine the pairs() plot to diagnose sampling problems

# Simplify: the variation among dates for phi will be constant between sites
# (phi_sigma_ds is a scalar now).

# Now it worked.
# 1290.21 / 60 = 21 min

# load model
glm1 <- readRDS(paste("fmc_models/model samples_", f, ".R", sep = ""))
sglm1 <- readRDS(paste("fmc_models/model summary_", f, ".R", sep = ""))

# Diagnosis with pairs plots ----------------------------------------------

# pairs(glm1, pars = c("mu_sigma_ds"))
# pairs(glm1, pars = c("mu_sigma_sites", "mu_sigma_points", "mu_lambda"))
# pairs(glm1, pars = c("sites_mu"))
# pairs(glm1, pars = c("phi_sigma_ds"))
# pairs(glm1, pars = c("phi_sigma_sites", "pred_phi"))
# pairs(glm1, pars = c("sites_phi")) # algunos valores muy extremos, y quizás tenga sense.

# Ningún patrón raro.
# t( apply(as.matrix(glm1, "sites_phi"), 2, quantile, 
#          probs = c(0.05, 0.98), method = 8) )

# Means by site and vpd  --------------------------------------------------

# compute site means
sites_mu <- as.matrix(glm1, pars = "sites_mu") %>% t
sites_phi <- as.matrix(glm1, pars = "sites_phi") %>% t

nsim <- ncol(sites_mu)

site_mean <- sites_mu # Gamma is parameterized directly with its mean, as
                      # opposed to truncated Normal.

site_summ <- apply(site_mean, 1, etimean, name = "fmc") %>% t
site_vpd_prediction <- cbind(fdata_site, site_summ)

# Get prediction of mean FMC as a function of vpd

vpd_prediction_mu <- as.matrix(glm1, pars = "pred_mu") %>% t
phi_mean <- as.matrix(glm1, pars = "pred_phi") %>% as.numeric

vpd_prediction_mean <- vpd_prediction_mu
# Gamma is parameterized directly with its mean, as
# opposed to truncated Normal.

vpd_prediction_summ <- apply(vpd_prediction_mean, 1, etimean, name = "fmc") %>% t %>% as.data.frame
vpd_temp <- data.frame(fuel_type = f, vpd = vpd_pred)
vpd_prediction <- cbind(vpd_temp, vpd_prediction_summ)

# plot to check

ggplot(site_vpd_prediction, 
       aes(x = vpd_max, #xmin = vpd_lower, xmax = vpd_upper,
           y = fmc_mean, ymin = fmc_lower, ymax = fmc_upper,
           color = Community)) +
  # sites estimates
  geom_point() +
  geom_errorbar() +
  #geom_errorbarh() +
  # prediction
  geom_line(data = vpd_prediction, 
            mapping = aes(x = vpd, y = fmc_mean),
            inherit.aes = FALSE) + 
  geom_ribbon(data = vpd_prediction, alpha = 0.3, 
              mapping = aes(x = vpd, y = fmc_mean, ymin = fmc_lower, 
                            ymax = fmc_upper),
              inherit.aes = FALSE) + 
  # observed data
  geom_point(mapping = aes(x = vpd_max, y = moisture), shape = 21,
             size = 2)


# Within transect differences -----------------------------------------------

# Compute difference as community_x - forest
diff_trans_raw <- do.call("rbind", lapply(levels(fdata_site$transect), function(tran) {
  #tran = "t1"
  coms <- unique(fdata_site[fdata_site$transect == tran, "Community"])
  coms <- coms[coms != "U forest"] %>% as.character
  
  # loop over communities inside transect
  one_transect <- do.call("rbind", lapply(coms, function(c) {
    #c = coms[1]
    filter1 <- which(fdata_site$transect == tran & fdata_site$Community == c)
    filter2 <- which(fdata_site$transect == tran & fdata_site$Community == "U forest")
    
    diff_post <- site_mean[filter1, ] - site_mean[filter2, ]
    
    diff_data <- data.frame(Community = c, transect = tran, diff = diff_post)
    return(diff_data)
  }))
  
  return(one_transect)
}))


# Function to compute summaries of distributions
# p is the probability of Y_{some community} - Y_{U forest} being <0
diff_p <- function(x) {
  summ <- etimean(x, name = "diff") 
  p <- length(x[x < 0]) / length(x)
  res <- c(summ, p = p)
  return(res)
}

# by site
diff_trans <- aggregate(diff ~ transect + Community, data = diff_trans_raw, 
                        FUN = diff_p)

# mean difference by community
diff_com_avg <- do.call("rbind", lapply(unique(diff_trans$Community), function(c) {
  #c = "B forest +reg"
  filter <- diff_trans_raw$Community == c
  diff_means <- apply(matrix(diff_trans_raw$diff[filter], 
                             ncol = length(unique(diff_trans_raw$transect[filter]))),
                      MARGIN = 1, mean) 
  summ <- diff_p(diff_means) %>% as.matrix() %>% t
  res <- cbind(data.frame(transect = "avg", Community = c), summ) %>% as.data.frame()
  return(res)
}))

# bring together both aggregations
diff_trans_2 <- do.call(data.frame, diff_trans)
names(diff_trans_2) <- c("transect", "Community", "diff_lower", "diff_mean", "diff_upper", "p")

(diff_table_avg <- rbind(diff_com_avg, diff_trans_2))

diff_table_avg$transect <- factor(diff_table_avg$transect, 
                                  levels = c("avg", unique(as.character(fdata$transect))))
diff_table_avg$Community <- factor(diff_table_avg$Community, 
                                   levels = levels(fdata$Community)[levels(fdata$Community) != "U forest"])


# diff_table_avg$p <- round(diff_table_avg$p, digits = 2)
ggplot(diff_table_avg, aes(x = Community, y = diff_mean, ymin = diff_lower,
                           ymax = diff_upper, colour = Community, 
                           shape = transect)) +
  geom_hline(yintercept = 0) + 
  geom_point(position = position_dodge(width = 0.8), size = 4) +
  geom_linerange(position = position_dodge(width = 0.8)) + 
  geom_text(mapping = aes(label = round(p, digits = 2), y = diff_lower), 
            position = position_dodge(width = 0.8),
            vjust = 1.5, color = "grey40")


# Data for appendix: FMC as a function of time (posterior check) ----------
# check coverage and fit

# Simulate from the fitted mu and scale parameters to get a predictive distribution
# and compute its HDI. At the same time, compute mean(x) with etruncnorm

# Create df without point_id. (ds for date-site)
fdata_ds <- fdata[unique_id, ]

# bring parameters by observation 
mu_hat <- as.matrix(glm1, par = "mu") %>% t
phi_hat <- as.matrix(glm1, par = "phi") %>% t

# get unique values
mu_uni <- mu_hat[unique_id, ]
phi_uni <- phi_hat[unique_id, ]

# Simulate and compute hdi by date-site

nsamp <- 50 # from each posterior sample and observation, take nsamp samples

y_hdi <- sapply(1:nrow(fdata_ds), function(row) {
  print(row)
  #row = 1
  filter <- which(fdata$time == fdata_ds$time[row] &
                  fdata$site == fdata_ds$site[row])
  
  if (length(filter) > 5 | length(filter) < 3) 
    stop(paste("Check error in filtering, row", row))
  
  mus <- as.numeric(mu_hat[filter, ])
  phis <- as.numeric(phi_hat[filter, ])
  
  y_hat <- rgamma(n = length(mus) * nsamp,
                  shape = 1 / phis,
                  rate = 1 / (phis * mus))
  
  return(hdint(y_hat, name = "y"))
}) %>% t %>% as.data.frame

# 2) Estimated mean FMC and CI by date-site

ds_means <- sglm1[grep("ds_means", rownames(sglm1)), 
                       c("2.5%", "mean", "97.5%")] %>% as.data.frame
rownames(ds_means) <- NULL
names(ds_means) <- c("mu_lower", "mu_mean", "mu_upper")

# merge
ds_predictions <- cbind(fdata_ds[, c("site", "transect", "Community", "time")],
                        y_hdi, ds_means)

# plot

ggplot(ds_predictions, aes(x = time, y = mu_mean, ymin = mu_lower, ymax = mu_upper)) +
  geom_line() +
  geom_ribbon(mapping = aes(ymin = y_lower, ymax = y_upper, y = mu_mean, x = time),
              color = NA, alpha = 0.20, inherit.aes = FALSE) + 
  geom_ribbon(color = NA, alpha = 0.40) +
  facet_grid(cols = vars(Community), rows = vars(transect)) +
  geom_point(fdata, mapping = aes(x = time, y = moisture), shape = 19,
             inherit.aes = FALSE, alpha = 0.7)

# Perfect


# P(ignition) mean by site ------------------------------------------------

# Bring data to compute ignition probability function

curve_data <- read.csv("data_Bianchi 2018 table 1_ignition freq and lfmc.csv")
curve_data$ig_prop[curve_data$ig_prop == 1] <- 1 - 1e-6
# remove non-native species from database
curve_data <- curve_data[curve_data$group == "native", ]

## simple model only to plot
# mtemp <- gam(ig_prop ~ lfmc, family = betar(), data = curve_data)
# nd <- data.frame(
#   lfmc = seq(min(curve_data$lfmc), max(curve_data$lfmc), length.out = 200)
# )
# nd$p <- predict(mtemp, nd, type = "response")
# plot(ig_prop ~ lfmc, data = curve_data, pch = 19)
# lines(p ~ lfmc, data = nd, col = "blue")

# Fit ignition probability model
ig_text <- "
data {
  int N;
  vector[N] y;
  vector[N] x;
}

parameters {
  real a; // intercept
  real b; // slope
  real<lower=0> phi; // dispersion parameter
}

transformed parameters {
  vector[N] mu = inv_logit(a + b * x);
  /*
  // In case the beta() function is used, then y ~ beta(alpha, beta)
  vector[N] alpha = mu * phi;
  vector[N] beta = (1 - mu) * phi;
  */
}

model {
  y ~ beta_proportion(mu, phi);
}
"
# compile
ig_code <- stan_model(model_code = ig_text, verbose = TRUE)

# sample
ig_model <- sampling(
  ig_code, data = list(N = nrow(curve_data), y = curve_data$ig_prop,
                       x = curve_data$lfmc), 
  seed = 564, refresh = 200,
  chains = nc, cores = nc, iter = ns + nw, warmup = nw,
  #chains = 1, cores = 1, iter = 20,
  control = list(adapt_delta = 0.9, max_treedepth = 15)
)

# get model parameters
a <- as.matrix(ig_model, "a")
b <- as.matrix(ig_model, "b")

# check:
# plot(ig_prop ~ lfmc, curve_data, ylim = c(0, 1))
# for(i in 1:2000) {
#   curve(plogis(a[i] + b[i] * x), add = TRUE, col = rgb(0, 0, 0, 0.01))
# }

# matrices to fill with summary of the posterior predictive distribution.
post_ig <- data.frame(matrix(NA, nrow(fdata_site), 4))
names(post_ig) <- paste("ig", c("lower", "mean", "upper", "obs"), sep = "_")

nsamp <- 10
N_samples <- ncol(mu_hat)
# Simulate predictions for observed sites, dates and sampling points.
# then, compute the mean.

for(i in 1:nrow(fdata_site)) {
  #i = 1
  print(i)
  
  filter <- which(as.character(fdata$transect) == as.character(fdata_site$transect[i]) & 
                    fdata$Community == fdata_site$Community[i])
  
  ig_means <- sapply(1:N_samples, function(s) {
    #s = 1
    mus <- as.numeric(mu_hat[filter, s])
    phis <- as.numeric(phi_hat[filter, s])
    
    # simulate fmc
    y_sim <- rgamma(n = length(mus) * nsamp, 
                    shape = 1 / phis, rate = 1 / (mus * phis))
    
    # compute mean ignition probability
    ig_mu <- plogis(a[s] + b[s] * y_sim)
    
    return(mean(ig_mu))
  })
  # summaryze ignition probability
  post_ig[i, 1:3] <- hdmean(ig_means, name = "ig")
  
  # the same for data, with posterior mean.
  moisture <- fdata$moisture[filter]
  # compute ignition probability distribution
  pmean <- sapply(1:N_samples, function(s) {
    mean(plogis(a[s] + b[s] * moisture))
  }) %>% mean
  
  post_ig[i, "ig_obs"] <- pmean
}

# merge predictions
ig_data <- cbind(site_vpd_prediction, post_ig)

# plot
alpha_bar <- 0.6
ggplot(ig_data, aes(x = fmc_mean, xmin = fmc_lower, xmax = fmc_upper,
                    y = ig_mean, ymin = ig_lower, ymax = ig_upper,
                    colour = Community, shape = Community)) +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  geom_point(size = 3, alpha = 1) + 
  geom_errorbar(alpha = alpha_bar) +
  geom_errorbarh(alpha = alpha_bar) +
  #ylim(0, 1) +
  geom_point(mapping = aes(y = ig_obs, x = moisture, fill = Community)) + 
  facet_wrap(vars(transect)) + 
  theme(panel.grid.minor = element_blank()) +
  ggtitle("live fuel mixture")

p_ig_data <- ig_data
p_ig_data$fuel_type <- f

# Number of days below thresholds (only observed) ------------------------

N_dates <- length(unique(fdata$date))
fdata_dates <- aggregate(moisture ~ date + time + transect + Community,
                         fdata, mean)
ddates_wide <- pivot_wider(fdata_dates, names_from = "Community", values_from = "moisture")
coms <- levels(fdata$Community)
# View(ddates_wide)

m <- as.matrix(ddates_wide[, names(ddates_wide) %in% coms])
th <- 110 
below <- (m < th) * 1

# this is part of the final table
percents_below <- aggregate(below ~ transect, ddates_wide, 
                            FUN = function(x) sum(x) / N_dates * 100,
                            na.action = "na.pass")

write.csv(percents_below, "fmc_models/table_live_thresholds.csv")

# Export ------------------------------------------------------------------

diff_table_avg$fuel_type <- f
site_vpd_prediction$fuel_type <- f
vpd_prediction$fuel_type <- f
ds_predictions$fuel_type <- f
p_ig_data$fuel_type <- f

export_list <- list(
  diff_table = diff_table_avg,
  
  # fmc ~ vpd by site
  site_vpd_prediction = site_vpd_prediction,
  vpd_prediction = vpd_prediction, # line and ribbon
  
  # for posterior predictive check (means as a function of time):
  fdata = fdata,
  ds_predictions = ds_predictions,
  
  # ignition probability predictions
  p_ig = p_ig_data,
  
  notes = paste("Created in \"", rstudioapi::getSourceEditorContext()$path, 
                "\". \nDate: ", Sys.Date(), sep = "")
)

saveRDS(export_list, 
        paste("fmc_models/model_predictions_r_object_", f, ".R", sep = ""))


# Ckeck fit to means by date-site and site ----------------------------------

# date-site means
fdata_ds_agg <- aggregate(moisture ~ time + site + Community + transect, fdata, mean)

means_check <- left_join(fdata_ds_agg, ds_predictions, 
                         by = c("time", "transect", "Community", "site"))

# mean fmc as a function of time
ggplot(means_check, aes(x = time, y = mu_mean, ymin = mu_lower, ymax = mu_upper)) +
  geom_line() +
  geom_ribbon(color = NA, alpha = 0.4) + 
  geom_point(mapping = aes(x = time, y = moisture), color = "red") +
  facet_grid(cols = vars(Community), rows = vars(transect)) 

# Means by site
sites_compare <- cbind(site_summ, fdata_site)

ggplot(sites_compare, aes(x = Community, y = fmc_mean, ymin = fmc_lower, 
                          ymax = fmc_upper)) +
  geom_point(size = 3) +
  geom_errorbar() + 
  geom_point(data = sites_compare, mapping = aes(x = Community, y = moisture), 
             color = "red", alpha = 1, size = 3) +
  facet_wrap(vars(transect))

# Perfect.