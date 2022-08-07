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
library(truncnorm) 
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

### Others (don't fit in this script)
# Live fuel mixture(150) 
# Caña(100)
# Coihue (150)
# Laura (150)

# For dead fuels, I run this script 3 times. 
# I don't loop because each model has to be checked and tuned 
# (stan parameters as max_treedepth and adapt_delta) carefully.

f <- "10 h fuel sticks"#"10 h fuel sticks"#
prior_mean <- 50
K <- 2 # basis functions for spline

# Data -------------------------------------------------------------------

setwd("/home/ivan/Insync/Humedad y EEA/fmc_alternative_states") # to run in raw R.

fdata <- readRDS("data_fmc_dataframe.R")
fdata <- fdata[fdata$fuel_type == f, ]
if(f == "Litter" & alternative) {
  fdata <- fdata[fdata$time != min(fdata$time), ]
}
#nrow(fdata) # 517 # or 468 (alternative for litter)

# Bring vapour pressure deficit data (vpd)
air_data <- readRDS("microclimate/predictions_vpd.R")
air_temp <- air_data$together
air_temp$site <- paste(air_temp$Community, air_temp$transect, sep = " / ")

# filter sites
sites_filter <- which(levels(fdata$site) %in% unique(fdata$site) &
                        levels(fdata$site) %in% unique(air_temp$site))

sites_present <- levels(fdata$site)[sites_filter]

fdata$site <- factor(fdata$site, levels = sites_present)

# remove rows from unused sites
fdata <- fdata[!is.na(fdata$site), ]

fdata <- fdata[order(fdata$site, fdata$time), ]
fdata$date_num <- factor(as.numeric(fdata$time), 
                         levels = as.character(unique(fdata$time)),
                         labels = as.character(1:length(unique(fdata$time))))

# datesite (ds)
fdata$ds <- paste(fdata$site, fdata$date_num, sep = " / ")
fdata$ds <- factor(fdata$ds, levels = unique(fdata$ds))

# table(fdata$site, fdata$site_num)

site_matrix <- model.matrix(moisture ~ site - 1, fdata)
#date_matrix <- model.matrix(moisture ~ date_num - 1, fdata)
ds_matrix <- model.matrix(moisture ~ ds - 1, fdata)

unique_id <- which(!duplicated(fdata[, c("site", "date_num", "ds")]))

# Merge fuel moisture and vpd data
fdata_site <- aggregate(moisture ~ transect + Community + site, fdata, mean)
fdata_site <- left_join(fdata_site, air_temp, by = c("transect", "Community", "site"))

# vpd sequence for prediction
N_pred <- 200
vpd_pred <- seq(min(fdata_site$vpd_lower), max(fdata_site$vpd_upper), 
                length.out = N_pred)

# scale vpd
# vpd_scaled <- scale(fdata_site$vpd_max) # "vpd_max" is observed mean by site
# vpd_scaled <- scale(fdata_site$vpd_mean) # posterior mean as point estimate# 
# data for prediction
# vpd_pred_scaled <- (vpd_pred - attr(vpd_scaled, "scaled:center")) / attr(vpd_scaled, "scaled:scale")

# Basis for spline of vpd

basis <- smoothCon(
  s(vpd_max, bs = "tp", k = K + 1), absorb.cons = TRUE, # vpd_max is observed mean by site (mean of the max)
  diagonal.penalty = T, data = fdata_site
)[[1]]

basis_fit <- basis$X
basis_pred <- PredictMat(basis, data.frame(vpd_max = vpd_pred))


# Prior predictive check --------------------------------------------------
N_sites <- 10
global_mean <- 10

prior_phi_int_mean = log(15)
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

# curve(dtruncnorm(x, 0, Inf, mean = global_mean, sd = phi_sample),
#       from = 0, to = 100, #add = TRUE,
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
  N_unique = length(unique_id),
  N_pred = N_pred,
  unique_id = unique_id,
  K = K, # basis functions for spline (without intercept)
  
  y = fdata$moisture,
  site_matrix = site_matrix,
  ds_matrix = ds_matrix,
  #vpd = vpd_scaled %>% as.numeric,
  #vpd_pred = vpd_pred_scaled,
  basis_fit = basis_fit,
  basis_pred = basis_pred,
  
  # Priors 
  # (log scale)
  
  # mu
  
  # fixed effects
  prior_mu_eta_mean = log(prior_mean),
  prior_mu_eta_sd = 1,
  prior_mu_lambda_sd = 2.5,
  
  # random effects
  prior_mu_sigma_sites_sd = 0.7,
  prior_mu_sigma_ds_sd = 0.7,
  # prior_mu_gamma_ds_sd = 0.7, 
  # prior_mu_tau_ds_sd = 0.7,    
   
  # phi
  
  # fixed effects
  prior_phi_eta_mean = log(15), # notice it's higher than for live fuels, 
                                # for this is a normal sd and live fuels have
                                # an inverse shape parameter for Gamma distribution.
  prior_phi_eta_sd = 1,
  
  # random effects
  prior_phi_sigma_sites_sd = 1,
  prior_phi_sigma_ds_sd = 1
  # prior_phi_gamma_ds_sd = 1, 
  # prior_phi_tau_ds_sd = 1    
)
# str(sdata)

# Sampling settings
nc <- 10    # chains and cores
ns <- 1000  # samples 
nw <- 1000  # warm up iterations

if(f == "1 h fuel sticks") {  # more samples needed to reach min(N_eff) ~ 1000
  nc <- 10    # chains and cores
  ns <- 2500  # samples
  nw <- 1000  # warm up iterations
}

N_samples <- nc * ns # to use later

# Model fit ----------------------------------------------------------------

# Compile, sample and save
stan_code <- stan_model("fmc_models/GLMM_Sticks.stan", verbose = TRUE)
glm1 <- sampling(
  stan_code, data = sdata, seed = 564, refresh = 10,
  chains = nc, cores = nc, iter = ns + nw, warmup = nw,
  #chains = 1, cores = 1, iter = 10,
  control = list(adapt_delta = 0.999, max_treedepth = 15),
  pars = c("mu_sigma_sites", "mu_sigma_ds",
           "mu_lambda",
           "phi_sigma_sites", "phi_sigma_ds",
           "sites_phi", "pred_phi",
           "phi",
           "sites_mu", "pred_mu",
           "mu")
) # mu are mu for truncnorm. They are not means due to truncation.
saveRDS(glm1, paste("fmc_models/model samples_", f, ".R", sep = ""))
sglm1 <- summary(glm1)[[1]]
saveRDS(sglm1, paste("fmc_models/model summary_", f, ".R", sep = ""))
head(sglm1)
min(sglm1[, "n_eff"]); max(sglm1[, "Rhat"])
# Done.

# 1 h stick:
# 1456.3 / 60 = 24 min;  # with 1000 iter min neff is 400. run more iters
#plot(sglm1[, "n_eff"], ylim = c(0, 1000))

# 10 h stick:
# 1142 / 60 =  19 min; all good

# load model
glm1 <- readRDS(paste("fmc_models/model samples_", f, ".R", sep = ""))
sglm1 <- readRDS(paste("fmc_models/model summary_", f, ".R", sep = ""))



# Diagnosis with pairs plots ----------------------------------------------

#   pars = c("mu_sigma_sites", "mu_sigma_ds", 
#            "mu_lambda",
#            "phi_sigma_sites", "phi_sigma_ds", 
#            "sites_phi", "pred_phi",
#            "phi",
#            "sites_mu", "pred_mu",
#            "mu")

# pairs(glm1, pars = c("mu_sigma_ds"))
# pairs(glm1, pars = c("mu_sigma_sites", "mu_lambda"))
# pairs(glm1, pars = c("sites_mu"))
# pairs(glm1, pars = c("phi_sigma_ds"))
# pairs(glm1, pars = c("phi_sigma_sites", "pred_phi"))
# pairs(glm1, pars = c("sites_phi")) # algunos valores muy extremos, y quizás tenga sense.

# Ningún patrón raro.
# t( apply(as.matrix(glm1, "sites_phi"), 2, quantile, 
#          probs = c(0.05, 0.98), method = 8) )

# Means by site and vpd  ----------------------------------------------------

# compute site means
sites_mu <- as.matrix(glm1, pars = "sites_mu") %>% t
sites_phi <- as.matrix(glm1, pars = "sites_phi") %>% t

nsim <- ncol(sites_mu)

site_mean <- sapply(1:nsim, function(i) {
  etruncnorm(mean = sites_mu[, i], sd = sites_phi[, i], a = 0)
})

site_summ <- apply(site_mean, 1, etimean, name = "fmc") %>% t
site_vpd_prediction <- cbind(fdata_site, site_summ)

# Get prediction of mean FMC as a function of vpd

vpd_prediction_mu <- as.matrix(glm1, pars = "pred_mu") %>% t
#str(vpd_prediction_mu)
phi_mean <- as.matrix(glm1, pars = "pred_phi") %>% as.numeric

vpd_prediction_mean <- sapply(1:nsim, function(i) {
  etruncnorm(a = 0, mean = vpd_prediction_mu[, i], sd = phi_mean[i])
})
str(vpd_prediction_mean)

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

# Notes for litter:
# It's good. The uncertain means are coherent with the variability that forests
# show. Let's see the posteriors:

# for(i in 1:N_sites) {
#   plot(density(site_mean[i, ], from = 0, to = 100), 
#        main = site_vpd_prediction$site[i],
#        xlim = c(0, 40))
# } # OK, very assymetric ones


# Within transect differences --------------------------------------------

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

# Sample from Truncnorm
full_summary <- sapply(1:nrow(fdata_ds), function(row) {
  print(row)
  
  mus <- as.numeric(mu_uni[row, ])
  phis <- as.numeric(phi_uni[row, ])
    
  y_hat <- rtruncnorm(n = length(mus) * nsamp, 
                      a = 0, mean = mus, sd = phis)
  
  y_ints <- hdint(y_hat, name = "y")
  mean_summ <- etimean(etruncnorm(a = 0, mean = mus, sd = phis), name = "mu")
  
  return(c(mean_summ, y_ints))
}) %>% t %>% as.data.frame

# merge
ds_predictions <- cbind(fdata_ds[, c("site", "transect", "Community", "time")],
                        full_summary)

# plot

ggplot(ds_predictions, aes(x = time, y = mu_mean, ymin = mu_lower, ymax = mu_upper)) +
  geom_line() +
  geom_ribbon(mapping = aes(ymin = y_lower, ymax = y_upper, y = mu_mean, x = time),
              color = NA, alpha = 0.20, inherit.aes = FALSE) + 
  geom_ribbon(color = NA, alpha = 0.40) +
  facet_grid(cols = vars(Community), rows = vars(transect)) +
  geom_point(fdata, mapping = aes(x = time, y = moisture), shape = 19,
             inherit.aes = FALSE, alpha = 0.7)

# Le emboca perfecto a los datos este último modelo. 


# Export ------------------------------------------------------------------

diff_table_avg$fuel_type <- f
site_vpd_prediction$fuel_type <- f
vpd_prediction$fuel_type <- f
ds_predictions$fuel_type <- f

export_list <- list(
  diff_table = diff_table_avg,
  
  # fmc ~ vpd by site
  site_vpd_prediction = site_vpd_prediction,
  vpd_prediction = vpd_prediction, # line and ribbon
  
  # for posterior predictive check (means as a function of time):
  fdata = fdata,
  ds_predictions = ds_predictions,
  
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
