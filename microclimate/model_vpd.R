# Check cholesky parameterization:

# cormat <- matrix(c(1, 0.5, 0.3,
#                    0.5, 1, 0.2,
#                    0.3, 0.2, 1),
#                  nrow = 3, byrow = T)
# sdd <- 2
# 
# cor_chol <- t(chol(cormat))
# 
# vvv <- sdd ^ 2 * cormat
# vvv_chol <- t(chol(vvv))
# 
# vvv_chol == cor_chol * sdd
# The left hand side will be used in stan

# Check mean computation for lognormal distribution with random effects.

# mu <- 1 # mean for 1 community at log scale
# mu_sd <- 0.5 # site random effects sd
# sigma <- 1 # unexplained variability at log scale
# mean1 <- exp(mu + 0.5 * (sigma ^ 2 + mu_sd ^ 2)) # the mean for each community is 
# # computed this way in stan
# 
# # verify with simulation
# nsim <- 1e7 # simulated sites
# sitt <- rnorm(nsim, mu, mu_sd)
# mean2 <- mean(exp(sitt + 0.5 * sigma ^ 2))
# 
# mean1; mean2 # perfect.


# Packages ----------------------------------------------------------------

# setwd("/home/ivan/Insync/Humedad y EEA/fmc_alternative_states") # necessary to run from raw R.

library(rstan)
library(tidyverse); theme_set(theme_bw(base_family = "Ubuntu"))
library(bayestestR)


# Functions ---------------------------------------------------------------

hdmed <- function(x, ci = 0.95, name = "mu") {
  ci <- hdi(x, ci = ci)
  result <- c(ci$CI_low, median(x), ci$CI_high)
  names(result) <- paste(rep(name, 3), c("lower", "median", "upper"), sep = "_")
  result
}

hdmean <- function(x, ci = 0.95, name = "mu") {
  ci <- hdi(x, ci = ci)
  result <- c(ci$CI_low, mean(x), ci$CI_high)
  names(result) <- paste(rep(name, 3), c("lower", "mean", "upper"), sep = "_")
  result
}


etimean <- function(x, ci = 0.95, name = "mu") {
  out <- 1 - ci
  q <- quantile(x, probs = c(out / 2, 1 - out / 2))
  result <- c(q[1], mean(x), q[2])
  names(result) <- paste(rep(name, 3), c("lower", "mean", "upper"), sep = "_")
  result
}

# Prepare data ------------------------------------------------------------

fdata_air <- readRDS("data_fmc_dataframe.R")

# remove t4, which did not have a datalogger
fdata_air <- fdata_air[fdata_air$transect != "t4", ]

# recreate factors levels 
fdata_air$transect <- factor(fdata_air$transect, 
                             levels = levels(fdata_air$transect)[levels(fdata_air$transect) %in% unique(fdata_air$transect)])
fdata_air$site <- factor(fdata_air$site, 
                         levels = levels(fdata_air$site)[levels(fdata_air$site) %in% unique(fdata_air$site)])
fdata_air$site_num <- as.numeric(fdata_air$site)

fdata_air$com_num <- as.numeric(fdata_air$Community)

# load air data
air <- readRDS("microclimate/data_microclimate_daily_dataframe.R")

# air por comunidad y transecta (raw)
air$Community <- factor(air$Community, levels = levels(fdata_air$Community))
air$transect <- factor(air$transect, levels = levels(fdata_air$transect))
air$site <- factor(air$site, levels = levels(fdata_air$site))
air$site_num <- as.numeric(air$site)
air$com_num <- as.numeric(air$Community)
with(air, table(com_num, Community))

# order by site and time
air <- air[air$fDate <= "2020-03-15", ]
air <- air[order(air$site_num, air$time), ]

# Data for Stan ------------------------------------------------------------

enes <- aggregate(vpd_max ~ site_num, air, length)
vpd_max <- matrix(air$vpd_max, nrow = enes[1, 2]) %>% log
colnames(vpd_max) <- enes$site_num

# comunity id for every column (site) in vpd_max matrix:
ttt <- table(air$com_num, air$site_num)
com <- sapply(1:ncol(ttt), function(x) which(ttt[, x] > 0)) %>% unname

d <- as.matrix(dist(air$time[air$site_num == 1]), upper = TRUE, diag = TRUE)
#View(d)

# prior predictive checks
# tau_mu <- 1.5
# meanlog <- log(2)
# tau_sample <- abs(rnorm(1, 0, tau_mu))
# curve(dlnorm(x, meanlog, tau_sample), from = 0, to = 10, ylim = c(0, 3))
# for(i in 1:500) {
#   tau_mu <- 1
#   tau_sample <- abs(rnorm(1, 0, tau_mu))
#   curve(dlnorm(x, meanlog, tau_sample), from = 0, to = 10, add = TRUE, 
#         col = rgb(0, 0, 0, 0.1))
# }
# 
# tau_sd <- 1
# tau_mean <- -1
# vpdml <- log(2)
# tau_sample <- exp(rnorm(1, tau_mean, tau_sd))
# curve(dlnorm(x, vpdml, tau_sample), from = 0, to = 6, ylim = c(0, 3))
# for(i in 1:500) {
#   tau_sample <- exp(rnorm(1, tau_mean, tau_sd))
#   curve(dlnorm(x, vpdml, tau_sample), add = TRUE, 
#         col = rgb(0, 0, 0, 0.1))
# }

sdata_vpd <- list(
  y = vpd_max,
  N = nrow(vpd_max),
  N_sites = ncol(vpd_max),
  N_com = length(unique(air$Community)),
  com = com,
  d = d,
  
  prior_eta_sd = 5,
  prior_sigma_sd = 1.5,
  prior_sigma_delta_sd = 1
  # # priors for {mu, sigma, rho}, at
  # # {identity, log, logit} scales.
  # # parameters affect the response at the log scale.
  # prior_eta_mean = c(0, -1, 0), # fixed effects prior means (normal)
  # prior_eta_sd = c(5, 1, 3),    # fixed effects prior sds (normal)
  # prior_tau_sd = c(1, 1, 1.5)   # random effects sd prior sds (half-normals)
)

# Model fit ---------------------------------------------------------------

# code <- stan_model("microclimate/model_vpd.stan", verbose = TRUE)
# m_vpd <- sampling(
#   code, data = sdata_vpd, seed = 123,
#   cores = 10, chains = 10, iter = 1500, warmup = 500, refresh = 10,
#   #cores = 1, chains = 1, iter = 5, refresh = 1,
#   control = list(adapt_delta = 0.99, max_treedepth = 15),
#   pars = c("eta", "delta", "sigma", "sigma_delta", "rho",
#            "mu_sites", "mu_exp", "mu_com")
# )
# saveRDS(m_vpd, "microclimate/model_vpd_samples.R")
# 1850.05 / 60 = 30 min

m_vpd <- readRDS("microclimate/model_vpd_samples.R")
sm_vpd <- summary(m_vpd)[[1]]
# sm_vpd2 <- summary(m_vpd2)[[1]]
#View(sm_vpd)
summary(sm_vpd) 
# max rhat 1.003
# min neff 2408

pairs(m_vpd, pars = c("rho", "sigma", "sigma_delta"))

# Check fit to the means

y_means <- aggregate(vpd_max ~ site_num, air, mean)
mu_exp <- as.matrix(m_vpd, pars = "mu_exp")

mu_exp_long <- matrix(mu_exp, ncol = 1)
mu_exp_df <- data.frame(
  site_num = rep(1:10, each = nrow(mu_exp)),
  mu_exp = mu_exp_long
)

mu_hat_means <- aggregate(mu_exp ~ site_num, mu_exp_df, mean)

ggplot(mu_exp_df, aes(x = mu_exp)) + 
  geom_density() + 
  facet_wrap(vars(site_num), scales = "free_y") +
  geom_vline(data = y_means, mapping = aes(xintercept = vpd_max), alpha = 0.7) +
  geom_vline(data = mu_hat_means, mapping = aes(xintercept = mu_exp), 
             colour = "red", alpha = 0.7) #+
  #xlim(0, 6)

plot(mu_hat_means$mu_exp ~ y_means$vpd_max); abline(0, 1)
# good 

# Predictions to export-----------------------------------------------------

# MEANS # 

# Means by site
air_agg <- aggregate(vpd_max ~ transect + Community + site + site_num, air, mean)
estimates <- apply(mu_exp, 2, etimean, name = "vpd") %>% t %>% as.data.frame
rownames(estimates) <- NULL
export_site <- cbind(air_agg, estimates)

# Means by community
com_agg <- aggregate(vpd_max ~ Community, air_agg, mean)
mu_com <- as.matrix(m_vpd, pars = "mu_com")
com_estimates <- apply(mu_com, 2, etimean, name = "vpd") %>% t %>% as.data.frame
rownames(com_estimates) <- NULL
export_com <- cbind(com_agg, com_estimates)
export_com$transect <- "average"
names(export_com) <- c("Community", "vpd_max", "vpd_lower", "vpd_mean", "vpd_upper", "transect")

avetran <- rbind(
  export_site[, c("transect", "Community", "vpd_max", "vpd_lower", "vpd_mean", "vpd_upper")],
  export_com[, c("transect", "Community", "vpd_max", "vpd_lower", "vpd_mean", "vpd_upper")]
)

avetran$type <- factor(rep(c("transects", "average"), 
                           c(nrow(export_site), nrow(export_com))),
                        levels = c("average", "transects"))
avetran$transect <- factor(avetran$transect, levels = c("average", "t1", "t2", "t3"))

# plot:
ggplot(avetran, aes(x = Community, y = vpd_mean, colour = transect, 
                    shape = transect)) + 
  geom_pointrange(position = position_dodge(width = 0.7),
                  mapping = aes(ymin = vpd_lower, ymax = vpd_upper), 
                  fatten = 4, size = 0.9) +
  scale_color_manual(values = c("black", "grey50", "grey50", "grey50")) +
  scale_shape_manual(values = c(19, 15, 17, 18)) + 
  theme(legend.title = element_blank()) +
  ylab("Vapour pressure deficit (Mpa)") + 
  ggtitle("Average vapour pressure deficit by community",
          subtitle = "Posterior medians and 95 % credible intervals")

# CONTRASTS (differences community_i - forest) # 

# Function to compute summaries of distributions
# p is the probability of Y_{some community} - Y_{U forest} being <0
diff_p <- function(x) {
  summ <- etimean(x, name = "diff") 
  p <- length(x[x > 0]) / length(x)
  res <- c(summ, p = p)
  return(res)
}

# Sites differences

site_samples <- mu_exp %>% t

# Compute and summarize difference as community_x - forest
diff_sites <- do.call("rbind", lapply(levels(air_agg$transect), function(tran) {
  #tran = "t1"
  coms <- unique(air_agg[air_agg$transect == tran, "Community"])
  coms <- coms[coms != "U forest"] %>% as.character
  
  # loop over communities inside transect
  one_transect <- do.call("rbind", lapply(coms, function(c) {
    #c = coms[1]
    filter1 <- which(air_agg$transect == tran & air_agg$Community == c)
    filter2 <- which(air_agg$transect == tran & air_agg$Community == "U forest")
    
    diff_post <- site_samples[filter1, ] - site_samples[filter2, ]
    d <- diff_p(diff_post)
    
    diff_data <- data.frame(Community = c, transect = tran, 
                            diff_lower = d["diff_lower"],
                            diff_mean = d["diff_mean"],
                            diff_upper = d["diff_upper"],
                            p = d["p"])
    rownames(diff_data) <- NULL
    return(diff_data)
  }))
  
  return(one_transect)
}))

# Difference by community
diff_com <- diff_sites[1:3, ]
diff_com$transect <- "average"

for(i in 1:nrow(diff_com)) {
  delta <- mu_com[, i+1] - mu_com[, 1]
  dsumm <- diff_p(delta)
  diff_com[i, names(dsumm)] <- dsumm
}

# merge sites and communities differences
diff_both <- rbind(diff_sites, diff_com)
diff_both$transect <- factor(diff_both$transect, 
                             levels = c("average", unique(as.character(air$transect))))
diff_both$Community <- factor(diff_both$Community, 
                              levels = levels(air$Community)[levels(air$Community) != "U forest"])

# plot
ggplot(diff_both, aes(x = Community, y = diff_mean, ymin = diff_lower,
                           ymax = diff_upper, colour = Community, 
                           shape = transect)) +
  geom_hline(yintercept = 0) + 
  geom_point(position = position_dodge(width = 0.8), size = 4) +
  geom_linerange(position = position_dodge(width = 0.8)) + 
  geom_text(mapping = aes(label = round(p, digits = 2), y = diff_lower), 
            position = position_dodge(width = 0.8),
            vjust = 1.5, color = "grey40")


# Export ------------------------------------------------------------------

export_all <- list(
  # means
  by_site = export_site,
  by_com = export_com,
  together = avetran,
  # differences
  diff = diff_both
)

saveRDS(export_all, "microclimate/predictions_vpd.R")

