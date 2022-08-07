# Informative prior for sd_sites in cania, based on posterior for coihue and laura.
library(tidyverse)
library(viridis)
setwd("/home/ivan/Insync/Humedad y EEA/fmc_alternative_states")

mcoi <- readRDS("fmc_models/model samples_Coihue.R")
mlau <- readRDS("fmc_models/model samples_Laura.R")

sigma_mu <- data.frame(
  coi = as.matrix(mcoi, "mu_sigma_sites") %>% as.numeric,
  lau = as.matrix(mlau, "mu_sigma_sites") %>% as.numeric
)

sigma_phi <- data.frame(
  coi = as.matrix(mcoi, "phi_sigma_sites") %>% as.numeric,
  lau = as.matrix(mlau, "phi_sigma_sites") %>% as.numeric
)

# mu
plot(density(sigma_mu$coi, from = 0), col = "blue", type = "l",
     ylim = c(0, 25))
lines(density(sigma_mu$lau, from = 0), col = "red")
curve(dnorm(x, 0, 0.1) * 2, from = 0, add = TRUE, col = "black",
      lwd = 2)

# phi
plot(density(sigma_phi$coi, from = 0), col = "blue", type = "l",
     ylim = c(0, 5))
lines(density(sigma_phi$lau, from = 0), col = "red")
curve(dnorm(x, 0, 0.5) * 2, from = 0, add = TRUE, col = "black",
      lwd = 2)

# For cania, use 
# mu_sigma_sites ~ normal(0, 0.1)T[0, ]
# phi_sigma_sites ~ normal(0, 0.5)T[0, ]


# plot --------------------------------------------------------------------

mcan <- readRDS("fmc_models/model samples_Cania.R")

sigma_mu$can <- as.matrix(mcan, "mu_sigma_sites") %>% as.numeric
sigma_phi$can <- as.matrix(mcan, "phi_sigma_sites") %>% as.numeric

# Repeat plots including cania

# mu
plot(density(sigma_mu$coi, from = 0), col = "blue", type = "l",
     ylim = c(0, 25))
lines(density(sigma_mu$lau, from = 0), col = "red")
lines(density(sigma_mu$can, from = 0), col = "green")
curve(dnorm(x, 0, 0.1) * 2, from = 0, add = TRUE, col = "black",
      lwd = 2)

# phi
plot(density(sigma_phi$coi, from = 0), col = "blue", type = "l",
     ylim = c(0, 5))
lines(density(sigma_phi$lau, from = 0), col = "red")
lines(density(sigma_phi$can, from = 0), col = "green")
curve(dnorm(x, 0, 0.5) * 2, from = 0, add = TRUE, col = "black",
      lwd = 2)

# Better plot.

sigma_mu$param <- "Sigma_mu"
sigma_phi$param <- "Sigma_phi"

# merge
sss <- rbind(sigma_mu, sigma_phi)
colnames(sss) <- c("N. dombeyi", "S. pagagonicus", "C. culeou")
slong <- pivot_longer(sss, 1:3, names_to = "Species", values_to = "y")

# densities better:
dd <- density(sigma_mu$coi)
dd$x
str(dd)

sp_names <- c("N. dombeyi", "S. pagagonicus", "C. culeou")
dens_mu <- do.call("rbind", lapply(1:3, function(i) {
  d <- density(sigma_mu[, i], from = 0)
  res <- data.frame(x = d$x, density = d$y, Species = sp_names[i],
                    param = "sigma_mu_sites")
  return(res)
}))

# compute prior density.
xp <- seq(min(dens_mu$x), max(dens_mu$x), length.out = nrow(dens_mu) / 3)
d_prior_mu <- dnorm(xp, 0, 0.1) * 2
df_prior_mu <- data.frame(x = xp, density = d_prior_mu, Species = "Prior for C. culeou",
                          param = "sigma_mu_sites")

dens_mu2 <- rbind(dens_mu, df_prior_mu)

# Repeat for phi
dens_phi <- do.call("rbind", lapply(1:3, function(i) {
d <- density(sigma_phi[, i], from = 0)
res <- data.frame(x = d$x, density = d$y, Species = sp_names[i],
                  param = "sigma_phi_sites")
return(res)
}))

# compute prior density.
xp <- seq(min(dens_phi$x), max(dens_phi$x), length.out = nrow(dens_phi) / 3)
d_prior_phi <- dnorm(xp, 0, 0.5) * 2
df_prior_phi <- data.frame(x = xp, density = d_prior_phi, Species = "Prior for C. culeou",
                          param = "sigma_phi_sites")

dens_phi2 <- rbind(dens_phi, df_prior_phi)


# merge and plot
dens <- rbind(dens_mu2, dens_phi2)
dens$Species <- factor(dens$Species, levels = c(
  "N. dombeyi", "S. pagagonicus", "C. culeou", "Prior for C. culeou"
))

colors <- viridis(4, option = "A", direction = -1, end = 0.8)
ggplot(dens, aes(x = x, y = density, colour = Species, linetype = Species,
                 size = Species)) + # , size = Species
  geom_line() +
  scale_size_manual(values = c(rep(0.7, 3), 0.5)) +
  scale_color_manual(values = colors) +
  scale_linetype_manual(values = c(1, 1, 1, 2)) +
  facet_wrap(vars(param), nrow = 1, scales = "free") +
  theme_bw() + 
  xlab("parameter") +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom")

ggsave("fmc_models/cania prior posterior sites.png",
       width = 17, height = 10, units = "cm")

