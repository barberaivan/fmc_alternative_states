library(tidyverse)
library(lme4)
f <- "Live fuel mixture"#"Litter"#

# Data -------------------------------------------------------------------

fdata <- readRDS("data_fmc_dataframe.R")
fdata <- fdata[fdata$fuel_type == f, ]
air_data <- readRDS("microclimate/predictions_vpd.R")
air_temp <- air_data$together

# Merge
fdata <- left_join(fdata, 
                   air_temp[, c("Community", "transect", "vpd_max")], 
                   by = c("transect", "Community"))
fdata$vpd <- scale(fdata$vpd_max) %>% as.numeric


m1 <- lm(moisture ~ vpd, data = fdata)
summary(m1)

m2 <- lm(moisture ~ Community + vpd, data = fdata)
summary(m2)

m3 <- lmer(moisture ~ vpd + (1 | site), data = fdata)
m4 <- lmer(moisture ~ Community + vpd + (1 | site), data = fdata)

fixef(m3)
fixef(m4)
cov2cor(vcov(m4))
cov2cor(vcov(m3))

cov2cor(vcov(m1))
cov2cor(vcov(m2))

# como comunidad y vpd están totalmente correlated, los parámetros también,
# dificultando terriblemente la interpretabilidad.
# esto es por culpa del N?

C <- 4
N <- 3
means <- seq(10, 40, length.out = C)
sdd <- 10
sigma <- 30
beta <- 10
x <- do.call("c", lapply(1:C, function(i) rnorm(N, means[i], sdd)))
y <- rnorm(N * C, beta * x, sigma)
d <- data.frame(y = y, x = x, com = rep(letters[1:C], each = N))

t1 <- lm(y ~ x, data = d)
t2 <- lm(y ~ com + x, data = d)

coef(t1); coef(t2)
cov2cor(vcov(t1))
cov2cor(vcov(t2))

# si solo hay efecto de vpd (x), no hay problema, el modelo se estima bien y
# no hay grandes correlaciones. Veamos qué pasa si hay eff de otra covariable
# muy correlated con com.

C <- 4
N <- 2
means <- seq(10, 40, length.out = C)
sdd <- 1
sigma <- 30
beta_x <- 10
beta_z <- 0
x <- do.call("c", lapply(1:C, function(i) rnorm(N, means[i], sdd)))
z <- do.call("c", lapply(1:C, function(i) rnorm(N, means[i], sdd)))
cor(x, z)
y <- rnorm(N * C, beta_x * x + beta_z * z, sigma)
d <- data.frame(y = y, x = x, com = rep(letters[1:C], each = N))

t1 <- lm(y ~ x, data = d)
t2 <- lm(y ~ com + x, data = d)

coef(t1); coef(t2)
cov2cor(vcov(t1))
cov2cor(vcov(t2))
anova(t2)

summary(t2)

ggplot(d, aes(x = com, y = x)) + geom_point()

# Empecé a poder reproducir lo que pasaba con los datos cuando aumenté la 
# correlación entre comunidad y vpd (x). Ahí el modelo que incluye comunidad
# tiene corr muy altas. (metiendo eff de z.)
# Sin meter eff de z pasa lo mismo.

# Si el N es muy bajo y la corr entre comunidad y vpd (x) es alta, meter
# comunidad compite con vpd, generando errores muy grandes en la estimación del
# beta_vpd.