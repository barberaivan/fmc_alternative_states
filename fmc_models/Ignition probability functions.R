# Regresión logística de FFMC = f(FMC) hecha por Bianchi y Defossé 2014
# para litter de cipresal y ñirantal.

# El objetivo de este script es obtener umbrales de FMC que estén asociados 
# a ciertas probabilidades de ignición o a cierto FFMC, que Bianchi describe
# como riesgo extremo y alto. 

# En realidad es más coherente pensar en umbrales de p(ig), pero quiero ver 
# a qué umbrales de p(ig) se corresponden las clases de FFMC.

# Bianchi hizo p(ig) = f(FMC) y luego FFMC = g(FMC).

# Umbrales de FFMC definidos por Bianchi:
# Low   //   High    //  Extreme
# <=85  // (85; 89]  //  >89

# Umbral de probabilidad:
# 0.5. A partir de ahí fires will start and spread. 

library(arm)  # logit
library(mgcv) # beta regression

setwd("/home/ivan/Insync/Humedad y EEA/Datos y scripts/Bianchi - ignition thresholds")
d_ig <- read.csv("Bianchi 2014 - p(ig) as f(MC) logistic.csv")
d_ffmc <- read.csv("Bianchi 2014 - FFMC as f(MC) model logistic 4 params.csv")

ffmc_model <- function(x, sp) {
  b1 <- d_ffmc[d_ffmc$sp == sp, "b1"]
  b2 <- d_ffmc[d_ffmc$sp == sp, "b2"]
  b3 <- d_ffmc[d_ffmc$sp == sp, "b3"]
  b4 <- d_ffmc[d_ffmc$sp == sp, "b4"]
  
  y <- b1 + (b2 - b1) / (1 + exp((b3 - x) / b4))
  return(y)
}

ig_model <- function(x, sp, d = d_ig_correct) {
  a <- d_ig[d_ig$sp == sp, "a"]
  b <- d_ig[d_ig$sp == sp, "b"]
  
  y <- plogis(a + b * x)
  return(y)
}

ig_model_inv <- function(p, sp) {
  a <- d_ig[d_ig$sp == sp, "a"]
  b <- d_ig[d_ig$sp == sp, "b"]
  
  p[p == 0] <- 1e-9
  p[p == 1] <- 1 - 1e-9
  
  # despejamos x
  # logit(p) <- a + b * x
  # (logit(p) - a) <- b * x
  # (logit(p) - a) / b <- x
  x <- (logit(p) - a) / b
  
  return(x)
}

curve(ffmc_model(x, "both"), to = 100)
curve(ig_model(x, "both"), to = 100)

# Bianchi dice que los MC para p = 0.5 son
# 25% para nire y 
# 21.5% para cip
ig_model(21.5, "nire")
ig_model(21.5, "cypress")

# nire
# dice 39.1%  <- pr(ig) = 0.01
ig_model(39.1, "nire") #3 me da 10%, no 0.1
ig_model_inv(0.5, "nire") # me da el umbral == 29.41 %

# cypress
# dice 21.5 % < pr(ig) = 0.5
ig_model_inv(0.5, "cypress") # 21.375
# Los parámetros de ñire están mal...
# El intercept que reporta debe estar más alto que el real.


#logit(p) <- a + b x
#logit(0.5) <- a + d_ig[d_ig$sp == "nire", "b"] * 25
a_nire_correct <- logit(0.5) - d_ig[d_ig$sp == "nire", "b"] * 25
a_nire_correct <- 5.47 # da 5.5, pero como decía 6.47, asumo que fue un error de tipeo.

# corrijo
d_ig[d_ig$sp == "nire", "a"] <- a_nire_correct
# probably both también tenga mal el intercept.

curve(ig_model(x, "nire"), to = 60)
abline(v = ig_model_inv(0.5, "nire"), lty = 2)
abline(h = 0)

ig_model_inv(0.5, "nire") # 24.86 %
ig_model_inv(0.8, "nire") # 18.56 %

ig_model_inv(0.95, "nire") # 11.48 %
ig_model_inv(0.05, "nire") # 38.25 %

ig_model_inv(0.99, "nire") # me da valores mucho más extremos que los que reporta en el paper
ig_model_inv(0.01, "nire")

ig_model(11, "nire") # 
ig_model(39.1, "nire") # 0.5

ffmc_model(11, "nire") # 87, pero el umbral que propone es 89
ffmc_model(5, "cypress") # 88.57, pero el umbral que propone es 89



# Conclu sobre modelos del paper ------------------------------------------
# Los params que reporta no coinciden con las curvas. Como las curvas son confiables, 
# reajustaré las curvas tomando "datos" de webplotdigitizer.

curve_data <- read.csv("ignition curve_web plot digitizer.csv")
curve_data$logit_p <- logit(curve_data$p)
nrow(curve_data)

# si quito los valores muy extremos, dará mejor el ajuste a su curva?
# curve_data <- curve_data[curve_data$p > 0.05 & curve_data$p < 0.95, ]
nrow(curve_data) # quité varios

m1_logit <- lm(logit_p ~ fmc * sp, curve_data)
m1 <- gam(p ~ fmc * sp, curve_data, method = "REML", family = betar())
coef(m1); coef(m1_logit)

plot(resid(m1) ~ predict(m1))
plot(resid(m1_logit) ~ predict(m1_logit))
# mejor m1, con regresión beta.

d_ig_correct <- d_ig
d_ig_correct[d_ig$sp == "nire", "a"] <- coef(m1)[1] + coef(m1)[3]
d_ig_correct[d_ig$sp == "nire", "b"] <- coef(m1)[2] + coef(m1)[4]

d_ig_correct[d_ig$sp == "cypress", "a"] <- coef(m1)[1]
d_ig_correct[d_ig$sp == "cypress", "b"] <- coef(m1)[2]

# recreo funciones
ig_model <- function(x, sp, d = d_ig_correct) {
  a <- d[d$sp == sp, "a"]
  b <- d[d$sp == sp, "b"]
  
  y <- plogis(a + b * x)
  return(y)
}

ig_model_inv <- function(p, sp, d = d_ig_correct) {
  a <- d[d$sp == sp, "a"]
  b <- d[d$sp == sp, "b"]
  
  p[p == 0] <- 1e-9
  p[p == 1] <- 1 - 1e-9
  
  # despejamos x
  # logit(p) <- a + b * x
  # (logit(p) - a) <- b * x
  # (logit(p) - a) / b <- x
  x <- (logit(p) - a) / b
  
  return(x)
}


ig_model(25, "nire") # close to 0.5: 0.5289062
ig_model(21.5, "cypress") # not so close to 0.5: 0.5472487


# Umbrales para usar:

ig_model_inv(0.5, "nire") # 25.3412
ig_model_inv(0.9, "nire") # 18.86454

ig_model_inv(0.5, "cypress") # 22.11533
ig_model_inv(0.9, "cypress") # 14.9829



# Con el evento extremo es fácil: mostrar la p(p(ig) >= 0.9),
# usando umbrales de cipresal para bosque y ñirantal para los demás.

# Con el evento de alto riesgo, qué conviene mostrar?
# p( 0.9 > p(ig) >= 0.5)
# or
# p(p(ig) >= 0.5) ?
# Me parece mejor la primera.

# O quizás hay alguna cuenta mejor. Por ejemplo, 
# dada una distribución estimada de FMC, 
# cuál es la esperanza de la p(ig)?

# Esas cuentas son más elegantes, pero van a ser menos entendidas. 
# mejor mostrar probabilidad de peligro alto y extremo.



# LFMC and ig freq (Bianchi 2018) -----------------------------------------

dlive <- read.csv("Bianchi 2018 table 1 - ig freq and lfmc.csv")
dlive$ig_prop[dlive$ig_prop == 1] <- 1 - 1e-6

plot(ig_prop ~ lfmc, dlive)

m2 <- gam(ig_prop ~ lfmc, data = dlive, method = "REML", 
          family = betar())
m3 <- gam(ig_prop ~ lfmc, data = dlive[dlive$group == "native", ], method = "REML", 
          family = betar())

# m2 <- gam(ig_prop ~ s(lfmc, k = 4), data = dlive, method = "REML", 
#           family = betar(), gamma = 0.6)
# no lineales no tienen mucho sense

x_new <- seq(60, 260, length.out = 200)
p <- predict(m2, data.frame(lfmc = x_new), type = "response")
p3 <- predict(m3, data.frame(lfmc = x_new), type = "response")
plot(ig_prop ~ lfmc, dlive[dlive$group == "native", ], xlim = c(60, 260), 
     ylim = c(0.3, 1), col = 2)
points(ig_prop ~ lfmc, dlive[dlive$group == "exotic", ], xlim = c(60, 260), col = 1)
lines(p ~ x_new, col = 1)
lines(p3 ~ x_new, col = 2)

coef(m2)

ig_live_inv <- function(p, model = m3) {
  a <- coef(model)[1]
  b <- coef(model)[2]
  
  p[p == 0] <- 1e-9
  p[p == 1] <- 1 - 1e-9
  
  # despejamos x
  # logit(p) <- a + b * x
  # (logit(p) - a) <- b * x
  # (logit(p) - a) / b <- x
  x <- (logit(p) - a) / b
  
  return(x)
}

ig_live_inv(0.5) # 209.8604 is lfmc with p(ig) = 0.5
ig_live_inv(0.9) # 124.5284 is lfmc with p(ig) = 0.9
