# Simulating data with temporal correlation and trying to recover the 
# real parameters. What happens with low N?

library(rstan)

model_code <- 
"
data{
  int N;
  vector[N] y;
  matrix[N, N] d;
  real prior_sd_mu;
  real prior_sd_sigma;
}

parameters {
  real mu;
  real<lower = -1, upper = 1> rho;
  real<lower = 0> sigma;
}

transformed parameters {
  matrix[N, N] Sigma; 
  for(i in 1:N) {
    for(j in 1:N) {
      Sigma[i, j] = rho ^ d[i, j] * sigma ^ 2;
    }
  }
}

model {
  mu ~ normal(0, prior_sd_mu);
  sigma ~ normal(0, prior_sd_sigma);
  y ~ multi_normal(rep_vector(mu, N), Sigma);
}
"

stan_code <- stan_model(model_code = model_code, verbose = TRUE)

mu <- 0 
sigma <- 1
rho <- 0.9
N <- 100   
time <- 1:N
d <- as.matrix(dist(time))
Sigma <- rho ^ d * sigma ^ 2
y <- mgcv::rmvn(1, rep(mu, N), Sigma)
plot(y ~ time, main = round(cor(y[-1], y[-N]), 3))
# curve(rho ^ x, from = 0, to = N, xlab = "temporal distance", 
#       ylab = "correlation")

stan_data <- list(y = y, N = N, d = d,
                  prior_sd_sigma = 10, 
                  prior_sd_mu = 10)

mod <- sampling(stan_code, data = stan_data, cores = 3, 
                chains = 3, iter = 2000)

pairs(mod, pars = c("mu", "sigma", "rho"))


# N = 20, rho = 0.9, sigma = 1, mu = 0, prior_sd = 10:
# I find the same patterns as with the data: terribly high correlation
# between rho and sigma. Rho can take huge values that are compensated by
# a huge variance and a very uncertain mean. It's similar with N = 9 and 
# rho = 0.4, but few data simulations were run.

# N = 100: the same happens.
# rho values close to 1 make possible extreme variances and mu. 
# Perhaps with lower rho or a regularizing prior on it, extreme variances 
# would show lower posterior density.

# N = 100, rho = 0.4: Now it works, although a high correlation between
# sigma and rho still remains.

# Conclusion: estimating correlation coefficients with very low N is
# challenging, requiring informative priors. It would even be necessary in 
# case of large N bur very high correlation (~ 0.9), where the prior for
# sigma would probably need an informative or at least regularizing prior.