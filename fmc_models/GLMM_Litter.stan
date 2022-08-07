/*
  mu_ for parameters affecting location,
  phi_for parameters affectinc dispersion.
  ds stands for date-site (date * site interaction).
  sigma is used for standard deviations of random effects.
  
  In this version, I remove the hierarchical variance for the date-site effects,
  making the variances fixed effects.
  If the ds variance is fixed among sites, the model fits poorly.
  
*/

data{
  int<lower = 1> N;
  int<lower = 1> N_dates;
  int<lower = 1> N_sites;
  int<lower = 1> N_unique;
  int<lower = 1> N_pred;
  int<lower = 2> K; // basis functions for spline (without intercept)

  vector[N] y;
  matrix[N, N_sites] site_matrix;
  matrix[N, N_sites * N_dates] ds_matrix;
  //vector[N_sites] vpd;
  //vector[N_pred] vpd_pred;
  matrix[N_sites, 2] basis_fit; // matrices for spline
  matrix[N_pred, 2] basis_pred;
  int unique_id[N_unique]; // indexes of unrepeated date-sites (eliminating point)

  // Priors --------------------------------------------------------------
  
  // mu
  
    // fixed effects
  real prior_mu_eta_mean;
  real prior_mu_eta_sd;
  real prior_mu_lambda_sd;
  
    // random effects
  real prior_mu_sigma_sites_sd;
  real prior_mu_sigma_ds_sd;

  // phi
  
    // fixed effects
  real prior_phi_eta_mean;
  real prior_phi_eta_sd;

    // random effects
  real prior_phi_sigma_sites_sd;
  real prior_phi_sigma_ds_sd;
}

parameters{
// for mu -----------------------------------------------------------------
  
  // fixed effects (intercept and VPD)
  real mu_eta_raw;
  vector[K] mu_lambda_raw;
  
  // site effect 
  vector[N_sites] mu_error_sites_raw;
  real<lower = 0> mu_sigma_sites_raw;  
  
  // date-site effect 
  matrix[N_dates, N_sites] mu_error_ds_raw; // ds for date - site
  vector<lower = 0>[N_sites] mu_sigma_ds_raw;

  // for phi ----------------------------------------------------------------
  
  // fixed effects (intercept)
  real phi_eta_raw;

  // site effect 
  vector[N_sites] phi_error_sites_raw;
  real<lower = 0> phi_sigma_sites_raw;
  
  // date-site effect 
  matrix[N_dates, N_sites] phi_error_ds_raw; // ds for date - site
  vector<lower = 0>[N_sites] phi_sigma_ds_raw;
}

transformed parameters{
  vector[N] mu;
  vector[N] phi;

  // for mu -----------------------------------------------------------------
  
  // fixed effects (intercept and VPD)
  real mu_eta = mu_eta_raw * prior_mu_eta_sd + prior_mu_eta_mean;
  vector[K] mu_lambda = mu_lambda_raw * prior_mu_lambda_sd;

  // site effect 
  real<lower = 0> mu_sigma_sites = mu_sigma_sites_raw * prior_mu_sigma_sites_sd;
  vector[N_sites] mu_error_sites = mu_error_sites_raw * mu_sigma_sites + 
                                   basis_fit * mu_lambda;
  
  // date * site effect 
  vector<lower = 0>[N_sites] mu_sigma_ds = mu_sigma_ds_raw * prior_mu_sigma_ds_sd; // variances between dates for every site
  matrix[N_dates, N_sites] mu_error_ds; // to be filled below
  
  // for phi ----------------------------------------------------------------
  
  // fixed effects (intercept)
  real phi_eta = phi_eta_raw * prior_phi_eta_sd + prior_phi_eta_mean;
  
  // site effect 
  real<lower = 0> phi_sigma_sites = phi_sigma_sites_raw * prior_phi_sigma_sites_sd;
  vector[N_sites] phi_error_sites = phi_error_sites_raw * phi_sigma_sites;
  
  // date * site effect 
  vector<lower = 0>[N_sites] phi_sigma_ds = phi_sigma_ds_raw * prior_phi_sigma_ds_sd; // variances between dates for every site
  matrix[N_dates, N_sites] phi_error_ds; // to be filled below
  
  
  // fill date-site random effects for mu and phi ---------------------------
  
  for(s in 1:N_sites) {
    mu_error_ds[, s] = mu_sigma_ds[s] * mu_error_ds_raw[, s];
    phi_error_ds[, s] = phi_sigma_ds[s] * phi_error_ds_raw[, s];
  }
  
  // Compute mu and phi vectors -----------------------------------------
  
  mu = exp(mu_eta + 
           site_matrix * mu_error_sites +
           ds_matrix * to_vector(mu_error_ds));
  
  phi = exp(phi_eta + 
            site_matrix * phi_error_sites +
            ds_matrix * to_vector(phi_error_ds)); 

}

model{
  // Priors ---------------------------------------------------------------
  
  // mu
  mu_eta_raw ~ std_normal();
  mu_lambda_raw ~ std_normal();

  mu_error_sites_raw ~ std_normal();
  to_vector(mu_error_ds_raw) ~ std_normal();
  
  mu_sigma_sites_raw ~ std_normal();
  mu_sigma_ds_raw ~ std_normal();

  // phi
  phi_eta_raw ~ std_normal();

  phi_error_sites_raw ~ std_normal();
  to_vector(phi_error_ds_raw) ~ std_normal();
  
  phi_sigma_sites_raw ~ std_normal();
  phi_sigma_ds_raw ~ std_normal();

  // Likelihood ---------------------------------------------------------
  
  for(n in 1:N)
    y[n] ~ normal(mu[n], phi[n])T[0, ];
}


generated quantities {
  
  /* 
    Compute site-level mu and phi marginal to the date random effects. 
    In the case of Truncated Normal likelihood, mu and phi will be used to 
    compute the mean.
  
    mus and phis follow a log-normal distribution, because there is a date*site
    normal random effect within sites, at the log scale.
    lognormal mean is = exp(mu + 0.5 sigma ^ 2)
  */
  
  vector[N_sites] sites_mu;  // to fill below
  vector[N_sites] sites_phi;
  
  // VPD prediction marginal to site
  
    // get variance between dates to compute their means (not perfect, but approximate)
  vector[N_sites] mu_sigma_ds2 = mu_sigma_ds .* mu_sigma_ds;
  vector[N_sites] phi_sigma_ds2 = phi_sigma_ds .* phi_sigma_ds;
  
    // compute prediction.
  vector[N_pred] pred_mu = exp(mu_eta + basis_pred * mu_lambda + 
                               0.5 * (mu_sigma_sites ^ 2 + mean(mu_sigma_ds2)));
  real pred_phi = exp(phi_eta + 
                      0.5 * (phi_sigma_sites ^ 2 + mean(phi_sigma_ds2)));
  
  // Fill mu and phi by site
  for(s in 1:N_sites) {
    sites_mu[s] = exp(mu_eta + mu_error_sites[s] + 0.5 * mu_sigma_ds[s] ^ 2);
    sites_phi[s] = exp(phi_eta + phi_error_sites[s] + 0.5 * phi_sigma_ds[s] ^ 2);
  }

  
}



