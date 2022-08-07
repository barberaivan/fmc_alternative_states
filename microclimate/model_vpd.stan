data {
  int N;       // days
  int N_sites; // sites
  int N_com;   // communities
  int com[N_sites];  // community id
  /* 
     Community id is only N_sites long because 
     I use a matrix to hold the response variable, where every column is a site
     and I need to identify the community for every site.
  */
  matrix[N, N_sites] y; // log VPD
  matrix[N, N] d; // tempral distance matrix (days)
  
  // priors
  real prior_eta_sd;
  real prior_sigma_sd;
  real prior_sigma_delta_sd;
  
}

parameters {
  vector[N_com] eta_raw;  
  vector[N_sites] delta_raw;
  vector[N_com] sigma_raw;
  real<lower = 0> sigma_delta_raw;
  real<lower = 0, upper = 1> rho;
}

transformed parameters {
  vector[N_com] eta = eta_raw * prior_eta_sd;  
  vector[N_com] sigma = sigma_raw * prior_sigma_sd;
  real<lower = 0> sigma_delta = sigma_delta_raw * prior_sigma_delta_sd;
  vector[N_sites] delta = delta_raw * sigma_delta;
  
  vector[N_sites] mu_sites; // site means at log scale
  vector[N_sites] mu_exp;   // site means at exp scale
  vector[N_com] mu_com;     // communities means at exp scale
  
  matrix[N, N] cor_mat; 
  matrix[N, N] cor_mat_chol;
  
  for(s in 1:N_sites) {
    mu_sites[s] = eta[com[s]] + delta[s];
    mu_exp[s] = exp(mu_sites[s] + 0.5 * sigma[com[s]] ^ 2);
  }
  
  for(c in 1:N_com) {
    real vari = sigma[c] ^ 2 + sigma_delta ^ 2;
    mu_com[c] = exp(eta[c] + 0.5 * vari);
  }
  
  for(i in 1:N) {
    for(j in 1:N) {
      cor_mat[i, j] = rho ^ d[i, j];
    }
  }
    
  cor_mat_chol = cholesky_decompose(cor_mat);
}

model {
  
  eta_raw ~ std_normal();
  delta_raw ~ std_normal();
  sigma_raw ~ std_normal();
  sigma_delta_raw ~ std_normal();
  
  for(s in 1:N_sites) {
    // compute variance-covariance matrix cholesky factor
    matrix[N, N] vcov_chol = sigma[com[s]] * cor_mat_chol;
    // likelihood
    y[, s] ~ multi_normal_cholesky(rep_vector(mu_sites[s], N), vcov_chol); 
  }
}  

