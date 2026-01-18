functions {
  //Eigenvalue function
  real lambda_1D(real L, int m) {
    return ((m*pi())/(2*L))^2;
  }
  vector lambda_nD(vector L, row_vector m, int D) {
    vector[D] lam;
    for(i in 1:D){
      lam[i] = ((m[i]*pi())/(2*L[i]))^2; }
    return lam;
  }
	
  //Square root of spectral densitiy function (for a squared exponential kernel)
  real sqrt_spd_1D(real gpscale, real lscale, real w) {
    return gpscale * sqrt(sqrt(2*pi()) * lscale) * exp(-.25*(lscale^2)*(w^2));
  }
  real sqrt_spd_2D(real gpscale, real lscale1, real lscale2, real w1, real w2) {
    return gpscale * sqrt(sqrt(2*pi())^2 * lscale1*lscale2) * exp(-.25*(lscale1^2*w1^2 + lscale2^2*w2^2));
  }
  real sqrt_spd_nD(real gpscale, vector lscale, vector w, int D) {
    return gpscale * sqrt(sqrt(2*pi())^D * prod(lscale)) * exp(-.25*((to_row_vector(lscale) .* to_row_vector(lscale)) * (w .* w)));
  }

  // vector diagSPD_Matern32(real alpha, real rho, real L, int M) {
  // return 2*alpha * (sqrt(3)/rho)^1.5 * inv((sqrt(3)/rho)^2 + ((pi()/2/L) * linspaced_vector(M, 1, M))^2);
  // }

  //Square root vector of spectral densities (for a squared exponential kernel)
  vector sqrt_diagSPD_1D(real gpscale, real lscale, real L, int M) {
    return gpscale * sqrt(sqrt(2*pi()) * lscale) * exp(-.25*(lscale*pi()/2/L)^2 * linspaced_vector(M, 1, M)^2);
  }
  // vector sqrt_diagSPD_nD(real gpscale, vector lscale, vector L, matrix indices, int D) {
  // return gpscale *  sqrt(sqrt(2*pi())^D * prod(lscale)) * exp(-.25 * (indices^2 * (lscale*pi() ./ (2*L))^2));
  // }
  vector sqrt_diagSPD_nD(real gpscale, real lscale, vector L, matrix indices, int D) {
    return gpscale *  sqrt(sqrt(2*pi())^D * lscale^2) * exp(-.25 * (indices^2 * (lscale*pi() ./ (2*L))^2));
  }
  
  //Eigenfunction
  vector phi_1D(real L, int m, vector x) {
    return 1/sqrt(L) * sin(m*pi()/(2*L) * (x+L));
  }
  vector phi_2D(real L1, real L2, int m1, int m2, vector x1, vector x2) {
    vector[rows(x1)] fi1;
    vector[rows(x1)] fi2;
    fi1 = 1/sqrt(L1)*sin(m1*pi()*(x1+L1)/(2*L1));
    fi2 = 1/sqrt(L2)*sin(m2*pi()*(x2+L2)/(2*L2));
    return fi1 .* fi2;
  }
  vector phi_nD(vector L, row_vector m, matrix x) {
    int c = cols(x);
    int r = rows(x);
    matrix[r,c] fi;
    vector[r] fi1;
    for (i in 1:c){
      fi[,i] = 1/sqrt(L[i])*sin(m[i]*pi()*(x[,i]+L[i])/(2*L[i]));
    }
    fi1 = fi[,1];
    for (i in 2:c){
      fi1 = fi1 .* fi[,i];
    }
    return fi1;
  }
	
  //Matrix of eigenfunction values
  matrix PHI_1D(int N, int M, real L, vector x) {
    matrix[N,M] PHI = sin(diag_post_multiply(rep_matrix(pi()/(2*L) * (x+L), M), linspaced_vector(M, 1, M)))/sqrt(L);
    return PHI;
  }
  matrix PHI_2D(int N, int M1, int M2, real L1, real L2, vector x1, vector x2) {
    matrix[N,M1*M2] PHI;
    matrix[N,M1] PHI_1 = sin(diag_post_multiply(rep_matrix(pi()/(2*L1) * (x1+L1), M1), linspaced_vector(M1, 1, M1)))/sqrt(L1);
    matrix[N,M2] PHI_2 = sin(diag_post_multiply(rep_matrix(pi()/(2*L2) * (x2+L2), M2), linspaced_vector(M2, 1, M2)))/sqrt(L2);
    PHI[,1:M2] = rep_matrix(PHI_1[,1], M2) .* PHI_2;
    for(i in 2:M1)
      PHI[,1:(M2*i)] = append_col(PHI[,1:(M2*(i-1))], rep_matrix(PHI_1[,i], M2) .* PHI_2);
    return PHI;
  }

  // distribution functions
  real hurdle_neg_binomial_log_logit_lpmf(int y, real eta, real phi, real hu) {
    if (y == 0) { 
      return bernoulli_logit_lpmf(1 | hu); 
    } else { 
      return bernoulli_logit_lpmf(0 | hu) +  
	neg_binomial_2_log_lpmf(y | eta, phi) - 
	log1m((phi / (exp(eta) + phi))^phi); 
    } 
  }
}

data {
  int<lower=0> N; // sample size
  int<lower=0> nAreas; // for the world: 6
  int<lower=0> y[N];      // outcomes
  int<lower=0> ngr[nAreas]; //number of obs in each area
  int ids_afr[ngr[1]];
  int ids_eur[ngr[2]];
  int ids_nam[ngr[3]];
  int ids_sam[ngr[4]];
  int ids_pap[ngr[5]];
  int ids_aus[ngr[6]];
}

parameters {
  // for the intercept
  real Intercept;
  real mu_area[nAreas];
  real<lower=0> sigma_intercept;
  // for the hurdle
  real hurdle_intercept;
  real hurdle_mu_area[nAreas];
  real hurdle_sigma_intercept;
  // for the neg_bin
  real<lower=0> shape;       // shape parameter
}

transformed parameters {

  vector[N] mu;
  vector[N] hurdle_mu;
  
  for (n in ids_afr) {
    mu[n] = Intercept + mu_area[1] * sigma_intercept; 
    hurdle_mu[n] = hurdle_intercept + hurdle_mu_area[1] * hurdle_sigma_intercept; 

  }
  for (n in ids_eur) {
    mu[n] = Intercept + mu_area[2] * sigma_intercept; 
    hurdle_mu[n] = hurdle_intercept + hurdle_mu_area[2] * hurdle_sigma_intercept; 

  }
  for (n in ids_nam) {
    mu[n] = Intercept + mu_area[3] * sigma_intercept; 
    hurdle_mu[n] = hurdle_intercept + hurdle_mu_area[3] * hurdle_sigma_intercept; 
  }
  for (n in ids_sam) {
    mu[n] = Intercept + mu_area[4] * sigma_intercept; 
    hurdle_mu[n] = hurdle_intercept + hurdle_mu_area[4] * hurdle_sigma_intercept; 
  }
  for (n in ids_pap) {
    mu[n] = Intercept + mu_area[5] * sigma_intercept; 
    hurdle_mu[n] = hurdle_intercept + hurdle_mu_area[5] * hurdle_sigma_intercept; 
  }
  for (n in ids_aus) {
    mu[n] = Intercept + mu_area[6] * sigma_intercept; 
    hurdle_mu[n] = hurdle_intercept + hurdle_mu_area[6] * hurdle_sigma_intercept; 
  }
}

model {

  // different intercept per area, with shrinkage
  Intercept         ~ normal(0, 2);
  sigma_intercept   ~ lognormal(0, 0.5);
  mu_area           ~ normal(0, 2);

  // Hurdle stuff

  hurdle_intercept       ~ normal(0, 2);
  shape                  ~ exponential(1); // lognormal(0, 5)
  hurdle_mu_area         ~ normal(0, 2);
  hurdle_sigma_intercept ~ lognormal(0, 0.5);
  

  // int y, real eta, real phi, real hu
  for (n in 1:N) {
    target += hurdle_neg_binomial_log_logit_lpmf(y[n] | mu[n], shape, hurdle_mu[n]);
  }

}

generated quantities{
  
  vector[N] log_lik;
  
  for (n in 1:N) {
    log_lik[n] = hurdle_neg_binomial_log_logit_lpmf(y[n] | mu[n], shape, hurdle_mu[n]);
  }
  
}

