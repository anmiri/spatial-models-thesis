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
  matrix[ngr[1],2] lonlat_afr; 
  matrix[ngr[2],2] lonlat_eur;
  matrix[ngr[3],2] lonlat_nam; 
  matrix[ngr[4],2] lonlat_sam;
  matrix[ngr[5],2] lonlat_pap;
  matrix[ngr[6],2] lonlat_aus;
  real<lower=0> max_lon_afr; //maximum of longitude and latitude
  real<lower=0> max_lat_afr;
  real<lower=0> max_lon_eur; 
  real<lower=0> max_lat_eur; 
  real<lower=0> max_lon_nam; 
  real<lower=0> max_lat_nam;
  real<lower=0> max_lon_sam;
  real<lower=0> max_lat_sam;
  real<lower=0> max_lon_pap; 
  real<lower=0> max_lat_pap;
  real<lower=0> max_lon_aus; 
  real<lower=0> max_lat_aus;
  // for the approx gp stuff
  int M[2];
  int<lower=1> M_nD;
  matrix[M_nD, 2] indices;
  real<lower=1> c1;
  real<lower=1> c2;
}

transformed data {
  vector[2] L_afr;
  vector[2] L_eur;
  vector[2] L_nam;
  vector[2] L_sam;
  vector[2] L_pap;
  vector[2] L_aus;

  // PHI_2D(N_pred, M[1], M[2], L[1], L[2], x_pred[,1], x_pred[,2]);
  
  L_afr[1] = c1 * max_lon_afr;
  L_afr[2] = c2 * max_lat_afr;
  matrix[ngr[1],M_nD] phi_afr =
    PHI_2D(ngr[1], M[1], M[2], L_afr[1], L_afr[2], lonlat_afr[,1], lonlat_afr[,2]);
  
  L_eur[1] = c1 * max_lon_eur;
  L_eur[2] = c2 * max_lat_eur;
  matrix[ngr[2],M_nD] phi_eur =
    PHI_2D(ngr[2], M[1], M[2], L_eur[1], L_eur[2], lonlat_eur[,1], lonlat_eur[,2]);
  
  L_nam[1] = c1 * max_lon_afr;
  L_nam[2] = c2 * max_lat_afr;
  matrix[ngr[3],M_nD] phi_nam =
    PHI_2D(ngr[3], M[1], M[2], L_nam[1], L_nam[2], lonlat_nam[,1], lonlat_nam[,2]);
  
  L_sam[1] = c1 * max_lon_sam;
  L_sam[2] = c2 * max_lat_sam;
  matrix[ngr[4],M_nD] phi_sam =
    PHI_2D(ngr[4], M[1], M[2], L_sam[1], L_sam[2], lonlat_sam[,1], lonlat_sam[,2]);
  
  L_pap[1] = c1 * max_lon_pap;
  L_pap[2] = c2 * max_lat_pap;
  matrix[ngr[5],M_nD] phi_pap =
    PHI_2D(ngr[5], M[1], M[2], L_pap[1], L_pap[2], lonlat_pap[,1], lonlat_pap[,2]);
  
  L_aus[1] = c1 * max_lon_aus;
  L_aus[2] = c2 * max_lat_aus;
  matrix[ngr[6],M_nD] phi_aus =
    PHI_2D(ngr[6], M[1], M[2], L_aus[1], L_aus[2], lonlat_aus[,1], lonlat_aus[,2]);
}

parameters {
  // for the intercept
  real Intercept;
  real mu_area[nAreas];
  real<lower=0> sigma_intercept;
  // for the hurdle
  real hurdle_intercept;
  real hurdle_mu_area[nAreas];
  real<lower=0> hurdle_sigma_intercept;
  // gp beta
  vector[M_nD] beta_afr;     // the basis functions coefficients 
  vector[M_nD] beta_eur;     // the basis functions coefficients 
  vector[M_nD] beta_nam;     // the basis functions coefficients 
  vector[M_nD] beta_sam;     // the basis functions coefficients 
  vector[M_nD] beta_pap;     // the basis functions coefficients
  vector[M_nD] beta_aus;     // the basis functions coefficients 
  vector[M_nD] hurdle_beta_afr;     // the basis functions coefficients 
  vector[M_nD] hurdle_beta_eur;     // the basis functions coefficients 
  vector[M_nD] hurdle_beta_nam;     // the basis functions coefficients 
  vector[M_nD] hurdle_beta_sam;     // the basis functions coefficients 
  vector[M_nD] hurdle_beta_pap;     // the basis functions coefficients
  vector[M_nD] hurdle_beta_aus;     // the basis functions coefficients 
  // for the neg_bin
  real<lower=0> shape;       // shape parameter
  // for the GP
  real<lower=0> lscale[nAreas];            // length scales vector[nGP]
  real<lower=0> sdgp[nAreas];              // sd of the GPs (or single GP)
  real<lower=0> hurdle_lscale[nAreas];            // length scales vector[nGP]
  real<lower=0> hurdle_sdgp[nAreas];              // sd of the GPs (or single GP)
}

transformed parameters {
  // spectral densities (for f1)
  // one phimat for each area; ids -> N
  vector[N] mu;
  vector[N] hurdle_mu;

  // GP basis stuff
  // 2 is the number of dimensions
  vector[M_nD] diagSPD_afr = sqrt_diagSPD_nD(sdgp[1], lscale[1], L_afr, indices, 2);
  vector[M_nD] diagSPD_eur = sqrt_diagSPD_nD(sdgp[2], lscale[2], L_eur, indices, 2);
  vector[M_nD] diagSPD_nam = sqrt_diagSPD_nD(sdgp[3], lscale[3], L_nam, indices, 2);
  vector[M_nD] diagSPD_sam = sqrt_diagSPD_nD(sdgp[4], lscale[4], L_sam, indices, 2);
  vector[M_nD] diagSPD_pap = sqrt_diagSPD_nD(sdgp[5], lscale[5], L_pap, indices, 2);
  vector[M_nD] diagSPD_aus = sqrt_diagSPD_nD(sdgp[6], lscale[6], L_aus, indices, 2);

  vector[M_nD] hurdle_diagSPD_afr = sqrt_diagSPD_nD(hurdle_sdgp[1], hurdle_lscale[1], L_afr, indices, 2);
  vector[M_nD] hurdle_diagSPD_eur = sqrt_diagSPD_nD(hurdle_sdgp[2], hurdle_lscale[2], L_eur, indices, 2);
  vector[M_nD] hurdle_diagSPD_nam = sqrt_diagSPD_nD(hurdle_sdgp[3], hurdle_lscale[3], L_nam, indices, 2);
  vector[M_nD] hurdle_diagSPD_sam = sqrt_diagSPD_nD(hurdle_sdgp[4], hurdle_lscale[4], L_sam, indices, 2);
  vector[M_nD] hurdle_diagSPD_pap = sqrt_diagSPD_nD(hurdle_sdgp[5], hurdle_lscale[5], L_pap, indices, 2);
  vector[M_nD] hurdle_diagSPD_aus = sqrt_diagSPD_nD(hurdle_sdgp[6], hurdle_lscale[6], L_aus, indices, 2);
  
  // intercept stuff

  vector[M_nD] SPD_afr = diagSPD_afr .* beta_afr;
  vector[M_nD] SPD_eur = diagSPD_eur .* beta_eur;
  vector[M_nD] SPD_nam = diagSPD_nam .* beta_nam;
  vector[M_nD] SPD_sam = diagSPD_sam .* beta_sam;
  vector[M_nD] SPD_pap = diagSPD_pap .* beta_pap;
  vector[M_nD] SPD_aus = diagSPD_aus .* beta_aus;

  vector[M_nD] hurdle_SPD_afr = hurdle_diagSPD_afr .* hurdle_beta_afr;
  vector[M_nD] hurdle_SPD_eur = hurdle_diagSPD_eur .* hurdle_beta_eur;
  vector[M_nD] hurdle_SPD_nam = hurdle_diagSPD_nam .* hurdle_beta_nam;
  vector[M_nD] hurdle_SPD_sam = hurdle_diagSPD_sam .* hurdle_beta_sam;
  vector[M_nD] hurdle_SPD_pap = hurdle_diagSPD_pap .* hurdle_beta_pap;
  vector[M_nD] hurdle_SPD_aus = hurdle_diagSPD_aus .* hurdle_beta_aus;
  
  // main mu
  
  mu[ids_afr] = Intercept + mu_area[1] * sigma_intercept + (phi_afr[,] * SPD_afr); 
  mu[ids_eur] = Intercept + mu_area[2] * sigma_intercept + (phi_eur[,] * SPD_eur); 
  mu[ids_nam] = Intercept + mu_area[3] * sigma_intercept + (phi_nam[,] * SPD_nam); 
  mu[ids_sam] = Intercept + mu_area[4] * sigma_intercept + (phi_sam[,] * SPD_sam); 
  mu[ids_pap] = Intercept + mu_area[5] * sigma_intercept + (phi_pap[,] * SPD_pap); 
  mu[ids_aus] = Intercept + mu_area[6] * sigma_intercept + (phi_aus[,] * SPD_aus);

  // hurdle mu

  hurdle_mu[ids_afr] = hurdle_intercept + hurdle_mu_area[1] * hurdle_sigma_intercept + (phi_afr[,] * hurdle_SPD_afr); 
  hurdle_mu[ids_eur] = hurdle_intercept + hurdle_mu_area[2] * hurdle_sigma_intercept + (phi_eur[,] * hurdle_SPD_eur); 
  hurdle_mu[ids_nam] = hurdle_intercept + hurdle_mu_area[3] * hurdle_sigma_intercept + (phi_nam[,] * hurdle_SPD_nam); 
  hurdle_mu[ids_sam] = hurdle_intercept + hurdle_mu_area[4] * hurdle_sigma_intercept + (phi_sam[,] * hurdle_SPD_sam); 
  hurdle_mu[ids_pap] = hurdle_intercept + hurdle_mu_area[5] * hurdle_sigma_intercept + (phi_pap[,] * hurdle_SPD_pap); 
  hurdle_mu[ids_aus] = hurdle_intercept + hurdle_mu_area[6] * hurdle_sigma_intercept + (phi_aus[,] * hurdle_SPD_aus);

}

model {

  // different intercept per area, with shrinkage
  Intercept         ~ normal(0, 2);
  sigma_intercept   ~ lognormal(0,0.5);
  mu_area           ~ normal(0, 2);
  
  beta_afr ~ std_normal();
  beta_eur ~ std_normal();
  beta_nam ~ std_normal();
  beta_sam ~ std_normal();
  beta_pap ~ std_normal();
  beta_aus ~ std_normal();

  // Hurdle stuff
  
  hurdle_intercept         ~ normal(0, 2);
  hurdle_sigma_intercept   ~ lognormal(0, 0.5);
  hurdle_mu_area           ~ normal(0, 2);
  shape                    ~ exponential(1);
  
  
  hurdle_beta_afr ~ std_normal();
  hurdle_beta_eur ~ std_normal();
  hurdle_beta_nam ~ std_normal();
  hurdle_beta_sam ~ std_normal();
  hurdle_beta_pap ~ std_normal();
  hurdle_beta_aus ~ std_normal();
  
  for (i in 1:nAreas){
    lscale[i] ~ inv_gamma(5, 5);
    sdgp[i] ~ std_normal(); 
    hurdle_lscale[i] ~ inv_gamma(5, 5);
    hurdle_sdgp[i] ~ std_normal(); 
  }

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
