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

data{
  int<lower=0> N;
  matrix[N,2] lonlat; 
  real<lower=0> max_lon; 
  real<lower=0> max_lat;
  int M[2];
  int<lower=1> M_nD;
  matrix[M_nD, 2] indices;
  real<lower=1> c1;
  real<lower=1> c2;
  real<lower=0> lscale;
  real<lower=0> sdgp; 
  real<lower=0> shape;// shape parameter for neg bin
  real hurdle_intercept;
  real<lower=0> hurdle_lscale;
  real<lower=0> hurdle_sdgp;
  real Intercept;
  real mu_area;
  real hurdle_mu_area;
  real<lower=0> sigma_intercept;
  real<lower=0> hurdle_sigma_intercept;
  vector[M_nD] beta;     // for one area
  vector[M_nD] hurdle_beta;
  vector[N] alpha; // phylo
}



generated quantities{
  vector[N] mu;
  vector[N] hurdle_mu;
  array[N] int preds;
  array[N] int hurdle_preds;
  vector[2] L;
  matrix[N,M_nD] phi;
  vector[M_nD] diagSPD;
  vector[M_nD] hurdle_diagSPD;
  
  L[1] = c1 * max_lon;
  L[2] = c2 * max_lat;
  
  phi = PHI_2D(N, M[1], M[2], L[1], L[2], lonlat[,1], lonlat[,2]);
  
  diagSPD = sqrt_diagSPD_nD(sdgp, lscale, L, indices, 2);
  hurdle_diagSPD = sqrt_diagSPD_nD(hurdle_sdgp, hurdle_lscale, L, indices, 2);

  vector[M_nD] SPD = diagSPD .* beta;
  vector[M_nD] hurdle_SPD = hurdle_diagSPD .* hurdle_beta;
  
  mu = Intercept + mu_area * sigma_intercept + alpha + (phi[,] * SPD);
  hurdle_mu = hurdle_intercept + hurdle_mu_area * hurdle_sigma_intercept + alpha + (phi[,] * hurdle_SPD);
  
  for (n in 1:N) {
    
    if (bernoulli_logit_rng(hurdle_mu[n])) {
      preds[n] = 0;
      hurdle_preds[n] = 0;
    } else {
      hurdle_preds[n] = 1;
      int w; // temporary variable
        w = neg_binomial_2_log_rng(mu[n], shape);
      while (w == 0) { 
        w = neg_binomial_2_log_rng(mu[n], shape);
      }
    preds[n] = w;
    }
  }
}
