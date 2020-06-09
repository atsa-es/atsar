data {
  int<lower=0> N; // number of data points
  int<lower=0> P; // number of time series of data
  int<lower=0> K; // number of trends
  int<lower=0> nZ; // number of unique z elements
  int<lower=0> row_indx[nZ];
  int<lower=0> col_indx[nZ];
  int<lower=0> nVariances;
  int<lower=0> varIndx[P];
  int<lower=0> nZero;
  int<lower=0> row_indx_z[nZero];
  int<lower=0> col_indx_z[nZero];
  int<lower=0> n_pos; // number of non-missing observations
  int<lower=0> row_indx_pos[n_pos]; // row indices of non-missing obs
  int<lower=0> col_indx_pos[n_pos]; // col indices of non-missing obs
  real y[n_pos]; // vectorized matrix of observations
  int<lower=0> n_na; // number of missing observations
  int<lower=0> row_indx_na[n_na]; // row indices of missing obs
  int<lower=0> col_indx_na[n_na]; // col indices of missing obs
  int<lower=0> num_covar; // number of unique covariates
  int<lower=0> num_unique_covar; // number of covar parameters to estimate
  matrix[num_covar,N] d_covar; // inputted covariate matrix
  int covar_indexing[P,num_covar]; // index of covariates to estimate
}
transformed data {
  int n_pcor; // dimension for cov matrix
  int n_loglik; // dimension for loglik calculation
  n_loglik = P;
  n_pcor = P;
  if(nVariances < 2) {
    n_pcor = 2;
  }
}
parameters {
  matrix[K,N] x; //vector[N] x[P]; // random walk-trends
  vector[nZ] z; // estimated loadings in vec form
  vector<lower=0>[K] zpos; // constrained positive values
  real<lower=0> sigma[nVariances];
  real ymiss[n_na];
  vector[num_unique_covar] D; // optional covariates to estimate
}
transformed parameters {
  matrix[P,N] pred; //vector[P] pred[N];
  matrix[P,K] Z;
  //vector[N] yall[P]; // combined vectors of missing and non-missing values
  matrix[P,N] yall;
  vector[P] sigma_vec;

  for(p in 1:P) {
    sigma_vec[p] = sigma[varIndx[p]]; // convert estimated sigmas to vec form
  }

  // Fill yall with non-missing values
  for(i in 1:n_pos) {
    yall[row_indx_pos[i], col_indx_pos[i]] = y[i];
  }
  // Include missing observations
  if(n_na > 0) {
    for(i in 1:n_na) {
      yall[row_indx_na[i], col_indx_na[i]] = ymiss[i];
    }
  }

  for(i in 1:nZ) {
    Z[row_indx[i],col_indx[i]] = z[i]; // convert z to from vec to matrix
  }
  // fill in zero elements
  if(nZero > 2) {
    for(i in 1:(nZero-2)) {
      Z[row_indx_z[i],col_indx_z[i]] = 0;
    }
  }

  for(k in 1:K) {
    Z[k,k] = zpos[k];// add constraint for Z diagonal
  }
  // N is sample size, P = time series, K = number trends
  // [PxN] = [PxK] * [KxN]
  pred = Z * x;

  // include covariates if specified
  if(num_covar > 0) {
    for(p in 1:P) {
      for(n in 1:N) {
        for(k in 1:num_covar) {
          pred[p,n] = pred[p,n] + d_covar[k,n] * D[covar_indexing[p,k]];
        }
      }
    }
  }
}
model {
  // initial state for each trend
  for(k in 1:K) {
    x[k,1] ~ cauchy(0,3);//normal(0, 1);
    for(t in 2:N) {
      x[k,t] ~ normal(x[k,t-1], 1);
    }
  }

  // prior on covar
  if(num_covar > 0) {
    D ~ normal(0,1);
  }
  // prior on loadings
  z ~ normal(0, 1);
  zpos ~ normal(0, 1);
  // observation variance
  sigma ~ student_t(3, 0, 2);

  // likelihood for independent
  for(i in 1:P){
    target += normal_lpdf(yall[i] | pred[i], sigma_vec[i]);
  }

}
generated quantities {
  vector[n_loglik] log_lik;
  //calculate looic based on regresssion example in loo() package
  for(i in 1:P) {
    log_lik[i] = normal_lpdf(yall[i] | pred[i], sigma_vec[i]);
  }
}
