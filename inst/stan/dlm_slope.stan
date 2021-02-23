data {
  int<lower=0> N;
  int<lower=0> n_pos;
  vector[n_pos] y;
  int<lower=0> K; // number of covariates
  matrix[N, K] x; // matrix of covariates
  int y_int[n_pos];
  int pos_indx[n_pos+2];
  int family; // 1 = normal, 2 = binomial, 3 = poisson, 4 = gamma, 5 = lognormal
}
parameters {
  real x0;
  vector[K] beta0;
  vector[K] pro_dev[N-1];
  real<lower=0> sigma_process[K];
  real<lower=0> sigma_obs;
}
transformed parameters {
  vector[N] pred;
  vector[K] beta[N]; // elements accessed [N,K]

  for(k in 1:K) {
   beta[1,k] = beta0[k];
   for(i in 2:N) {
    beta[i,k] = beta[i-1,k] + sigma_process[k]*pro_dev[i-1,k];
   }
  }
  for(i in 1:N) {
    pred[i] = x[i] * beta[i] + x0;
  }

}
model {
  x0 ~ normal(0,10);
  sigma_obs ~ student_t(3,0,2);
  for(k in 1:K) {
    beta0[k] ~ normal(0,1);
    sigma_process[k] ~ student_t(3,0,2);
    pro_dev[k] ~ std_normal();//normal(0, sigma_process[k]);
  }
  if(family==1) {
    for(i in 1:(n_pos)) {
      //y ~ normal(pred, sigma_obs);
      y[i] ~ normal(pred[pos_indx[i]], sigma_obs);
    }
  }
  if(family==2) {
    for(i in 1:(n_pos)) {
      y_int[i] ~ bernoulli_logit(pred[pos_indx[i]]);
    }
  }
  if(family==3) {
    for(i in 1:(n_pos)) {
      y_int[i] ~ poisson_log(pred[pos_indx[i]]);
    }
  }
  if(family==4) {
    for(i in 1:(n_pos)) {
      y[i] ~ gamma(sigma_obs, sigma_obs ./ exp(pred[pos_indx[i]]));
    }
  }
  if(family==5) {
    for(i in 1:(n_pos)) {
      y[i] ~ lognormal(pred[pos_indx[i]], sigma_obs);
    }
  }
}
generated quantities {
  vector[n_pos] log_lik;
  // regression example in loo() package
  if(family==1) for (n in 1:n_pos) log_lik[n] = normal_lpdf(y[n] | pred[pos_indx[n]], sigma_obs);
  if(family==2) for (n in 1:n_pos) log_lik[n] = bernoulli_lpmf(y_int[n] | inv_logit(pred[pos_indx[n]]));
  if(family==3) for (n in 1:n_pos) log_lik[n] = poisson_lpmf(y_int[n] | exp(pred[pos_indx[n]]));
  if(family==4) for (n in 1:n_pos) log_lik[n] = gamma_lpdf(y[n] | sigma_obs, sigma_obs ./ exp(pred[pos_indx[n]]));
  if(family==5) for (n in 1:n_pos) log_lik[n] = lognormal_lpdf(y[n] | pred[pos_indx[n]], sigma_obs);
}
