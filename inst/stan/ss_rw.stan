data {
  int<lower=0> N;
  int<lower=0> n_pos;
  vector[n_pos] y;
  int pos_indx[n_pos+2];
}
parameters {
  real x0;
  vector[N-1] pro_dev;
  real<lower=0> sigma_process;
  real<lower=0> sigma_obs;
}
transformed parameters {
  vector[N] pred;
  pred[1] = x0;
  for(i in 2:N) {
    pred[i] = pred[i-1] + pro_dev[i-1];
  }
}
model {
  x0 ~ normal(0,10);
  sigma_process ~ cauchy(0,5);
  sigma_obs ~ cauchy(0,5);
  pro_dev ~ normal(0, sigma_process);
  for(i in 1:n_pos) {
    y[i] ~ normal(pred[pos_indx[i]], sigma_obs);
  }
}
generated quantities {
  vector[n_pos] log_lik;
  // regresssion example in loo() package
  for (n in 1:n_pos) log_lik[n] = normal_lpdf(y[n] | pred[pos_indx[n]], sigma_obs);
}
