data {
  int<lower=0> N;
  vector[N] y;
}
parameters {
  real x0;
  real<lower=-0.999,upper=0.999> phi;
  vector[N-1] pro_dev;
  real<lower=0> sigma_process;
  real<lower=0> sigma_obs;
}
transformed parameters {
  vector[N] pred;
  pred[1] = x0;
  for(i in 2:N) {
    pred[i] = phi*pred[i-1] + pro_dev[i-1];
  }
}
model {
  x0 ~ normal(0,10);
  phi ~ normal(0,1);
  sigma_process ~ cauchy(0,5);
  sigma_obs ~ cauchy(0,5);
  pro_dev ~ normal(0, sigma_process);
  y ~ normal(pred, sigma_obs);
}
generated quantities {
  vector[N] log_lik;
  # regresssion example in loo() package 
  for (n in 1:N) log_lik[n] = normal_lpdf(y[n] | pred[n], sigma_obs);
}
