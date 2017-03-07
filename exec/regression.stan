data {
  int<lower=0> N;
  int<lower=0> K;
  matrix[N, K] x;
  vector[N] y;
  int y_int[N];
  int family; # 1 = normal, 2 = binomial, 3 = poisson, 4 = gamma, 5 = lognormal
}
parameters {
  vector[K] beta;
  real<lower=0> sigma;
}
transformed parameters {
  vector[N] pred;
  pred = x * beta;
}
model {
  beta ~ normal(0,2);
  sigma ~ cauchy(0, 5);
  if(family==1) y ~ normal(pred, sigma);
  if(family==2) y_int ~ bernoulli_logit(pred);
  if(family==3) y_int ~ poisson_log(pred);
  if(family==4) y ~ gamma(sigma, sigma ./ exp(pred));
  if(family==5) y ~ lognormal(pred, sigma);  
}
generated quantities {
  vector[N] log_lik;
  # regresssion example in loo() package 
  if(family==1) for (n in 1:N) log_lik[n] = normal_lpdf(y[n] | pred[n], sigma);
  if(family==2) for (n in 1:N) log_lik[n] = bernoulli_lpmf(y_int[n] | inv_logit(pred[n]));
  if(family==3) for (n in 1:N) log_lik[n] = poisson_lpmf(y_int[n] | exp(pred[n]));
  if(family==4) for (n in 1:N) log_lik[n] = gamma_lpdf(y[n] | sigma, sigma ./ exp(pred[n]));
  if(family==5) for (n in 1:N) log_lik[n] = lognormal_lpdf(y[n] | pred[n], sigma);  
}
