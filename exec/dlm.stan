data {
  int<lower=0> N;
  vector[N] y;
  int<lower=0> K; # number of covariates
  matrix[N, K] x; # matrix of covariates
  int y_int[N];
  int family; # 1 = normal, 2 = binomial, 3 = poisson, 4 = gamma, 5 = lognormal  
}
parameters {
  #real x0;
  vector[K] beta0;  
  vector[K] pro_dev[N-1];
  real<lower=0> sigma_process[K];
  real<lower=0> sigma_obs;
}
transformed parameters {
  vector[N] pred;
  vector[K] beta[N]; # elements accessed [N,K]
  
  for(k in 1:K) {
   beta[1,k] = beta0[k];
   for(i in 2:N) {
    beta[i,k] = beta[i-1,k] + pro_dev[i-1,k];
   }
  }
  for(i in 1:N) {
    pred[i] = x[i] * beta[i];# + x0;   
  }
 
}
model {
  #x0 ~ normal(0,10);
  sigma_obs ~ cauchy(0,5);
  for(k in 1:K) {
    beta0[k] ~ normal(0,1);
    sigma_process[k] ~ cauchy(0,5);
    pro_dev[k] ~ normal(0, sigma_process[k]);
  }
  if(family==1) y ~ normal(pred, sigma_obs);
  if(family==2) y_int ~ bernoulli_logit(pred);
  if(family==3) y_int ~ poisson_log(pred);
  if(family==4) y ~ gamma(sigma_obs, sigma_obs ./ exp(pred));
  if(family==5) y ~ lognormal(pred, sigma_obs); 
}
generated quantities {
  vector[N] log_lik;
  # regresssion example in loo() package 
  # regresssion example in loo() package 
  if(family==1) for (n in 1:N) log_lik[n] = normal_lpdf(y[n] | pred[n], sigma_obs);
  if(family==2) for (n in 1:N) log_lik[n] = bernoulli_lpmf(y_int[n] | inv_logit(pred[n]));
  if(family==3) for (n in 1:N) log_lik[n] = poisson_lpmf(y_int[n] | exp(pred[n]));
  if(family==4) for (n in 1:N) log_lik[n] = gamma_lpdf(y[n] | sigma_obs, sigma_obs ./ exp(pred[n]));
  if(family==5) for (n in 1:N) log_lik[n] = lognormal_lpdf(y[n] | pred[n], sigma_obs); 
}
