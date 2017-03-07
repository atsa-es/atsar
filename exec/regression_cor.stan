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
  real<lower=-1,upper=1> phi;
}
transformed parameters {
  vector[N] pred;
  vector[N] epsilon;
  real sigma_cor;
  pred[1] = x[1] * beta;
  epsilon[1] = y[1]-pred[1];
  for(i in 2:N) {
  pred[i] = x[i] * beta;
  epsilon[i] = (y[i] - pred[i]) - phi*epsilon[i-1];
  }
  sigma_cor = sqrt(sigma*sigma * (1-phi*phi)); # Var = sigma2 * (1-rho^2)
}
model {
  phi ~ normal(0,1);
  beta ~ normal(0,2);
  sigma ~ cauchy(0,5);
  if(family==1) y ~ normal(pred, sigma_cor);
  if(family==2) y_int ~ bernoulli_logit(pred);
  if(family==3) y_int ~ poisson_log(pred);
  if(family==4) y ~ gamma(sigma_cor, sigma_cor ./ exp(pred));
  if(family==5) y ~ lognormal(pred, sigma_cor);
}
generated quantities {
  vector[N] log_lik;
  # regresssion example in loo() package 
  if(family==1) for (n in 1:N) log_lik[n] = normal_lpdf(y[n] | pred[n], sigma_cor);
  if(family==2) for (n in 1:N) log_lik[n] = bernoulli_lpmf(y_int[n] | inv_logit(pred[n]));
  if(family==3) for (n in 1:N) log_lik[n] = poisson_lpmf(y_int[n] | exp(pred[n]));
  if(family==4) for (n in 1:N) log_lik[n] = gamma_lpdf(y[n] | sigma_cor, sigma_cor ./ exp(pred[n]));
  if(family==5) for (n in 1:N) log_lik[n] = lognormal_lpdf(y[n] | pred[n], sigma_cor); 
}
