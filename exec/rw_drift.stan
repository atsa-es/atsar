data {
  int<lower=0> N;
  real y[N];
}
parameters {
  real<lower=0> sigma;  // outcome noise
  real mu;
}
transformed parameters {
  real pred[N];
  pred[1] = 0;
  for(i in 2:N) {
    pred[i] = y[i-1] + mu;
  }
}
model {
  y[2:N] ~ normal(pred[2:N], sigma);
  sigma ~ cauchy(0, 5);
  mu ~ normal(0,1);
}
generated quantities {
  vector[N-1] log_lik;
  # regresssion example in loo() package 
  for (n in 2:N) log_lik[n-1] = normal_lpdf(y[n] | pred[n], sigma);
}
