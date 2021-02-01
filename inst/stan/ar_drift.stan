data {
  int<lower=0> N;
  real y[N];
  int<lower=1> P;
}
parameters {
  real<lower=0> sigma;  // outcome noise
  real<lower=-0.999,upper=0.999> phi[P];
  real mu;
}
transformed parameters {
  real pred[N];
  for(i in 1:P) {
    pred[i] = 0;
  }
  for(i in (P+1):N) {
    pred[i] = mu;
    for(j in 1:P) {
      pred[i] = pred[i] + phi[j]*y[i-j];
    }
  }
}
model {
  y[P:N] ~ normal(pred[P:N], sigma);
  sigma ~ cauchy(0, 5);
  phi ~ normal(0,1);
  mu ~ normal(0,1);
}

