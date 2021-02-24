data {
  int<lower=0> N;
  real y[N];
  int<lower=0> est_drift;
  int<lower=0> est_nu;    
}
parameters {
  real<lower=0> sigma;  // outcome noise
  real mu[est_drift];
  real<lower=2> nu[est_nu];  
}
transformed parameters {
  real pred[N];
  real temp;
  pred[1] = 0;
  temp = 0;
  if(est_drift == 1) {
    temp = mu[1];
  }
  for(i in 2:N) {
    pred[i] = y[i-1] + temp;
  }
}
model {
  if(est_nu==1) {
    nu ~ student_t(3,2,2);
  }   
  if(est_nu==0) {
    y[2:N] ~ normal(pred[2:N], sigma);
  } else {
    y[2:N] ~ student_t(nu, pred[2:N], sigma);
  }
  sigma ~ student_t(3,0,2);
  if(est_drift==1) {
    mu ~ normal(0,1);
  }
}
generated quantities {
  vector[N-1] log_lik;
  // regresssion example in loo() package
  for (n in 2:N) log_lik[n-1] = normal_lpdf(y[n] | pred[n], sigma);
}
