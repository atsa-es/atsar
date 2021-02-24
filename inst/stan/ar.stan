data {
  int<lower=0> N;
  real y[N];
  int<lower=1> P;
  int<lower=0> est_drift;
  int<lower=0> est_nu;
}
parameters {
  real<lower=0> sigma;  // outcome noise
  real<lower=-0.999,upper=0.999> phi[P];
  real mu[est_drift];
  real<lower=2> nu[est_nu];
}
transformed parameters {
  real pred[N];
  for(i in 1:P) {
    pred[i] = 0;
  }
  for(i in (P+1):N) {
    if(est_drift==1) {
      pred[i] = mu[1];
    } else {
      pred[i] = 0;
    }
    for(j in 1:P) {
      pred[i] = pred[i] + phi[j]*y[i-j];
    }
  }
}
model {
  if(est_drift==1) {
    mu ~ normal(0,1);
  }
  if(est_nu == 0) {
    y[P:N] ~ normal(pred[P:N], sigma);
  } else {
    y[P:N] ~ student_t(nu, pred[P:N], sigma);
  }
  sigma ~ student_t(3,0,2);
  phi ~ normal(0,1);
  if(est_nu==1) {
    nu ~ student_t(3,2,2);
  }
}

