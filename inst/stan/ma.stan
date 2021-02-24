data {
  int<lower=0> Q;  // num previous noise terms
  int<lower=3> N;  // num observations
  vector[N] y;     // observation at time t
  int<lower=0> est_nu; 
}
parameters {
  real mu;              // mean
  real<lower=0> sigma;  // error scale
  vector[Q] theta;      // error coeff, lag -t
  real<lower=2> nu[est_nu];    
}
transformed parameters {
  vector[N] epsilon;    // error term at time t
  vector[N] pred;
  for (t in 1:N) {
    epsilon[t] = y[t] - mu;
    for (q in 1:min(t - 1, Q))
      epsilon[t] = epsilon[t] - theta[q] * epsilon[t - q];
  }
  for (t in 1:N) {
    pred[t] = mu;
    for (q in 1:min(t - 1, Q))
      pred[t] = pred[t] + theta[q] * epsilon[t - q];
  }
}
model {
  if(est_nu==1) {
    nu ~ student_t(3,2,2);
  }    
  mu ~ student_t(3,0,2);
  theta ~ student_t(3,0,2);
  sigma ~ student_t(3,0,2);
  if(est_nu==0) {
    y ~ normal(pred, sigma);
  } else {
    y ~ student_t(nu, pred, sigma);
  }
}
