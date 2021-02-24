data {
  int<lower=2> N;  // number of observations
  vector[N] y;     // observation at time T
  int<lower=0> est_nu; 
}
parameters {
  real mu;              // mean
  real<lower=0> sigma;  // error scale
  real theta;      // lag coefficients
  real<lower=2> nu[est_nu];   
}
transformed parameters {
  vector[N] epsilon;    // error terms
  vector[N] pred;    // predictions
  epsilon[1] = y[1] - mu;
  pred[1] = y[1];
  for (t in 2:N) {
    epsilon[t] = (y[t] - mu - theta * epsilon[t - 1]);
    pred[t] = mu + theta * epsilon[t - 1];
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
    for (t in 2:N) {
      y[t] ~ normal(pred[t],sigma);
    }
  } else {
    for (t in 2:N) {
      y[t] ~ student_t(nu, pred[t],sigma);
    } 
  }
}
