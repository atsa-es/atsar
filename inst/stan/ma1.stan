data {
  int<lower=2> N;  // number of observations
  vector[N] y;     // observation at time T
}
parameters {
  real mu;              // mean
  real<lower=0> sigma;  // error scale
  real theta;      // lag coefficients
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
  mu ~ student_t(3,0,2);
  theta ~ student_t(3,0,2);
  sigma ~ student_t(3,0,2);
  for (t in 2:N) {
    y[t] ~ normal(pred[t],sigma);
  }
}
