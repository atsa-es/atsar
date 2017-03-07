data {
  int<lower=0> Q;  // num previous noise terms
  int<lower=3> N;  // num observations
  vector[N] y;     // observation at time t
}
parameters {
  real mu;              // mean
  real<lower=0> sigma;  // error scale
  vector[Q] theta;      // error coeff, lag -t
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
  mu ~ cauchy(0, 2.5);
  theta ~ cauchy(0, 2.5);
  sigma ~ cauchy(0, 2.5);
  y ~ normal(pred, sigma);
}
