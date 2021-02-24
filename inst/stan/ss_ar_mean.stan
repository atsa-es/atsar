data {
  int<lower=0> N;
  int<lower=0> n_pos;
  vector[n_pos] y;
  int y_int[n_pos];  
  int pos_indx[n_pos+2];
  int<lower=0> est_nu;  
  int family; // 1 = normal, 2 = binomial, 3 = poisson, 4 = gamma, 5 = lognormal  
}
parameters {
  real x0;
  real<lower=-0.999,upper=0.999> phi;
  real mu;
  vector[N-1] pro_dev;
  real<lower=0> sigma_process;
  real<lower=0> sigma_obs;
  real<lower=2> nu[est_nu];   
}
transformed parameters {
  vector[N] pred;
  pred[1] = x0;
  for(i in 2:N) {
    pred[i] = mu + phi*(pred[i-1] - mu) + sigma_process*pro_dev[i-1];
  }
}
model {
  x0 ~ normal(0,10);
  phi ~ normal(0,10);
  mu ~ normal(0,10);
  if(est_nu==1) {
    nu ~ student_t(3,2,2);
  }      
  sigma_process ~ student_t(3,0,2);
  sigma_obs ~ student_t(3,0,2);
  if(est_nu==0) {
    pro_dev ~ std_normal();//normal(0, sigma_process);
  } else {
    pro_dev ~ student_t(nu,0,1);
  }

  if(family==1) {
    for(i in 1:(n_pos)) {
      //y ~ normal(pred, sigma_obs);
      y[i] ~ normal(pred[pos_indx[i]], sigma_obs);
    }
  }
  if(family==2) {
    for(i in 1:(n_pos)) {
      y_int[i] ~ bernoulli_logit(pred[pos_indx[i]]);
    }
  }
  if(family==3) {
    for(i in 1:(n_pos)) {
      y_int[i] ~ poisson_log(pred[pos_indx[i]]);
    }
  }
  if(family==4) {
    for(i in 1:(n_pos)) {
      y[i] ~ gamma(sigma_obs, sigma_obs ./ exp(pred[pos_indx[i]]));
    }
  }
  if(family==5) {
    for(i in 1:(n_pos)) {
      y[i] ~ lognormal(pred[pos_indx[i]], sigma_obs);
    }
  }
}
generated quantities {
  vector[n_pos] log_lik;
  // regresssion example in loo() package
  if(family==1) for (n in 1:n_pos) log_lik[n] = normal_lpdf(y[n] | pred[pos_indx[n]], sigma_obs);
  if(family==2) for (n in 1:n_pos) log_lik[n] = bernoulli_lpmf(y_int[n] | inv_logit(pred[pos_indx[n]]));
  if(family==3) for (n in 1:n_pos) log_lik[n] = poisson_lpmf(y_int[n] | exp(pred[pos_indx[n]]));
  if(family==4) for (n in 1:n_pos) log_lik[n] = gamma_lpdf(y[n] | sigma_obs, sigma_obs ./ exp(pred[pos_indx[n]]));
  if(family==5) for (n in 1:n_pos) log_lik[n] = lognormal_lpdf(y[n] | pred[pos_indx[n]], sigma_obs);
}
