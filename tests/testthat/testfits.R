test_that("fitmodels", {
  
  set.seed(123)
  y <- cumsum(rnorm(30))
  
  set.seed(123)
  mod <- fit_stan(y, model_name = "ar", est_drift = TRUE, mcmc_list=list(n_chain = 1, n_mcmc=300, n_burn=100, n_thin=1) )
  pars <- rstan::extract(mod)
  expect_equal(mean(pars$sigma), 1.028974, tolerance = 0.01)
  
  set.seed(123)
  mod <- fit_stan(y, model_name = "ar", est_drift = FALSE, mcmc_list=list(n_chain = 1, n_mcmc=300, n_burn=100, n_thin=1) )
  pars <- rstan::extract(mod)
  expect_equal(mean(pars$sigma), 1.011565, tolerance = 0.01)
  
  set.seed(123)
  mod <- fit_stan(y, model_name = "rw", est_drift = TRUE, mcmc_list=list(n_chain = 1, n_mcmc=300, n_burn=100, n_thin=1) )
  pars <- rstan::extract(mod)
  expect_equal(mean(pars$sigma), 1.032227, tolerance = 0.01)
  
  set.seed(123)
  mod <- fit_stan(y, model_name = "rw", est_drift = FALSE, mcmc_list=list(n_chain = 1, n_mcmc=300, n_burn=100, n_thin=1) )
  pars <- rstan::extract(mod)
  expect_equal(mean(pars$sigma), 1.03083, tolerance = 0.01)
  
  set.seed(123)
  mod <- fit_stan(y, x = 1:30, model_name = "regression", est_drift = FALSE, mcmc_list=list(n_chain = 1, n_mcmc=300, n_burn=100, n_thin=1) )
  pars <- rstan::extract(mod)
  expect_equal(mean(pars$sigma), 2.18601, tolerance = 0.01)
  
  set.seed(123)
  mod <- fit_stan(y, model_name = "ss_ar", est_drift = TRUE, mcmc_list=list(n_chain = 1, n_mcmc=300, n_burn=100, n_thin=1) )
  pars <- rstan::extract(mod)
  expect_equal(mean(pars$sigma_obs), 0.2898811, tolerance = 0.01)
  
  set.seed(123)
  mod <- fit_stan(y, model_name = "ss_ar", est_drift = FALSE, mcmc_list=list(n_chain = 1, n_mcmc=300, n_burn=100, n_thin=1) )
  pars <- rstan::extract(mod)
  expect_equal(mean(pars$sigma_obs), 0.2503624, tolerance = 0.01)
  
  set.seed(123)
  mod <- fit_stan(y, model_name = "ss_rw", est_drift = TRUE, mcmc_list=list(n_chain = 1, n_mcmc=300, n_burn=100, n_thin=1) )
  pars <- rstan::extract(mod)
  expect_equal(mean(pars$sigma_obs), 0.3714019, tolerance = 0.01)
  
  set.seed(123)
  mod <- fit_stan(y, model_name = "ss_rw", est_drift = FALSE, mcmc_list=list(n_chain = 1, n_mcmc=300, n_burn=100, n_thin=1) )
  pars <- rstan::extract(mod)
  expect_equal(mean(pars$sigma_obs), 0.4040785, tolerance = 0.01)

})
