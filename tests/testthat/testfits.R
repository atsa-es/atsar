test_that("fitmodels", {
  
  y = c(1, 5, 3, 12, 5, 6, 6, 2, 8, 14)
  set.seed(123)
  mod = lm(y~1)
  #mod = statss::fit_stan(y, model_name="ar", est_drift = FALSE, mcmc_list=list(n_chain = 1, n_mcmc=300, n_burn=100, n_thin=1) )
  #testthat::expect_equal(mean(rstan::extract(mod)$lp__), -23.75921, tol=10)

  #mod = statss::fit_stan(y, model="ar", est_drift = TRUE, mcmc_list=list(n_chain = 1, n_mcmc=100, n_burn=50, n_thin=1) )
  #testthat::expect_equal(mean(rstan::extract(mod)$lp__), -23.18948, tol=0.1)
  
  #mod = statss::fit_stan(y, model="rw", mcmc_list=list(n_chain = 1, n_mcmc=100, n_burn=50, n_thin=1) )
  #testthat::expect_equal(mean(rstan::extract(mod)$lp__), -18.68706, tol=0.1)
  
  #mod = statss::fit_stan(y, model="rw", est_drift = TRUE, mcmc_list=list(n_chain = 1, n_mcmc=100, n_burn=50, n_thin=1) )
  #testthat::expect_equal(mean(rstan::extract(mod)$lp__), -19.05754, tol=0.1)
  
})