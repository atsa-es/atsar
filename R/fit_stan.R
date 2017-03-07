#' fit_stan is the primary function which calls pre-written stan scripts for time series data.
#'
#' @param y The response variable (numeric)
#' @param x The predictors, either a vector or matrix
#' @param model_name The specific name of the model to be fitted. Currently supported are 'regression', 'ar', 'rw', 'ma', 'ss_ar' (state space univariate AR), or 'ss_rw' (state space univariate random walk).
#' @param est_drift Whether or not to estimate a drift parameter (default = FALSE). Only applicable to the rw and ar models.
#' @param est_mean Whether to estimate a mean or not (for state space autoregressive model only)
#' @param P The order of the ar model, with minimum value = 1 (default).
#' @param Q The order of the ma model, with minimum value = 1 (default).
#' @param mcmc_list A list of MCMC control parameters. These include the number of 'iterations' (default = 1000), burn in or warmup (default = 500), chains (default = 3), and thinning (default = 1)
#' @param family A named distribution for the observation model, defaults to gaussian
#' @param marss A named list containing the following elements for specifying marss models: (states=NULL, obsVariances=NULL, proVariances=NULL, trends=NULL)
#' 
#' @return an object of class 'rstan'
#' @export
#'
fit_stan <- function(y, x=NA, model_name = NA, est_drift = FALSE, est_mean = FALSE, P = 1, Q = 1, mcmc_list = list(n_mcmc = 1000, n_burn = 500, n_chain = 3, n_thin = 1), family="gaussian", marss = list(states=NULL, obsVariances=NULL, proVariances=NULL, trends=NULL)) {
  stan_dir = find.package("statss")
  
  dist = c("gaussian", "binomial", "poisson", "gamma", "lognormal")
  family = which(dist==family)
  
  if(model_name == "regression") {
    if(class(x)!="matrix") x = matrix(x,ncol=1)
    mod = rstan::stan(paste0(stan_dir, "/exec/regression.stan"), data = list("N"=length(y),"K"=dim(x)[2],"x"=x,"y"=y,"y_int"=round(y), "family"=family),
      pars = c("beta","sigma","pred","log_lik"), chains = mcmc_list$n_chain, iter = mcmc_list$n_mcmc, thin = mcmc_list$n_thin)
  }
  if(model_name == "regression_cor") {
    if(class(x)!="matrix") x = matrix(x,ncol=1)
    mod = rstan::stan(paste0(stan_dir,"/exec/regression_cor.stan"), data = list("N"=length(y),"K"=dim(x)[2],"x"=x,"y"=y,"y_int"=round(y), "family"=family),
      pars = c("beta","sigma","pred","phi","sigma_cor","log_lik"), chains = mcmc_list$n_chain, iter = mcmc_list$n_mcmc, thin = mcmc_list$n_thin)
  }
  if(model_name == "rw" & est_drift == FALSE) {
    mod = rstan::stan(paste0(stan_dir,"/exec/rw.stan"), data = list("y"=y,"N"=length(y)), pars = c("sigma","pred"),
      chains = mcmc_list$n_chain, iter = mcmc_list$n_mcmc, thin = mcmc_list$n_thin)
  }
  if(model_name == "rw" & est_drift == TRUE) {
    mod = rstan::stan(paste0(stan_dir,"/exec/rw_drift.stan"), data = list("y"=y,"N"=length(y)), pars = c("sigma","pred","mu"),
      chains = mcmc_list$n_chain, iter = mcmc_list$n_mcmc, thin = mcmc_list$n_thin)
  }
  if(model_name == "ar" & est_drift == FALSE) {
    mod = rstan::stan(paste0(stan_dir,"/exec/ar.stan"), data = list("y"=y,"N"=length(y)), pars = c("sigma","pred","phi"),
      chains = mcmc_list$n_chain, iter = mcmc_list$n_mcmc, thin = mcmc_list$n_thin)
  }
  if(model_name == "ar" & est_drift == TRUE) {
    mod = rstan::stan(paste0(stan_dir,"/exec/ar_drift.stan"), data = list("y"=y,"N"=length(y)), pars = c("sigma","pred","mu","phi"),
      chains = mcmc_list$n_chain, iter = mcmc_list$n_mcmc, thin = mcmc_list$n_thin)
  }
  if(model_name == "ma" & Q == 1) {
    mod = rstan::stan(paste0(stan_dir,"/exec/ma1.stan"), data = list("y"=y,"N"=length(y)), pars = c("sigma","pred","mu","theta"),
      chains = mcmc_list$n_chain, iter = mcmc_list$n_mcmc, thin = mcmc_list$n_thin)
  }
  if(model_name == "ma" & Q > 1) {
    mod = rstan::stan(paste0(stan_dir,"/exec/ma.stan"), data = list("Q"=Q,"y"=y,"N"=length(y)), pars = c("sigma","pred","mu","theta"),
      chains = mcmc_list$n_chain, iter = mcmc_list$n_mcmc, thin = mcmc_list$n_thin)
  }
  if(model_name == "ss_rw" & est_drift == FALSE) {
    mod = rstan::stan(paste0(stan_dir,"/exec/ss_rw.stan"), data = list("y"=y,"N"=length(y)), pars = c("sigma_process","pred", "sigma_obs"),
      chains = mcmc_list$n_chain, iter = mcmc_list$n_mcmc, thin = mcmc_list$n_thin)
  }
  if(model_name == "ss_rw" & est_drift == TRUE) {
    mod = rstan::stan(paste0(stan_dir,"/exec/ss_rw_drift.stan"), data = list("y"=y,"N"=length(y)), pars = c("sigma_process","pred", "sigma_obs", "mu"),
      chains = mcmc_list$n_chain, iter = mcmc_list$n_mcmc, thin = mcmc_list$n_thin)
  }
  if(model_name == "ss_ar" & est_drift == FALSE) {
    mod = rstan::stan(paste0(stan_dir,"/exec/ss_ar.stan"), data = list("y"=y,"N"=length(y)), pars = c("sigma_process","pred", "sigma_obs", "phi"),
      chains = mcmc_list$n_chain, iter = mcmc_list$n_mcmc, thin = mcmc_list$n_thin)
  }
  if(model_name == "ss_ar" & est_drift == TRUE) {
    mod = rstan::stan(paste0(stan_dir,"/exec/ss_ar_drift.stan"), data = list("y"=y,"N"=length(y)), pars = c("sigma_process","pred", "sigma_obs", "mu", "phi"),
      chains = mcmc_list$n_chain, iter = mcmc_list$n_mcmc, thin = mcmc_list$n_thin)
  }
  if(model_name == "ss_ar" & est_mean == TRUE) {
    mod = rstan::stan(paste0(stan_dir,"/exec/ss_ar_mean.stan"), data = list("y"=y,"N"=length(y)), pars = c("sigma_process","pred", "sigma_obs", "mu", "phi"),
      chains = mcmc_list$n_chain, iter = mcmc_list$n_mcmc, thin = mcmc_list$n_thin)
  }
  if(model_name == "arma11") {
    mod = rstan::stan(paste0(stan_dir,"/exec/arma11.stan"), data = list("y"=y,"N"=length(y)), pars = c("sigma", "theta", "mu", "phi"),
      chains = mcmc_list$n_chain, iter = mcmc_list$n_mcmc, thin = mcmc_list$n_thin)
  }
  if(model_name == "dlm-intercept") {
    # constant slope, and time -varying intercept model
    if(is.na(x)) {
      x = matrix(0, nrow=length(y), ncol=1)
    } 
    if(class(x)!="matrix") x = matrix(x,ncol=1)

    mod = rstan::stan(paste0(stan_dir, "/exec/dlm_int.stan"), data = list("N"=length(y),"K"=dim(x)[2],"x"=x,"y"=y,"y_int"=round(y), "family"=family),
      pars = c("beta","sigma_obs","sigma_process","pred","intercept","log_lik"), chains = mcmc_list$n_chain, iter = mcmc_list$n_mcmc, thin = mcmc_list$n_thin)
  }  
  if(model_name == "dlm-slope") {
    # constant estimated intercept, and time varying slopes
    if(class(x)!="matrix") x = matrix(x,ncol=1)
    
    mod = rstan::stan(paste0(stan_dir, "/exec/dlm_slope.stan"), data = list("N"=length(y),"K"=dim(x)[2],"x"=x,"y"=y,"y_int"=round(y), "family"=family),
      pars = c("beta","sigma_obs","sigma_process","pred","log_lik"), chains = mcmc_list$n_chain, iter = mcmc_list$n_mcmc, thin = mcmc_list$n_thin)
  }    
  if(model_name == "dlm") {
    # this is just a time-varying model with time varying intercept and slopes
    if(class(x)!="matrix") x = matrix(x,ncol=1)
    
    mod = rstan::stan(paste0(stan_dir, "/exec/dlm.stan"), data = list("N"=length(y),"K"=dim(x)[2],"x"=x,"y"=y,"y_int"=round(y), "family"=family),
      pars = c("beta","sigma_obs","sigma_process","pred","log_lik"), chains = mcmc_list$n_chain, iter = mcmc_list$n_mcmc, thin = mcmc_list$n_thin)
  }     
  if(model_name == "marss") {
    if(is.null(marss$states)) states = rep(1, nrow(y))
    if(is.null(marss$obsVariances)) obsVariances = rep(1, nrow(y))
    if(is.null(marss$proVariances)) proVariances = 1
    if(is.null(marss$trends)) trends = 1
    proVariances = c(proVariances, 0) # to keep types in stan constant
    trends = c(trends, 0)     # to keep types in stan constant
    N = ncol(y)
    M = nrow(y)
    row_indx_pos = matrix((rep(1:M, N)), M, N)[which(!is.na(y))]
    col_indx_pos = matrix(sort(rep(1:N, M)), M, N)[which(!is.na(y))]
    n_pos = length(row_indx_pos)
    y = y[which(!is.na(y))]
    
    mod = rstan::stan(paste0(stan_dir, "/exec/marss.stan"), 
      data = list("N"=nrow(y),"M"=ncol(y), "y"=y,
      "states"=states, "S" = max(states), "obsVariances"=obsVariances,
      "n_obsvar" = max(obsVariances), "proVariances" = proVariances, 
        "n_provar" = max(proVariances),
      "trends"=trends, "n_trends" = max(trends),
        "n_pos" = n_pos,
        "col_indx_pos" = col_indx_pos,
        "row_indx_pos" = row_indx_pos,
        "y_int"=round(y), 
        "family"=family),
      pars = c("pred","log_lik"), chains = mcmc_list$n_chain, 
      iter = mcmc_list$n_mcmc, thin = mcmc_list$n_thin)
  }      
  return(mod)
}
