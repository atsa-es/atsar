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
#' @param est_nu Boolean, whether to model process deviations as Student-t or not (default).
#' @param marss A named list containing the following elements for specifying marss models: (states=NULL, obsVariances=NULL, proVariances=NULL, trends=NULL
#' @param map_estimation Whether to do maximum a posteriori estimation via [rstan::optimizing()] (defualts to FALSE)
#' @param hessian Whether to return hessian if map_estimation is TRUE via [rstan::optimizing()]
#' @param ... Any other arguments passed to [rstan::sampling()].
#' @return an object of class 'rstan'
#' @importFrom rstan sampling
#' @export
#'
fit_stan <- function(y, x = NA, model_name = NA,
                     est_drift = FALSE,
                     est_mean = FALSE,
                     P = 1,
                     Q = 1,
                     mcmc_list = list(n_mcmc = 1000, n_burn = 500, n_chain = 3, n_thin = 1),
                     family = "gaussian",
                     est_nu = FALSE,
                     marss = list(states = NULL, obsVariances = NULL, proVariances = NULL, trends = NULL),
                     map_estimation = FALSE,
                     hessian = FALSE, ...) {
  dist <- c("gaussian", "binomial", "poisson", "gamma", "lognormal")
  family <- which(dist == family)

  # process potential NAs in data
  if (!is.matrix(y)) {
    N <- length(y)
    pos_indx <- which(!is.na(y))
    y <- y[pos_indx]
    n_pos <- length(pos_indx)
    # catch case where -- needs to be 2 elements for stan vec
    if (length(pos_indx) == 0) {
      pos_indx <- rep(0, 2)
    } else {
      pos_indx <- c(pos_indx, 0, 0)
    }

    if (length(pos_indx) != (length(y) + 2)) {
      # include any other errors?
      if (model_name %in% c("regression", "regression_cor", "ar", "rw", "ma", "arma11")) {
        stop("Error: data cannot contain NAs for the specified model")
      }
    }
  }

  data <- NA
  if (model_name == "regression") {
    if (is.matrix(x) == FALSE) x <- matrix(x, ncol = 1)
    object <- stanmodels$regression
    data <- list("N" = length(y), "K" = dim(x)[2], "x" = x, "y" = y, "y_int" = round(y), "family" = family)
    pars <- c("beta", "sigma", "pred", "log_lik")
  }
  if (model_name == "regression_cor") {
    if (is.matrix(x) == FALSE) x <- matrix(x, ncol = 1)
    object <- stanmodels$regression_cor
    data <- list("N" = length(y), "K" = dim(x)[2], "x" = x, "y" = y, "y_int" = round(y), "family" = family)
    pars <- c("beta", "sigma", "pred", "phi", "sigma_cor", "log_lik")
  }
  if (model_name == "rw") {
    object <- stanmodels$rw
    data <- list(
      "y" = y, "N" = length(y), est_drift = ifelse(est_drift == FALSE, 0, 1),
      est_nu = ifelse(est_nu == FALSE, 0, 1)
    )
    pars <- c("sigma", "pred")
    if (est_drift == TRUE) pars <- c(pars, "mu")
    if (est_nu == TRUE) pars <- c(pars, "nu")
  }
  if (model_name == "ar") {
    object <- stanmodels$ar
    data <- list(
      "y" = y, "N" = length(y), est_drift = ifelse(est_drift == FALSE, 0, 1),
      est_nu = ifelse(est_nu == FALSE, 0, 1)
    )
    pars <- c("sigma", "pred", "phi")
    if (est_drift == TRUE) pars <- c(pars, "mu")
    if (est_nu == TRUE) pars <- c(pars, "nu")
  }
  if (model_name == "ma" & Q == 1) {
    object <- stanmodels$ma1
    data <- list("y" = y, "N" = length(y), est_nu = ifelse(est_nu == FALSE, 0, 1))
    pars <- c("sigma", "pred", "mu", "theta")
    if (est_nu == TRUE) pars <- c(pars, "nu")
  }
  if (model_name == "ma" & Q > 1) {
    object <- stanmodels$ma
    data <- list("Q" = Q, "y" = y, "N" = length(y), est_nu = ifelse(est_nu == FALSE, 0, 1))
    pars <- c("sigma", "pred", "mu", "theta")
    if (est_nu == TRUE) pars <- c(pars, "nu")
  }
  if (model_name == "ss_rw") {
    object <- stanmodels$ss_rw
    data <- list(
      "y" = y, "N" = N, "n_pos" = n_pos, "pos_indx" = pos_indx, "y_int" = round(y),
      est_drift = ifelse(est_drift == FALSE, 0, 1),
      "family" = family,
      est_nu = ifelse(est_nu == FALSE, 0, 1)
    )
    pars <- c("sigma_process", "pred", "sigma_obs")
    if (est_drift == TRUE) pars <- c(pars, "mu")
    if (est_nu == TRUE) pars <- c(pars, "nu")
  }
  if (model_name == "ss_ar" & est_mean == FALSE) {
    object <- stanmodels$ss_ar
    data <- list(
      "y" = y, "N" = N, "n_pos" = n_pos, "pos_indx" = pos_indx, "y_int" = round(y),
      est_drift = ifelse(est_drift == FALSE, 0, 1),
      "family" = family,
      est_nu = ifelse(est_nu == FALSE, 0, 1)
    )
    pars <- c("sigma_process", "pred", "sigma_obs", "phi")
    if (est_drift == TRUE) pars <- c(pars, "mu")
    if (est_nu == TRUE) pars <- c(pars, "nu")
  }
  if (model_name == "ss_ar" & est_mean == TRUE) {
    object <- stanmodels$ss_ar_mean
    data <- list(
      "y" = y, "N" = N, "n_pos" = n_pos, "pos_indx" = pos_indx, "y_int" = round(y),
      "family" = family,
      est_nu = ifelse(est_nu == FALSE, 0, 1)
    )
    pars <- c("sigma_process", "pred", "sigma_obs", "mu", "phi")
    if (est_nu == TRUE) pars <- c(pars, "nu")
  }
  # if(model_name == "arma11") {
  #   object <- stanmodels$arma11
  #   data <- list("y"=y,"N"=length(y))
  #   pars <- c("sigma", "theta", "mu", "phi")
  # }
  if (model_name == "dlm-intercept") {
    object <- stanmodels$dlm_int
    # constant slope, and time -varying intercept model
    if (is.na(x)) {
      x <- matrix(0, nrow = length(y), ncol = 1)
    }
    if (is.matrix(x) == FALSE) x <- matrix(x, ncol = 1)
    data <- list("N" = N, "K" = dim(x)[2], "x" = x, "y" = y, "y_int" = round(y), "family" = family, "n_pos" = n_pos, "pos_indx" = pos_indx)
    pars <- c("beta", "sigma_obs", "sigma_process", "pred", "intercept", "log_lik")
  }
  if (model_name == "dlm-slope") {
    object <- stanmodels$dlm_slope
    # constant estimated intercept, and time varying slopes
    if (is.matrix(x) == FALSE) x <- matrix(x, ncol = 1)
    data <- list("N" = N, "K" = dim(x)[2], "x" = x, "y" = y, "y_int" = round(y), "family" = family, "n_pos" = n_pos, "pos_indx" = pos_indx)
    pars <- c("beta", "sigma_obs", "sigma_process", "pred", "log_lik")
  }
  if (model_name == "dlm") {
    object <- stanmodels$dlm
    # this is just a time-varying model with time varying intercept and slopes
    if (is.matrix(x) == FALSE) x <- matrix(x, ncol = 1)
    data <- list("N" = N, "K" = dim(x)[2], "x" = x, "y" = y, "y_int" = round(y), "family" = family, "n_pos" = n_pos, "pos_indx" = pos_indx)
    pars <- c("beta", "sigma_obs", "sigma_process", "pred", "log_lik")
  }
  if (model_name == "marss") {
    if (is.null(marss$states)) states <- rep(1, nrow(y))
    if (is.null(marss$obsVariances)) obsVariances <- rep(1, nrow(y))
    if (is.null(marss$proVariances)) proVariances <- rep(1, max(states))
    if (is.null(marss$trends)) trends <- rep(1, max(states))
    proVariances <- c(proVariances, 0) # to keep types in stan constant
    trends <- c(trends, 0) # to keep types in stan constant
    N <- ncol(y)
    M <- nrow(y)
    row_indx_pos <- matrix((rep(1:M, N)), M, N)[which(!is.na(y))]
    col_indx_pos <- matrix(sort(rep(1:N, M)), M, N)[which(!is.na(y))]
    n_pos <- length(row_indx_pos)
    y <- y[which(!is.na(y))]

    # mod <- rstan::stan(paste0(stan_dir, "/exec/marss.stan"),
    #   data = list("N"=nrow(y),"M"=ncol(y), "y"=y,
    #   "states"=states, "S" = max(states), "obsVariances"=obsVariances,
    #   "n_obsvar" = max(obsVariances), "proVariances" = proVariances,
    #     "n_provar" = max(proVariances),
    #   "trends"=trends, "n_trends" = max(trends),
    #     "n_pos" = n_pos,
    #     "col_indx_pos" = col_indx_pos,
    #     "row_indx_pos" = row_indx_pos,
    #     "y_int"=round(y),
    #     "family"=family),
    #   pars = c("pred","log_lik"), chains = mcmc_list$n_chain,
    #   iter = mcmc_list$n_mcmc, thin = mcmc_list$n_thin)
  }
  if (map_estimation == FALSE) {
    out <- rstan::sampling(
      object = object,
      data = data,
      pars = pars,
      control = list(adapt_delta = 0.8, max_treedepth = 20),
      warmup = mcmc_list$n_burn,
      iter = mcmc_list$n_mcmc,
      thin = mcmc_list$n_thin,
      chains = mcmc_list$n_chain, ...
    )
  } else {
    out <- rstan::optimizing(
      object = object,
      data = data,
      hessian = hessian, ...
    )
  }

  return(out)
}
