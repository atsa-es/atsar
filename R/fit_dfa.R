#' fit_dfa is the primary function for fitting dynamic factor analysis to time series data.
#'
#' @param y The response variable (numeric matrix)
#' @param covar The optional matrix of covariates (time on columns)
#' @param covar_index The matrix of index values, indicating which covariate effects on time series are shared. 
#' @param num_trends The number of DFA trends (random walks) to estimate. Process variance fixed at 1.
#' @param varIndx The vector indicating which time series share variances. By default, (1,2,3,...,n)
#' @param zscore Whether or not to z-score data. Defaults to true.
#' @param iter Number of MCMC iterations, defaults to 4000.
#' @param chains Number of MCMC chains, defaults to 1 (for DFA this is the best choice)
#' @param control Default control list for stanfit objects.
#'
#' @return an object of class 'rstan'
#' @export
#'
fit_dfa <- function(y = y, 
  covar=NULL,
  covar_index=NULL,
  num_trends = 2,
  varIndx = NULL,
  zscore = TRUE,
  iter = 4000,
  chains = 1,
  control = list(adapt_delta = 0.99)) {
  
  stan_dir = find.package("statss")
  model = paste0(stan_dir, "/exec/dfa.stan")
  
  # parameters for DFA
  N = ncol(y)
  P = nrow(y)
  K = num_trends # number of dfa trends
  nZ = P * K - sum(1:K) + K
  d_covar = covar;
  num_covar = nrow(d_covar)
  covar_indexing = covar_index
  if(!is.null(d_covar) & !is.null(covar_indexing)) {
    num_unique_covar = max(covar_indexing)
  } 
  if(!is.null(d_covar) & is.null(covar_indexing)) {
    # covariates included but index matrix not, assume independent for all elements
    covar_indexing = matrix(seq(1,num_covar*P),P,num_covar)
    num_unique_covar = max(covar_indexing)
  }
  if(is.null(d_covar)) {
    covar_indexing = matrix(0,P,0)
    d_covar = matrix(0,0,N)
    num_covar = 0
    num_unique_covar = 0
  }
  
  if (zscore == TRUE) {
    for (i in 1:P) {
      y[i, ] = scale(y[i, ], center = TRUE, scale = TRUE)
    }
  }
  
  mat_indx = matrix(0, P, K)
  start = 1
  for (k in 1:K) {
    if (k == 1)
      mat_indx[, k] = (1:nZ)[start:(start + P - k)]
    if (k > 1)
      mat_indx[-c(0:(k - 1)), k] = (1:nZ)[start:(start + P - k)]
    start = start + (P - k + 1)
  }
  
  row_indx = matrix((rep(1:P, K)), P, K)[which(mat_indx > 0)]
  col_indx = rep(1:K, times = P:(P - K + 1))
  row_indx_z = matrix((rep(1:P, K)), P, K)[which(mat_indx == 0)]
  col_indx_z = matrix(sort(rep(1:K, P)), P, K)[which(mat_indx == 0)]
  row_indx_z = c(row_indx_z, 0, 0)# +2 zeros for making stan ok with data types
  col_indx_z = c(col_indx_z, 0, 0)# +2 zeros for making stan ok with data types
  nZero = length(row_indx_z)
  
  # set the model up to have shared variances between first two time series,
  # third is different
  if (is.null(varIndx))
    varIndx = rep(1, P)
  nVariances = length(unique(varIndx))
  
  # indices of positive values - stan can't handle NAs
  row_indx_pos = matrix((rep(1:P, N)), P, N)[which(!is.na(y))]
  col_indx_pos = matrix(sort(rep(1:N, P)), P, N)[which(!is.na(y))]
  n_pos = length(row_indx_pos)
  y = y[which(!is.na(y))]
  
  data_list = list(
    N,
    P,
    K,
    nZ,
    y,
    row_indx,
    col_indx,
    nZero,
    row_indx_z,
    col_indx_z,
    nZero,
    row_indx_z,
    col_indx_z,
    row_indx_pos,
    col_indx_pos,
    n_pos,
    d_covar,
    num_covar,
    covar_indexing,
    num_unique_covar
  )
  pars <- c("x", "Z", "sigma", "log_lik", "pred")
  if(!is.null(covar)) pars = c(pars, "D")
  mod = rstan::stan(
    data = data_list,
    pars = pars,
    file = model[[1]],
    chains = chains,
    iter = iter,
    control = control
  )
  return(mod)
}
