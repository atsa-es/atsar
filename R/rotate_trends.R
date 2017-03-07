#' rotate_trends is the primary function which calls pre-written stan scripts for time series data.
#'
#' @param fitted_model A fitted stanfit object
#' 
#' @return a list object, with the rotated trends for each MCMC chain (trends), rotated Z for every MCMC chain (Z_rot), mean Z (Z_rot_mean), mean trends (trends_mean), lower 2.5% interval on trends (trends_lower), upper 97.5% interval for trends (trends_upper)
#' @export
#'
rotate_trends = function(fitted_model) {
  # Illustrate how to get the trends out of the model
  # get the inverse of the rotation matrix
  n_mcmc = dim(rstan::extract(fitted_model)$Z)[1]
  Z = rstan::extract(fitted_model)$Z
  x = rstan::extract(fitted_model)$x
  n_ts = dim(Z)[2]
  n_trends = dim(x)[2]
  n_years = dim(x)[3]
  mcmc_trends_rot = array(0, dim = c(n_mcmc, n_trends, n_years))
  mcmc_Z_rot = array(0, dim = c(n_mcmc, n_ts, n_trends))
  for(i in 1:n_mcmc) {
    Zest = Z[i,,]
    H.inv = stats::varimax(Zest)$rotmat
  
    # rotate factor loadings
    Z.rot = Zest %*% H.inv
    mcmc_Z_rot[i,,] = Z.rot
    
    # rotate trends
    states = x[i,,]
    trends.rot = solve(H.inv) %*% states
    mcmc_trends_rot[i,,] = trends.rot
  }
  return(list("Z_rot"=mcmc_Z_rot, "trends"=mcmc_trends_rot,
    "Z_rot_mean" = apply(mcmc_Z_rot,c(2,3),mean),
    "trends_mean" = apply(mcmc_trends_rot,c(2,3),mean),
  "trends_lower" = apply(mcmc_trends_rot,c(2,3),stats::quantile,0.025),
"trends_upper" = apply(mcmc_trends_rot,c(2,3),stats::quantile,0.975)))
}