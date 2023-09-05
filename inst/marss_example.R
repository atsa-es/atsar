y = matrix(0, 4, 30)
for(i in 1:nrow(y)) y[i,] = cumsum(rnorm(ncol(y)))


mcmc_list = list(n_mcmc = 100, n_burn = 50, n_chain = 3, n_thin = 1)
         
marss = list(states = c(1,2,3,4), obsVariances = c(1,2,2,1), 
             proVariances = c(1,2,3,3), trends = c(1,2,3,4),
             est_trend = FALSE, est_B = TRUE)

fit_stan(y = y, model_name = "marss", marss = marss,
         mcmc_list = mcmc_list)


# Or if you want to run the code yourself you could run this part below
marss = list(states = c(1,2,3,3), obsVariances = c(1,2,2,1), 
             proVariances = c(1,2,3), trends = c(1,2,3),
             est_trend = TRUE, est_B = TRUE)

# Don't change anything below here
if (is.null(marss$states)) marss$states <- rep(1, nrow(y))
if(length(marss$states) != nrow(y)) stop("Error: state vector must be same length as number of time series in y")
if (is.null(marss$obsVariances)) marss$obsVariances <- rep(1, nrow(y))
if(length(marss$obsVariances) != nrow(y)) stop("Error: vector of observation error variances must be same length as number of time series in y")
if (is.null(marss$proVariances)) marss$proVariances <- rep(1, max(marss$states))
if(length(marss$proVariances) < max(marss$states)) stop("Error: vector of process error variances is fewer than the number of states")
if(length(marss$proVariances) > max(marss$states)) stop("Error: vector of process error variances is larger than the number of states")
if (is.null(marss$trends)) marss$trends <- rep(1, max(marss$states))
if(length(marss$trends) < max(marss$states)) stop("Error: vector of trends is fewer than the number of states")
if(length(marss$trends) > max(marss$states)) stop("Error: vector of trends is larger than the number of states")
if (marss$est_trend == FALSE) est_trend = FALSE
if (marss$est_B == FALSE) est_B = FALSE
proVariances <- c(marss$proVariances, 0) # to keep types in stan constant
trends <- c(marss$trends, 0) # to keep types in stan constant
N <- ncol(y)
M <- nrow(y)
row_indx_pos <- matrix((rep(1:M, N)), M, N)[which(!is.na(y))]
col_indx_pos <- matrix(sort(rep(1:N, M)), M, N)[which(!is.na(y))]
n_pos <- length(row_indx_pos)
y <- y[which(!is.na(y))]

data_list = list("N"=N,"M"=M, "y"=y,
            "states"=states, "S" = max(marss$states), "obsVariances"=marss$obsVariances,
            "n_obsvar" = max(marss$obsVariances), "proVariances" = proVariances,
            "n_provar" = max(proVariances),
            "trends"=trends, "n_trends" = max(trends),
            "n_pos" = n_pos,
            "col_indx_pos" = col_indx_pos,
            "row_indx_pos" = row_indx_pos,
            "y_int"=round(y),
            "family"=1,
            "est_trend" = as.numeric(est_trend),
            "est_B" = as.numeric(est_B))

mod <- rstan::stan(paste0("inst/stan/marss.stan"),
data = data_list, chains = mcmc_list$n_chain,
iter = mcmc_list$n_mcmc, thin = mcmc_list$n_thin)

pars = rstan::extract(mod)