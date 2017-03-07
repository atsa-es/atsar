<!-- README.md is generated from README.Rmd. Please edit that file -->
statss (Stan for Analyzing Time Series)
=======================================

The statss R package implements Bayesian time series models, primarily for illustrative purposes and teaching (University of Washington's Fish 507, Winter quarter 2017). You can cite the package as:

#### Citation: Ward, E.J., M.D. Scheuerell, and E.E. Holmes. 2017. 'statss': Stan for Analyzing Time Series: an introduction to time series analysis for ecological and fisheries data. [![DOI](https://zenodo.org/badge/78160922.svg)](https://zenodo.org/badge/latestdoi/78160922)

You can install the development version of the package with:

``` r
# install.packages("devtools")
devtools::install_github("eric-ward/safs-timeseries/statss")
```

An example model
----------------

Simulate data:

``` r
library(rstan)
library(statss)
#> 
#> Attaching package: 'statss'
#> The following object is masked from 'package:STATS':
#> 
#>     fit_stan
set.seed(123)
s = cumsum(rnorm(50))
```

``` r
plot(s)
```

![](README-figs/plot-1.png)

Fit several models to this data:

``` r
# Regression, no slope
regression_model = fit_stan(y = s, x = model.matrix(lm(s~1)), model_name="regression")

# Regression, with slope
regression_model = fit_stan(y = s, x = model.matrix(lm(s~seq(1,length(s)))), model_name="regression")

# AR(1) time series model
ar1_model = fit_stan(y = s, est_drift=FALSE, P = 1, model_name = "ar")

# ARMA(1,1) time series model
arma1_model = fit_stan(y = s, model_name = "arma11")

# univariate ss model -- without drift but mean reversion estimated
ss_model = fit_stan(y = s, model_name = "ss_ar", est_drift=FALSE)
```

References
==========

[Fish 507 class website](https://catalyst.uw.edu/workspace/fish203/35553/243766) ...
