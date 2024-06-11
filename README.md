<!-- badges: start --> [![atsar status
badge](https://atsa-es.r-universe.dev/badges/atsar)](https://atsa-es.r-universe.dev)[![R-CMD-check](https://github.com/nwfsc-timeseries/atsar/workflows/R-CMD-check/badge.svg)](https://github.com/nwfsc-timeseries/atsar/actions)
<!-- badges: end -->

<!-- badges: start --> [![Codecov test
coverage](https://codecov.io/gh/atsa-es/atsar/branch/master/graph/badge.svg)](https://app.codecov.io/gh/atsa-es/atsar?branch=master)
<!-- badges: end -->

<style>
.nav{
    border:1px solid #ccc;
    border-width:1px 0;
    list-style:none;
    margin:0;
    padding:0;
    text-align:center;
}
.nav li{
    display:inline-block;
}
.nav a{
    display:inline-block;
    padding:5px;
}
</style>
<ul class="nav">
<li>
<a href="#install">Install</a>
</li>
<li>
<a href="#documentation">Documentation</a>
</li>
<li>
<a href="#example">Example</a>
</li>
<li>
<a href="#citation">Citation</a>
</li>
<li>
<a href="#license">License</a>
</li>
<li>
<a href="https://github.com/nwfsc-timeseries/atsar">GitHub</a>
</li>
</ul>

The atsar R package implements Bayesian time series models using Stan,
primarily for illustrative purposes and teaching (University of
Washington’s Fish 507, Winter quarters). The Stan webpage, and
appropriate citation guidelines are [here](http://mc-stan.org/).

### INSTALL

# Install atsar from the atsa-es r-universe repository (binaries for Windows and Mac-Intel)

    install.packages('atsar', repos = c('https://atsa-es.r-universe.dev', 'https://cloud.r-project.org'))

You can build the development version of the package from the source
here. Note you need to use this if you have a M1/M2 Mac.

``` r
# install.packages("remotes")
remotes::install_github("nwfsc-timeseries/atsar")
```

### EXAMPLE

Simulate data:

``` r
library(rstan)
#> Warning: package 'rstan' was built under R version 4.3.2
#> Loading required package: StanHeaders
#> Warning: package 'StanHeaders' was built under R version 4.3.2
#> 
#> rstan version 2.32.6 (Stan version 2.32.2)
#> For execution on a local, multicore CPU with excess RAM we recommend calling
#> options(mc.cores = parallel::detectCores()).
#> To avoid recompilation of unchanged Stan programs, we recommend calling
#> rstan_options(auto_write = TRUE)
#> For within-chain threading using `reduce_sum()` or `map_rect()` Stan functions,
#> change `threads_per_chain` option:
#> rstan_options(threads_per_chain = 1)
library(atsar)
set.seed(123)
s = cumsum(rnorm(50))
```

``` r
plot(s)
```

![](README-figs/plot-1.png)<!-- -->

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

To see the Stan mode code behind each of these, look in the `inst/stan`
folder on the GitHub repository. Note that `fit_stan.R` does some data
preparation to deal with Stan not accepting NAs in the data.

### DOCUMENTATION

- [ATSA lab book](https://nwfsc-timeseries.github.io/atsa-labs/) - Many
  applications are covered in our Applied Time Series Analysis book
  developed from the labs in our course.
- [ATSA course website](https://nwfsc-timeseries.github.io/atsa/) - We
  have lectures and all material from our course on our course website.
- Additional information can be found on the NWFSC time series page
  which includes several additional books and packages, [NWFSC time
  series page](https://nwfsc-timeseries.github.io/)

### CITATION

Ward, E.J., M.D. Scheuerell, and E.E. Holmes. 2018. ‘atsar’: Applied
Time Series Analysis in R: an introduction to time series analysis for
ecological and fisheries data with Stan.
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1158021.svg)](https://doi.org/10.5281/zenodo.1158021)

### NOAA Disclaimer

This repository is a scientific product and is not official
communication of the National Oceanic and Atmospheric Administration, or
the United States Department of Commerce. All NOAA GitHub project code
is provided on an ‘as is’ basis and the user assumes responsibility for
its use. Any claims against the Department of Commerce or Department of
Commerce bureaus stemming from the use of this GitHub project will be
governed by all applicable Federal law. Any reference to specific
commercial products, processes, or services by service mark, trademark,
manufacturer, or otherwise, does not constitute or imply their
endorsement, recommendation or favoring by the Department of Commerce.
The Department of Commerce seal and logo, or the seal and logo of a DOC
bureau, shall not be used in any manner to imply endorsement of any
commercial product or activity by DOC or the United States Government.

<img src="https://raw.githubusercontent.com/nmfs-general-modeling-tools/nmfspalette/main/man/figures/noaa-fisheries-rgb-2line-horizontal-small.png" height="75" alt="NOAA Fisheries">

[U.S. Department of Commerce](https://www.commerce.gov/) \| [National
Oceanographic and Atmospheric Administration](https://www.noaa.gov) \|
[NOAA Fisheries](https://www.fisheries.noaa.gov/)
