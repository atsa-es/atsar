## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----install, eval=TRUE, warning=FALSE, message=FALSE, results='hide'---------
library(rstan)
library(devtools)
library(atsar)
# for optimizing stan on your machine,
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

## -----------------------------------------------------------------------------
library(MARSS)
data(SalmonSurvCUI)

## ---- message=FALSE, warning=FALSE, results='hide'----------------------------
mod = fit_stan(y = SalmonSurvCUI$logit.s, model_name="dlm-intercept")

