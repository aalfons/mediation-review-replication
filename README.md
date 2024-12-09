# Robust mediation analysis: What we talk about when we talk about robustness

## About
Mediation analysis allows empirical researchers to study how an exposure
variable affects an outcome variable through one or more intervening variables.
Over the years, mediation analysis has become one of the most widely used
statistical techniques in the social, behavioral, and medical sciences. Yet
the most popular techniques for mediation analysis rely strongly on normality
assumptions and are therefore susceptible to the influence of outliers or
nonnormal distributions. We review common mediation models and discuss various
approaches for estimation and inference, including their implementation in
software packages. Moreover, we present robust alternatives, thereby clarifying
different notions of robustness. Finally, we consider a setting where the 
mediation model holds in a latent space but where measurement issues create
deviations from normality assumptions in the observed variable space, which
is a setting not commonly considered in the literature on robust mediation
analysis, and we obtain preleminary results via a simulation study.

More information can be found in our paper:

A. Alfons and D.R. Schley (2024).
Robust mediation analysis: What we talk about when we talk about robustness. 
PsyArXiv. doi: [10.31234/osf.io/2hqdy](https://doi.org/10.31234/osf.io/2hqdy).


## Reproduce results
This repository provides a collection of [R](https://CRAN.R-project.org/) 
scripts to reproduce the simulations and figures in our article.

The easiest way to reproduce the results is to clone this repository with 
[RStudio](https://rstudio.com/products/rstudio/download/).  Running the 
scripts within the resulting RStudio project ensures that there are no issues 
with file paths for storing or reading results, or for producing files 
containing plots.

Please note that this repository is rather large because it also contains an R 
data file with the simulation results.  This way, if you only want to quickly
reproduce the figures with simulation results, you do not actually need to run 
the simulations first.
