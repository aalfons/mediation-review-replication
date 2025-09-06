# ************************************
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ************************************


# load packages
library("lavaan")    # for SEM
library("parallel")  # for parallel computing

# load script with additional utility functions for simulation
source("simulations/utils.R")


## control parameters for simulation
seed <- 20241026  # seed for the random number generator
K <- 1000         # number of replications


## control parameters for data generation
# parameters to loop over
settings <- data.frame(
  # whether X should be binary or continuous
  exposure = "continuous",
  # effect of M on Y
  b = rep.int(c(0, 0.4), times = 7),
  # noise-to-signal ratio when adding random measurement error
  noise_ratio = rep(c(0, 0.25, 0, 0, 0, 0, 0), each = 2),
  # transformation to skewed or heavy-tailed distributions
  transformation = rep(c("identity", "identity", "skewness", "heavy tails",
                         "identity", "identity", "identity"), each = 2),
  # bunching by rounding observation to nearest integer with given probability
  rounding_probability = rep(c(0, 0, 0, 0, 0.2, 0, 0), each = 2),
  # replace observation by outlier with given probability
  outlier_probability = rep(c(0, 0, 0, 0, 0, 0.02, 0), each = 2),
  # probability mass in left tail to be censored (for a normal distribution,
  # we censor at one standard deviation from the mean)
  censoring_fraction = rep(c(0, 0, 0, 0, 0, 0, pnorm(-1)), each = 2),
  # other arguments
  stringsAsFactors = FALSE
)
# parameters to remain fixed
n <- 100                      # number of observations
a <- 0.4                      # effect of X on M
c <- 0.4                      # direct effect of X on Y
sigma_e_m_latent <- 1         # error scale in (latent) regression of M on X
sigma_e_y_latent <- 1         # error scale in (latent) regression of Y on M and X
transformation_lambda <- 1/3  # parameter of inverse Yeo-Johnson transformation
outlier_level_m <- 0.001      # probability for quantile for shift in m
outlier_level_y <- 0.999      # probability for quantile for shift in y
outlier_multiplier <- 0.1     # scaling factor of outliers before shift


## control parameters for parallel computing
# number of CPU cores to be used
nb_cores <- if (.Platform$OS.type == "windows") 1 else 2
# list for splitting settings
# (each list element gives settings to parallelize over)
nb_settings <- nrow(settings)
settings_index_list <- split(seq_len(nb_settings),
                             seq(from = 0, to = nb_settings - 1) %/% nb_cores)


## input for methods
model <- ' # direct effect
           y ~ c*x
           # mediator
           m ~ a*x
           y ~ b*m
           # indirect effect
           ab := a*b
           # total effect
           total := (a*b) + c
         '
level <- 0.95       # confidence level
alpha <- 1 - level  # significance level


## run simulation
set.seed(seed)
cat(format(Sys.time()), ": starting ...\n")
results_list <- lapply(seq_len(K), function(k) {

  # print information on replication
  cat(format(Sys.time()), sprintf(":   replication = %d\n", k))

  # generate a uniform variable for X (which will be transformed to different
  # distributions) as well as normal regression and measurement error terms
  x_uniform <- runif(n)
  e_m_latent <- rnorm(n)
  e_y_latent <- rnorm(n)
  e_m_measurement <- rnorm(n)
  e_y_measurement <- rnorm(n)

  # generate uniform random variables: if lower than a given probability, the
  # corresponding observation is replaced with rounded value or outlier
  u_rounding <- runif(n)
  u_outlier <- runif(n)

  # loop over list of settings for parallelization
  results_k <- lapply(settings_index_list, function(settings_indices) {

    # parallelize over different settings
    results_parallel <- mclapply(settings_indices, function(l) {

      # print information on current setting
      cat(format(Sys.time()), sprintf(":     setting = %d\n", l))

      # extract information from data frame
      exposure <- settings[l, "exposure"]
      b <- settings[l, "b"]
      noise_ratio <- settings[l, "noise_ratio"]
      transformation <- settings[l, "transformation"]
      rounding_probability <- settings[l, "rounding_probability"]
      outlier_probability <- settings[l, "outlier_probability"]
      censoring_fraction  <- settings[l, "censoring_fraction"]

      # generate exposure variable (mean 0 and variance 1)
      x <- switch(exposure,
                  binary = ifelse(x_uniform < 0.5, -1, 1),
                  continuous = qnorm(x_uniform))

      # generate latent variables
      m_latent <- a * x + sigma_e_m_latent * e_m_latent
      y_latent <- b * m_latent + c * x + sigma_e_y_latent * e_y_latent

      # compute means and standard deviations of latent variables
      mu_m_latent <- 0
      mu_y_latent <- 0
      sigma_m_latent <- sqrt(a^2 + sigma_e_m_latent^2)
      sigma_y_latent <- sqrt(b^2 * (a^2 + sigma_e_m_latent^2) + c^2 +
                               2 * a * b * c + sigma_e_y_latent^2)

      # generate measured variables
      if (noise_ratio > 0) {
        # compute standard deviation of measurement error
        sigma_e_m_measurement <- noise_ratio * sigma_m_latent
        sigma_e_y_measurement <- noise_ratio * sigma_y_latent
        # add measurement error to latent variables variables
        m <- m_latent + sigma_e_m_measurement * e_m_measurement
        y <- y_latent + sigma_e_y_measurement * e_y_measurement
        # compute marginal standard deviations of measured variables
        sigma_m <- sqrt(1 + noise_ratio^2) * sigma_m_latent
        sigma_y <- sqrt(1 + noise_ratio^2) * sigma_y_latent
      } else {
        # no measurement error
        m <- m_latent
        y <- y_latent
        sigma_m <- sigma_m_latent
        sigma_y <- sigma_y_latent
      }
      # compute marginal means of measured variables
      mu_m <- 0
      mu_y <- 0

      # transform to skewed or heavy-tailed distributions
      if (transformation != "identity") {
        # transform the variables
        transformation_antisymmetric <- transformation == "heavy tails"
        m <- transform(m, lambda = transformation_lambda,
                       antisymmetric = transformation_antisymmetric)
        y <- transform(y, lambda = transformation_lambda,
                       antisymmetric = transformation_antisymmetric)
      }

      # replace observations by value rounded to nearest integer with given
      # probability
      if (rounding_probability > 0) {
        m <- replace_with_rounding(m, u = u_rounding,
                                   probability = rounding_probability)
        y <- replace_with_rounding(y, u = u_rounding,
                                   probability = rounding_probability)
      }


      # prepare list of parameters describing the distributions of M and Y
      if (outlier_probability > 0 || censoring_fraction > 0) {
        parameters_m <- list(mu = mu_m, sigma = sigma_m)
        parameters_y <- list(mu = mu_y, sigma = sigma_y)
        if (transformation != "identity") {
          parameters_m$lambda <- transformation_lambda
          parameters_m$antisymmetric <- transformation_antisymmetric
          parameters_y$lambda <- transformation_lambda
          parameters_y$antisymmetric <- transformation_antisymmetric
        }
      }

      # replace observations by outliers with given probability
      if (outlier_probability > 0) {
        # generate outliers
        m <- replace_with_outliers(m, u = u_outlier,
                                   probability = outlier_probability,
                                   level = outlier_level_m,
                                   multiplier = outlier_multiplier,
                                   parameters = parameters_m)
        y <- replace_with_outliers(y, u = u_outlier,
                                   probability = outlier_probability,
                                   level = outlier_level_y,
                                   multiplier = outlier_multiplier,
                                   parameters = parameters_y)
      }

      # censor the measured variables
      if (censoring_fraction > 0) {
        m <- left_censor(m, fraction = censoring_fraction,
                         parameters = parameters_m)
        y <- left_censor(y, fraction = censoring_fraction,
                         parameters = parameters_y)
      }

      # construct data frame of measured variables
      measured_data <- data.frame(x = x, m = m, y = y)

      # information on current setting to be stored with results
      info <- data.frame(replication = k, exposure, a, b, noise_ratio,
                         transformation, rounding_probability,
                         outlier_probability, censoring_fraction)

      # ML estimation with bootstrap inference
      df_ml_boot <- tryCatch({
        # fit SEM and perform bootstrap
        suppressWarnings(
          fit <- sem(model, data = measured_data, se = "bootstrap",
                     bootstrap = 1000)
        )
        coefficients <- summary(fit)$pe
        # extract relevant indices
        which_a <- which(coefficients$label == "a")
        which_b <- which(coefficients$label == "b")
        which_indirect <- which(coefficients$label == "ab")
        # extract bootstrap replicates
        boot <- fit@boot$coef
        indirect <- boot[, which_a] * boot[, which_b]
        ci <- quantile(indirect, probs = c(alpha/2,  1-alpha/2), na.rm = TRUE)
        # bootstrap test for the indirect effect ab
        reject_index <- prod(ci) > 0
        # joint bootstrap z tests for a and b
        p_values_joint <- coefficients[c(which_a, which_b), "pvalue"]
        reject_joint <- all(p_values_joint < alpha)
        # add to data frame
        cbind(info, Method = "ml_boot",
              estimate_boot = mean(indirect, na.rm = TRUE),
              estimate_data = coefficients[which_indirect, "est"],
              reject_index, reject_joint)
      }, error = function(condition) NULL)

      # ML estimation with robust inference
      df_ml_robust <- tryCatch({
        # fit SEM and perform bootstrap
        fit <- sem(model, data = measured_data, estimator = "MLM")
        coefficients <- summary(fit)$pe
        # extract relevant indices
        which_a <- which(coefficients$label == "a")
        which_b <- which(coefficients$label == "b")
        which_indirect <- which(coefficients$label == "ab")
        # z test for the indirect effect ab
        p_value_index <- coefficients[which_indirect, "pvalue"]
        reject_index <- p_value_index < alpha
        # joint z tests for a and b
        p_values_joint <- coefficients[c(which_a, which_b), "pvalue"]
        reject_joint <- all(p_values_joint < alpha)
        # add to data frame
        cbind(info, Method = "ml_robust", estimate_boot = NA_real_,
              estimate_data = coefficients[which_indirect, "est"],
              reject_index, reject_joint)
      }, error = function(condition) NULL)

      # combine results
      rbind(df_ml_boot, df_ml_robust)

    },
    # Argument 'mc.set.seed = FALSE' is set such that seed of the random
    # number generator is passed on from main thread to the worker threads.
    # Together with splitting the indices of settings into a list ensures
    # that the same random numbers are used in each setting for the current
    # replication. However, lavaan will draw different bootstrap samples in
    # for different settings of the same replication. Variability in the
    # simulated data across the replications is guaranteed by performing the
    # random draws outside the parallelized loop.
    mc.set.seed = FALSE,
    mc.cores = length(settings_indices))

    # combine results from parallelization into data frame
    do.call(rbind, results_parallel)

  })

  # combine results for current replication into data frame
  do.call(rbind, results_k)

})

## combine results into data frame
results <- do.call(rbind, results_list)
row.names(results) <- NULL
cat(format(Sys.time()), ": finished.\n")

## store results
session_info <- sessionInfo()
file <- "simulations/results/results_lavaan.RData"
save(results, n, c, level, sigma_e_m_latent, sigma_e_y_latent,
     transformation_lambda, outlier_level_m, outlier_level_y,
     outlier_multiplier, session_info, file = file)
