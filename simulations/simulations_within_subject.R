# **************************************
# Authors: Andreas Alfons and Dan Schley
#          Erasmus University Rotterdam
# **************************************


# load packages
library("robmed")      # regression-based mediation analysis
library("mediation")   # causal mediation analysis
library("quantreg")    # for median regression
library("parallel")    # for parallel computing

# load script with additional utility functions for simulation
source("simulations/utils.R")


## control parameters for simulation
seed <- 20241026  # seed for the random number generator
K <- 1000         # number of replications


## control parameters for data generation
# parameters to loop over
settings <- expand.grid(
  # whether X should be binary or continuous
  exposure = c("binary", "continuous"),
  # effect of M on Y
  b = c(0, 0.4),
  # noise-to-signal ratio when adding random measurement error
  noise_ratio = c(0, 0.25),
  # transformation to skewed or heavy-tailed distributions
  transformation = c("identity", "skewness", "heavy tails"),
  # bunching by rounding observation to nearest integer with given probability
  rounding_probability = c(0, 0.2),
  # replace observation by outlier with given probability
  outlier_probability = c(0, 0.02),
  # probability mass in left tail to be censored (for a normal distribution,
  # we censor at one standard deviation from the mean)
  censoring_fraction = c(0, pnorm(-1)),
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
nb_cores <- if (.Platform$OS.type == "windows") 1 else 12
# list for splitting settings
# (each list element gives settings to parallelize over)
nb_settings <- nrow(settings)
settings_index_list <- split(seq_len(nb_settings),
                             seq(from = 0, to = nb_settings - 1) %/% nb_cores)


## control parameters for methods
R <- 5000                                               # bootstrap samples
level <- 0.95                                           # confidence level
lmrob_control <- MM_reg_control(max_iterations = 5000)  # MM-regression
median_control <- median_reg_control(algorithm = "fn")  # median regression


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

  # ensure that the same bootstrap samples are used for all parameter
  # settings and all bootstrap tests for maximum comparability
  boot_indices <- boot_samples(n, R = R)

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

      # OLS bootstrap
      df_ols_boot <- tryCatch({
        # perform bootstrap test
        test <- test_mediation(measured_data, x = "x", y = "y", m = "m",
                               test = "boot", R = R, level = level,
                               type = "perc", indices = boot_indices,
                               method = "regression", robust = FALSE,
                               family = "gaussian")
        # bootstrap test for the indirect effect ab
        reject_index <- prod(test$ci) > 0
        # joint normal-theory t tests for a and b
        p_values_joint_t <- p_value(test, parm = c("a", "b"), type = "data")
        reject_joint_t <- all(p_values_joint_t < 1 - level)
        # joint bootstrap z tests for a and b
        p_values_joint_z <- p_value(test, parm = c("a", "b"), type = "boot")
        reject_joint_z <- all(p_values_joint_z < 1 - level)
        # add to data frame
        cbind(info, Method = "ols_boot", estimate_boot = test$indirect,
              estimate_data = test$fit$indirect, reject_index,
              reject_joint_t, reject_joint_z)
      }, error = function(condition) NULL)

      # winsorized bootstrap
      df_winsorized_boot <- tryCatch({
        # perform bootstrap test
        test <- test_mediation(measured_data, x = "x", y = "y", m = "m",
                               test = "boot", R = R, level = level,
                               type = "perc", indices = boot_indices,
                               method = "covariance", robust = TRUE)
        # bootstrap test for the indirect effect ab
        reject_index <- prod(test$ci) > 0
        # joint normal-theory t tests for a and b
        p_values_joint_t <- p_value(test, parm = c("a", "b"), type = "data")
        reject_joint_t <- all(p_values_joint_t < 1 - level)
        # joint bootstrap z tests for a and b
        p_values_joint_z <- p_value(test, parm = c("a", "b"), type = "boot")
        reject_joint_z <- all(p_values_joint_z < 1 - level)
        # add to data frame
        cbind(info, Method = "winsorized_boot", estimate_boot = test$indirect,
              estimate_data = test$fit$indirect, reject_index,
              reject_joint_t, reject_joint_z)
      }, error = function(condition) NULL)

      # median bootstrap
      df_median_boot <- tryCatch({
        # perform bootstrap test
        # (many warnings due to nonunique solutions on bootstrap samples)
        test <- suppressWarnings(
          test_mediation(measured_data, x = "x", y = "y", m = "m",
                         test = "boot", R = R, level = level,
                         type = "perc", indices = boot_indices,
                         method = "regression", robust = "median",
                         control = median_control)
        )
        # bootstrap test for the indirect effect ab
        reject_index <- prod(test$ci) > 0
        # joint normal-theory t tests for a and b
        p_values_joint_t <- p_value(test, parm = c("a", "b"), type = "data")
        reject_joint_t <- all(p_values_joint_t < 1 - level)
        # joint bootstrap z tests for a and b
        p_values_joint_z <- p_value(test, parm = c("a", "b"), type = "boot")
        reject_joint_z <- all(p_values_joint_z < 1 - level)
        # add to data frame
        cbind(info, Method = "median_boot", estimate_boot = test$indirect,
              estimate_data = test$fit$indirect, reject_index,
              reject_joint_t, reject_joint_z)
      }, error = function(condition) NULL)

      # ROBMED
      df_ROBMED <- tryCatch({
        # perform bootstrap test
        test <- test_mediation(measured_data, x = "x", y = "y", m = "m",
                               test = "boot", R = R, level = level,
                               type = "perc", indices = boot_indices,
                               method = "regression", robust = TRUE,
                               control = lmrob_control)
        # bootstrap test for the indirect effect ab
        reject_index <- prod(test$ci) > 0
        # joint normal-theory t tests for a and b
        p_values_joint_t <- p_value(test, parm = c("a", "b"), type = "data")
        reject_joint_t <- all(p_values_joint_t < 1 - level)
        # joint bootstrap z tests for a and b
        p_values_joint_z <- p_value(test, parm = c("a", "b"), type = "boot")
        reject_joint_z <- all(p_values_joint_z < 1 - level)
        # add to data frame
        cbind(info, Method = "ROBMED", estimate_boot = test$indirect,
              estimate_data = test$fit$indirect, reject_index,
              reject_joint_t, reject_joint_z)
      }, error = function(condition) NULL)

      # for binary exposure variable, apply causal mediation analysis
      # (implementation is rather slow; also note that unlike for the methods
      # above, we cannot supply the same bootstrap samples as an argument)
      if (exposure == "binary") {

        # causal mediation analysis based on OLS regressions
        df_ols_causal <- tryCatch({
          # perform causal mediation analysis
          suppressMessages({
            # apply OLS regressions
            fit_m <- lm(m ~ x, data = measured_data)
            fit_y <- lm(y ~ m + x, data = measured_data)
            # apply causal bootstrap test
            test <- mediate(fit_m, fit_y, sims = 1000, boot = TRUE,
                            boot.ci.type = "perc", treat = "x",
                            mediator = "m", conf.level = level)
          })
          # point estimates of the average conditional mediation effect
          acme_boot <- mean(test$d.avg.sims, na.rm = TRUE)
          acme_data <- test$d.avg
          # bootstrap test for the average conditional mediation effect
          reject_index <- prod(test$d.avg.ci) > 0
          # add to data frame
          cbind(info, Method = "ols_causal", estimate_boot = acme_boot,
                estimate_data = acme_data, reject_index,
                reject_joint_t = NA, reject_joint_z = NA)
        }, error = function(condition) NULL)

        # causal mediation analysis based on median regressions
        df_median_causal <- tryCatch({
          # perform causal mediation analysis
          suppressMessages({
            suppressWarnings({
              # apply median regressions
              fit_m <- rq(m ~ x, tau = 0.5, data = measured_data,
                          method = median_control$method)
              fit_y <- rq(y ~ m + x, tau = 0.5, data = measured_data,
                          method = median_control$method)
              # apply causal bootstrap test
              test <- mediate(fit_m, fit_y, sims = 1000, boot = TRUE,
                              boot.ci.type = "perc", treat = "x",
                              mediator = "m", conf.level = level)
            })
          })
          # point estimates of the average conditional mediation effect
          acme_boot <- mean(test$d.avg.sims, na.rm = TRUE)
          acme_data <- test$d.avg
          # bootstrap test for the average conditional mediation effect
          reject_index <- prod(test$d.avg.ci) > 0
          # add to data frame
          cbind(info, Method = "median_causal", estimate_boot = acme_boot,
                estimate_data = acme_data, reject_index,
                reject_joint_t = NA, reject_joint_z = NA)
        }, error = function(condition) NULL)

      } else {
        df_ols_causal <- NULL
        df_median_causal <- NULL
      }

      # combine results
      rbind(df_ols_boot, df_winsorized_boot, df_median_boot, df_ROBMED,
            df_ols_causal, df_median_causal)

    },
    # Argument 'mc.set.seed = FALSE' is set such that seed of the random
    # number generator is passed on from main thread to the worker threads.
    # Together with splitting the indices of settings into a list ensures
    # that the same random numbers are used in each setting for the current
    # replication. The only place inside the parallelized loop where random
    # numbers are needed is the subsampling algorithm of the MM-estimator.
    # Variability in the simulated data across the replications is guaranteed
    # by performing the random draws outside the parallelized loop.
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
cat(format(Sys.time()), ": finished.\n")

## store results
session_info <- sessionInfo()
file <- "simulations/results/results_within_subject.RData"
save(results, n, c, level, sigma_e_m_latent, sigma_e_y_latent,
     transformation_lambda, outlier_level_m, outlier_level_y,
     outlier_multiplier, session_info, file = file)
