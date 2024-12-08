# ************************************
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ************************************


# apply the inverse of the Yeo-Johnson transformation or an antisymmetrized
# version to obtain a skewed or heavy tailed distribution
# Note: The argument 'inverse' corresponds to the desired transformation for
# our purposes.  That is, 'inverse = FALSE' gives the inverse (antisymmetrized)
# Yeo-Johnson transformation that yields a skewed or heavy-tailed distribution,
# while 'inverse = TRUE' gives the (antisymmetrized) Yeo-Johnson transformation
# that transforms back to a normal distribution.
transform <- function(x, lambda, antisymmetric = FALSE, inverse = FALSE) {
  if (antisymmetric) {
    sign(x) * VGAM::yeo.johnson(abs(x), lambda = lambda, inverse = !inverse)
  } else {
    VGAM::yeo.johnson(x, lambda = lambda, inverse = !inverse)
  }
}


# replace observation by value rounded to nearest integer with given probability
replace_with_rounding <- function(x, u = runif(x), probability = 0.2) {
  # replace selected observations with rounded value
  ifelse(u < probability, round(x), x)
}


# replace observation by outlier with given probability: the function assumes
# that the two supplied variables have the same type of distribution and are
# positively correlated so that it places the outliers in the top left or
# bottom right corner
# level ........ probability for quantile to be used in outlier shift
# multiplier ... multiplication factor to be applied to the original value to
#                bring it closer to zero before shift (so we still have some
#                variability in the outliers, but less than for the majority)
# parameters ... list of parameters describing the distribution
# Notes to Dan:
# 1) If we want harmful outliers with a high leverage effect to influence the a
#    path, we would also need to include outliers in X, which we wanted to keep
#    perfectly measured. So I suggest to fix the a path to a certain value, and
#    let b vary to generate settings with and without mediation. We can generate
#    outliers with high leverage for the b path by placing outliers in M and Y.
# 2) For now, I implemented this such that we do the same as in my ORM paper
#    in the case of a normal distribution, but implemented more generally with
#    the outlier shift based on a quantile rather than standard deviations.
# 3) In the simulations below, I call this function with a negative shift in M
#    and a positive shift in Y to go against the positive relationship for the
#    majority of the observations.
replace_with_outliers <- function(x, u = runif(x), probability = 0.05,
                                  level = 0.001, multiplier = 0.1,
                                  parameters = list(mu = 0, sigma = 1)) {
  # check if any observations should be replaced with outliers
  which <- which(u < probability)
  if (length(which) == 0L) return(x)
  # obtain value of outlier shift for normal distribution
  shift <- qnorm(level, mean = parameters$mu, sd = parameters$sigma)
  if (!is.null(parameters$lambda)) {
    # apply transformation to obtain quantile of transformed distribution
    shift <- transform(shift, lambda = parameters$lambda,
                       antisymmetric = parameters$antisymmetric)
  }
  # compute values of outliers from original values
  outliers <- x[which] * multiplier + shift
  # replace observations with outliers by alternating the candidate values
  x[which] <- outliers
  x
}


# left-censor a variable
left_censor <- function(x, fraction = pnorm(-1),
                        parameters = list(mu = 0, sigma = 1)) {
  # obtain lower bound for censoring
  lower <- qnorm(fraction, mean = parameters$mu, sd = parameters$sigma)
  if (!is.null(parameters$lambda)) {
    # apply transformation to obtain quantile of transformed distribution
    lower <- transform(lower, lambda = parameters$lambda,
                       antisymmetric = parameters$antisymmetric)
  }
  # apply left-censoring
  x[x < lower] <- lower
  x
}
