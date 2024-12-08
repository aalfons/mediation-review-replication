# ************************************
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ************************************


# load packages
library("robmed")
library("scales")
library("sn")


## control parameters for simulation
seed <- 20210617  # seed for the random number generator
n <- 200          # number of observations
K <- 1000         # number of simulation runs
R <- 5000         # number of bootstrap samples
# coefficients in mediation model
a <- b <- c <- 0.4
# error scales in mediation model
sigma_m <- sqrt(1 - a^2)
sigma_y <- sqrt(1 - b^2 - c^2 - 2*a*b*c)

## parameter combinations and fancy labels for plots
alpha <- c(0, 0, -Inf, Inf)
nu <- c(Inf, 2, Inf, Inf)
labels <- c("Normal", "Heavy tails", "Left-skewed", "Right-skewed")

## control parameters for methods
level <- 0.95                                        # confidence level
lmrob_control <- reg_control(max_iterations = 5000)  # MM-estimator

# function to compute mean of skew-t distribution
get_mean <- function(ksi = 0, omega = 1, alpha = 0, nu = Inf) {
  if (is.infinite(alpha)) delta <- sign(alpha)
  else delta <- alpha / sqrt(1 + alpha^2)
  if (is.infinite(nu)) b_nu <- sqrt(2 / pi)
  else b_nu <- sqrt(nu) * beta((nu-1)/2, 0.5) / (sqrt(pi) * gamma(0.5))
  # compute mean
  ksi + omega * delta * b_nu
}

## set seed of the random number generator for reproducibility
set.seed(seed)

## generate independent variable and error terms such that results are
## comparable across different values of the parameters
X <- rnorm(n)
# First uniformly distributed values are drawn, which are later transformed
# to the given error distribution by applying the inverse CDF.
u_m <- runif(n)
u_y <- runif(n)

## loop over parameter settings
df_weights <- NULL
for (i in seq_along(labels)) {

  # transform errors
  if (is.infinite(nu[i])) {
    # normal or skew-normal distribution
    if (is.finite(alpha[i]) && alpha[i] == 0) {
      # normal errors
      e_m <- qnorm(u_m)
      e_y <- qnorm(u_y)
    } else {
      # skew-normal errors
      e_m <- qsn(u_m, alpha = alpha[i])
      e_y <- qsn(u_y, alpha = alpha[i])
    }
  } else {
    # t or skew-t distribution
    if (is.finite(alpha[i]) && alpha[i] == 0) {
      # t-distributed errors
      e_m <- qt(u_m, df = nu[i])
      e_y <- qt(u_y, df = nu[i])
    } else {
      # skew-t distributed errors
      e_m <- qst(u_m, alpha = alpha[i], nu = nu[i])
      e_y <- qst(u_y, alpha = alpha[i], nu = nu[i])
    }
  }

  # compute mean of error terms such that it can be subtracted
  mu_m <- get_mean(omega = sigma_m, alpha = alpha[i], nu = nu[i])
  mu_y <- get_mean(omega = sigma_y, alpha = alpha[i], nu = nu[i])

  # transform error terms (scale and center)
  e_m_transformed <- sigma_m * e_m - mu_m
  e_y_transformed <- sigma_y * e_y - mu_y

  # generate hypothesized mediator and dependent variable
  M <- a * X + e_m_transformed
  Y <- b * M + c * X + e_y_transformed
  simulated_data <- data.frame(X, Y, M)

  # fit mediation model
  robust_fit <- fit_mediation(simulated_data, "X", "Y", "M",
                              control = lmrob_control)

  # obtain relevant information for diagnostic plot (only keep one of the
  # two plots since we illustrate different deviations from normality)
  setup <- setup_weight_plot(robust_fit, outcome = "M")

  # extract relevant data frame
  tmp <- data.frame(Distribution = factor(labels[i], levels = labels),
                    setup$data)
  df_weights <- rbind(df_weights, tmp)

}


# create plot
xlab <- "Weight threshold"
ylab <- "Percentage of observations with weight lower than threshold"
plt <- ggplot() +
  geom_line(aes(x = Threshold, y = Percentage, color = Weights),
            data = df_weights) +
  scale_color_manual("", values = c("black", "#00BFC4")) +
  facet_grid(Distribution ~ Tail) +
  labs(x = xlab, y = ylab) +
  theme_bw() +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 13),
        legend.direction = "vertical",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "line"),
        strip.text = element_text(size = 12))
class(plt) <- c("gg_weight_plot", class(plt))

# save plot to file
pdf("illustration/diagnostic_plot.pdf", width = 6.875, height = 7.5)
print(plt)
dev.off()
