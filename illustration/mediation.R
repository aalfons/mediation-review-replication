# ************************************
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ************************************


## load packages
library("colorspace")
library("ggh4x")
library("mvtnorm")
library("robmed")
library("viridis")

## function to let color fade
fade_color <- function(color, alpha = 1) {
  # 'alpha; behaves similarly to transparancy in ggplot2, but colors are faded
  # to white and not transparant (the default is to use fully opaque colors)
  rgb <- col2rgb(color) / 255
  faded <- mixcolor(1 - alpha, RGB(rgb[1, ], rgb[2, ], rgb[3, ]), RGB(1, 1, 1))
  rgb(faded@coords[, 1], faded@coords[, 2], faded@coords[, 3])
}

## function to convert colors to grayscale
col2gray <- function (color, method = c("luminosity", "average")) {
  # initializations
  method <- match.arg(method)
  # convert to RGB values
  rgb <- col2rgb(color)
  # compute weighted combination of RGB values
  weights <- switch(method, luminosity = c(0.2126, 0.7152, 0.0722),
                    average = rep.int(1/3, 3))
  gray <- crossprod(weights, rgb)
  # convert to text string
  rgb(gray, gray, gray, maxColorValue = 255)
}


## control parameters for data generation
n <- 100          # number of observations
p_w <- 5          # number of items within a scale
sigma <- 0.8      # standard deviations of continuous variables
rho_w <- 0.5      # correlation of items within a scale
rho_b <- 0.3      # correlation of items between scales
seed <- 20240528  # seed for the random number generator


## control parameters for discretization
nb_cat <- 9                                             # number of categories
values <- 1:nb_cat                                      # rating-scale values
midpoint <- (1 + nb_cat) / 2                            # rating-scale midpoint
breaks <- c(-Inf, (values[-1]+values[-nb_cat])/2, Inf)  # cut points


## construct covariance matrix for generation of scales
# loadings onto latent constructs
alpha <- rep.int(1/p_w, times = p_w)
# generate diagonal block for scale measuring latent construct behind M
Sigma_MM <- matrix(rho_w, nrow = p_w, ncol = p_w)
diag(Sigma_MM) <- sigma^2
# compute variance of latent construct behind M
sigma_M <- sqrt(drop(t(alpha) %*% Sigma_MM %*% alpha))
# generate diagonal block for scale measuring latent construct behind Y
Sigma_YY <- matrix(rho_w, nrow = p_w, ncol = p_w)
diag(Sigma_YY) <- sigma^2
# compute variance of latent construct behind Y
sigma_Y <- sqrt(drop(t(alpha) %*% Sigma_YY %*% alpha))
# covariance matrix between scales
Sigma_YM <- matrix(rho_b * sigma^2, nrow = p_w, ncol = p_w)
sigma_YM <- drop(t(alpha) %*% Sigma_YM %*% alpha)
# put blockdiagnonal matrix together for data generation
Sigma <- cbind(
  rbind(Sigma_MM, Sigma_YM),
  rbind(t(Sigma_YM), Sigma_YY)
)

## true coefficients in mediation model
a <- 1.4
b <- sigma_YM / sigma_M^2
c <- 1.4

## means of M and Y
mu0 <- c(M = 0, Y = 0)
mu1 <- c(M = a, Y = b*a + c)

## generate data with multiple continuous variables for M and Y
set.seed(seed)
X <- sort(sample(0:1, n, replace = TRUE))
n0 <- sum(X == 0)
n1 <- n - n0
MY0 <- rmvnorm(n0, mean = rep(mu0, each = p_w), sigma = Sigma)
MY1 <- rmvnorm(n1, mean = rep(mu1, each = p_w), sigma = Sigma)
MY <- rbind(MY0, MY1)

## discretize the continuous variables for M and Y
shift <- c(M = 2.4, Y = 5.2)
MY <- sweep(MY, 2, rep(shift, each = p_w), "+")
MY <- apply(MY, 2, function(x) as.numeric(cut(x, breaks = breaks)))
M <- rowMeans(MY[, 1:p_w])
Y <- rowMeans(MY[, p_w + (1:p_w)])

## different situations for outlier
mu <- c(X = 1, M = 9, Y = 1)    # means for contaminated data


## apply OLS and robust regression to data with and without outlier
# create data frames
data_without <- data.frame(X, M, Y)
data_with <- rbind(data_without, mu)
# apply regression
fit_OLS_without <- fit_mediation(data_without, "X", "Y", "M", robust = FALSE)
fit_ROBMED_without <- fit_mediation(data_without, "X", "Y", "M", robust = TRUE)
fit_OLS_with <- fit_mediation(data_with, "X", "Y", "M", robust = FALSE)
fit_ROBMED_with <- fit_mediation(data_with, "X", "Y", "M", robust = TRUE)
# combine eveything into a list of list for loop to create plot
scenarios <- c("Without outlier", "Including outlier")
methods <- c("OLS", "ROBMED")
settings <- data.frame(scenario = rep(scenarios, each = 2),
                       method = rep(methods, times = 2))
fit_list <- list(fit_OLS_without, fit_ROBMED_without,
                 fit_OLS_with, fit_ROBMED_with)


## graphical parameters for pointers to outlier
offset <- c(X = NA, M = -0.75, Y = 1.5)
label <- "Outlier"
col_pointer <- "darkgray"
pointer <- arrow(length = unit(0.15, "cm"), ends = "last", type = "closed")

## additional pointer to indirect effect that is pulled flat by outlier
affected <- arrow(angle = 0, length = unit(0, "cm"))

## control parameters for plots
max_Y <- max(Y)
col_lines <- viridis_pal(begin = 0.3, end = 0.7, direction = -1)(2)
# col_points <- fade_color(col_lines, alpha = 0.8)  # faded colors for points
col_points <- col_lines                           # no color fade for points
line_types <- c("dotdash", "solid", "dashed")
line_size <- 2/3
effect <- arrow(length = unit(0.2, "cm"), ends = "both", type = "open")
arrow_size <- 2/3


## loop over scenario and method
df_points <- NULL
df_MX <- df_YX <- df_YMX <- NULL
df_a <- df_c <- df_c_prime <- NULL
df_ab <- list()
df_label_top <- df_label_left <- df_label_right <- df_label_arrow <- NULL
df_outlier_arrow <- df_outlier_label <- NULL
for (k in seq_along(fit_list)) {

  ## extract relevant information
  scenario <- settings[k, "scenario"]
  method <- settings[k, "method"]
  fit <- fit_list[[k]]
  data <- if (scenario == scenarios[1]) data_without else data_with

  ## observations to be plotted
  data <- data.frame(data, Scenario = factor(scenario, levels = scenarios),
                     Method = factor(method, levels = methods))
  data$X <- as.factor(data$X)
  df_points <- rbind(df_points, data)

  ## extract coefficients
  coef_MX <- coef(fit$fit_mx)
  coef_YMX <- coef(fit$fit_ymx)
  if (method == methods[1]) {
    coef_YX <- coef(fit$fit_yx)
  } else if (method == methods[2]) {
    i2 <- coef_YMX[1] + coef_YMX[2] * coef_MX[1]
    coef_YX <- c(i2, fit$total)
  }

  ## equations to be plotted
  tmp_MX <- data.frame(Equation = c("i1", "i1 + a"),
                       Value = c(unname(coef_MX)[1], sum(coef_MX)),
                       Scenario = factor(scenario, levels = scenarios),
                       Method = factor(method, levels = methods))
  df_MX <- rbind(df_MX, tmp_MX)
  tmp_YMX <- data.frame(Equation = c("i2 + b*M", "i2 + b*M + c"),
                        Intercept = c(unname(coef_YMX)[1], sum(coef_YMX[-2])),
                        Slope = rep.int(unname(coef_YMX)[2], 2),
                        Scenario = factor(scenario, levels = scenarios),
                        Method = factor(method, levels = methods))
  df_YMX <- rbind(df_YMX, tmp_YMX)
  tmp_YX <- data.frame(Equation = c("i3", "i3 + c'"),
                       Value = c(unname(coef_YX)[1], sum(coef_YX)),
                       Scenario = factor(scenario, levels = scenarios),
                       Method = factor(method, levels = methods))
  df_YX <- rbind(df_YX, tmp_YX)

  ## effects to be plotted
  tmp_a <- data.frame(x = tmp_MX$Value, y = min(values),
                      Scenario = factor(scenario, levels = scenarios),
                      Method = factor(method, levels = methods))
  df_a <- rbind(df_a, tmp_a)
  tmp_c_prime <- data.frame(x = max(values), y = tmp_YX$Value,
                            Scenario = factor(scenario, levels = scenarios),
                            Method = factor(method, levels = methods))
  df_c_prime <- rbind(df_c_prime, tmp_c_prime)
  intersection_Y <- tmp_YMX[1, "Intercept"] + tmp_YMX[1, "Slope"] * tmp_MX[2, "Value"]
  tmp_c <- data.frame(x = tmp_MX[2, "Value"],
                      y = c(intersection_Y, tmp_YX[2, "Value"]),
                      Scenario = factor(scenario, levels = scenarios),
                      Method = factor(method, levels = methods))
  df_c <- rbind(df_c, tmp_c)
  tmp_ab <- data.frame(x = tmp_MX[2, "Value"],
                       y = c(tmp_YX[1, "Value"], intersection_Y),
                       Scenario = factor(scenario, levels = scenarios),
                       Method = factor(method, levels = methods))
  df_ab[[paste(scenario, method, sep = "_")]] <- tmp_ab

  ## effect labels to be plotted
  # labels plotted on top of arrows
  tmp <- data.frame(x = mean(tmp_a$x),
                    y = mean(tmp_a$y),
                    Label = "italic(hat(a))",
                    Scenario = factor(scenario, levels = scenarios),
                    Method = factor(method, levels = methods))
  df_label_top <- rbind(df_label_top, tmp)
  # labels plotted to the left of arrows
  tmp <- data.frame(x = mean(tmp_c_prime$x),
                    y = mean(tmp_c_prime$y),
                    Label = "italic(paste(hat(c), \"'\"))",
                    Scenario = factor(scenario, levels = scenarios),
                    Method = factor(method, levels = methods))
  df_label_left <- rbind(df_label_left, tmp)
  # labels plotted to the right of arrows
  labels_right <- c("italic(hat(c))", "italic(widehat(ab))")
  if(scenario == scenarios[2] && method == methods[1]) {
    # no space for arrow of standard method when the outlier is included,
    # so label should be offset a bit
    offset_ab <- c(x = 0.4, y = -0.4)
    tmp <- data.frame(x = c(mean(tmp_c$x),
                            mean(tmp_ab$x) + offset_ab["x"]),
                      y = c(mean(tmp_c$y),
                            mean(tmp_ab$y) + offset_ab["y"]),
                      Label = labels_right,
                      Scenario = factor(scenario, levels = scenarios),
                      Method = factor(method, levels = methods))
    df_label_right <- rbind(df_label_right, tmp)
    # add an arrow to connect the label to the effect
    start <- c(mean(tmp_ab$x) + offset_ab["x"],
               mean(tmp_ab$y) + offset_ab["y"])
    end <- c(x = mean(tmp_ab$x), y = mean(tmp_ab$y))
    cut_start <- start - 0.1 * (end - start)
    cut_end <- start + 0.95 * (end - start)
    tmp <- data.frame(x = cut_start["x"], y = cut_start["y"],
                      x_end = cut_end["x"], y_end = cut_end["y"],
                      Scenario = factor(scenario, levels = scenarios),
                      Method = factor(method, levels = methods))
    df_label_arrow <- rbind(df_label_arrow, tmp)
  } else {
    # enough space for all the arrows and labels
    tmp <- data.frame(x = c(mean(tmp_c$x), mean(tmp_ab$x)),
                      y = c(mean(tmp_c$y), mean(tmp_ab$y)),
                      Label = labels_right,
                      Scenario = factor(scenario, levels = scenarios),
                      Method = factor(method, levels = methods))
    df_label_right <- rbind(df_label_right, tmp)
  }

  ## extra information on outlier
  if(scenario == scenarios[2]) {
    ## arrow to outlier to be plotted
    start <- c(mu["M"] + offset["M"], mu["Y"] + offset["Y"])
    end <- c(mu["M"], mu["Y"])
    names(start) <- names(end) <- c("x", "y")
    cut_start <- start + 0.1 * (end - start)
    cut_end <- start + 0.8 * (end - start)
    tmp <- data.frame(x = cut_start["x"], y = cut_start["y"],
                      x_end = cut_end["x"], y_end = cut_end["y"],
                      Scenario = factor(scenario, levels = scenarios),
                      Method = factor(method, levels = methods))
    df_outlier_arrow <- rbind(df_outlier_arrow, tmp)

    ## outlier label to be plotted
    tmp <- data.frame(x = mu["M"] + offset["M"],
                      y = mu["Y"] + offset["Y"],
                      Label = label,
                      Scenario = factor(scenario, levels = scenarios),
                      Method = factor(method, levels = methods))
    df_outlier_label <- rbind(df_outlier_label, tmp)
  }
}

# fix data frame for arrows for indirect effect
df_ab <- list(arrow = rbind(df_ab[[1]], df_ab[[2]], df_ab[[4]]),
              line = df_ab[[3]])

# equation labels for legend
equations <- expression(italic(hat(M) == hat(i)[1]),
                        italic(paste(hat(M) == hat(i)[1] + hat(a), "  ")),
                        italic(hat(Y) == hat(i)[2] + hat(b) %.% M),
                        italic(hat(Y) == hat(i)[2] + hat(b) %.% M + hat(c)),
                        italic(hat(Y) == hat(i)[3]),
                        italic(paste(hat(Y) == hat(i)[3] + hat(c), "'  ")))


## plot data and mediation effects
# create plot object
plt_estimation <- ggplot() +
  geom_point(aes(x = M, y = Y, fill = X), data = df_points,
             color = "transparent", shape = 21, size = 2.5) +
  geom_vline(aes(xintercept = Value, linetype = Equation, color = Equation),
             data = df_MX, linewidth = line_size) +
  geom_hline(aes(yintercept = Value, linetype = Equation, color = Equation),
             data = df_YX, linewidth = line_size) +
  geom_abline(aes(intercept = Intercept, slope = Slope, linetype = Equation,
                  color = Equation), data = df_YMX, linewidth = line_size) +
  geom_segment(aes(x = x, y = y, xend = x_end, yend = y_end),
               data = df_label_arrow, color = col_pointer,
               arrow = pointer) +
  geom_line(aes(x = x, y = y), data = df_a, arrow = effect,
            linewidth = arrow_size) +
  geom_line(aes(x = x, y = y), data = df_c, arrow = effect,
            linewidth = arrow_size) +
  geom_line(aes(x = x, y = y), data = df_c_prime, arrow = effect,
            linewidth = arrow_size) +
  geom_line(aes(x = x, y = y), data = df_ab[["arrow"]], arrow = effect,
            linewidth = arrow_size) +
  geom_line(aes(x = x, y = y), data = df_ab[["line"]], arrow = affected,
            linewidth = arrow_size) +
  geom_text(aes(x = x, y = y, label = Label), data = df_label_top, size = 4,
            hjust = "center", vjust = "bottom", parse = TRUE, nudge_y = 0.2) +
  geom_text(aes(x = x, y = y, label = Label), data = df_label_left, size = 4,
            hjust = "right", vjust = "middle", parse = TRUE, nudge_x = -0.25) +
  geom_text(aes(x = x, y = y, label = Label), data = df_label_right, size = 4,
            hjust = "left", vjust = "middle", parse = TRUE, nudge_x = 0.25) +
  geom_segment(aes(x = x, y = y, xend = x_end, yend = y_end),
               data = df_outlier_arrow, color = col_pointer,
               arrow = pointer) +
  geom_text(aes(x = x, y = y, label = Label), data = df_outlier_label, size = 4,
            hjust = "center", vjust = "bottom") +
  scale_x_continuous(breaks = values, minor_breaks = NULL,
                     limits = range(values)) +
  scale_y_continuous(breaks = values, minor_breaks = NULL,
                     limits = range(values)) +
  scale_fill_manual(values = col_points) +
  scale_color_manual(values = rep.int(col_lines, 3), labels = equations) +
  scale_linetype_manual(values = rep(line_types, each = 2), labels = equations) +
  facet_grid(Method~Scenario) +
  theme_bw() +
  coord_fixed() +
  theme(axis.text = element_blank(),
        axis.title = element_text(size = 13, face = "italic"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        plot.margin = unit(c(0, 0.25, 0.25, 0.25), "line"),
        strip.text = element_text(size = 12)) +
  guides(color = guide_legend(keywidth = unit(1.5, "line")),
         linetype = guide_legend(keywidth = unit(1.5, "line")),
         fill = guide_legend(title.theme = element_text(face = "italic",
                                                        angle = 0)))
# plot to file
pdf("illustration/mediation_estimation.pdf", width = 6.875, height = 5)
print(plt_estimation)
dev.off()


## graphical parameters for plot of bootstrap distributions
expand_mult <- 0.05

## apply bootstrap tests
seed <- 20240927
R <- 10000
df_density <- df_ci <- widths <- NULL
for (scenario in scenarios) {

  # ensure that the same bootstrap samples are used for all parameter
  # settings and all bootstrap tests for maximum comparability
  set.seed(seed)
  n <- if (scenario == scenarios[1]) nrow(data_without) else nrow(data_with)
  indices <- boot_samples(n, R = R)

  # find relevant regression fits
  if (scenario == scenarios[1]) {
    fit_OLS <- fit_OLS_without
    fit_ROBMED <- fit_ROBMED_without
  } else {
    fit_OLS <- fit_OLS_with
    fit_ROBMED <- fit_ROBMED_with
  }

  # apply OLS bootstrap and ROBMED
  test_OLS <- test_mediation(fit_OLS, type = "perc", indices = indices)
  test_ROBMED <- test_mediation(fit_ROBMED, type = "perc", indices = indices)

  # obtain relevant information for density plot
  setup <- setup_density_plot(list(OLS = test_OLS, ROBMED = test_ROBMED))

  # extract data frame for density estimate
  tmp <- data.frame(Scenario = factor(scenario, levels = scenarios),
                    setup$density)
  df_density <- rbind(df_density, tmp)

  # compute range of x-axis limits (to be used for panel sizes)
  x_range <- range(tmp$Indirect)
  x_limits <- x_range + c(-1, 1) * expand_mult * diff(x_range)
  widths <- c(widths, diff(x_limits))

  # extract data frame for confidence interval
  tmp <- data.frame(Scenario = factor(scenario, levels = scenarios),
                    setup$ci)
  df_ci <- rbind(df_ci, tmp)

}

# create plot of bootstrap distribution
plt_boot <- ggplot() +
  geom_density(aes(x = Indirect, y = Density, color = Method),
               data = df_density, stat = "identity") +
  geom_vline(aes(xintercept = Estimate, color = Method),
             data = df_ci) +
  geom_rect(aes(xmin = Lower, xmax = Upper, ymin = -Inf, ymax = Inf,
                fill = Method),
            data = df_ci, color = NA, alpha = 0.2) +
  facet_grid(. ~ Scenario, scale = "free_x") +
  facetted_pos_scales(x = list(
    scale_x_continuous(breaks = c(0.25, 0.75, 1.25),
                       expand = expansion(mult = expand_mult)),
    scale_x_continuous(expand = expansion(mult = expand_mult))
  )) +
  force_panelsizes(cols = widths) +
  labs(title = NULL, x = "Indirect effect") +
  theme_bw() +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 13),
        legend.direction = "vertical",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "line"),
        strip.text = element_text(size = 12))

# plot to file
pdf("illustration/mediation_bootstrap.pdf", width = 6.875, height = 3)
print(plt_boot)
dev.off()
