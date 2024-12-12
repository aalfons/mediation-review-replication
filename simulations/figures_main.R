# *************************************
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# *************************************

## load required packages
library("dplyr")
library("ggplot2")
library("grid")
library("scales")
library("stringr")
library("tidyr")

## load file containing results
file_results <- "simulations/results/results_within_subject.RData"
load(file_results)

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

## methods to include in paper
methods <- c("ols_boot", "ols_causal", "median_boot", "median_causal",
             "winsorized_boot", "ROBMED")

## nice labels to be used for plots
method_labels <- c(ols_boot = "Reg-OLS",
                   ols_causal = "CMA-OLS",
                   median_boot = "Reg-Median",
                   median_causal = "CMA-Median",
                   winsorized_boot = "Winsorized",
                   ROBMED = "ROBMED")
setting_labels <- c(latent = "\"Latent\"", noise = "\"Noise\"",
                    skewness = "\"Skewness\"", heavy_tails = "\"Heavy tails\"",
                    rounding = "\"Rounding\"", outliers = "\"Outliers\"",
                    censoring = "\"Censoring\"")
parameter_labels <- c(ab = "paste(\"Indirect effect \", italic(ab))",
                      correct = "\"Rate of rejection\nwith correct sign\"",
                      reject = "\"Rejection rate\"")

## filter results for selected methods
## (we only compare bootstrap methods, so we always take the bootstrap point
## estimates here)
results <- results %>%
  filter(Method %in% methods)


## loop over different settings
settings <- results %>% distinct(exposure, a, b)
for (i in 1:nrow(settings)) {

  # filter results for current mediation setting
  current_results <- results %>%
    filter(exposure == settings[i, "exposure"],
           a == settings[i, "a"],
           b == settings[i, "b"]) %>%
    mutate(is_latent = (noise_ratio == 0 &
                          transformation == "identity" &
                          rounding_probability == 0 &
                          outlier_probability == 0 &
                          censoring_fraction == 0),
           have_random_noise = (noise_ratio != 0 &
                                  transformation == "identity" &
                                  rounding_probability == 0 &
                                  outlier_probability == 0 &
                                  censoring_fraction == 0),
           have_censoring = (noise_ratio == 0 &
                               transformation == "identity" &
                               rounding_probability == 0 &
                               outlier_probability == 0 &
                               censoring_fraction != 0),
           have_heavy_tails = (noise_ratio == 0 &
                                 transformation == "heavy tails" &
                                 rounding_probability == 0 &
                                 outlier_probability == 0 &
                                 censoring_fraction == 0),
           have_skewness = (noise_ratio == 0 &
                              transformation == "skewness" &
                              rounding_probability == 0 &
                              outlier_probability == 0 &
                              censoring_fraction == 0),
           have_rounding = (noise_ratio == 0 &
                              transformation == "identity" &
                              rounding_probability > 0 &
                              outlier_probability == 0 &
                              censoring_fraction == 0),
           have_outliers = (noise_ratio == 0 &
                              transformation == "identity" &
                              rounding_probability == 0 &
                              outlier_probability > 0 &
                              censoring_fraction == 0),
           Setting = ifelse(is_latent, "latent",
                            ifelse(have_random_noise, "noise",
                                   ifelse(have_heavy_tails, "heavy_tails",
                                          ifelse(have_skewness, "skewness",
                                                 ifelse(have_rounding, "rounding",
                                                        ifelse(have_outliers, "outliers",
                                                               ifelse(have_censoring, "censoring",
                                                                      NA_character_)))))))) %>%
    filter(!is.na(Setting))

  # slightly different plots for mediation and nonmediation
  # rejection rate (with correct sign) is plotted in the bottom row
  if (settings[i, "a"] * settings[i, "b"] == 0) {
    reject <- "reject"
    # aggregate results
    df_reject <- current_results %>%
      group_by(exposure, a, b, Setting, Method) %>%
      summarize(Value = mean(reject_index), .groups = "drop") %>%
      mutate(Parameter = reject)
  } else {
    reject <- "correct"
    # for mediation, compute rate of rejection with correct sign
    current_results <- current_results %>%
      mutate(correct = reject_index & sign(estimate_boot) == sign(a * b))
    # aggregate results
    df_reject <- current_results %>%
      group_by(exposure, a, b, Setting, Method) %>%
      summarize(Value = mean(correct), .groups = "drop") %>%
      mutate(Parameter = reject)
  }

  # boxplots of the point estimates are plotted in the top row
  df_ab <- current_results %>%
    group_by(exposure, a, b, Setting, Method) %>%
    summarize(Min = boxplot.stats(estimate_boot)$stats[1],
              Lower = boxplot.stats(estimate_boot)$stats[2],
              Middle = boxplot.stats(estimate_boot)$stats[3],
              Upper = boxplot.stats(estimate_boot)$stats[4],
              Max = boxplot.stats(estimate_boot)$stats[5],
              .groups = "drop") %>%
    mutate(Parameter = "ab")

  # nice labels for plotting
  df_ab <- df_ab %>%
    mutate(Method = factor(method_labels[Method],
                           levels = method_labels),
           Setting = factor(setting_labels[Setting],
                            levels = setting_labels),
           Parameter = factor(parameter_labels[Parameter],
                              levels = parameter_labels))
  df_reject <- df_reject %>%
    mutate(Method = factor(method_labels[Method],
                           levels = method_labels),
           Setting = factor(setting_labels[Setting],
                            levels = setting_labels),
           Parameter = factor(parameter_labels[Parameter],
                              levels = parameter_labels))

  # data frame for adding true value
  df_true <- data.frame(Parameter = parameter_labels["ab"],
                        Value = settings$a[i] * settings$b[i])
  if (settings[i, "a"] * settings[i, "b"] == 0) {
    df_alpha <- data.frame(Parameter = parameter_labels[reject],
                           Value = 1 - level)
    df_true <- rbind(df_true, df_alpha)
  }
  df_true$Parameter <- factor(df_true$Parameter, levels = parameter_labels)

  # expand axis limits of power facet with blank layer
  df_ylim <- data.frame(Parameter = factor(unname(parameter_labels[reject]),
                                           levels = parameter_labels),
                        Method = unname(method_labels[1]), Value = 0:1)

  # colors and plot symbols
  if (settings[i, "exposure"] == "binary") {
    colors <- c("#F8766D", "#F8766D", "#619CFF", "#619CFF", "#F564E3", "#00BFC4")
    symbols <- c(21, 22, 21, 22, 21, 21)
  } else {
    colors <- c("#F8766D", "#619CFF", "#F564E3", "#00BFC4")
    symbols <- rep(21, 4)
  }

  # plot results
  p <- ggplot() +
    geom_boxplot(mapping = aes(x = Method, ymin = Min, lower = Lower,
                               middle = Middle, upper = Upper, ymax = Max,
                               fill = Method),
                 data = df_ab, stat = "identity", show.legend = FALSE) +
    geom_hline(aes(yintercept = Value), data = df_true) +
    geom_blank(mapping = aes(x = Method, y = Value), data = df_ylim) +
    geom_point(mapping = aes(x = Method, y = Value, fill = Method,
                             shape = Method),
               data = df_reject, size = 3, show.legend = FALSE) +
    facet_grid(Parameter ~ Setting, scales = "free_y",
               labeller = "label_parsed") +
    scale_shape_manual(values = symbols) +
    labs(x = NULL, y = NULL) + theme_bw() +
    theme(axis.text.x = element_text(size = 11, angle = 90,
                                     hjust = 1, vjust = 0.5),
          axis.text.y = element_text(size = 11),
          strip.text = element_text(size = 12))

  # save plot to file
  suffix <- paste(settings[i, "exposure"],
                   if (settings[i, "b"] == 0) "nonmediation" else "mediation",
                   sep = "_")
  if (settings[i, "b"] == 0) {
    file_plot <- "simulations/figures/figure_%s.pdf"
    # file_plot_bw <- "simulations/figures/figure_bw_%s.pdf"
  } else {
    file_plot <- "simulations/figures/figure_%s.pdf"
    # file_plot_bw <- "simulations/figures/figure_bw_%s.pdf"
  }
  # in color
  pdf(file = sprintf(file_plot, suffix),
      width = 8.5, height = 4.25)
  print(p + scale_fill_manual(values = colors))
  dev.off()
  # # in grayscale
  # pdf(file = sprintf(file_plot_bw, suffix),
  #     width = 8.5, height = 4.25)
  # print(p + scale_fill_manual(values = col2gray(colors)))
  # dev.off()

}
