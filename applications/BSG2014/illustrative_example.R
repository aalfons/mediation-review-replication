# ************************************
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ************************************


# load packages and data
library("cellWise")
library("dplyr")
library("ggplot2")
library("mediation")
library("robmed")
library("robmedExtra")
data("BSG2014")

# illustrative application
x <- "ValueDiversity"
y <- "TeamPerformance"
m <- "TaskConflict"

# seed tp be used for the random number generator
seed <- 20250423

# OLS bootstrap
set.seed(seed)
ols_boot <-  test_mediation(BSG2014, x = x, y = y, m = m, test = "boot",
                            method = "regression", robust = FALSE)
summary(ols_boot)
p_value(ols_boot, parm = "indirect")

# ROBMED
set.seed(seed)
robust_boot <- test_mediation(BSG2014, x = x, y = y, m = m, test = "boot",
                              method = "regression", robust = TRUE)
summary(robust_boot, plot = FALSE)
p_value(robust_boot, parm = "indirect")

# create LaTeX table
to_latex(list("Reg-OLS" = ols_boot, "ROBMED" = robust_boot))

# generate diagnostic plot (hack for facet labels)
setup_weight <- setup_weight_plot(robust_boot)
setup_weight$data <- setup_weight$data %>%
  mutate(Outcome = recode_factor(Outcome,
                                 "TaskConflict" = "Task conflict",
                                 "TeamPerformance" = "Team performance"))
plt_weight <- weight_plot(setup_weight) +
  scale_color_manual("", values = c("black", "#00BFC4")) +
  labs(y = "       Percentage of observations with weight lower than threshold") +
  theme_bw() +
  theme(axis.text = element_text(size = 11),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 12),
        legend.position = "top",
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 12))

# save plot to file
pdf("applications/BSG2014/diagnostic_plot_BSG2014.pdf",
    width = 4.5375, height = 4.5375)
print(plt_weight)
dev.off()


## plot of standardized residuals against exposure

# create data frame for points
fit <- robust_boot$fit$fit_m
# fit <- robust_boot$fit$fit_m[[1]]
df_plot <- data.frame(
  X = BSG2014[, x, drop = TRUE],
  StdResidual = residuals(fit) / fit$scale,
  weight = weights(fit, type = "robustness")
)

# create data frame for horizontal reference lines indicating outlying residuals
df_hline <- data.frame(
  yintercept = c(-3, 3)
)

# create data frame for vertical reference lines of quantiles
df_vline <- data.frame(
  xintercept = quantile(BSG2014[, x, drop = TRUE], probs = c(0.1, 0.9))
)

# generate plot
plt_leverage <- ggplot() +
  geom_hline(aes(yintercept = yintercept), data = df_hline,
             color = "darkgray") +
  geom_vline(aes(xintercept = xintercept), data = df_vline,
             color = "darkgray") +
  geom_point(aes(x = X, y = StdResidual, fill = weight),
             data = df_plot, pch = 21) +
  scale_fill_gradient("Weight    ", low = "white", high = "black",
                      limits = c(0, 1)) +
  labs(x = "Value diversity", y = "Standardized residual") +
  theme_bw() +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 13),
        legend.position = "top",
        legend.key.width = unit(25, "pt"),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 13),
        strip.text = element_text(size = 12))

# save plot to file
pdf("applications/BSG2014/leverage_plot_BSG2014.pdf",
    width = 4.5375, height = 4.5375)
print(plt_leverage)
dev.off()

# identify leverage points
leverage <- which(df_plot$X < df_vline[1, "xintercept"] &
                    df_plot$StdResidual > df_hline[2, "yintercept"])


## robustness checks

# Box-Cox transformation of M followed by OLS bootstrap
transformation <- transfo(BSG2014[, m, drop = FALSE], type = "BC",
                          robust = FALSE, standardize = FALSE,
                          checkPars = list(coreOnly = TRUE))
BSG2014_bc <- cbind(BSG2014[, c(x, y)], transformation$Y)
set.seed(seed)
bc_boot <- test_mediation(BSG2014_bc, x = x, y = y, m = m, test = "boot",
                          method = "regression", robust = FALSE)
summary(bc_boot)
p_value(bc_boot, parm = "indirect")

# robustness check: OLS bootstrap without leverage points
set.seed(seed)
ols_boot_cleaned <-  test_mediation(BSG2014[-leverage, ], x = x, y = y, m = m,
                                    test = "boot", method = "regression",
                                    robust = FALSE)
summary(ols_boot_cleaned)
p_value(ols_boot_cleaned, parm = "indirect")

# create LaTeX table: needs a hack to make it work
summary_bc_boot <- summary(bc_boot)
summary_ols_boot_cleaned <- summary(ols_boot_cleaned)
summary_ols_boot_cleaned$summary$n <- summary_bc_boot$summary$n
to_latex(list("Box-Cox transformation" = summary_bc_boot,
              "Without leverage points" = summary_ols_boot_cleaned))


## potential outcome approach

# fit OLS regressions with all data points
fit_m <- lm(TaskConflict ~ ValueDiversity, data = BSG2014)
fit_y <- lm(TeamCommitment ~ TaskConflict + ValueDiversity, data = BSG2014)
# apply bootstrap test
set.seed(seed)
po_boot <- mediate(fit_m, fit_y, sims = 5000, boot = TRUE,
                   treat = x, mediator = m)
summary(po_boot)
# conduct sensitivity analysis on sequential ignorability assumption
po_sens <- medsens(po_boot)
summary(po_sens)

# fit OLS regressions without leverage points
fit_m_cleaned <- lm(TaskConflict ~ ValueDiversity,
                    data = BSG2014[-leverage, ])
fit_y_cleaned <- lm(TeamCommitment ~ TaskConflict + ValueDiversity,
                    data = BSG2014[-leverage, ])
# apply bootstrap test
set.seed(seed)
po_boot_cleaned <- mediate(fit_m_cleaned, fit_y_cleaned, sims = 5000,
                   boot = TRUE, treat = x, mediator = m)
summary(po_boot_cleaned)
# conduct sensitivity analysis on sequential ignorability assumption
po_sens_cleaned <- medsens(po_boot_cleaned)
summary(po_sens_cleaned)

# extract relevant information as in plot.medsens()
df_sens <- rbind(
  data.frame(Data = "Full sample", rho = po_sens$rho, ACME = po_sens$d0,
             Lower = po_sens$lower.d0, Upper = po_sens$upper.d0),
  data.frame(Data = "Without leverage points", rho = po_sens_cleaned$rho,
             ACME = po_sens_cleaned$d0, Lower = po_sens_cleaned$lower.d0,
             Upper = po_sens_cleaned$upper.d0)
)

df_est <- df_sens %>% filter(rho == 0)

# generate plot
plt_sens <- ggplot() +
  geom_hline(yintercept = 0, color = "black") +
  geom_vline(xintercept = 0, color = "black") +
  geom_hline(aes(yintercept = ACME), data = filter(df_sens, rho == 0),
             color = "darkgray") +
  geom_ribbon(aes(x = rho, ymin = Lower, ymax = Upper), data = df_sens,
              alpha = 0.3) +
  geom_line(aes(x = rho, y = ACME), data = df_sens) +
  facet_grid(. ~ Data) +
  labs(x = expression("Sensitivity parameter"~italic(rho))) +
  theme_bw() +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 13),
        strip.text = element_text(size = 12))

# save plot to file
pdf("applications/BSG2014/sensitivity_plot_BSG2014.pdf",
    width = 9.167, height = 3.875)
print(plt_sens)
dev.off()

# check that the plot is indeed a nicer version of plot.medsens()
plot(po_sens, ylim = c(-1.5, 1.1))
plot(po_sens_cleaned, ylim = c(-1.5, 1.1))
