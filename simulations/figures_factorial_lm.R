# *************************************
# Author: Dan Schley
#         Erasmus University Rotterdam
# *************************************

# Load necessary libraries
library("tidyverse")
library("broom")
library("stringr")
library("purrr")


# Important, clear the workspace before running the loops below otherwise rows will stack duplicates.

# Load data
load("simulations/results/results_within_subject.RData")

# Assuming 'results' is your dataframe and necessary packages are loaded

results <- results %>%
  mutate(
    heavy_tails = ifelse(transformation == "heavy tails", 1, 0),
    skewness = ifelse(transformation == "skewness", 1, 0))

# Define moderators and set up combinations
moderators <- c("b", "noise_ratio", "heavy_tails", "skewness",
                "rounding_probability", "outlier_probability",
                "censoring_fraction", "Method")

results <- results %>%
  mutate(across(all_of(moderators), as.factor))

fixed_moderators <- c("b", "Method")
other_moderators <- setdiff(moderators, fixed_moderators)

# Create 'estimate_error' which is just the estimated indirect effect minus true indirect effect
results <- results %>%
  mutate(estimate_error = ifelse(b == 0,
                                 estimate_boot - 0,
                                 estimate_boot - 0.16))

# Get all unique combinations of 'Method' and 'b' as characters
method_b_combinations <- results %>%
  distinct(Method, b) %>%
  mutate(Method = as.character(Method),
         b = as.character(b)) %>%
  arrange(Method, b)

# Initialize a list to store the results
results_list <- list()

max_interaction_level <- 5 # Even though there are 6 variables a 5-way interaction is all that is possible since we exclude interactions involving both 'heavy_tails' and 'skewness' as they are not actually crossed in the design (i.e., these are dummies of a variable with 3 levels)

# Loop through each combination of 'Method' and 'b'
for (k in seq_len(nrow(method_b_combinations))) {
  current_Method <- method_b_combinations$Method[k]
  current_b <- method_b_combinations$b[k]

  # Subset data for current 'Method' and 'b'
  data_sub <- results %>%
    filter(as.character(Method) == current_Method,
           as.character(b) == current_b)

  # Drop unused factor levels
  data_sub <- droplevels(data_sub)

  # Check if there are enough observations
  if (nrow(data_sub) < 2) {
    next
  }

  # This code was here to deal with the fact that the 'exposure' variable had no variability for the Imai simulations. Because we are not modeling the exposure variable, this code is not necessary.
  usable_moderators <- other_moderators[sapply(data_sub[other_moderators], function(x) nlevels(x) >=2)]

  if (length(usable_moderators) == 0) {
    # No usable moderators, skip this iteration
    next
  }

  # Loop over interaction levels from 1 to max_interaction_level
  for (interaction_level in 1:max_interaction_level) {
    # Generate all combinations of 'usable_moderators' up to the current interaction level
    all_combinations <- unlist(lapply(1:interaction_level, function(comb_size) {
      if (length(usable_moderators) >= comb_size) {
        combn(usable_moderators, comb_size, simplify = FALSE)
      } else {
        list()
      }
    }), recursive = FALSE)

    if (length(all_combinations) == 0) {
      # No combinations possible at this interaction level
      next
    }

    # Exclude combinations involving both 'heavy_tails' and 'skewness'
    all_combinations <- all_combinations[!sapply(all_combinations, function(terms) {
      all(c('heavy_tails', 'skewness') %in% terms)
    })]

    # Generate terms
    all_terms <- sapply(all_combinations, function(terms) paste(terms, collapse = ":"))

    if (length(all_terms) == 0) {
      # No terms to include in the model
      next
    }

    # Build formula string
    formula_str <- paste("estimate_error ~", paste(all_terms, collapse = " + "))

    # Fit the model
    model <- lm(as.formula(formula_str), data = data_sub)

    # Extract coefficients, confidence intervals, and t-values
    tidy_model <- broom::tidy(model, conf.int = TRUE)

    # Compute interaction order for each term
    tidy_model <- tidy_model %>%
      mutate(
        interaction_order = sapply(term, function(x) {
          if (x == "(Intercept)") {
            return(0)
          } else {
            return(length(unlist(strsplit(x, ":"))))
          }
        }),
        focal_level = ifelse(interaction_order == interaction_level, 1, 0)
      )

    # Add 'Method', 'b', and 'interaction_level' to the results
    tidy_model <- tidy_model %>%
      mutate(Method = current_Method,
             b = current_b,
             interaction_level = interaction_level)

    # Store the result in the list
    results_list[[length(results_list) + 1]] <- tidy_model
  }
}

# Combine all the results into a single dataframe
final_results <- bind_rows(results_list)




# Define replacements and the primary custom order for variable names
combo_replacements <- c(
  "noise_ratio" = "Noise",
  "skewness" = "Skewness",
  "heavy_tails" = "Heavy tails",
  "rounding_probability" = "Rounding",
  "outlier_probability" = "Outliers",
  "censoring_fraction" = "Censoring"
)

# Apply replacements iteratively to account for patterns with numeric suffixes
plot_results <- final_results %>%
  filter(focal_level == 1) %>%
  mutate(
    term_relabel = term  # Start with the original term
  )

# Function to relabel terms, including interactions
relabel_term <- function(term) {
  # Split the term into individual components (split by ":")
  components <- str_split(term, ":", simplify = TRUE)

  # Replace components based on combo_replacements, handling unmatched terms
  relabeled_components <- map_chr(components, ~ {
    if (.x %in% names(combo_replacements)) {
      combo_replacements[.x]
    } else {
      .x  # Keep the original term if no replacement is found
    }
  })

  # Combine relabeled components with " X " for interactions
  paste(relabeled_components, collapse = " X ")
}



# Create new labels for terms (including interactions)
plot_results <- plot_results %>%
  mutate(
    term_relabel = map_chr(term, relabel_term),
    b = case_when(
      b == 0 ~ "Without mediation",
      b == 0.4 ~ "With mediation",
      TRUE ~ as.character(b)
    ),
    b = factor(b, levels = c("Without mediation", "With mediation"))
  )

# Convert to factor with custom order
custom_order <- c("Treatment Variable", "Noise", "Skewness", "Heavy tails", "Rounding", "Outliers", "Censoring")
plot_results <- plot_results %>%
  mutate(
    term_relabel = factor(term_relabel, levels = unique(term_relabel))
  )


## plot only methods without power transformations

# Reorder Method with custom labels
method_order <- c("ols_boot", "ols_causal", "median_boot", "median_causal",
                  "winsorized_boot", "ROBMED")
method_labels <- c("Reg-OLS", "PO-OLS", "Reg-Median", "PO-Median",
                   "Winsorized", "ROBMED")
plot_results_main <- plot_results %>%
  filter(Method %in% method_order) %>%
  mutate(Method = factor(Method, levels = method_order, labels = method_labels))


# Group by interaction_level to calculate boundaries and label positions
interaction_levels <- plot_results_main %>%
  group_by(interaction_level) %>%
  summarize(
    x_start = min(as.numeric(term_relabel)),
    x_end = max(as.numeric(term_relabel)),
    x_position = mean(as.numeric(term_relabel))
  )

# Create the plot with updated labels and aesthetics
interaction_plot <-
  ggplot(plot_results_main,
         aes(x = term_relabel, y = estimate, ymin = conf.low,
             ymax = conf.high)) +
  geom_vline(xintercept = interaction_levels$x_end + 0.5,
             linetype = "dashed", color = "gray")+
  geom_point(aes(color = grepl("outlier_probability0.02", term)), size = 1) + # any term with outlier_probability0.02 is colored red
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"),
                     guide = "none") + # Map TRUE to red, FALSE to black
  facet_grid(b ~ Method) +                   # Facet by relabeled "b" and Method
  labs(x = NULL,
       y = expression("Bias in estimated indirect effect"~italic(ab))) +  # Relabel y-axis
  theme_bw() +                          # Use a minimal theme
  theme(
    axis.title = element_text(size = 12),   # Customize y-axis title size
    axis.text.x = element_blank(),          # Remove x-axis labels
    axis.text.y = element_text(size = 11),  # Customize y-axis text size
    axis.ticks.x = element_blank(),         # Remove tick labels
    panel.grid.major.x = element_blank(),   # Remove major gridlines on x-axis
    panel.grid.minor = element_blank(),     # Remove minor gridlines
    strip.text = element_text(size = 12)    # Customize facet labels
  )

# Print the plot
print(interaction_plot)

# save plot to file
file_plot <- "simulations/figures/figure_factorial_lm.pdf"

# in color
pdf(file = file_plot, width = 8.5, height = 4)
print(interaction_plot)
dev.off()


## plot OLS bootstrap and ROBMED without and with transformations

# Reorder Method with custom labels
method_order <- c("ols_boot", "yj_boot", "ROBMED", "ryj_boot")
method_labels <- c("Reg-OLS", "YJ-Reg-OLS", "ROBMED", "rYJ-ROBMED")
plot_results_transformations <- plot_results %>%
  filter(Method %in% method_order) %>%
  mutate(Method = factor(Method, levels = method_order, labels = method_labels))


# Group by interaction_level to calculate boundaries and label positions
interaction_levels <- plot_results_transformations %>%
  group_by(interaction_level) %>%
  summarize(
    x_start = min(as.numeric(term_relabel)),
    x_end = max(as.numeric(term_relabel)),
    x_position = mean(as.numeric(term_relabel))
  )

# Create the plot with updated labels and aesthetics
interaction_plot_transformations <-
  ggplot(plot_results_transformations,
         aes(x = term_relabel, y = estimate, ymin = conf.low,
             ymax = conf.high)) +
  geom_vline(xintercept = interaction_levels$x_end + 0.5,
             linetype = "dashed", color = "gray")+
  geom_point(aes(color = grepl("outlier_probability0.02", term)), size = 1) + # any term with outlier_probability0.02 is colored red
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"),
                     guide = "none") + # Map TRUE to red, FALSE to black
  facet_grid(b ~ Method) +                   # Facet by relabeled "b" and Method
  labs(x = NULL,
       y = expression("Bias in estimated indirect effect"~italic(ab))) +  # Relabel y-axis
  theme_bw() +                          # Use a minimal theme
  theme(
    axis.title = element_text(size = 12),   # Customize y-axis title size
    axis.text.x = element_blank(),          # Remove x-axis labels
    axis.text.y = element_text(size = 11),  # Customize y-axis text size
    axis.ticks.x = element_blank(),         # Remove tick labels
    panel.grid.major.x = element_blank(),   # Remove major gridlines on x-axis
    panel.grid.minor = element_blank(),     # Remove minor gridlines
    strip.text = element_text(size = 12)    # Customize facet labels
  )

# Print the plot
print(interaction_plot_transformations)

# save plot to file
file_plot <- "simulations/figures/figure_factorial_lm_transformations.pdf"

# in color
pdf(file = file_plot, width = 8.5, height = 4)
print(interaction_plot_transformations)
dev.off()
