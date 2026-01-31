#' ---
#' title: "Diameter Distribution Analysis for Triplochiton scleroxylon"
#' author: "Ogundipe Oluwadamilola"
#' description: "Weibull distribution fitting for forest inventory data"
#' location: "Gambari Forest Reserve, Southwest Nigeria"
#' ---

# =============================================================================
# DIAMETER DISTRIBUTION ANALYSIS USING WEIBULL FUNCTION
# Species: Triplochiton scleroxylon (Obeche)
# Location: Gambari Forest Reserve, Nigeria
# =============================================================================

# Load required packages
if (!require("fitdistrplus")) install.packages("fitdistrplus")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")

library(fitdistrplus)
library(ggplot2)
library(dplyr)

# -----------------------------------------------------------------------------
# 1. DATA LOADING
# -----------------------------------------------------------------------------

# Read diameter data
dbh_data <- read.csv("data/obeche_dbh_gambari.csv", comment.char = "#")

cat("=== OBECHE DIAMETER DATA SUMMARY ===\n")
cat(sprintf("Number of trees: %d\n", nrow(dbh_data)))
cat(sprintf("DBH range: %.1f - %.1f cm\n", min(dbh_data$dbh_cm), max(dbh_data$dbh_cm)))
cat(sprintf("Mean DBH: %.2f cm\n", mean(dbh_data$dbh_cm)))
cat(sprintf("Std Dev: %.2f cm\n", sd(dbh_data$dbh_cm)))

# -----------------------------------------------------------------------------
# 2. WEIBULL DISTRIBUTION FITTING (3-PARAMETER)
# -----------------------------------------------------------------------------

#' Fit 3-parameter Weibull distribution
#' The 3-parameter Weibull is: f(x) = (k/λ) * ((x-θ)/λ)^(k-1) * exp(-((x-θ)/λ)^k)
#' where:
#'   k = shape parameter
#'   λ = scale parameter  
#'   θ = location parameter (minimum threshold)

# Shift data for 3-parameter Weibull (location parameter)
dbh_values <- dbh_data$dbh_cm
theta <- min(dbh_values) - 0.1  # Location parameter (threshold)
dbh_shifted <- dbh_values - theta

# Fit 2-parameter Weibull to shifted data
weibull_fit <- fitdist(dbh_shifted, "weibull", method = "mle")

cat("\n=== WEIBULL PARAMETERS (3-PARAMETER) ===\n")
cat(sprintf("Shape (k): %.4f\n", weibull_fit$estimate["shape"]))
cat(sprintf("Scale (λ): %.4f\n", weibull_fit$estimate["scale"]))
cat(sprintf("Location (θ): %.4f cm\n", theta))

# Goodness of fit statistics
gof <- gofstat(weibull_fit)
cat("\n=== GOODNESS OF FIT ===\n")
cat(sprintf("Kolmogorov-Smirnov statistic: %.4f\n", gof$ks))
cat(sprintf("AIC: %.2f\n", weibull_fit$aic))
cat(sprintf("BIC: %.2f\n", weibull_fit$bic))

# -----------------------------------------------------------------------------
# 3. CREATE DIAMETER CLASS DISTRIBUTION
# -----------------------------------------------------------------------------

# Define diameter classes (30cm intervals as in original study)
breaks <- seq(30, 300, by = 30)
dbh_data$class <- cut(dbh_data$dbh_cm, breaks = breaks, right = FALSE,
                      labels = paste0(head(breaks, -1), "-", tail(breaks, -1)))

# Count observed frequencies
observed_freq <- dbh_data %>%
  group_by(class) %>%
  summarise(
    observed = n(),
    midpoint = mean(as.numeric(gsub("-.*", "", as.character(class)))) + 15
  ) %>%
  filter(!is.na(class))

# Calculate expected frequencies from Weibull
shape <- weibull_fit$estimate["shape"]
scale <- weibull_fit$estimate["scale"]
n_trees <- nrow(dbh_data)

observed_freq$expected <- sapply(1:nrow(observed_freq), function(i) {
  lower <- as.numeric(gsub("-.*", "", as.character(observed_freq$class[i]))) - theta
  upper <- lower + 30
  prob <- pweibull(upper, shape, scale) - pweibull(lower, shape, scale)
  return(round(prob * n_trees, 2))
})

cat("\n=== DIAMETER CLASS DISTRIBUTION ===\n")
print(observed_freq, n = 20)

# -----------------------------------------------------------------------------
# 4. VISUALIZATION
# -----------------------------------------------------------------------------

# Create comparison plot
p1 <- ggplot(observed_freq, aes(x = midpoint)) +
  geom_bar(aes(y = observed, fill = "Observed"), stat = "identity", 
           alpha = 0.7, width = 25) +
  geom_line(aes(y = expected, color = "Weibull Expected"), 
            linewidth = 1.2) +
  geom_point(aes(y = expected, color = "Weibull Expected"), size = 3) +
  scale_fill_manual(values = c("Observed" = "steelblue")) +
  scale_color_manual(values = c("Weibull Expected" = "darkred")) +
  labs(
    title = "Diameter Distribution of Triplochiton scleroxylon",
    subtitle = "Gambari Forest Reserve, Nigeria (n = 499 trees)",
    x = "Diameter at Breast Height (cm)",
    y = "Frequency",
    fill = "", color = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "gray40")
  )

# Save plot
ggsave("figures/diameter_distribution_weibull.png", p1, 
       width = 10, height = 7, dpi = 300)

cat("\n✓ Plot saved to figures/diameter_distribution_weibull.png\n")

# -----------------------------------------------------------------------------
# 5. CHI-SQUARE GOODNESS OF FIT TEST
# -----------------------------------------------------------------------------

# Combine classes with expected < 5 for chi-square test
chi_data <- observed_freq %>%
  mutate(
    group = cumsum(expected >= 5 | lag(expected, default = 0) < 5)
  )

chi_result <- chisq.test(observed_freq$observed, p = observed_freq$expected/sum(observed_freq$expected))

cat("\n=== CHI-SQUARE GOODNESS OF FIT TEST ===\n")
cat(sprintf("Chi-square statistic: %.4f\n", chi_result$statistic))
cat(sprintf("Degrees of freedom: %d\n", chi_result$parameter))
cat(sprintf("P-value: %.4f\n", chi_result$p.value))

if (chi_result$p.value > 0.05) {
  cat("→ Weibull distribution provides adequate fit (p > 0.05)\n")
} else {
  cat("→ Significant deviation from Weibull distribution (p < 0.05)\n")
}

# -----------------------------------------------------------------------------
# 6. EXPORT RESULTS
# -----------------------------------------------------------------------------

# Save fitted parameters
params <- data.frame(
  parameter = c("shape_k", "scale_lambda", "location_theta", "n_trees", "aic", "bic"),
  value = c(shape, scale, theta, n_trees, weibull_fit$aic, weibull_fit$bic)
)
write.csv(params, "data/weibull_parameters.csv", row.names = FALSE)

# Save diameter class comparison
write.csv(observed_freq, "data/diameter_class_distribution.csv", row.names = FALSE)

cat("\n✓ Results exported to data/ folder\n")
cat("\n=== ANALYSIS COMPLETE ===\n")
