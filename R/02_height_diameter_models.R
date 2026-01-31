#' ---
#' title: "Height-Diameter Relationship Models"
#' author: "Ogundipe Oluwadamilola"
#' description: "Non-linear regression models for tree height prediction"
#' species: "Triplochiton scleroxylon (Obeche)"
#' ---

# =============================================================================
# HEIGHT-DIAMETER MODELS FOR TRIPLOCHITON SCLEROXYLON
# Non-linear regression approach for height prediction from DBH
# =============================================================================

if (!require("minpack.lm")) install.packages("minpack.lm")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")

library(minpack.lm)
library(ggplot2)
library(dplyr)

# -----------------------------------------------------------------------------
# 1. SAMPLE HEIGHT-DIAMETER DATA
# -----------------------------------------------------------------------------
# Note: This is example data structure. Replace with your actual measurements.
# Typically collected from a subsample of trees in inventory plots.

# Generate realistic H-D data based on tropical hardwood relationships
set.seed(42)
n <- 150  # Number of trees with height measurements

# Simulate data following typical obeche growth patterns
dbh <- runif(n, 30, 150)  # DBH range 30-150 cm
# True relationship with some noise
height <- 1.3 + (35 - 1.3) * (1 - exp(-0.025 * dbh)) + rnorm(n, 0, 2)

hd_data <- data.frame(
  tree_id = 1:n,
  dbh_cm = round(dbh, 1),
  height_m = round(pmax(height, 5), 1)  # Minimum height 5m
)

cat("=== HEIGHT-DIAMETER DATA SUMMARY ===\n")
cat(sprintf("Number of trees: %d\n", nrow(hd_data)))
cat(sprintf("DBH range: %.1f - %.1f cm\n", min(hd_data$dbh_cm), max(hd_data$dbh_cm)))
cat(sprintf("Height range: %.1f - %.1f m\n", min(hd_data$height_m), max(hd_data$height_m)))

# -----------------------------------------------------------------------------
# 2. FIT NON-LINEAR HEIGHT-DIAMETER MODELS
# -----------------------------------------------------------------------------

#' Model 1: Inverse/Hyperbolic Function
#' H = a + b/D
#' Simple, often works well for tropical species

model1 <- tryCatch({
  nlsLM(height_m ~ a + b/dbh_cm, 
        data = hd_data,
        start = list(a = 30, b = -100))
}, error = function(e) NULL)

#' Model 2: Chapman-Richards Function
#' H = a * (1 - exp(-b * D))^c
#' Flexible, commonly used in forestry

model2 <- tryCatch({
  nlsLM(height_m ~ a * (1 - exp(-b * dbh_cm))^c,
        data = hd_data,
        start = list(a = 35, b = 0.02, c = 1))
}, error = function(e) NULL)

#' Model 3: Wykoff (Logistic-type)
#' H = 1.3 + exp(a + b/(D + 1))
#' Good asymptotic behavior

model3 <- tryCatch({
  nlsLM(height_m ~ 1.3 + exp(a + b/(dbh_cm + 1)),
        data = hd_data,
        start = list(a = 3.5, b = -50))
}, error = function(e) NULL)

# -----------------------------------------------------------------------------
# 3. MODEL EVALUATION
# -----------------------------------------------------------------------------

#' Calculate model fit statistics
#' @param model nls model object
#' @param data original data frame
#' @param model_name name for display
evaluate_model <- function(model, data, model_name) {
  if (is.null(model)) {
    return(data.frame(
      model = model_name,
      RMSE = NA, AIC = NA, BIC = NA, R2 = NA
    ))
  }
  
  pred <- predict(model)
  obs <- data$height_m
  n <- length(obs)
  k <- length(coef(model))
  
  # Root Mean Square Error
  rmse <- sqrt(mean((obs - pred)^2))
  
  # R-squared
  ss_res <- sum((obs - pred)^2)
  ss_tot <- sum((obs - mean(obs))^2)
  r2 <- 1 - ss_res/ss_tot
  
  # AIC and BIC
  aic <- AIC(model)
  bic <- BIC(model)
  
  data.frame(
    model = model_name,
    RMSE = round(rmse, 3),
    AIC = round(aic, 2),
    BIC = round(bic, 2),
    R2 = round(r2, 4)
  )
}

# Evaluate all models
results <- rbind(
  evaluate_model(model1, hd_data, "Inverse (H = a + b/D)"),
  evaluate_model(model2, hd_data, "Chapman-Richards"),
  evaluate_model(model3, hd_data, "Wykoff")
)

cat("\n=== MODEL COMPARISON ===\n")
print(results, row.names = FALSE)

# Best model selection
best_idx <- which.min(results$AIC)
cat(sprintf("\n→ Best model by AIC: %s\n", results$model[best_idx]))

# -----------------------------------------------------------------------------
# 4. PARAMETER ESTIMATES
# -----------------------------------------------------------------------------

cat("\n=== PARAMETER ESTIMATES ===\n")

if (!is.null(model1)) {
  cat("\nModel 1 (Inverse): H = a + b/D\n")
  print(summary(model1)$coefficients)
}

if (!is.null(model2)) {
  cat("\nModel 2 (Chapman-Richards): H = a * (1 - exp(-b*D))^c\n")
  print(summary(model2)$coefficients)
}

if (!is.null(model3)) {
  cat("\nModel 3 (Wykoff): H = 1.3 + exp(a + b/(D+1))\n")
  print(summary(model3)$coefficients)
}

# -----------------------------------------------------------------------------
# 5. VISUALIZATION
# -----------------------------------------------------------------------------

# Create prediction data
pred_dbh <- seq(min(hd_data$dbh_cm), max(hd_data$dbh_cm), length.out = 100)
pred_df <- data.frame(dbh_cm = pred_dbh)

# Add predictions from each model
if (!is.null(model1)) pred_df$Inverse <- predict(model1, pred_df)
if (!is.null(model2)) pred_df$ChapmanRichards <- predict(model2, pred_df)
if (!is.null(model3)) pred_df$Wykoff <- predict(model3, pred_df)

# Reshape for plotting
library(tidyr)
pred_long <- pred_df %>%
  pivot_longer(cols = -dbh_cm, names_to = "Model", values_to = "height_m")

# Plot
p <- ggplot() +
  geom_point(data = hd_data, aes(x = dbh_cm, y = height_m),
             alpha = 0.5, color = "gray40") +
  geom_line(data = pred_long, aes(x = dbh_cm, y = height_m, color = Model),
            linewidth = 1) +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Height-Diameter Models for Triplochiton scleroxylon",
    subtitle = "Comparison of three non-linear regression models",
    x = "Diameter at Breast Height (cm)",
    y = "Total Height (m)",
    color = "Model"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14)
  )

ggsave("figures/height_diameter_models.png", p, 
       width = 10, height = 7, dpi = 300)

cat("\n✓ Plot saved to figures/height_diameter_models.png\n")

# -----------------------------------------------------------------------------
# 6. RESIDUAL DIAGNOSTICS
# -----------------------------------------------------------------------------

if (!is.null(model2)) {  # Using Chapman-Richards as example
  hd_data$predicted <- predict(model2)
  hd_data$residual <- hd_data$height_m - hd_data$predicted
  
  p_resid <- ggplot(hd_data, aes(x = predicted, y = residual)) +
    geom_point(alpha = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    geom_smooth(method = "loess", se = FALSE, color = "blue") +
    labs(
      title = "Residual Plot - Chapman-Richards Model",
      x = "Predicted Height (m)",
      y = "Residual (m)"
    ) +
    theme_minimal()
  
  ggsave("figures/residual_plot_height.png", p_resid,
         width = 8, height = 6, dpi = 300)
  
  cat("✓ Residual plot saved\n")
}

# -----------------------------------------------------------------------------
# 7. EXPORT RESULTS
# -----------------------------------------------------------------------------

write.csv(results, "data/height_model_comparison.csv", row.names = FALSE)
write.csv(hd_data, "data/height_diameter_data.csv", row.names = FALSE)

cat("\n=== ANALYSIS COMPLETE ===\n")

# -----------------------------------------------------------------------------
# 8. PREDICTION FUNCTION FOR PRACTICAL USE
# -----------------------------------------------------------------------------

#' Predict tree height from DBH
#' @param dbh Diameter at breast height in cm
#' @param model Choice of model: "inverse", "chapman", "wykoff"
#' @return Predicted height in meters
predict_height <- function(dbh, model = "chapman") {
  if (model == "inverse" && !is.null(model1)) {
    coef1 <- coef(model1)
    return(coef1["a"] + coef1["b"]/dbh)
  } else if (model == "chapman" && !is.null(model2)) {
    coef2 <- coef(model2)
    return(coef2["a"] * (1 - exp(-coef2["b"] * dbh))^coef2["c"])
  } else if (model == "wykoff" && !is.null(model3)) {
    coef3 <- coef(model3)
    return(1.3 + exp(coef3["a"] + coef3["b"]/(dbh + 1)))
  }
  return(NA)
}

# Example usage
cat("\nExample predictions:\n")
test_dbh <- c(40, 60, 80, 100)
for (d in test_dbh) {
  cat(sprintf("DBH = %d cm → Height = %.1f m\n", d, predict_height(d)))
}
