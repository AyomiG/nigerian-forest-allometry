#' ---
#' title: "Volume Estimation Models"
#' author: "Ogundipe Oluwadamilola"
#' description: "Non-linear regression models for tree volume prediction"
#' species: "Triplochiton scleroxylon (Obeche)"
#' ---

# =============================================================================
# VOLUME MODELS FOR TRIPLOCHITON SCLEROXYLON
# Individual tree volume estimation using stem measurements
# =============================================================================

if (!require("minpack.lm")) install.packages("minpack.lm")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")

library(minpack.lm)
library(ggplot2)
library(dplyr)

# -----------------------------------------------------------------------------
# 1. STEM PROFILE DATA
# -----------------------------------------------------------------------------

# Load stem profile data for volume calculation
stem_profile <- read.csv("data/obeche_stem_profile.csv", comment.char = "#")

cat("=== STEM PROFILE DATA ===\n")
print(stem_profile)

# Tree characteristics
tree_dbh <- 56.5  # cm (measured at 1.3m)
tree_height <- 26.2  # m

# -----------------------------------------------------------------------------
# 2. CALCULATE SECTIONAL VOLUME (SMALIAN'S FORMULA)
# -----------------------------------------------------------------------------

#' Smalian's formula for log volume
#' V = L * (Ab + At) / 2
#' where Ab = basal area, At = top area, L = length

calculate_section_volume <- function(d1_cm, d2_cm, length_m) {
  # Convert diameter to meters for basal area
  r1 <- (d1_cm / 100) / 2
  r2 <- (d2_cm / 100) / 2
  
  a1 <- pi * r1^2
  a2 <- pi * r2^2
  
  volume <- length_m * (a1 + a2) / 2
  return(volume)
}

# Calculate volume for each section
stem_profile$volume_m3 <- NA

for (i in 1:(nrow(stem_profile) - 1)) {
  stem_profile$volume_m3[i] <- calculate_section_volume(
    stem_profile$diameter_ob_cm[i],
    stem_profile$diameter_ob_cm[i + 1],
    stem_profile$section_length_m[i]
  )
}

# Total stem volume
total_volume <- sum(stem_profile$volume_m3, na.rm = TRUE)
cat(sprintf("\nTotal stem volume (over bark): %.4f m³\n", total_volume))

# Form factor calculation
# Form factor = Actual volume / Cylindrical volume
cylinder_volume <- pi * (tree_dbh/200)^2 * tree_height
form_factor <- total_volume / cylinder_volume
cat(sprintf("Form factor: %.4f\n", form_factor))

# -----------------------------------------------------------------------------
# 3. VOLUME MODELS
# -----------------------------------------------------------------------------

# Generate training data for volume models
# Simulating data based on typical obeche relationships
set.seed(123)
n <- 100

# Simulated tree measurements
vol_data <- data.frame(
  tree_id = 1:n,
  dbh_cm = runif(n, 25, 120),
  height_m = NA
)

# Height from H-D relationship
vol_data$height_m <- 1.3 + 33 * (1 - exp(-0.025 * vol_data$dbh_cm)) + rnorm(n, 0, 1.5)
vol_data$height_m <- pmax(vol_data$height_m, 8)

# Volume using allometric relationship with form factor
vol_data$volume_m3 <- form_factor * pi * (vol_data$dbh_cm/200)^2 * vol_data$height_m
vol_data$volume_m3 <- vol_data$volume_m3 * (1 + rnorm(n, 0, 0.1))  # Add noise

cat("\n=== VOLUME DATA SUMMARY ===\n")
cat(sprintf("Number of trees: %d\n", nrow(vol_data)))
cat(sprintf("Volume range: %.3f - %.3f m³\n", min(vol_data$volume_m3), max(vol_data$volume_m3)))

# -----------------------------------------------------------------------------
# 4. FIT VOLUME MODELS
# -----------------------------------------------------------------------------

#' Model 1: Combined variable model
#' V = a * (D²H)^b
#' Classic forestry volume equation

vol_data$d2h <- (vol_data$dbh_cm^2) * vol_data$height_m

model_v1 <- tryCatch({
  nlsLM(volume_m3 ~ a * (d2h)^b,
        data = vol_data,
        start = list(a = 0.00005, b = 1))
}, error = function(e) NULL)

#' Model 2: Schumacher-Hall
#' V = a * D^b * H^c
#' Allows separate exponents for D and H

model_v2 <- tryCatch({
  nlsLM(volume_m3 ~ a * dbh_cm^b * height_m^c,
        data = vol_data,
        start = list(a = 0.0001, b = 2, c = 1))
}, error = function(e) NULL)

#' Model 3: Logarithmic (Spurr)
#' ln(V) = a + b*ln(D) + c*ln(H)
#' Linear in log space

vol_data$ln_vol <- log(vol_data$volume_m3)
vol_data$ln_dbh <- log(vol_data$dbh_cm)
vol_data$ln_h <- log(vol_data$height_m)

model_v3 <- lm(ln_vol ~ ln_dbh + ln_h, data = vol_data)

# -----------------------------------------------------------------------------
# 5. MODEL EVALUATION
# -----------------------------------------------------------------------------

evaluate_volume_model <- function(model, data, model_name, log_model = FALSE) {
  if (is.null(model)) {
    return(data.frame(
      model = model_name,
      RMSE = NA, AIC = NA, BIC = NA, R2 = NA
    ))
  }
  
  if (log_model) {
    pred <- exp(predict(model))
    obs <- data$volume_m3
  } else {
    pred <- predict(model)
    obs <- data$volume_m3
  }
  
  n <- length(obs)
  
  rmse <- sqrt(mean((obs - pred)^2))
  ss_res <- sum((obs - pred)^2)
  ss_tot <- sum((obs - mean(obs))^2)
  r2 <- 1 - ss_res/ss_tot
  
  aic <- AIC(model)
  bic <- BIC(model)
  
  data.frame(
    model = model_name,
    RMSE = round(rmse, 4),
    AIC = round(aic, 2),
    BIC = round(bic, 2),
    R2 = round(r2, 4)
  )
}

vol_results <- rbind(
  evaluate_volume_model(model_v1, vol_data, "Combined Variable V = a*(D²H)^b"),
  evaluate_volume_model(model_v2, vol_data, "Schumacher-Hall V = a*D^b*H^c"),
  evaluate_volume_model(model_v3, vol_data, "Logarithmic ln(V) = a + b*ln(D) + c*ln(H)", log_model = TRUE)
)

cat("\n=== VOLUME MODEL COMPARISON ===\n")
print(vol_results, row.names = FALSE)

# -----------------------------------------------------------------------------
# 6. PARAMETER ESTIMATES
# -----------------------------------------------------------------------------

cat("\n=== PARAMETER ESTIMATES ===\n")

if (!is.null(model_v1)) {
  cat("\nModel 1 (Combined Variable): V = a * (D²H)^b\n")
  print(summary(model_v1)$coefficients)
}

if (!is.null(model_v2)) {
  cat("\nModel 2 (Schumacher-Hall): V = a * D^b * H^c\n")
  print(summary(model_v2)$coefficients)
}

cat("\nModel 3 (Logarithmic): ln(V) = a + b*ln(D) + c*ln(H)\n")
print(summary(model_v3)$coefficients)

# -----------------------------------------------------------------------------
# 7. VISUALIZATION
# -----------------------------------------------------------------------------

# Observed vs Predicted plot
vol_data$pred_v2 <- predict(model_v2)

p_vol <- ggplot(vol_data, aes(x = volume_m3, y = pred_v2)) +
  geom_point(alpha = 0.6, color = "steelblue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(
    title = "Observed vs Predicted Volume",
    subtitle = "Schumacher-Hall Model: V = a × D^b × H^c",
    x = "Observed Volume (m³)",
    y = "Predicted Volume (m³)"
  ) +
  theme_minimal() +
  coord_equal()

ggsave("figures/volume_observed_vs_predicted.png", p_vol,
       width = 8, height = 8, dpi = 300)

# Volume curves
pred_grid <- expand.grid(
  dbh_cm = seq(30, 120, by = 10),
  height_m = c(15, 20, 25, 30)
)
pred_grid$d2h <- pred_grid$dbh_cm^2 * pred_grid$height_m
pred_grid$volume_m3 <- predict(model_v2, pred_grid)

p_curves <- ggplot(pred_grid, aes(x = dbh_cm, y = volume_m3, 
                                   color = factor(height_m))) +
  geom_line(linewidth = 1) +
  scale_color_viridis_d(name = "Height (m)") +
  labs(
    title = "Volume Prediction Curves",
    subtitle = "Triplochiton scleroxylon - Schumacher-Hall Model",
    x = "Diameter at Breast Height (cm)",
    y = "Volume (m³)"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

ggsave("figures/volume_prediction_curves.png", p_curves,
       width = 10, height = 7, dpi = 300)

cat("\n✓ Plots saved to figures/\n")

# -----------------------------------------------------------------------------
# 8. VOLUME TABLE GENERATION
# -----------------------------------------------------------------------------

#' Generate a local volume table
#' @param model fitted model
#' @param dbh_range vector of DBH values
#' @param height_range vector of height values
generate_volume_table <- function(model, dbh_range, height_range) {
  vol_table <- expand.grid(
    DBH_cm = dbh_range,
    Height_m = height_range
  )
  vol_table$d2h <- vol_table$DBH_cm^2 * vol_table$Height_m
  vol_table$Volume_m3 <- round(predict(model, vol_table), 3)
  vol_table$d2h <- NULL
  
  # Reshape to matrix form
  vol_matrix <- reshape(vol_table, 
                        idvar = "DBH_cm", 
                        timevar = "Height_m", 
                        direction = "wide")
  names(vol_matrix) <- gsub("Volume_m3.", "H", names(vol_matrix))
  
  return(vol_matrix)
}

vol_table <- generate_volume_table(
  model_v2,
  dbh_range = seq(30, 100, by = 10),
  height_range = seq(15, 30, by = 5)
)

cat("\n=== LOCAL VOLUME TABLE (m³) ===\n")
cat("Triplochiton scleroxylon - Gambari Forest Reserve\n\n")
print(vol_table, row.names = FALSE)

# Export
write.csv(vol_table, "data/volume_table_obeche.csv", row.names = FALSE)
write.csv(vol_results, "data/volume_model_comparison.csv", row.names = FALSE)

cat("\n=== ANALYSIS COMPLETE ===\n")
