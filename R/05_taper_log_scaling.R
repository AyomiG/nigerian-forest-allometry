#' ---
#' title: "Taper Equations and Log Scaling for Sawmill Applications"
#' author: "Ogundipe Oluwadamilola"
#' description: "Stem taper models, log volume scaling, and bucking optimization"
#' application: "Sawmill industry, timber trade, log inventory"
#' ---

# =============================================================================
# TAPER EQUATIONS AND LOG SCALING
# Practical tools for sawmill biometrics and timber industry
# =============================================================================

if (!require("minpack.lm")) install.packages("minpack.lm")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")

library(minpack.lm)
library(ggplot2)
library(dplyr)

# -----------------------------------------------------------------------------
# 1. LOAD STEM PROFILE DATA
# -----------------------------------------------------------------------------

stem_data <- read.csv("data/obeche_stem_profile.csv", comment.char = "#")

# Tree specifications
DBH <- 56.5      # cm (diameter at breast height, 1.3m)
TOTAL_HT <- 26.2 # m (total tree height)

cat("=============================================================\n")
cat("TAPER EQUATIONS AND LOG SCALING\n")
cat("Species: Triplochiton scleroxylon (Obeche)\n")
cat("Application: Sawmill Industry Biometrics\n")
cat("=============================================================\n\n")

cat("--- TREE SPECIFICATIONS ---\n")
cat(sprintf("DBH: %.1f cm\n", DBH))
cat(sprintf("Total Height: %.1f m\n", TOTAL_HT))
cat(sprintf("Sections measured: %d\n", nrow(stem_data)))

# Preview data
cat("\n--- STEM PROFILE DATA ---\n")
print(head(stem_data, 10))

# -----------------------------------------------------------------------------
# 2. TAPER EQUATIONS
# -----------------------------------------------------------------------------

# Prepare data for taper modeling
taper_data <- stem_data %>%
  mutate(
    d = diameter_ob_cm,           # diameter at height h
    h = height_m,                 # height from ground
    H = TOTAL_HT,                 # total tree height
    D = DBH,                      # DBH
    rel_h = h / H,                # relative height (0-1)
    rel_d = d / D                 # relative diameter
  )

cat("\n--- FITTING TAPER EQUATIONS ---\n")

#' =========================================================================
#' TAPER MODEL 1: Kozak (1988) Variable-Exponent Taper
#' =========================================================================
#' d/D = a0 * X^(a1 * z^4 + a2 * (exp(-D/H)) + a3 * X^0.1 + a4 * (1/exp(D/H)) + a5 * h^(1/3) + a6 * X)
#' where X = (1 - z^0.5) / (1 - p^0.5), z = h/H, p = 1.3/H
#' 
#' Simplified version for our data:

# Model 1: Simple polynomial taper
# d = a + b1*(h/H) + b2*(h/H)^2
model_poly <- tryCatch({
  nlsLM(d ~ D * (a + b1 * rel_h + b2 * rel_h^2),
        data = taper_data,
        start = list(a = 1.2, b1 = -0.5, b2 = -0.3))
}, error = function(e) NULL)

#' =========================================================================
#' TAPER MODEL 2: Max & Burkhart (1976) Segmented Polynomial
#' =========================================================================
#' Classic model used in many timber industries

# Model 2: Ormerod (1973) simple ratio
# d/D = (H-h)/(H-1.3))^b
model_ormerod <- tryCatch({
  nlsLM(d ~ D * ((H - h) / (H - 1.3))^b,
        data = taper_data %>% filter(h >= 1.3),
        start = list(b = 0.8))
}, error = function(e) NULL)

#' =========================================================================
#' TAPER MODEL 3: Kozak (2004) - Simplified Variable Exponent
#' =========================================================================

model_kozak <- tryCatch({
  nlsLM(d ~ a0 * D^a1 * H^a2 * ((1 - rel_h^0.5)/(1 - 0.22^0.5))^(b1 + b2/exp(D/H)),
        data = taper_data,
        start = list(a0 = 1.0, a1 = 1.0, a2 = -0.1, b1 = 0.5, b2 = 0.3),
        control = nls.lm.control(maxiter = 200))
}, error = function(e) NULL)

# -----------------------------------------------------------------------------
# 3. EVALUATE TAPER MODELS
# -----------------------------------------------------------------------------

evaluate_taper <- function(model, data, model_name) {
  if (is.null(model)) {
    return(data.frame(model = model_name, RMSE = NA, R2 = NA, Bias = NA))
  }
  
  pred <- predict(model)
  obs <- data$d[1:length(pred)]
  
  rmse <- sqrt(mean((obs - pred)^2))
  bias <- mean(pred - obs)
  ss_res <- sum((obs - pred)^2)
  ss_tot <- sum((obs - mean(obs))^2)
  r2 <- 1 - ss_res/ss_tot
  
  data.frame(
    model = model_name,
    RMSE = round(rmse, 3),
    R2 = round(r2, 4),
    Bias = round(bias, 3)
  )
}

taper_results <- rbind(
  evaluate_taper(model_poly, taper_data, "Polynomial (d = D*(a + b1*z + b2*z²))"),
  evaluate_taper(model_ormerod, taper_data %>% filter(h >= 1.3), "Ormerod ((H-h)/(H-1.3))^b"),
  evaluate_taper(model_kozak, taper_data, "Kozak Variable-Exponent")
)

cat("\n--- TAPER MODEL COMPARISON ---\n")
print(taper_results, row.names = FALSE)

# Best model
best_taper <- taper_results %>% filter(!is.na(R2)) %>% slice_max(R2)
cat(sprintf("\n→ Best taper model: %s (R² = %.4f)\n", best_taper$model, best_taper$R2))

# Print coefficients
if (!is.null(model_poly)) {
  cat("\nPolynomial Taper Coefficients:\n")
  print(coef(model_poly))
}

# -----------------------------------------------------------------------------
# 4. PREDICT DIAMETER AT ANY HEIGHT
# -----------------------------------------------------------------------------

#' Function to predict diameter at any height using fitted taper model
#' @param h Height from ground (m)
#' @param D DBH (cm)
#' @param H Total height (m)
#' @param model Fitted taper model
predict_diameter <- function(h, D, H, model = model_poly) {
  if (is.null(model)) return(NA)
  
  newdata <- data.frame(
    h = h, H = H, D = D,
    rel_h = h / H, d = NA
  )
  predict(model, newdata)
}

# Example: predict diameter at different heights
cat("\n--- DIAMETER PREDICTIONS ALONG STEM ---\n")
heights <- seq(0.5, 20, by = 2)
predictions <- data.frame(
  height_m = heights,
  predicted_diameter_cm = sapply(heights, function(h) {
    round(predict_diameter(h, DBH, TOTAL_HT), 1)
  })
)
print(predictions)

# -----------------------------------------------------------------------------
# 5. LOG SCALING RULES
# -----------------------------------------------------------------------------

cat("\n=============================================================\n")
cat("LOG VOLUME SCALING RULES\n")
cat("=============================================================\n")

#' =========================================================================
#' LOG SCALING FORMULAS (Volume in cubic meters)
#' =========================================================================

#' Smalian's Formula (most common in sawmills)
#' V = L * (Ab + At) / 2
#' where Ab = basal area, At = top area, L = length
smalian_volume <- function(d1_cm, d2_cm, length_m) {
  r1 <- (d1_cm / 100) / 2
  r2 <- (d2_cm / 100) / 2
  a1 <- pi * r1^2
  a2 <- pi * r2^2
  length_m * (a1 + a2) / 2
}

#' Huber's Formula
#' V = L * Am
#' where Am = cross-sectional area at mid-point
huber_volume <- function(d_mid_cm, length_m) {
  r <- (d_mid_cm / 100) / 2
  am <- pi * r^2
  length_m * am
}

#' Newton's Formula (most accurate)
#' V = L * (Ab + 4*Am + At) / 6
newton_volume <- function(d1_cm, d_mid_cm, d2_cm, length_m) {
  r1 <- (d1_cm / 100) / 2
  rm <- (d_mid_cm / 100) / 2
  r2 <- (d2_cm / 100) / 2
  a1 <- pi * r1^2
  am <- pi * rm^2
  a2 <- pi * r2^2
  length_m * (a1 + 4*am + a2) / 6
}

#' =========================================================================
#' BOARD FOOT RULES (US/International lumber trade)
#' =========================================================================

#' Doyle Rule (conservative, favors buyer)
#' BF = (D - 4)^2 * L / 16
#' where D = small-end diameter (inches), L = length (feet)
doyle_bf <- function(d_small_cm, length_m) {
  d_inch <- d_small_cm / 2.54
  l_feet <- length_m * 3.281
  ((d_inch - 4)^2 * l_feet) / 16
}

#' Scribner Rule (moderate)
#' Uses log tables - simplified approximation
scribner_bf <- function(d_small_cm, length_m) {
  d_inch <- d_small_cm / 2.54
  l_feet <- length_m * 3.281
  (0.79 * d_inch^2 - 2 * d_inch - 4) * (l_feet / 16)
}

#' International 1/4-inch Rule (most accurate for board feet)
#' BF = 0.905 * (0.22*D^2 - 0.71*D) * (L/4)
international_bf <- function(d_small_cm, length_m) {
  d_inch <- d_small_cm / 2.54
  l_feet <- length_m * 3.281
  0.905 * (0.22 * d_inch^2 - 0.71 * d_inch) * (l_feet / 4)
}

# -----------------------------------------------------------------------------
# 6. CALCULATE LOG VOLUMES FROM STEM PROFILE
# -----------------------------------------------------------------------------

cat("\n--- LOG VOLUME CALCULATIONS ---\n")

# Calculate section volumes using different methods
stem_data$smalian_vol <- NA
stem_data$huber_vol <- NA
stem_data$newton_vol <- NA

for (i in 1:(nrow(stem_data) - 1)) {
  d1 <- stem_data$diameter_ob_cm[i]
  d2 <- stem_data$diameter_ob_cm[i + 1]
  d_mid <- (d1 + d2) / 2
  len <- stem_data$height_m[i + 1] - stem_data$height_m[i]
  
  stem_data$smalian_vol[i] <- smalian_volume(d1, d2, len)
  stem_data$huber_vol[i] <- huber_volume(d_mid, len)
  stem_data$newton_vol[i] <- newton_volume(d1, d_mid, d2, len)
}

# Total volumes by method
cat("\nTotal Stem Volume by Scaling Method:\n")
cat(sprintf("  Smalian's Formula: %.4f m³\n", sum(stem_data$smalian_vol, na.rm = TRUE)))
cat(sprintf("  Huber's Formula:   %.4f m³\n", sum(stem_data$huber_vol, na.rm = TRUE)))
cat(sprintf("  Newton's Formula:  %.4f m³\n", sum(stem_data$newton_vol, na.rm = TRUE)))

# -----------------------------------------------------------------------------
# 7. LOG BUCKING OPTIMIZATION
# -----------------------------------------------------------------------------

cat("\n=============================================================\n")
cat("LOG BUCKING OPTIMIZATION\n")
cat("=============================================================\n")

#' Simulate different bucking scenarios to maximize value
#' Sawmills want logs of specific lengths (e.g., 2.5m, 3m, 4m, 5m)

# Define log length options and minimum diameters
LOG_LENGTHS <- c(2.5, 3.0, 4.0, 5.0)  # meters
MIN_SMALL_END_DIAM <- 20  # cm (minimum merchantable diameter)
STUMP_HEIGHT <- 0.3  # m

#' Function to buck a stem into optimal logs
#' @param stem_profile Data frame with height_m and diameter_ob_cm
#' @param log_lengths Vector of acceptable log lengths
#' @param min_diam Minimum small-end diameter
buck_stem <- function(stem_profile, log_lengths = LOG_LENGTHS, min_diam = MIN_SMALL_END_DIAM) {
  
  # Start from stump
  current_h <- STUMP_HEIGHT
  logs <- list()
  log_num <- 1
  
  # Find merchantable height (where diameter = min_diam)
  merch_idx <- which(stem_profile$diameter_ob_cm < min_diam)[1]
  if (is.na(merch_idx)) merch_idx <- nrow(stem_profile)
  merch_height <- stem_profile$height_m[merch_idx]
  
  while (current_h < merch_height - min(log_lengths)) {
    best_value <- 0
    best_length <- NULL
    
    for (len in sort(log_lengths, decreasing = TRUE)) {
      end_h <- current_h + len
      if (end_h > merch_height) next
      
      # Get diameters at start and end of potential log
      d_large <- approx(stem_profile$height_m, stem_profile$diameter_ob_cm, current_h)$y
      d_small <- approx(stem_profile$height_m, stem_profile$diameter_ob_cm, end_h)$y
      
      if (is.na(d_small) || d_small < min_diam) next
      
      # Calculate value (volume * price factor based on diameter)
      vol <- smalian_volume(d_large, d_small, len)
      # Larger diameter logs command premium price
      price_factor <- ifelse(d_small >= 40, 1.5, ifelse(d_small >= 30, 1.2, 1.0))
      value <- vol * price_factor
      
      if (value > best_value) {
        best_value <- value
        best_length <- len
        best_d_large <- d_large
        best_d_small <- d_small
        best_vol <- vol
      }
    }
    
    if (is.null(best_length)) break
    
    logs[[log_num]] <- data.frame(
      log_number = log_num,
      start_height_m = current_h,
      end_height_m = current_h + best_length,
      length_m = best_length,
      large_end_diam_cm = round(best_d_large, 1),
      small_end_diam_cm = round(best_d_small, 1),
      volume_m3 = round(best_vol, 4)
    )
    
    current_h <- current_h + best_length
    log_num <- log_num + 1
  }
  
  do.call(rbind, logs)
}

# Perform bucking optimization
bucking_result <- buck_stem(stem_data)

cat("\n--- OPTIMAL BUCKING PATTERN ---\n")
print(bucking_result, row.names = FALSE)

cat(sprintf("\nTotal logs: %d\n", nrow(bucking_result)))
cat(sprintf("Total merchantable volume: %.4f m³\n", sum(bucking_result$volume_m3)))
cat(sprintf("Merchantable length: %.1f m\n", sum(bucking_result$length_m)))

# -----------------------------------------------------------------------------
# 8. LOG TABLE GENERATION
# -----------------------------------------------------------------------------

cat("\n=============================================================\n")
cat("LOG VOLUME TABLE\n")
cat("=============================================================\n")

#' Generate a log volume table for field use
#' Shows volume for different diameter and length combinations

generate_log_table <- function(diameters, lengths) {
  table_data <- expand.grid(
    small_end_diam_cm = diameters,
    length_m = lengths
  )
  
  # Assume 10% taper per meter for large end estimation
  table_data$large_end_diam_cm <- table_data$small_end_diam_cm + 
                                   (table_data$length_m * 1.0)  # 1 cm per meter taper
  
  table_data$volume_m3 <- mapply(smalian_volume, 
                                  table_data$large_end_diam_cm,
                                  table_data$small_end_diam_cm,
                                  table_data$length_m)
  
  table_data$volume_m3 <- round(table_data$volume_m3, 4)
  
  # Reshape to matrix
  vol_matrix <- reshape(table_data[, c("small_end_diam_cm", "length_m", "volume_m3")],
                        idvar = "small_end_diam_cm",
                        timevar = "length_m",
                        direction = "wide")
  names(vol_matrix) <- gsub("volume_m3.", "L", names(vol_matrix))
  
  vol_matrix
}

# Create log table
log_table <- generate_log_table(
  diameters = seq(20, 60, by = 5),
  lengths = c(2.5, 3.0, 4.0, 5.0)
)

cat("\nLog Volume Table (m³) - Smalian's Formula\n")
cat("Rows = Small-end diameter (cm), Columns = Log length (m)\n\n")
print(log_table, row.names = FALSE)

# -----------------------------------------------------------------------------
# 9. BOARD FEET CONVERSION TABLE
# -----------------------------------------------------------------------------

cat("\n--- BOARD FEET SCALING ---\n")
cat("(For international timber trade)\n\n")

# Example: convert our logs to board feet
if (nrow(bucking_result) > 0) {
  bucking_result$doyle_bf <- mapply(doyle_bf, 
                                     bucking_result$small_end_diam_cm,
                                     bucking_result$length_m)
  bucking_result$scribner_bf <- mapply(scribner_bf,
                                        bucking_result$small_end_diam_cm,
                                        bucking_result$length_m)
  bucking_result$intl_bf <- mapply(international_bf,
                                    bucking_result$small_end_diam_cm,
                                    bucking_result$length_m)
  
  cat("Board Feet by Scaling Rule:\n")
  bf_summary <- bucking_result %>%
    summarise(
      Doyle = round(sum(doyle_bf)),
      Scribner = round(sum(scribner_bf)),
      International = round(sum(intl_bf))
    )
  print(bf_summary)
}

# -----------------------------------------------------------------------------
# 10. VISUALIZATIONS
# -----------------------------------------------------------------------------

# Stem taper profile
p_taper <- ggplot(taper_data, aes(x = d, y = h)) +
  geom_point(size = 3, color = "darkgreen") +
  geom_path(color = "darkgreen", linewidth = 1) +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "red", alpha = 0.7) +
  annotate("text", x = max(taper_data$d) * 0.9, y = 1.5, 
           label = "DBH (1.3m)", color = "red", size = 3) +
  coord_flip() +
  labs(
    title = "Stem Taper Profile - Triplochiton scleroxylon",
    subtitle = sprintf("DBH = %.1f cm, Height = %.1f m", DBH, TOTAL_HT),
    x = "Diameter (cm)",
    y = "Height (m)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

ggsave("figures/stem_taper_profile.png", p_taper, width = 8, height = 10, dpi = 300)

# Log bucking diagram
if (nrow(bucking_result) > 0) {
  p_bucking <- ggplot(bucking_result) +
    geom_rect(aes(xmin = start_height_m, xmax = end_height_m,
                  ymin = -small_end_diam_cm/2, ymax = small_end_diam_cm/2,
                  fill = factor(log_number)),
              alpha = 0.7, color = "black") +
    geom_text(aes(x = (start_height_m + end_height_m)/2, y = 0,
                  label = paste0("Log ", log_number, "\n", 
                                 length_m, "m\n",
                                 round(volume_m3, 3), " m³")),
              size = 3) +
    scale_fill_brewer(palette = "Set3", name = "Log #") +
    labs(
      title = "Optimal Log Bucking Pattern",
      subtitle = sprintf("Total: %d logs, %.3f m³ merchantable volume",
                         nrow(bucking_result), sum(bucking_result$volume_m3)),
      x = "Distance from ground (m)",
      y = "Diameter (cm)"
    ) +
    theme_minimal() +
    theme(legend.position = "none")
  
  ggsave("figures/log_bucking_pattern.png", p_bucking, width = 12, height = 6, dpi = 300)
}

cat("\n✓ Plots saved to figures/\n")

# -----------------------------------------------------------------------------
# 11. EXPORT RESULTS
# -----------------------------------------------------------------------------

# Save log table
write.csv(log_table, "data/log_volume_table.csv", row.names = FALSE)

# Save bucking results
write.csv(bucking_result, "data/bucking_optimization.csv", row.names = FALSE)

# Save taper model coefficients
if (!is.null(model_poly)) {
  taper_coefs <- data.frame(
    model = "Polynomial",
    parameter = names(coef(model_poly)),
    value = as.numeric(coef(model_poly))
  )
  write.csv(taper_coefs, "data/taper_coefficients.csv", row.names = FALSE)
}

cat("\n=============================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("=============================================================\n")
cat("\nOutputs generated:\n")
cat("  - figures/stem_taper_profile.png\n")
cat("  - figures/log_bucking_pattern.png\n")
cat("  - data/log_volume_table.csv\n")
cat("  - data/bucking_optimization.csv\n")
cat("  - data/taper_coefficients.csv\n")

cat("\n--- KEY TAKEAWAYS FOR SAWMILL APPLICATIONS ---\n")
cat("1. Taper equation allows diameter prediction at any height\n")
cat("2. Use Smalian's formula for quick field scaling\n")
cat("3. Newton's formula is most accurate for payment/inventory\n")
cat("4. Optimal bucking maximizes value, not just volume\n")
cat("5. Board foot rules vary - know your market!\n")
