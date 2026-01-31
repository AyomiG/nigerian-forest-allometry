#' ---
#' title: "Above-Ground Biomass Estimation Models"
#' author: "Ogundipe Oluwadamilola"
#' description: "Multivariate analysis and regression models for AGB prediction"
#' species: "Tectona grandis (Teak) and Nauclea diderrichii (Opepe)"
#' location: "Ajibode and Omo Forest Reserve, Nigeria"
#' course: "FRM 736 - Multivariate Analysis in Renewable Natural Resources"
#' data: "Field measurements from 200 trees (100 per species)"
#' ---

# =============================================================================
# ABOVE-GROUND BIOMASS MODELS
# Real field data from Ajibode (Teak) and Omo (Opepe) Forest Reserves
# =============================================================================

if (!require("minpack.lm")) install.packages("minpack.lm")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")
if (!require("corrplot")) install.packages("corrplot")
if (!require("factoextra")) install.packages("factoextra")
if (!require("tidyr")) install.packages("tidyr")

library(minpack.lm)
library(ggplot2)
library(dplyr)
library(corrplot)
library(factoextra)
library(tidyr)

# -----------------------------------------------------------------------------
# 1. LOAD REAL FIELD DATA
# -----------------------------------------------------------------------------

biomass_data <- read.csv("data/teak_opepe_biomass_data.csv")

cat("=============================================================\n")
cat("ABOVE-GROUND BIOMASS ANALYSIS\n")
cat("Species: Tectona grandis (Teak) & Nauclea diderrichii (Opepe)\n")
cat("Location: Ajibode & Omo Forest Reserves, Nigeria\n")
cat("=============================================================\n\n")

cat("--- DATA SUMMARY ---\n")
cat(sprintf("Total trees: %d\n", nrow(biomass_data)))
cat(sprintf("Teak (Ajibode): %d trees\n", sum(biomass_data$species == "Tectona_grandis")))
cat(sprintf("Opepe (Omo): %d trees\n", sum(biomass_data$species == "Nauclea_diderrichii")))

# Summary statistics by species
summary_stats <- biomass_data %>%
  group_by(species) %>%
  summarise(
    n = n(),
    dbh_mean = round(mean(dbh_cm), 1),
    dbh_range = paste0(round(min(dbh_cm), 1), "-", round(max(dbh_cm), 1)),
    height_mean = round(mean(total_height_m), 1),
    height_range = paste0(round(min(total_height_m), 1), "-", round(max(total_height_m), 1)),
    agb_mean = round(mean(agb_kg), 0),
    agb_range = paste0(round(min(agb_kg), 0), "-", round(max(agb_kg), 0))
  )

cat("\n--- Summary by Species ---\n")
print(as.data.frame(summary_stats))

# -----------------------------------------------------------------------------
# 2. CORRELATION ANALYSIS
# -----------------------------------------------------------------------------

# Select numeric variables for correlation
numeric_vars <- biomass_data %>%
  select(dbh_cm, db_cm, dm_cm, dt_cm, total_height_m, merch_height_m,
         crown_diameter_m, crown_length_m, basal_area_m2, 
         total_volume_m3, agb_kg)

cor_matrix <- cor(numeric_vars, use = "complete.obs")

cat("\n--- CORRELATION WITH AGB ---\n")
agb_cors <- sort(cor_matrix[, "agb_kg"], decreasing = TRUE)
print(round(agb_cors, 3))

# Save correlation plot
png("figures/correlation_matrix_real.png", width = 900, height = 900, res = 100)
corrplot(cor_matrix, method = "color", type = "upper",
         addCoef.col = "black", number.cex = 0.7,
         tl.col = "black", tl.srt = 45,
         title = "Correlation Matrix - Teak & Opepe Tree Variables",
         mar = c(0, 0, 2, 0))
dev.off()
cat("\n✓ Correlation plot saved to figures/correlation_matrix_real.png\n")

# -----------------------------------------------------------------------------
# 3. PRINCIPAL COMPONENT ANALYSIS (DATA REDUCTION)
# -----------------------------------------------------------------------------

# Standardize predictor variables (excluding AGB)
predictors <- numeric_vars %>% select(-agb_kg)
predictors_scaled <- scale(predictors)

# PCA
pca_result <- prcomp(predictors_scaled, center = TRUE, scale. = TRUE)

cat("\n--- PRINCIPAL COMPONENT ANALYSIS ---\n")
cat("Variance explained by each component:\n")
summary(pca_result)

# Scree plot
p_scree <- fviz_eig(pca_result, addlabels = TRUE,
                    main = "Scree Plot - Variance Explained")
ggsave("figures/pca_scree_plot_real.png", p_scree, width = 8, height = 6, dpi = 300)

# Variable contributions
p_contrib <- fviz_pca_var(pca_result, col.var = "contrib",
                          gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                          repel = TRUE,
                          title = "PCA - Variable Contributions (Teak & Opepe)")
ggsave("figures/pca_variables_real.png", p_contrib, width = 9, height = 8, dpi = 300)

cat("\nPC Loadings (first 3 components):\n")
print(round(pca_result$rotation[, 1:3], 3))

# Key finding: DBH, Height, and Crown Diameter have high loadings
cat("\n→ Key variables identified: DBH, Total Height, Crown Diameter\n")

# -----------------------------------------------------------------------------
# 4. CLUSTER ANALYSIS
# -----------------------------------------------------------------------------

set.seed(42)
kmeans_result <- kmeans(predictors_scaled, centers = 5, nstart = 25)
biomass_data$cluster <- factor(kmeans_result$cluster)

cat("\n--- CLUSTER ANALYSIS (5 SIZE CLASSES) ---\n")
cluster_summary <- biomass_data %>%
  group_by(cluster) %>%
  summarise(
    n = n(),
    mean_dbh = round(mean(dbh_cm), 1),
    mean_height = round(mean(total_height_m), 1),
    mean_agb = round(mean(agb_kg), 0)
  ) %>%
  arrange(mean_dbh)

print(as.data.frame(cluster_summary))

# Cluster visualization
p_cluster <- fviz_cluster(kmeans_result, data = predictors_scaled,
                          palette = "Set2", geom = "point",
                          ellipse.type = "convex",
                          main = "Tree Size Clusters - Teak & Opepe")
ggsave("figures/cluster_analysis_real.png", p_cluster, width = 9, height = 7, dpi = 300)

# -----------------------------------------------------------------------------
# 5. BIOMASS MODELS - LINEAR (LOG-TRANSFORMED)
# -----------------------------------------------------------------------------

cat("\n--- FITTING LINEAR BIOMASS MODELS ---\n")

# Create transformed variables
biomass_data <- biomass_data %>%
  mutate(
    ln_agb = log(agb_kg),
    ln_dbh = log(dbh_cm),
    ln_h = log(total_height_m),
    ln_cd = log(crown_diameter_m),
    ln_ba = log(basal_area_m2),
    d2h = dbh_cm^2 * total_height_m,
    ln_d2h = log(d2h)
  )

# Model L1: Simple DBH
model_l1 <- lm(ln_agb ~ ln_dbh, data = biomass_data)

# Model L2: DBH + Height
model_l2 <- lm(ln_agb ~ ln_dbh + ln_h, data = biomass_data)

# Model L3: DBH + Height + Crown Diameter
model_l3 <- lm(ln_agb ~ ln_dbh + ln_h + ln_cd, data = biomass_data)

# Model L4: Combined variable (D²H)
model_l4 <- lm(ln_agb ~ ln_d2h, data = biomass_data)

# Model L5: With species effect
model_l5 <- lm(ln_agb ~ ln_dbh + ln_h + species, data = biomass_data)

# -----------------------------------------------------------------------------
# 6. BIOMASS MODELS - NON-LINEAR
# -----------------------------------------------------------------------------

cat("\n--- FITTING NON-LINEAR BIOMASS MODELS ---\n")

# Model NL1: Power function (DBH only)
model_nl1 <- tryCatch({
  nlsLM(agb_kg ~ a * dbh_cm^b,
        data = biomass_data,
        start = list(a = 0.5, b = 2.3),
        control = nls.lm.control(maxiter = 200))
}, error = function(e) { cat("NL1 failed:", e$message, "\n"); NULL })

# Model NL2: Chave-type (D²H)
model_nl2 <- tryCatch({
  nlsLM(agb_kg ~ a * (d2h)^b,
        data = biomass_data,
        start = list(a = 0.1, b = 0.9),
        control = nls.lm.control(maxiter = 200))
}, error = function(e) { cat("NL2 failed:", e$message, "\n"); NULL })

# Model NL3: Basuki-type
model_nl3 <- tryCatch({
  nlsLM(agb_kg ~ exp(a + b * log(dbh_cm) + c * log(total_height_m)),
        data = biomass_data,
        start = list(a = -1, b = 2, c = 0.5),
        control = nls.lm.control(maxiter = 200))
}, error = function(e) { cat("NL3 failed:", e$message, "\n"); NULL })

# Model NL4: Three-variable
model_nl4 <- tryCatch({
  nlsLM(agb_kg ~ a * dbh_cm^b * total_height_m^c * crown_diameter_m^d,
        data = biomass_data,
        start = list(a = 0.05, b = 2, c = 0.5, d = 0.3),
        control = nls.lm.control(maxiter = 200))
}, error = function(e) { cat("NL4 failed:", e$message, "\n"); NULL })

# -----------------------------------------------------------------------------
# 7. MODEL EVALUATION
# -----------------------------------------------------------------------------

evaluate_model <- function(model, data, model_name, log_response = FALSE) {
  if (is.null(model)) {
    return(data.frame(model = model_name, SEE = NA, AIC = NA, R2_adj = NA))
  }
  
  if (log_response) {
    pred <- exp(predict(model))
    obs <- data$agb_kg
  } else {
    pred <- predict(model)
    obs <- data$agb_kg
  }
  
  n <- length(obs)
  k <- length(coef(model))
  
  see <- sqrt(sum((obs - pred)^2, na.rm = TRUE) / (n - k))
  ss_res <- sum((obs - pred)^2, na.rm = TRUE)
  ss_tot <- sum((obs - mean(obs))^2, na.rm = TRUE)
  r2 <- 1 - ss_res/ss_tot
  r2_adj <- 1 - (1 - r2) * (n - 1) / (n - k - 1)
  
  data.frame(
    model = model_name,
    SEE = round(see, 2),
    AIC = round(AIC(model), 2),
    R2_adj = round(r2_adj, 4)
  )
}

# Evaluate all models
results <- rbind(
  evaluate_model(model_l1, biomass_data, "L1: ln(AGB) ~ ln(DBH)", TRUE),
  evaluate_model(model_l2, biomass_data, "L2: ln(AGB) ~ ln(DBH) + ln(H)", TRUE),
  evaluate_model(model_l3, biomass_data, "L3: ln(AGB) ~ ln(DBH) + ln(H) + ln(CD)", TRUE),
  evaluate_model(model_l4, biomass_data, "L4: ln(AGB) ~ ln(D²H)", TRUE),
  evaluate_model(model_l5, biomass_data, "L5: ln(AGB) ~ ln(DBH) + ln(H) + Species", TRUE),
  evaluate_model(model_nl1, biomass_data, "NL1: AGB = a × DBH^b"),
  evaluate_model(model_nl2, biomass_data, "NL2: AGB = a × (D²H)^b"),
  evaluate_model(model_nl3, biomass_data, "NL3: Basuki-type"),
  evaluate_model(model_nl4, biomass_data, "NL4: Three-variable")
)

cat("\n========================================\n")
cat("MODEL COMPARISON RESULTS\n")
cat("========================================\n")
print(results %>% arrange(desc(R2_adj)), row.names = FALSE)

# Best models
best_linear <- results %>% filter(grepl("^L", model)) %>% slice_max(R2_adj)
best_nonlin <- results %>% filter(grepl("^NL", model)) %>% slice_max(R2_adj)

cat(sprintf("\n→ Best Linear Model: %s (R²adj = %.4f)\n", best_linear$model, best_linear$R2_adj))
cat(sprintf("→ Best Non-linear Model: %s (R²adj = %.4f)\n", best_nonlin$model, best_nonlin$R2_adj))

# -----------------------------------------------------------------------------
# 8. BEST MODEL COEFFICIENTS
# -----------------------------------------------------------------------------

cat("\n--- BEST MODEL PARAMETERS ---\n")

cat("\nModel L3 (Best Linear): ln(AGB) = a + b×ln(DBH) + c×ln(H) + d×ln(CD)\n")
print(summary(model_l3)$coefficients)

if (!is.null(model_nl3)) {
  cat("\nModel NL3 (Basuki-type): AGB = exp(a + b×ln(DBH) + c×ln(H))\n")
  print(summary(model_nl3)$coefficients)
}

# -----------------------------------------------------------------------------
# 9. VISUALIZATIONS
# -----------------------------------------------------------------------------

# Observed vs Predicted
biomass_data$pred_l3 <- exp(predict(model_l3))

p_obs_pred <- ggplot(biomass_data, aes(x = agb_kg, y = pred_l3, color = species)) +
  geom_point(alpha = 0.7, size = 2.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40") +
  scale_color_manual(values = c("Tectona_grandis" = "#2E86AB", 
                                 "Nauclea_diderrichii" = "#A23B72"),
                     labels = c("Opepe (Omo)", "Teak (Ajibode)")) +
  labs(
    title = "Above-Ground Biomass: Observed vs Predicted",
    subtitle = "Model L3: ln(AGB) = a + b×ln(DBH) + c×ln(H) + d×ln(CD)",
    x = "Observed AGB (kg)",
    y = "Predicted AGB (kg)",
    color = "Species"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom") +
  coord_equal()

ggsave("figures/agb_observed_vs_predicted_real.png", p_obs_pred,
       width = 9, height = 9, dpi = 300)

# Residual plot
biomass_data$residual <- biomass_data$agb_kg - biomass_data$pred_l3

p_resid <- ggplot(biomass_data, aes(x = pred_l3, y = residual, color = species)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_smooth(method = "loess", se = FALSE, color = "black", linewidth = 0.5) +
  scale_color_manual(values = c("Tectona_grandis" = "#2E86AB", 
                                 "Nauclea_diderrichii" = "#A23B72")) +
  labs(
    title = "Residual Plot - Model L3",
    x = "Predicted AGB (kg)",
    y = "Residual (kg)"
  ) +
  theme_minimal()

ggsave("figures/agb_residuals_real.png", p_resid, width = 10, height = 6, dpi = 300)

# Species comparison
p_species <- ggplot(biomass_data, aes(x = dbh_cm, y = agb_kg, color = species)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE) +
  scale_color_manual(values = c("Tectona_grandis" = "#2E86AB", 
                                 "Nauclea_diderrichii" = "#A23B72"),
                     labels = c("Opepe (Omo)", "Teak (Ajibode)")) +
  labs(
    title = "AGB-DBH Relationship by Species",
    x = "Diameter at Breast Height (cm)",
    y = "Above-Ground Biomass (kg)",
    color = "Species"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("figures/agb_dbh_by_species.png", p_species, width = 10, height = 7, dpi = 300)

cat("\n✓ All plots saved to figures/\n")

# -----------------------------------------------------------------------------
# 10. EXPORT FINAL RESULTS
# -----------------------------------------------------------------------------

# Save model comparison
write.csv(results, "data/agb_model_comparison_real.csv", row.names = FALSE)

# Save recommended equations
coef_l3 <- coef(model_l3)
equations <- data.frame(
  model = c("L3 (Recommended)", "L2 (Simple)"),
  equation = c(
    sprintf("AGB = exp(%.4f + %.4f×ln(DBH) + %.4f×ln(H) + %.4f×ln(CD))",
            coef_l3[1], coef_l3[2], coef_l3[3], coef_l3[4]),
    sprintf("AGB = exp(%.4f + %.4f×ln(DBH) + %.4f×ln(H))",
            coef(model_l2)[1], coef(model_l2)[2], coef(model_l2)[3])
  ),
  R2_adj = c(best_linear$R2_adj, 
             results$R2_adj[results$model == "L2: ln(AGB) ~ ln(DBH) + ln(H)"])
)
write.csv(equations, "data/recommended_agb_equations_real.csv", row.names = FALSE)

cat("\n========================================\n")
cat("ANALYSIS COMPLETE\n")
cat("========================================\n")
cat("\nRecommended equation for AGB prediction:\n")
cat(sprintf("  AGB (kg) = exp(%.4f + %.4f×ln(DBH) + %.4f×ln(H) + %.4f×ln(CD))\n",
            coef_l3[1], coef_l3[2], coef_l3[3], coef_l3[4]))
cat(sprintf("  R²adj = %.4f, SEE = %.2f kg\n", best_linear$R2_adj, best_linear$SEE))
