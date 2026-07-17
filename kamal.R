# ======================================================
# FULL SCRIPT: Water Quality Analysis (15 wells)
# Extended with normality, Spearman, regression,
# scatter plots, group comparison, and distance check
# ======================================================

# Load required packages
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(writexl)

# add Files dir.
# Read
kamal <- read_excel("kamal.xlsx")
View(kamal)

out_dir <- "output"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


# Clean the data

# Remove first empty column (if present)
kamal <- kamal[, -1]

# Assign column names
# NOTE: Turbidity values in data are ~0.005 NTU (very low; all below WHO limit of 5)
names(kamal) <- c("Sample", "pH", "TDS_mgL", "EC_uScm",
                  "Turbidity_NTU", "TCC_10x3_cfu_ml",
                  "Nitrate_mgL", "Ecoli_cfu_100ml", "Distance_m")

# Remove rows with missing Sample ID
kamal <- kamal %>% filter(!is.na(Sample))


# Descriptive statistics

quality_vars <- kamal %>% select(where(is.numeric))

descriptives_raw <- quality_vars %>%
  summarise(across(everything(),
                   list(mean = ~mean(., na.rm = TRUE),
                        min  = ~min(.,  na.rm = TRUE),
                        max  = ~max(.,  na.rm = TRUE),
                        sd   = ~sd(.,   na.rm = TRUE))))

descriptives <- descriptives_raw %>%
  pivot_longer(everything(), names_to = "var", values_to = "value") %>%
  mutate(
    Parameter = sub("_(mean|min|max|sd)$", "", var),
    Statistic = sub(".*_(mean|min|max|sd)$", "\\1", var)
  ) %>%
  select(-var) %>%
  pivot_wider(names_from = Statistic, values_from = value) %>%
  select(Parameter, mean, min, max, sd) %>%
  arrange(Parameter)

print(descriptives)
write_xlsx(descriptives, file.path(out_dir, "descriptive_statistics.xlsx"))
cat("Saved: descriptive_statistics.xlsx\n")


# Compliance with WHO / NSDWQ standards
kamal <- kamal %>%
  rowwise() %>%
  mutate(
    pH_complies        = between(pH, 6.5, 8.5),
    TDS_complies       = TDS_mgL <= 500,
    Nitrate_complies   = Nitrate_mgL <= 50,
    Ecoli_complies     = Ecoli_cfu_100ml == 0,
    Turbidity_complies = Turbidity_NTU <= 5,
    TCC_complies       = TCC_10x3_cfu_ml == 0
  ) %>%
  ungroup()

comparison <- kamal %>% select(Sample, ends_with("complies"))
write_xlsx(comparison, file.path(out_dir, "compliance_WHO_NSDWQ.xlsx"))
cat("Saved: compliance_WHO_NSDWQ.xlsx\n")


# Relationship: distance vs. contamination
relationship_table <- kamal %>%
  select(Sample, Distance_m, Nitrate_mgL, Ecoli_cfu_100ml, TCC_10x3_cfu_ml)
write_xlsx(relationship_table, file.path(out_dir, "distance_contamination.xlsx"))
cat("Saved: distance_contamination.xlsx\n")


# Graph: ALL parameters in ONE combined plot (% of WHO limit)
who_upper <- c(
  pH              = 8.5,
  TDS_mgL         = 500,
  Nitrate_mgL     = 50,
  Ecoli_cfu_100ml = 1,     # limit is 0; use 1 to avoid /0
  Turbidity_NTU   = 5
)

param_labels <- c(
  pH              = "pH",
  TDS_mgL         = "TDS (mg/L)",
  Nitrate_mgL     = "Nitrate (mg/L)",
  Ecoli_cfu_100ml = "E. coli (cfu/100 mL)",
  Turbidity_NTU   = "Turbidity (NTU)"
)

plot_data <- kamal %>%
  select(Sample, names(who_upper)) %>%
  pivot_longer(cols = -Sample, names_to = "Parameter", values_to = "Value") %>%
  mutate(
    WHO_limit       = who_upper[Parameter],
    Pct_of_limit    = (Value / WHO_limit) * 100,
    Parameter_label = factor(param_labels[Parameter],
                             levels = param_labels)
  )

p_who <- ggplot(plot_data, aes(x = Parameter_label, y = Pct_of_limit,
                               color = Parameter_label)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = 100,
           fill = "#d4edda", alpha = 0.35) +
  geom_hline(yintercept = 100, linetype = "dashed",
             colour = "#c0392b", linewidth = 1) +
  annotate("text", x = 0.55, y = 103, label = "WHO limit (100 %)",
           hjust = 0, size = 3.2, colour = "#c0392b") +
  geom_jitter(width = 0.18, size = 3, alpha = 0.75) +
  stat_summary(fun = median, geom = "crossbar",
               width = 0.45, linewidth = 0.6, colour = "grey20") +
  scale_color_manual(values = c("#1f4e79","#2980b9","#27ae60",
                                "#e67e22","#8e44ad")) +
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     breaks = seq(0, max(plot_data$Pct_of_limit, na.rm = TRUE) + 50,
                                  by = 50)) +
  scale_x_discrete() +   # <--- FIX: explicitly declare discrete x-axis
  labs(
    title    = "Water Quality Parameters vs. WHO Standards",
    subtitle = "Values expressed as % of WHO upper limit  |  n = 15 wells  |  crossbar = median",
    x        = NULL,
    y        = "% of WHO limit",
    caption  = "Note: E. coli WHO limit = 0 cfu/100 mL; plotted against reference unit of 1."
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title       = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle    = element_text(hjust = 0.5, size = 10, colour = "grey40"),
    axis.text.x      = element_text(angle = 25, hjust = 1, size = 11),
    legend.position  = "none",
    panel.grid.minor = element_blank()
  )

print(p_who)
ggsave(file.path(out_dir, "plot_who_comparison.png"), p_who,
       width = 10, height = 6, dpi = 300)
cat("Saved: plot_who_comparison.png\n")


# INFERENTIAL STATISTICS

# Pearson Correlation matrix 
cat("\n--- Pearson Correlation (numeric parameters) ---\n")

num_vars <- kamal %>%
  select(pH, TDS_mgL, EC_uScm, Turbidity_NTU,
         TCC_10x3_cfu_ml, Nitrate_mgL, Ecoli_cfu_100ml, Distance_m)

cor_matrix <- cor(num_vars, use = "pairwise.complete.obs", method = "pearson")
print(round(cor_matrix, 3))

cor_pairs <- combn(names(num_vars), 2, simplify = FALSE)

cor_results <- lapply(cor_pairs, function(pair) {
  test <- cor.test(num_vars[[pair[1]]], num_vars[[pair[2]]], method = "pearson")
  data.frame(
    Var1    = pair[1],
    Var2    = pair[2],
    r       = round(test$estimate, 4),
    t_stat  = round(test$statistic, 4),
    df      = test$parameter,
    p_value = round(test$p.value, 4),
    Sig     = ifelse(test$p.value < 0.001, "***",
                     ifelse(test$p.value < 0.01,  "**",
                            ifelse(test$p.value < 0.05,  "*", "ns")))
  )
})

cor_table <- bind_rows(cor_results)
print(cor_table)
write_xlsx(cor_table, file.path(out_dir, "inferential_correlation.xlsx"))
cat("Saved: inferential_correlation.xlsx\n")

# ── 6b. One-sample t-test (sample mean vs. WHO upper limit) ────
cat("\n--- One-sample t-test (observed mean vs. WHO upper limit) ---\n")

who_ttest <- list(
  pH            = 8.5,
  TDS_mgL       = 500,
  Nitrate_mgL   = 50,
  Turbidity_NTU = 5
)

ttest_results <- lapply(names(who_ttest), function(param) {
  vals <- kamal[[param]]
  mu   <- who_ttest[[param]]
  test <- t.test(vals, mu = mu, alternative = "two.sided")
  data.frame(
    Parameter   = param,
    n           = length(vals),
    Sample_mean = round(mean(vals, na.rm = TRUE), 4),
    WHO_limit   = mu,
    t_stat      = round(test$statistic, 4),
    df          = test$parameter,
    p_value     = round(test$p.value, 4),
    CI_lower    = round(test$conf.int[1], 4),
    CI_upper    = round(test$conf.int[2], 4),
    Sig         = ifelse(test$p.value < 0.001, "***",
                         ifelse(test$p.value < 0.01,  "**",
                                ifelse(test$p.value < 0.05,  "*", "ns")))
  )
})

ttest_table <- bind_rows(ttest_results)
print(ttest_table)
write_xlsx(ttest_table,  file.path(out_dir, "inferential_ttest.xlsx"))
cat("Saved: inferential_ttest.xlsx\n")

# ── 6c. Kruskal–Wallis test ─────────────────────────────────────
cat("\n--- Kruskal-Wallis test (parameter ~ Sample) ---\n")

kw_params <- c("pH", "TDS_mgL", "EC_uScm", "Turbidity_NTU",
               "TCC_10x3_cfu_ml", "Nitrate_mgL", "Ecoli_cfu_100ml")

kw_results <- lapply(kw_params, function(param) {
  test <- kruskal.test(kamal[[param]] ~ factor(kamal$Sample))
  data.frame(
    Parameter = param,
    chi_sq    = round(test$statistic, 4),
    df        = test$parameter,
    p_value   = round(test$p.value, 4),
    Sig       = ifelse(test$p.value < 0.001, "***",
                       ifelse(test$p.value < 0.01,  "**",
                              ifelse(test$p.value < 0.05,  "*", "ns")))
  )
})

kw_table <- bind_rows(kw_results)
print(kw_table)
write_xlsx(kw_table, file.path(out_dir, "inferential_kruskal_wallis.xlsx"))
cat("Saved: inferential_kruskal_wallis.xlsx\n")

cat("\nSignificance codes:  *** p<0.001  ** p<0.01  * p<0.05  ns = not significant\n")


# ======================================================
# 7. NEW: NORMALITY TESTING (Shapiro-Wilk)
# ======================================================
# Shapiro-Wilk is appropriate for small samples (n=15).
# Null hypothesis: data are normally distributed.
# p < 0.05 → reject normality → use non-parametric tests.
# This justifies using Spearman correlation and Mann-Whitney U
# for bacterial counts (E. coli, TCC) which are typically non-normal.

cat("\n--- Shapiro-Wilk Normality Test ---\n")

shapiro_params <- c("pH", "TDS_mgL", "EC_uScm", "Turbidity_NTU",
                    "TCC_10x3_cfu_ml", "Nitrate_mgL",
                    "Ecoli_cfu_100ml", "Distance_m")

shapiro_results <- lapply(shapiro_params, function(param) {
  vals <- na.omit(kamal[[param]])
  test <- shapiro.test(vals)
  data.frame(
    Parameter  = param,
    n          = length(vals),
    W_statistic = round(test$statistic, 4),
    p_value    = round(test$p.value, 4),
    Normal     = ifelse(test$p.value >= 0.05, "YES", "NO"),
    Recommended_test = ifelse(test$p.value >= 0.05,
                              "Parametric (t-test / Pearson)",
                              "Non-parametric (Mann-Whitney / Spearman)")
  )
})

shapiro_table <- bind_rows(shapiro_results)
print(shapiro_table)
write_xlsx(shapiro_table, file.path(out_dir, "normality_shapiro_wilk.xlsx"))
cat("Saved: normality_shapiro_wilk.xlsx\n")



# SPEARMAN RANK CORRELATION
# Spearman is robust to non-normality and outliers —
# important here since E. coli and TCC are likely non-normal.
# We compute both Pearson and Spearman for all pairs and
# present a combined comparison table.

cat("\n--- Spearman Rank Correlation ---\n")

spearman_results <- lapply(cor_pairs, function(pair) {
  test <- cor.test(num_vars[[pair[1]]], num_vars[[pair[2]]], method = "spearman",
                   exact = FALSE)   # exact=FALSE avoids warning with ties
  data.frame(
    Var1      = pair[1],
    Var2      = pair[2],
    rho       = round(test$estimate, 4),
    S_stat    = round(test$statistic, 4),
    p_value   = round(test$p.value, 4),
    Sig       = ifelse(test$p.value < 0.001, "***",
                       ifelse(test$p.value < 0.01,  "**",
                              ifelse(test$p.value < 0.05,  "*", "ns")))
  )
})

spearman_table <- bind_rows(spearman_results)
print(spearman_table)

# Combined Pearson + Spearman comparison table
combined_cor <- cor_table %>%
  select(Var1, Var2, Pearson_r = r, Pearson_p = p_value, Pearson_Sig = Sig) %>%
  left_join(
    spearman_table %>%
      select(Var1, Var2, Spearman_rho = rho, Spearman_p = p_value, Spearman_Sig = Sig),
    by = c("Var1", "Var2")
  )

print(combined_cor)
write_xlsx(
  list(Pearson   = cor_table,
       Spearman  = spearman_table,
       Combined  = combined_cor), 
  file.path(
  out_dir, "correlation_pearson_spearman.xlsx"
))
cat("Saved: correlation_pearson_spearman.xlsx (3 sheets)\n")


# SIMPLE LINEAR REGRESSION (Distance vs Contaminants)
# For each key contaminant, fit: Contaminant ~ Distance_m
# Output: slope, 95% CI, R², adjusted R², p-value, interpretation.
# NOTE: With n=15, power is limited — results are indicative.

cat("\n--- Simple Linear Regression: Distance vs Contaminants ---\n")

reg_targets <- c("Ecoli_cfu_100ml", "TCC_10x3_cfu_ml", "Nitrate_mgL")

reg_results <- lapply(reg_targets, function(outcome) {
  formula  <- as.formula(paste(outcome, "~ Distance_m"))
  model    <- lm(formula, data = kamal)
  s        <- summary(model)
  ci       <- confint(model)
  slope    <- coef(model)["Distance_m"]
  slope_ci <- ci["Distance_m", ]
  p_slope  <- coef(s)["Distance_m", "Pr(>|t|)"]
  
  # Human-readable interpretation
  direction <- ifelse(slope < 0, "decreases", "increases")
  interp <- sprintf(
    "Each additional 1 m from septic tank %s %s by %.2f units (95%% CI: %.2f to %.2f).",
    direction, outcome, abs(slope), slope_ci[1], slope_ci[2]
  )
  
  cat("\n>>>", outcome, "~ Distance_m\n")
  print(s)
  cat("Interpretation:", interp, "\n")
  
  data.frame(
    Outcome       = outcome,
    Intercept     = round(coef(model)[1], 3),
    Slope         = round(slope, 3),
    CI_lower_95   = round(slope_ci[1], 3),
    CI_upper_95   = round(slope_ci[2], 3),
    R_squared     = round(s$r.squared, 4),
    Adj_R_squared = round(s$adj.r.squared, 4),
    F_stat        = round(s$fstatistic[1], 4),
    p_value_slope = round(p_slope, 4),
    Sig           = ifelse(p_slope < 0.001, "***",
                           ifelse(p_slope < 0.01,  "**",
                                  ifelse(p_slope < 0.05,  "*", "ns"))),
    Interpretation = interp
  )
})

reg_table <- bind_rows(reg_results)
print(reg_table %>% select(-Interpretation))
write_xlsx(reg_table, file.path( out_dir, "regression_distance_contaminants.xlsx"))
cat("Saved: regression_distance_contaminants.xlsx\n")



# 10. NEW: SCATTER PLOTS WITH REGRESSION LINES
# Three publication-ready plots: Distance vs E. coli, Nitrate, TCC.
# Grey band = 95% confidence interval around the regression line.

scatter_specs <- list(
  list(y = "Ecoli_cfu_100ml",   ylab = "E. coli (cfu/100 mL)",
       file = "project res/scatter_distance_ecoli.png"),
  list(y = "Nitrate_mgL",       ylab = "Nitrate (mg/L)",
       file = "project res/scatter_distance_nitrate.png"),
  list(y = "TCC_10x3_cfu_ml",   ylab = "Total Coliform (×10³ cfu/mL)",
       file = "project res/scatter_distance_tcc.png")
)

for (spec in scatter_specs) {
  y_col  <- spec$y
  r2_val <- reg_table %>% filter(Outcome == y_col) %>% pull(R_squared)
  pval   <- reg_table %>% filter(Outcome == y_col) %>% pull(p_value_slope)
  sig    <- reg_table %>% filter(Outcome == y_col) %>% pull(Sig)
  
  p_scatter <- ggplot(kamal, aes_string(x = "Distance_m", y = y_col)) +
    geom_point(size = 3.5, colour = "#1f4e79", alpha = 0.85) +
    geom_smooth(method = "lm", se = TRUE,
                colour = "#c0392b", fill = "#f1948a", linewidth = 1.1) +
    geom_text(aes(label = Sample), vjust = -0.8, size = 3, colour = "grey40") +
    # Annotate R² and p-value in top-right corner
    annotate("text",
             x = max(kamal$Distance_m, na.rm = TRUE) * 0.72,
             y = max(kamal[[y_col]], na.rm = TRUE) * 0.95,
             label = paste0("R\u00b2 = ", round(r2_val, 3),
                            "\np = ", round(pval, 4), " ", sig),
             size = 3.8, hjust = 0, colour = "#2c3e50") +
    labs(
      title   = paste("Distance to Septic Tank vs.", spec$ylab),
      x       = "Distance to Septic Tank (m)",
      y       = spec$ylab,
      caption = "Red line = OLS regression; shaded band = 95% CI"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title   = element_text(hjust = 0.5, face = "bold"),
      plot.caption = element_text(size = 9, colour = "grey50")
    )
  
  print(p_scatter)
  ggsave(file.path(out_dir, plot = p_scatter),
         width = 8, height = 5, dpi = 300)
  cat("Saved:", spec$file, "\n")
}



# GROUP COMPARISON – CLOSE vs FAR WELLS
# Split at median distance OR 15 m (Nigerian/WHO practical guideline).
# Data: Distance range 8–25 m, median computed below.
# We use 15 m as the threshold — it is the commonly cited minimum
# safe separation in Nigerian rural water guidelines (SON 2015),
# and it falls near the median of this dataset.

med_dist <- median(kamal$Distance_m, na.rm = TRUE)
cat(sprintf("\nMedian distance = %.1f m  |  Using threshold = 15 m for grouping\n",
            med_dist))

kamal <- kamal %>%
  mutate(
    Dist_Group = ifelse(Distance_m < 15, "Close (<15 m)", "Far (\u226515 m)"),
    Dist_Group = factor(Dist_Group, levels = c("Close (<15 m)", "Far (\u226515 m)"))
  )

# Count wells per group
cat("\nWells per distance group:\n")
print(table(kamal$Dist_Group))

# Mann-Whitney U test for each contaminant between the two groups
# (preferred given non-normality of bacterial counts)
group_params <- c("Ecoli_cfu_100ml", "TCC_10x3_cfu_ml", "Nitrate_mgL")

mw_group_results <- lapply(group_params, function(param) {
  close_vals <- kamal %>% filter(Dist_Group == "Close (<15 m)") %>% pull(!!sym(param))
  far_vals   <- kamal %>% filter(Dist_Group == "Far (\u226515 m)")   %>% pull(!!sym(param))
  
  # Mann-Whitney U (Wilcoxon rank-sum)
  mw <- wilcox.test(close_vals, far_vals, exact = FALSE)
  
  data.frame(
    Parameter     = param,
    Close_median  = round(median(close_vals, na.rm = TRUE), 3),
    Far_median    = round(median(far_vals,   na.rm = TRUE), 3),
    Close_n       = sum(!is.na(close_vals)),
    Far_n         = sum(!is.na(far_vals)),
    W_statistic   = round(mw$statistic, 3),
    p_value       = round(mw$p.value, 4),
    Sig           = ifelse(mw$p.value < 0.001, "***",
                           ifelse(mw$p.value < 0.01,  "**",
                                  ifelse(mw$p.value < 0.05,  "*", "ns"))),
    Interpretation = ifelse(
      mw$p.value < 0.05,
      paste(param, "significantly higher in close wells (p<0.05)"),
      paste(param, "not significantly different between groups (p\u22650.05)")
    )
  )
})

mw_group_table <- bind_rows(mw_group_results)
cat("\n--- Mann-Whitney U: Close vs Far Wells ---\n")
print(mw_group_table)
write_xlsx(mw_group_table, file.path(out_dir,  "group_comparison_close_far.xlsx"))
cat("Saved: group_comparison_close_far.xlsx\n")

# Boxplot: E. coli by distance group
p_box_group <- ggplot(kamal, aes(x = Dist_Group, y = Ecoli_cfu_100ml,
                                 fill = Dist_Group)) +
  geom_boxplot(alpha = 0.75, outlier.shape = 21,
               outlier.fill = "white", outlier.size = 2) +
  geom_jitter(width = 0.12, size = 2.8, alpha = 0.8, colour = "grey25") +
  geom_text(aes(label = Sample), vjust = -0.9, size = 3, colour = "grey35") +
  scale_fill_manual(
    values = c("Close (<15 m)" = "#e74c3c", "Far (\u226515 m)" = "#2980b9"),
    name   = "Distance Group"
  ) +
  labs(
    title    = "E. coli Contamination by Distance to Septic Tank",
    subtitle = paste0("Threshold = 15 m  |  Mann-Whitney p = ",
                      mw_group_table %>%
                        filter(Parameter == "Ecoli_cfu_100ml") %>%
                        pull(p_value)),
    x        = "Distance Group",
    y        = "E. coli (cfu/100 mL)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10, colour = "grey40"),
    legend.position = "none"
  )

print(p_box_group)
ggsave(file.path(out_dir, "plot_who_comparison.png"), p_who,
       width = 10, height = 6, dpi = 300)
cat("Saved: plot_who_comparison.png\n")

ggsave(file.path(out_dir, "boxplot_ecoli_distance_group.png"),
       plot = p_box_group, width = 7, height = 5, dpi = 300)
cat("Saved: boxplot_ecoli_distance_group.png\n")

# ================================================================
# NEW: MINIMUM SAFE DISTANCE CHECK
# Compare observed well-to-septic-tank distances against
# published guidelines:
#   - WHO / FAO recommendation: ≥ 30 m
#   - Nigerian SON 2015 / FMEnv guideline: ≥ 15 m
#   - UNICEF rural water standard: ≥ 15–30 m (context-dependent)

cat("\n--- Minimum Safe Distance Assessment ---\n")

who_min_dist  <- 30   # WHO/FAO recommended minimum (metres)
son_min_dist  <- 15   # Nigerian SON 2015 practical minimum (metres)

distance_check <- kamal %>%
  select(Sample, Distance_m) %>%
  mutate(
    Meets_SON_15m  = ifelse(Distance_m >= son_min_dist, "YES", "NO"),
    Meets_WHO_30m  = ifelse(Distance_m >= who_min_dist, "YES", "NO")
  )

cat("\nDistance compliance summary:\n")
print(distance_check)

cat(sprintf(
  "\nDistance range: %.0f m – %.0f m (mean = %.1f m, median = %.1f m)\n",
  min(kamal$Distance_m), max(kamal$Distance_m),
  mean(kamal$Distance_m), median(kamal$Distance_m)
))

n_son <- sum(distance_check$Meets_SON_15m == "YES")
n_who <- sum(distance_check$Meets_WHO_30m == "YES")

cat(sprintf("  Wells meeting Nigerian SON ≥15 m guideline: %d / %d (%.0f%%)\n",
            n_son, nrow(kamal), 100 * n_son / nrow(kamal)))
cat(sprintf("  Wells meeting WHO ≥30 m guideline:          %d / %d (%.0f%%)\n",
            n_who, nrow(kamal), 100 * n_who / nrow(kamal)))


write_xlsx(distance_check, file.path(
  out_dir, "minimum_safe_distance_check.xlsx"))
cat("Saved: minimum_safe_distance_check.xlsx\n")


# ======================================================
# COMPLETION SUMMARY
cat("\n====================================================\n")
cat("  All exports completed. Files saved in 'project res/':\n")
cat("  EXISTING:\n")
cat("    - descriptive_statistics.xlsx\n")
cat("    - compliance_WHO_NSDWQ.xlsx\n")
cat("    - distance_contamination.xlsx\n")
cat("    - plot_who_comparison.png\n")
cat("    - inferential_correlation.xlsx\n")
cat("    - inferential_ttest.xlsx\n")
cat("    - inferential_kruskal_wallis.xlsx\n")
cat("  NEW:\n")
cat("    - normality_shapiro_wilk.xlsx\n")
cat("    - correlation_pearson_spearman.xlsx  (3 sheets)\n")
cat("    - regression_distance_contaminants.xlsx\n")
cat("    - scatter_distance_ecoli.png\n")
cat("    - scatter_distance_nitrate.png\n")
cat("    - scatter_distance_tcc.png\n")
cat("    - group_comparison_close_far.xlsx\n")
cat("    - boxplot_ecoli_distance_group.png\n")
cat("    - minimum_safe_distance_check.xlsx\n")
cat("====================================================\n")

