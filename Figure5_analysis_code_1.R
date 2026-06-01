################################################################################
# Title:   Statistical Analysis and Visualization Code for Figure 5
# Journal: Nature Communications
# Description: This script reproduces all panels of Figure 5, including:
#              (a) Univariate Cox proportional hazards forest plot
#              (b) Kaplan-Meier survival curve (generated separately in GraphPad)
#              (c) Multiple ROC curves for immune correlates of protection
#              (d) hSBA titer grouped by infection clearance day (generated separately)
#              (e) Per-mouse infection clearance dot plot
#              (Supplementary) Spearman correlation heatmap
#
# Input data: Tab-delimited .txt files containing per-animal immunological
#             measurements and infection outcome data from the murine
#             gonococcal challenge model.
#
# Variables:
#   hSBA_Titre   - Human serum bactericidal assay (hSBA) titer (pre-challenge)
#   Lysozyme     - Cervicovaginal lysozyme activity
#   Serum_IgA    - Serum IgA endpoint titer
#   Serum_IgG    - Total serum IgG endpoint titer
#   Serum_IgG1   - Serum IgG1 subclass titer
#   Serum_IgG2c  - Serum IgG2c subclass titer
#   Serum_IgG3   - Serum IgG3 subclass titer
#   Vaginal_IgA  - Vaginal IgA endpoint titer
#   Vaginal_IgG  - Vaginal IgG endpoint titer
#   status       - Binary infection outcome: 1 = cleared, 0 = not cleared
#   clearance    - Day of infection clearance (Day3, Day5, Day7, or "never")
#
# Export format: PDF (vector) — preferred by Nature Communications for line art
#                and ggplot-based figures. Resolution-independent and fully
#                editable by journal production staff.
#                TIFF at 600 DPI also exported as fallback raster format.
#
# Nature Communications figure width specifications:
#   Single column : 88 mm  (3.46 in)
#   1.5 column    : 120 mm (4.72 in)
#   Double column : 180 mm (7.09 in)
#
# R version: 4.3.x or later recommended
# Required packages: survival, broom, purrr, ggplot2, dplyr, forcats,
#                    rstatix, pROC, Hmisc
#
# Authors: [Author names]
# Date:    [Submission date]
################################################################################


# ==============================================================================
# 0. PACKAGE LOADING
# ==============================================================================

library(survival)    # Cox proportional hazards models
library(broom)       # Tidying model outputs (tidy())
library(purrr)       # Functional programming tools (map_dfr())
library(ggplot2)     # Data visualization
library(dplyr)       # Data manipulation
library(forcats)     # Factor reordering (fct_inorder())
library(rstatix)     # p_format() for significance annotation
library(pROC)        # ROC curve generation and AUC calculation
library(Hmisc)       # rcorr() for Spearman correlation with p-values


# ==============================================================================
# EXPORT SETTINGS — Nature Communications specifications
# ==============================================================================
# Adjust width per panel based on how many columns the figure occupies in print.
# All ggplot panels are saved as both:
#   (1) PDF  — vector format, preferred for line art (resolution-independent)
#   (2) TIFF — raster fallback at 600 DPI (Nature Communications requirement
#              for line art submitted as bitmap)
#
# Create an output folder for all exported figures
if (!dir.exists("figures")) dir.create("figures")

# Helper function: saves a ggplot object as both PDF and TIFF
# Arguments:
#   plot_obj  - a ggplot object
#   filename  - base name without extension (e.g. "Figure5a")
#   width_in  - width in inches (use NComms column widths below)
#   height_in - height in inches
save_figure <- function(plot_obj, filename, width_in, height_in) {
  
  # PDF — vector, preferred submission format
  ggsave(
    filename  = file.path("figures", paste0(filename, ".pdf")),
    plot      = plot_obj,
    width     = width_in,
    height    = height_in,
    units     = "in",
    device    = "pdf"
  )
  
  # TIFF — 600 DPI raster fallback for line art per NComms guidelines
  ggsave(
    filename  = file.path("figures", paste0(filename, ".tiff")),
    plot      = plot_obj,
    width     = width_in,
    height    = height_in,
    units     = "in",
    dpi       = 600,
    device    = "tiff",
    compression = "lzw"   # LZW compression keeps TIFF file size manageable
  )
  
  message("Saved: ", filename, ".pdf and ", filename, ".tiff")
}

# Nature Communications column width presets (in inches)
W_SINGLE <- 3.46   # 88 mm  — single column
W_HALF   <- 4.72   # 120 mm — 1.5 column
W_DOUBLE <- 7.09   # 180 mm — double column


# ==============================================================================
# 1. DATA LOADING
# ==============================================================================
# Replace the file path below with the path to your local copy of the data file.
# The input file is a tab-delimited text file with one row per animal and
# columns for each immune measurement and outcome variable.
#
# Expected column structure:
#   mouse | hSBA_Titre | Lysozyme | Serum_IgA | Serum_IgG | Serum_IgG1 |
#   Serum_IgG2c | Serum_IgG3 | Vaginal_IgA | Vaginal_IgG | status | clearance
# ------------------------------------------------------------------------------

dd <- read.table("your_data_file.txt", sep = "\t", header = TRUE, check.names = FALSE)

# ------------------------------------------------------------------------------
# Note: 'status' must be binary (1 = infection cleared, 0 = not cleared).
# 'hSBA_Titre' and all antibody/enzyme variables should be numeric.
# ------------------------------------------------------------------------------


# ==============================================================================
# FIGURE 5a: UNIVARIATE COX PROPORTIONAL HAZARDS FOREST PLOT
# ==============================================================================
# Purpose: Assess whether each pre-challenge immune measurement independently
#          predicts time to infection clearance using Cox regression.
#          Hazard ratios (HR) > 1 indicate faster clearance with higher values.
# ------------------------------------------------------------------------------

# --- Step 1: Fit univariate Cox models for each predictor -------------------

# Define predictors to test
predictors <- c("hSBA_Titre", "Lysozyme", "Serum_IgA", "Serum_IgG",
                "Serum_IgG1", "Serum_IgG2c", "Serum_IgG3",
                "Vaginal_IgA", "Vaginal_IgG")

# Fit one Cox model per predictor using a named list
# 'status' is the event indicator; time-to-event is encoded in the data structure
univ_models <- lapply(predictors, function(var) {
  formula <- as.formula(paste("Surv(time, status) ~", var))
  coxph(formula, data = dd)
})
names(univ_models) <- predictors

# Tidy all models into a single data frame with exponentiated coefficients (HR)
# and 95% confidence intervals
results <- map_dfr(univ_models, tidy,
                   exponentiate = TRUE,
                   conf.int     = TRUE,
                   .id          = "model_id")
rr <- data.frame(results)

# Preview of expected output (actual values will depend on your dataset):
#   model_id        term    estimate    std.error  statistic      p.value    conf.low    conf.high
# 1        1  hSBA_Titre   1.0002561 6.660053e-05  3.8447820 0.0001206596 1.000125538     1.000387
# 2        2    Lysozyme   1.0019360 8.272869e-03  0.2337905 0.8151475982 0.985821074     1.018314
# 3        3   Serum_IgA   1.9149864 2.889146e-01  2.2487979 0.0245253547 1.087029062     3.373574
# ...

# --- Step 2: Add significance labels based on p-value thresholds ------------

r1 <- results
r1$sig <- NA
r1$sig[r1$p.value < 0.001]                        <- "***"
r1$sig[r1$p.value < 0.01  & r1$p.value >= 0.001]  <- "**"
r1$sig[r1$p.value < 0.05  & r1$p.value >= 0.01]   <- "*"

# Sort by descending p-value so the most significant predictor appears at the top
r1 <- r1 %>% arrange(desc(p.value))

# --- Step 3: Build and export forest plot -----------------------------------
# Width: double column (7.09 in) — needed to accommodate the annotation columns
# Height: 5 in — one row per predictor (9 rows) with comfortable spacing

fig5a <- ggplot(r1, aes(x = estimate, y = fct_inorder(term), color = estimate > 1)) +

  # Log10 x-axis to handle the wide range of HRs (from ~0.94 to ~131)
  scale_x_log10(limits = c(0.001, 1e9),
                breaks  = c(0.001, 0.1, 1, 10, 100, 1000, 10000, 100000)) +

  # Vertical dashed reference line at HR = 1 (null hypothesis)
  geom_vline(xintercept = 1, linetype = "dashed", color = "brown") +

  # Point-range geom: square point (pch=15) with 95% CI whiskers
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high),
                  size      = 1.2,
                  pch       = 15,
                  linewidth = 1.2) +

  # Right-hand annotation column 1: HR [95% CI] values
  # Placed at x = 1e6 (beyond the largest CI upper bound of ~63,454)
  geom_text(aes(label = paste0(round(estimate, 2), " (",
                               round(conf.low, 2), "-",
                               round(conf.high, 2), ")")),
            x     = 1e6,
            hjust = 0,
            size  = 4,
            color = "black") +

  # Right-hand annotation column 2: formatted p-values with significance stars
  # Placed at x = 3e7
  geom_text(aes(label = paste0(p_format(p.value, accuracy = 0.001, add_p = TRUE),
                               " ", ifelse(is.na(sig), "", sig))),
            x     = 3e7,
            hjust = 0,
            size  = 4,
            color = "black") +

  # Color: red = HR > 1 (protective direction), green = HR <= 1
  scale_color_manual(values = c("green4", "red3"), guide = "none") +

  # Title uses spacing to align "HR [95% CI]" and "p-value" over their columns
  labs(
    title = paste0(paste(rep(" ", 18), collapse = ""), "Hazard Ratio Plot",
                   paste(rep(" ", 22), collapse = ""), "HR [95% CI]",
                   paste(rep(" ", 10), collapse = ""), "p-value"),
    x = "Hazard Ratio -- 95% CI (Log Scale)",
    y = ""
  ) +

  theme_minimal() +
  theme(
    text             = element_text(size = 15, face = "bold"),
    axis.text.x      = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position  = "none",
    axis.title.x     = element_text(color = "dimgray", size = 14, face = "bold",
                                    margin = margin(t = 20)),
    axis.title.y     = element_text(color = "dimgray", size = 14, face = "bold"),
    panel.grid.minor = element_blank()
  )

# Export Figure 5a
# Double-column width to fit the HR + p-value annotation columns on the right
save_figure(fig5a, filename = "Figure5a_forest_plot", width_in = W_DOUBLE, height_in = 5)


# ==============================================================================
# FIGURE 5b: KAPLAN-MEIER SURVIVAL CURVE
# ==============================================================================
# Note: The Kaplan-Meier curve comparing mice with hSBA >= 2048 vs < 2048
#       (log-rank test, p = 0.0061) was generated and exported using
#       GraphPad Prism (version 10). The survival object can also be constructed
#       in R as shown below for reproducibility reference.
# ------------------------------------------------------------------------------

# Dichotomize hSBA titer at the threshold of 2048
dd$hSBA_group <- ifelse(dd$hSBA_Titre >= 2048, "hSBA >= 2048", "hSBA < 2048")

# Fit Kaplan-Meier survival curve
km_fit <- survfit(Surv(time, status) ~ hSBA_group, data = dd)

# Log-rank test
survdiff(Surv(time, status) ~ hSBA_group, data = dd)

# Note: To export the KM plot from R if needed (e.g., using survminer):
# library(survminer)
# fig5b <- ggsurvplot(km_fit, data = dd, pval = TRUE, conf.int = FALSE,
#                     palette = c("orange", "darkblue"))$plot
# save_figure(fig5b, filename = "Figure5b_KM_curve",
#             width_in = W_SINGLE, height_in = 3.5)


# ==============================================================================
# FIGURE 5c: MULTIPLE ROC CURVES FOR ALL IMMUNE PREDICTORS
# ==============================================================================
# Purpose: Evaluate the discriminative ability of each immune measurement
#          for predicting infection clearance (AUC = area under the ROC curve).
#          'status' is used as the binary outcome (1 = cleared, 0 = not cleared).
# ------------------------------------------------------------------------------

# Build ROC objects for each predictor
roc_list <- list(
  hSBA        = roc(dd$status, dd$hSBA_Titre),
  Lysozyme    = roc(dd$status, dd$Lysozyme),
  Serum_IgA   = roc(dd$status, dd$Serum_IgA),
  Serum_IgG   = roc(dd$status, dd$Serum_IgG),
  Serum_IgG1  = roc(dd$status, dd$Serum_IgG1),
  Serum_IgG2c = roc(dd$status, dd$Serum_IgG2c),
  Serum_IgG3  = roc(dd$status, dd$Serum_IgG3),
  Vaginal_IgA = roc(dd$status, dd$Vaginal_IgA),
  Vaginal_IgG = roc(dd$status, dd$Vaginal_IgG)
)

# Generate legend labels with AUC values embedded
auc_labels <- sapply(names(roc_list), function(name) {
  auc_val <- round(auc(roc_list[[name]]), 2)
  paste0(name, " (AUC = ", auc_val, ")")
})

# Plot all ROC curves on a single panel
# legacy.axes = TRUE converts x-axis to 1-Specificity (FPR) convention
# hSBA receives a thicker line (linewidth = 2) as the primary predictor
fig5c <- ggroc(roc_list, legacy.axes = TRUE) +
  aes(linewidth = name) +

  scale_x_continuous(breaks = seq(0, 1, 0.2),
                     labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0")) +

  # Thicker line for hSBA (index 1), standard width for all others
  scale_linewidth_manual(values = c(2, rep(1, 8)), guide = "none") +

  # Custom color palette: brown for hSBA, distinct palette for remaining curves
  scale_color_manual(
    values = c("Brown", "#F3DF6C", "#CEAB07", "#D5D5D3", "#24281A",
               "#85D4E3", "#F4B5BD", "#9C964A", "#CDC08C", "#FAD77B"),
    labels = auc_labels
  ) +

  # Diagonal reference line (random classifier, AUC = 0.5)
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
               color     = "lightgrey",
               linetype  = "dashed",
               linewidth = 0.5) +

  labs(x = "1 - Specificity (FPR)", y = "Sensitivity (TPR)") +

  theme(
    text             = element_text(size = 15, face = "bold"),
    axis.text.x      = element_text(angle = 0, hjust = 0.1, vjust = 0.2),
    legend.position  = "right",
    legend.title     = element_blank(),
    axis.title.x     = element_text(color = "dimgray", size = 14, face = "bold",
                                    margin = margin(t = 20, r = 20, b = 0, l = 0)),
    axis.title.y     = element_text(color = "dimgray", size = 14, face = "bold")
  )

# Export Figure 5c
# Double-column width to accommodate the legend with AUC values on the right
save_figure(fig5c, filename = "Figure5c_ROC_curves", width_in = W_DOUBLE, height_in = 5)


# ==============================================================================
# FIGURE 5d: hSBA TITER GROUPED BY INFECTION CLEARANCE DAY
# ==============================================================================
# Note: Box plots with individual data points comparing hSBA titers across
#       clearance groups (Day 3, Day 5, Day 7, Never) with Kruskal-Wallis
#       and Dunn post-hoc tests (p = 0.05 and p = 0.0023) were generated
#       and exported using GraphPad Prism (version 10).
# ------------------------------------------------------------------------------


# ==============================================================================
# FIGURE 5e: PER-MOUSE INFECTION CLEARANCE DOT PLOT
# ==============================================================================
# Purpose: Visualize individual mouse infection status at each time point
#          post-challenge. Each row is one mouse, ordered by descending hSBA
#          titer. Point size encodes bacterial burden; color encodes the day
#          of clearance. The hSBA titer label is printed to the right.
#
# Required data frame columns:
#   mouse      - Mouse ID (numeric or character)
#   Day        - Day of post-challenge assessment (1, 3, 5, 7)
#   value      - Bacterial CFU count or equivalent burden measure (numeric)
#   hSBA_Titre - Pre-challenge hSBA titer (used for y-axis ordering)
#   clearance  - Day of clearance or "never" (character/factor)
# ------------------------------------------------------------------------------

fig5e <- ggplot(df, aes(x     = Day,
                        y     = reorder(mouse, as.numeric(hSBA_Titre)),
                        size  = as.numeric(value),
                        color = as.character(clearance))) +
  geom_point() +

  # X-axis shows only assessment days; extended to 8.5 to accommodate text labels
  scale_x_continuous(limits = c(1, 8.5),
                     breaks = c(1, 3, 5, 7),
                     labels = c("Day1", "Day3", "Day5", "Day7")) +

  # Print hSBA titer values to the right of each row for direct reference
  geom_text(aes(label = hSBA_Titre, x = 8, hjust = 0, size = 6),
            color = "dimgray") +

  # Color encodes day of clearance; orange = never cleared
  scale_color_manual(values = c("Day3"  = "darkblue",
                                "Day5"  = "skyblue3",
                                "Day7"  = "skyblue1",
                                "never" = "darkorange3")) +

  # Minimalist theme: remove all background elements
  theme(
    panel.background    = element_blank(),
    plot.background     = element_blank(),
    panel.grid.major    = element_blank(),
    panel.grid.minor    = element_blank(),
    panel.border        = element_blank(),
    text                = element_text(size = 15, face = "bold"),
    axis.text.x         = element_text(angle = 0, hjust = 0.1, vjust = 0.2),
    legend.position     = "right",
    legend.title        = element_blank(),
    plot.title.position = "plot",
    plot.title          = element_text(size = 15),
    axis.title.x        = element_text(color = "dimgray", size = 14, face = "bold",
                                       margin = margin(t = 20, r = 20, b = 0, l = 0)),
    axis.title.y        = element_text(color = "dimgray", size = 14, face = "bold")
  ) +

  xlab("Days post-challenge") +
  ylab("") +

  # Title spacing aligns "Mouse" and "hSBA titre" over their respective columns
  labs(title = paste0("Mouse", paste(rep(" ", 135), collapse = ""), "hSBA titre"))

# Export Figure 5e
# Double-column width: 19 mouse rows need full width; height scales with n rows
save_figure(fig5e, filename = "Figure5e_per_mouse_clearance", width_in = W_DOUBLE, height_in = 6.5)


# ==============================================================================
# SUPPLEMENTARY FIGURE: SPEARMAN CORRELATION HEATMAP
# ==============================================================================
# Purpose: Pairwise Spearman correlations between hSBA titer and all measured
#          immune parameters. Tiles are color-coded by correlation coefficient
#          (R); significance stars are overlaid (*** p<0.001, ** p<0.01, * p<0.05).
#
# Input: Tab-delimited file with animals as rows, immune measurements as columns.
#        Columns 2 and 6-14 are selected for correlation (adjust indices if needed).
# ------------------------------------------------------------------------------

# Load supplementary heatmap data
# Replace the path below with the path to your local copy of the data file
df_heatmap <- read.table("your_heatmap_data_file.txt",
                         sep         = "\t",
                         header      = TRUE,
                         check.names = FALSE)

# Compute pairwise Spearman correlations and p-values using Hmisc::rcorr()
# Columns 2 and 6-14 correspond to hSBA_Titre and all immune variables
cor_result <- rcorr(as.matrix(df_heatmap[, c(2, 6:14)]))

# Reshape from matrix to long-format data frame using all unique variable pairs
nm <- rownames(cor_result$r)
m  <- t(combn(nm, 2))
d  <- cbind(data.frame(m),
            R = cor_result$r[m],
            P = cor_result$P[m])

# Add significance annotations to correlation labels
d$label <- round(d$R, 2)
d$label[d$P < 0.001]                <- paste0(d$label[d$P < 0.001], "\n***")
d$label[d$P < 0.01 & d$P >= 0.001] <- paste0(d$label[d$P < 0.01 & d$P >= 0.001], "\n**")
d$label[d$P < 0.05 & d$P >= 0.01]  <- paste0(d$label[d$P < 0.05 & d$P >= 0.01], "\n*")

# Set factor levels to maintain variable order on both axes
d$X1 <- factor(d$X1, nm)
d$X2 <- factor(d$X2, rev(nm))  # Reversed for upper-triangle display

# Close any open graphics devices before plotting
graphics.off()

fig_supp_heatmap <- ggplot(d, aes(X1, X2, fill = R, label = label)) +

  # Tiles filled by Spearman R; white grid lines separate tiles
  geom_tile(color = "white") +

  # Diverging color scale: blue = negative correlation, brown = positive
  scale_fill_gradient2(low      = "skyblue4",
                       high     = "brown",
                       mid      = "white",
                       midpoint = 0,
                       limits   = c(-1, 1)) +

  # White text on strongly positive tiles (R >= 0.5) for legibility
  geom_text(color = ifelse(d$R < 0.5, "black", "white")) +

  theme_bw() +
  coord_equal() +

  theme(
    text             = element_text(size = 15, face = "bold"),
    axis.text.x      = element_text(angle = 90, hjust = 0, vjust = 0.05),
    legend.position  = "right",
    legend.title     = element_blank(),
    axis.title.x     = element_text(color = "dimgray", size = 14, face = "bold",
                                    margin = margin(t = 20, r = 20, b = 0, l = 0)),
    axis.title.y     = element_text(color = "dimgray", size = 14, face = "bold")
  ) +

  xlab("") +
  ylab("")

# Export Supplementary heatmap
# Square dimensions for the symmetric heatmap; double-column for label readability
save_figure(fig_supp_heatmap, filename = "Supplementary_correlation_heatmap",
            width_in = W_DOUBLE, height_in = 7)


# ==============================================================================
# EXPORT SUMMARY
# ==============================================================================
# After running this script, the figures/ folder will contain:
#
#   figures/Figure5a_forest_plot.pdf            — Figure 5a (vector)
#   figures/Figure5a_forest_plot.tiff           — Figure 5a (600 DPI raster)
#   figures/Figure5c_ROC_curves.pdf             — Figure 5c (vector)
#   figures/Figure5c_ROC_curves.tiff            — Figure 5c (600 DPI raster)
#   figures/Figure5e_per_mouse_clearance.pdf    — Figure 5e (vector)
#   figures/Figure5e_per_mouse_clearance.tiff   — Figure 5e (600 DPI raster)
#   figures/Supplementary_correlation_heatmap.pdf  — Supp. heatmap (vector)
#   figures/Supplementary_correlation_heatmap.tiff — Supp. heatmap (600 DPI)
#
# Panels 5b and 5d were generated in GraphPad Prism (version 10) and should
# be exported from Prism directly as PDF or TIFF at 600 DPI.
#
# For Nature Communications submission, submit PDF files as the primary format.
# TIFF files are provided as backup if the editorial system requires raster input.
# ==============================================================================


# ==============================================================================
# SESSION INFORMATION
# Paste the output of sessionInfo() as a comment block here before depositing
# to GitHub/Zenodo, so reviewers can reproduce the exact package environment.
# ==============================================================================
sessionInfo()
