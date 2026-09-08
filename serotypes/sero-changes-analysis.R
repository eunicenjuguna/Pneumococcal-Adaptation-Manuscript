
# ============================================================
# 6A-REFERENCED FISHER'S EXACT TEST
# + 0.5 HALDANE-ANScombe CORRECTED OR / 95% CI
# + PUBLICATION-STYLE PREVALENCE + FOREST PLOTS
#
# ANALYSES:
#   1. Invasive disease
#   2. Carriage
#
# REFERENCE SEROTYPE:
#   6A
#
# OR INTERPRETATION:
#   OR > 1  = increased relative to 6A
#   OR < 1  = decreased relative to 6A
#   OR = 1  = no difference
#
# SIGNIFICANCE:
#   BH-adjusted Fisher's exact p-value < 0.05
#
# EFFECT SIZE:
#   0.5 Haldane-Anscombe corrected OR and 95% CI
#
# IMPORTANT:
#   Fisher's exact test uses RAW counts.
#   The 0.5 correction is applied ONLY to OR and 95% CI.
# ============================================================


# ============================================================
# 1. LOAD PACKAGES
# ============================================================

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(patchwork)
library(writexl)


# ============================================================
# 2. LOAD RAW DATA
# ============================================================

data <- read_excel(
  "./KLF_population_ex2011.xlsx"
)


# ============================================================
# 3. CLEAN VARIABLES
# ============================================================

data <- data %>%
  mutate(
    Serotype_merged = trimws(
      as.character(Serotype_merged)
    ),
    
    disease_carriage = trimws(
      as.character(disease_carriage)
    ),
    
    vaccine_period = trimws(
      as.character(vaccine_period)
    )
  )


# ============================================================
# 4. CHECK RAW DATA
# ============================================================

cat("\n============================================\n")
cat("RAW DATA CHECK\n")
cat("============================================\n\n")

cat("Dataset dimensions:\n")
print(dim(data))

cat("\nDisease/carriage:\n")
print(
  table(
    data$disease_carriage,
    useNA = "always"
  )
)

cat("\nVaccine period:\n")
print(
  table(
    data$vaccine_period,
    useNA = "always"
  )
)

cat("\nSerotype 6A:\n")
print(
  table(
    data$Serotype_merged == "6A",
    useNA = "always"
  )
)


# ============================================================
# 5. DEFINE VACCINE TYPE
#
# 6A is the REFERENCE serotype.
#
# PCV10 vaccine-type serotypes:
#   1, 4, 5, 6B, 9V, 14, 18C, 19F, 23F
#
# Everything else = NVT
#
# 6A is deliberately NOT included in VT here because it is
# treated separately as the reference serotype.
# ============================================================

VT_serotypes <- c(
  "1",
  "4",
  "5",
  "6B",
  "9V",
  "14",
  "18C",
  "19F",
  "23F"
)


data <- data %>%
  mutate(
    vaccine_type = case_when(
      
      Serotype_merged %in% VT_serotypes ~ "VT",
      
      TRUE ~ "NVT"
    )
  )


data$vaccine_type <- factor(
  data$vaccine_type,
  levels = c("VT", "NVT")
)


# ============================================================
# 6. CHECK 6A COUNTS
# ============================================================

cat("\n============================================\n")
cat("6A COUNTS\n")
cat("============================================\n\n")

sixA_counts <- data %>%
  filter(
    Serotype_merged == "6A"
  ) %>%
  count(
    disease_carriage,
    vaccine_period
  )

print(sixA_counts)


# ============================================================
# 7. FISHER'S EXACT TEST FUNCTION
#
# Contingency table:
#
#                 PRE     POST
# Serotype         a        b
# 6A               c        d
#
# OR = (b*c)/(a*d)
#
# This estimates the POST vs PRE odds ratio for the serotype
# relative to the corresponding change in 6A.
# ============================================================

run_raw_fisher <- function(
    pre_serotype,
    post_serotype,
    pre_6A,
    post_6A
) {
  
  mat <- matrix(
    c(
      pre_serotype,
      post_serotype,
      pre_6A,
      post_6A
    ),
    nrow = 2,
    byrow = TRUE
  )
  
  rownames(mat) <- c(
    "Serotype",
    "6A"
  )
  
  colnames(mat) <- c(
    "Pre",
    "Post"
  )
  
  ft <- fisher.test(mat)
  
  tibble(
    fisher_OR = unname(
      ft$estimate
    ),
    
    fisher_CI_lower = ft$conf.int[1],
    
    fisher_CI_upper = ft$conf.int[2],
    
    fisher_p = ft$p.value
  )
}


# ============================================================
# 8. 0.5 HALDANE-ANScombe CORRECTION FUNCTION
#
# Adds 0.5 to ALL FOUR CELLS.
#
# Corrected OR:
#
#             b*c
# OR = -----------------
#             a*d
#
# Corrected 95% CI:
#
# log(OR) +/- 1.96 * SE
# ============================================================

calculate_0.5_corrected <- function(
    pre_serotype,
    post_serotype,
    pre_6A,
    post_6A
) {
  
  # ----------------------------------------------------------
  # Haldane-Anscombe correction
  # ----------------------------------------------------------
  
  a <- pre_serotype + 0.5
  
  b <- post_serotype + 0.5
  
  c <- pre_6A + 0.5
  
  d <- post_6A + 0.5
  
  
  # ----------------------------------------------------------
  # Corrected OR
  # ----------------------------------------------------------
  
  OR <- (
    b * c
  ) / (
    a * d
  )
  
  
  # ----------------------------------------------------------
  # Standard error
  # ----------------------------------------------------------
  
  SE <- sqrt(
    1 / a +
      1 / b +
      1 / c +
      1 / d
  )
  
  
  # ----------------------------------------------------------
  # Log OR
  # ----------------------------------------------------------
  
  log_OR <- log(OR)
  
  
  # ----------------------------------------------------------
  # 95% CI
  # ----------------------------------------------------------
  
  CI_lower <- exp(
    log_OR - 1.96 * SE
  )
  
  CI_upper <- exp(
    log_OR + 1.96 * SE
  )
  
  
  tibble(
    corrected_OR = OR,
    
    corrected_CI_lower = CI_lower,
    
    corrected_CI_upper = CI_upper
  )
}


# ============================================================
# 9. ANALYSIS FUNCTION
#
# This function performs the complete analysis for either:
#   "disease"
# or
#   "carriage"
#
# It:
#   1. Gets total serotype counts
#   2. Applies eligibility threshold
#   3. Gets pre/post counts
#   4. Extracts 6A reference counts
#   5. Runs Fisher's exact tests
#   6. Applies BH correction
#   7. Calculates 0.5-corrected OR/CI
#   8. Returns the final results table
# ============================================================

run_6A_analysis <- function(
    dataset,
    group_name,
    minimum_total
) {
  
  cat("\n\n============================================\n")
  cat(
    toupper(group_name),
    "ANALYSIS\n"
  )
  cat("============================================\n\n")
  
  
  # ----------------------------------------------------------
  # 9A. TOTAL SEROTYPE COUNTS
  # ----------------------------------------------------------
  
  totals <- dataset %>%
    filter(
      disease_carriage == group_name
    ) %>%
    count(
      Serotype_merged,
      name = "Total"
    )
  
  
  # ----------------------------------------------------------
  # 9B. IDENTIFY ELIGIBLE SEROTYPES
  #
  # 6A is always retained as the reference.
  # ----------------------------------------------------------
  
  eligible <- totals %>%
    filter(
      Total >= minimum_total |
        Serotype_merged == "6A"
    )
  
  
  cat(
    "Eligibility threshold:",
    minimum_total,
    "isolates\n\n"
  )
  
  cat(
    "Eligible",
    group_name,
    "serotypes:\n"
  )
  
  print(
    eligible,
    n = Inf
  )
  
  
  # ----------------------------------------------------------
  # 9C. PRE/POST COUNTS
  # ----------------------------------------------------------
  
  counts <- dataset %>%
    filter(
      disease_carriage == group_name,
      
      Serotype_merged %in%
        eligible$Serotype_merged,
      
      vaccine_period %in%
        c("pre", "post")
    ) %>%
    count(
      Serotype_merged,
      vaccine_period
    ) %>%
    pivot_wider(
      names_from = vaccine_period,
      
      values_from = n,
      
      values_fill = 0
    )
  
  
  # Make sure both columns exist
  if (!"pre" %in% names(counts)) {
    counts$pre <- 0
  }
  
  if (!"post" %in% names(counts)) {
    counts$post <- 0
  }
  
  
  # ----------------------------------------------------------
  # 9D. EXTRACT 6A REFERENCE COUNTS
  # ----------------------------------------------------------
  
  sixA <- counts %>%
    filter(
      Serotype_merged == "6A"
    )
  
  
  # If 6A is missing from the eligible dataset, obtain its
  # counts directly from the raw data.
  
  if (nrow(sixA) == 0) {
    
    sixA <- dataset %>%
      filter(
        disease_carriage == group_name,
        
        Serotype_merged == "6A",
        
        vaccine_period %in%
          c("pre", "post")
      ) %>%
      count(
        Serotype_merged,
        vaccine_period
      ) %>%
      pivot_wider(
        names_from = vaccine_period,
        
        values_from = n,
        
        values_fill = 0
      )
    
  }
  
  
  # Make sure both reference columns exist
  if (!"pre" %in% names(sixA)) {
    sixA$pre <- 0
  }
  
  if (!"post" %in% names(sixA)) {
    sixA$post <- 0
  }
  
  
  if (nrow(sixA) != 1) {
    stop(
      paste0(
        "6A reference counts could not be identified correctly for ",
        group_name,
        "."
      )
    )
  }
  
  
  sixA_pre <- sixA$pre
  
  sixA_post <- sixA$post
  
  
  cat("\n6A reference counts:\n")
  
  cat(
    "Pre  =",
    sixA_pre,
    "\n"
  )
  
  cat(
    "Post =",
    sixA_post,
    "\n"
  )
  
  
  # ----------------------------------------------------------
  # 9E. RUN RAW FISHER'S EXACT TEST
  # ----------------------------------------------------------
  
  results <- counts %>%
    filter(
      Serotype_merged != "6A"
    ) %>%
    rowwise() %>%
    mutate(
      
      fisher = list(
        run_raw_fisher(
          pre_serotype = pre,
          post_serotype = post,
          pre_6A = sixA_pre,
          post_6A = sixA_post
        )
      ),
      
      fisher_OR =
        fisher$fisher_OR,
      
      fisher_CI_lower =
        fisher$fisher_CI_lower,
      
      fisher_CI_upper =
        fisher$fisher_CI_upper,
      
      fisher_p =
        fisher$fisher_p
      
    ) %>%
    ungroup()
  
  
  # ----------------------------------------------------------
  # 9F. ADD TOTAL COUNTS
  # ----------------------------------------------------------
  
  results <- results %>%
    left_join(
      totals,
      by = "Serotype_merged"
    )
  
  
  # ----------------------------------------------------------
  # 9G. RENAME SEROTYPE
  # ----------------------------------------------------------
  
  results <- results %>%
    rename(
      Serotype = Serotype_merged
    )
  
  
  # ----------------------------------------------------------
  # 9H. ADD VACCINE TYPE
  # ----------------------------------------------------------
  
  vaccine_lookup <- dataset %>%
    distinct(
      Serotype_merged,
      vaccine_type
    ) %>%
    rename(
      Serotype = Serotype_merged
    )
  
  
  results <- results %>%
    left_join(
      vaccine_lookup,
      by = "Serotype"
    )
  
  
  # ----------------------------------------------------------
  # 9I. BH MULTIPLE-TESTING CORRECTION
  # ----------------------------------------------------------
  
  results <- results %>%
    mutate(
      
      p_adj = p.adjust(
        fisher_p,
        method = "BH"
      ),
      
      Significant =
        p_adj < 0.05,
      
      stars = case_when(
        
        p_adj < 0.001 ~ "***",
        
        p_adj < 0.01 ~ "**",
        
        p_adj < 0.05 ~ "*",
        
        TRUE ~ ""
      )
    )
  
  
  # ----------------------------------------------------------
  # 9J. CALCULATE 0.5-CORRECTED OR AND CI
  # ----------------------------------------------------------
  
  results <- results %>%
    rowwise() %>%
    mutate(
      
      corrected = list(
        calculate_0.5_corrected(
          pre_serotype = pre,
          post_serotype = post,
          pre_6A = sixA_pre,
          post_6A = sixA_post
        )
      ),
      
      corrected_OR =
        corrected$corrected_OR,
      
      corrected_CI_lower =
        corrected$corrected_CI_lower,
      
      corrected_CI_upper =
        corrected$corrected_CI_upper
      
    ) %>%
    ungroup()
  
  
  # ----------------------------------------------------------
  # 9K. CHECK WHETHER CORRECTED CI EXCLUDES 1
  # ----------------------------------------------------------
  
  results <- results %>%
    mutate(
      
      Corrected_CI_excludes_1 =
        corrected_CI_lower > 1 |
        corrected_CI_upper < 1
    )
  
  
  # ----------------------------------------------------------
  # 9L. ADD DIRECTION OF CHANGE
  # ----------------------------------------------------------
  
  results <- results %>%
    mutate(
      
      Direction = case_when(
        
        corrected_OR > 1 ~ "Increased",
        
        corrected_OR < 1 ~ "Decreased",
        
        TRUE ~ "No change"
      )
    )
  
  
  # ----------------------------------------------------------
  # 9M. FINAL COLUMN ORDER
  # ----------------------------------------------------------
  
  results <- results %>%
    select(
      Serotype,
      vaccine_type,
      pre,
      post,
      Total,
      
      fisher_OR,
      fisher_CI_lower,
      fisher_CI_upper,
      fisher_p,
      
      p_adj,
      
      corrected_OR,
      corrected_CI_lower,
      corrected_CI_upper,
      
      Significant,
      stars,
      
      Corrected_CI_excludes_1,
      
      Direction
    ) %>%
    rename(
      Pre = pre,
      Post = post
    ) %>%
    arrange(
      fisher_p
    )
  
  
  return(
    list(
      results = results,
      eligible = eligible,
      counts = counts,
      sixA_pre = sixA_pre,
      sixA_post = sixA_post
    )
  )
}


# ============================================================
# 10. RUN INVASIVE DISEASE ANALYSIS
#
# Eligibility:
#   >=5 isolates
# ============================================================

disease_analysis <- run_6A_analysis(
  dataset = data,
  group_name = "disease",
  minimum_total = 5
)


disease_results <- disease_analysis$results

disease_eligible <- disease_analysis$eligible

disease_counts <- disease_analysis$counts


# ============================================================
# 11. RUN CARRIAGE ANALYSIS
#
# Eligibility:
#   >=10 isolates
# ============================================================

carriage_analysis <- run_6A_analysis(
  dataset = data,
  group_name = "carriage",
  minimum_total = 10
)


carriage_results <- carriage_analysis$results

carriage_eligible <- carriage_analysis$eligible

carriage_counts <- carriage_analysis$counts


# ============================================================
# 12. PRINT FINAL RESULTS
# ============================================================

cat("\n\n============================================\n")
cat("INVASIVE DISEASE — RESULTS VS 6A\n")
cat("============================================\n\n")

print(
  disease_results,
  n = Inf
)


cat("\n\n============================================\n")
cat("CARRIAGE — RESULTS VS 6A\n")
cat("============================================\n\n")

print(
  carriage_results,
  n = Inf
)


# ============================================================
# 13. PRINT SIGNIFICANT RESULTS
# ============================================================

cat("\n\n============================================\n")
cat("SIGNIFICANT INVASIVE DISEASE RESULTS\n")
cat("============================================\n\n")

print(
  disease_results %>%
    filter(
      Significant == TRUE
    ),
  n = Inf
)


cat("\n\n============================================\n")
cat("SIGNIFICANT CARRIAGE RESULTS\n")
cat("============================================\n\n")

print(
  carriage_results %>%
    filter(
      Significant == TRUE
    ),
  n = Inf
)


# ============================================================
# 14. SAVE STATISTICAL RESULTS
# ============================================================

write_xlsx(
  disease_results,
  "disease_Fisher_vs_6A_raw_and_0.5_corrected.xlsx"
)

write_xlsx(
  carriage_results,
  "carriage_Fisher_vs_6A_raw_and_0.5_corrected.xlsx"
)

write_xlsx(
  disease_eligible,
  "eligible_disease_serotypes.xlsx"
)

write_xlsx(
  carriage_eligible,
  "eligible_carriage_serotypes.xlsx"
)


# ============================================================
# 15. FUNCTION TO PREPARE PLOT DATA
#
# This function creates:
#   - prevalence/bar data
#   - forest data
#   - star data
#
# It ensures that ALL THREE panels use exactly the same
# serotype order.
# ============================================================

prepare_plot_data <- function(
    results,
    dataset,
    group_name
) {
  
  # ----------------------------------------------------------
  # 15A. Add direction
  # ----------------------------------------------------------
  
  results_plot <- results %>%
    mutate(
      Direction = case_when(
        
        corrected_OR > 1 ~ "Increased",
        
        corrected_OR < 1 ~ "Decreased",
        
        TRUE ~ "No change"
      )
    )
  
  
  # ----------------------------------------------------------
  # 15B. Create common serotype order
  #
  # Order by corrected OR.
  # 6A is always placed at the end.
  # ----------------------------------------------------------
  
  other_serotypes <- results_plot %>%
    filter(
      Serotype != "6A"
    ) %>%
    arrange(
      corrected_OR
    ) %>%
    pull(
      Serotype
    )
  
  
  serotype_order <- c(
    other_serotypes,
    "6A"
  )
  
  
  serotype_order <- unique(
    serotype_order
  )
  
  
  cat(
    "\n\n",
    toupper(group_name),
    "SEROTYPE PLOTTING ORDER:\n",
    sep = ""
  )
  
  print(
    serotype_order
  )
  
  
  # ----------------------------------------------------------
  # 15C. Prepare bar-plot data
  # ----------------------------------------------------------
  
  counts <- results_plot %>%
    select(
      Serotype,
      vaccine_type,
      Pre,
      Post
    ) %>%
    pivot_longer(
      cols = c(
        Pre,
        Post
      ),
      names_to = "Period",
      values_to = "Count"
    ) %>%
    mutate(
      Group = paste0(
        vaccine_type,
        "_",
        Period
      )
    )
  
  
  # ----------------------------------------------------------
  # 15D. Add 6A reference counts
  # ----------------------------------------------------------
  
  ref_counts <- dataset %>%
    filter(
      Serotype_merged == "6A",
      
      disease_carriage == group_name,
      
      vaccine_period %in%
        c("pre", "post")
    ) %>%
    count(
      vaccine_period
    ) %>%
    mutate(
      
      Period = case_when(
        
        vaccine_period == "pre" ~ "Pre",
        
        vaccine_period == "post" ~ "Post"
      )
    ) %>%
    select(
      Period,
      Count = n
    ) %>%
    mutate(
      
      Serotype = "6A",
      
      # 6A is displayed using the NVT colour scheme
      vaccine_type = "NVT",
      
      Group = paste0(
        vaccine_type,
        "_",
        Period
      )
    )
  
  
  # ----------------------------------------------------------
  # 15E. Combine all serotypes + 6A
  # ----------------------------------------------------------
  
  counts_all <- bind_rows(
    counts,
    ref_counts
  )
  
  
  # ----------------------------------------------------------
  # 15F. Calculate prevalence
  # ----------------------------------------------------------
  
  totals <- counts_all %>%
    group_by(
      Period
    ) %>%
    summarise(
      Total = sum(Count),
      .groups = "drop"
    )
  
  
  counts_all <- counts_all %>%
    left_join(
      totals,
      by = "Period"
    ) %>%
    mutate(
      
      Percentage =
        (Count / Total) * 100,
      
      Serotype = factor(
        Serotype,
        levels = serotype_order
      ),
      
      Group = factor(
        Group,
        levels = c(
          "VT_Pre",
          "VT_Post",
          "NVT_Pre",
          "NVT_Post"
        )
      )
    )
  
  
  # ----------------------------------------------------------
  # 15G. Prepare forest data
  # ----------------------------------------------------------
  
  forest_data <- results_plot %>%
    filter(
      Serotype != "6A"
    ) %>%
    mutate(
      Serotype = factor(
        Serotype,
        levels = serotype_order
      )
    )
  
  
  # ----------------------------------------------------------
  # 15H. Add 6A reference row
  # ----------------------------------------------------------
  
  reference_6A <- tibble(
    
    Serotype = factor(
      "6A",
      levels = serotype_order
    ),
    
    vaccine_type = "NVT",
    
    corrected_OR = 1,
    
    corrected_CI_lower = 1,
    
    corrected_CI_upper = 1,
    
    stars = "",
    
    Significant = FALSE,
    
    Corrected_CI_excludes_1 = FALSE,
    
    Direction = "Reference"
  )
  
  
  forest_data <- forest_data %>%
    bind_rows(
      reference_6A
    ) %>%
    mutate(
      Serotype = factor(
        Serotype,
        levels = serotype_order
      )
    )
  
  
  # ----------------------------------------------------------
  # 15I. Prepare star data
  # ----------------------------------------------------------
  
  star_data <- forest_data %>%
    filter(
      stars != ""
    ) %>%
    mutate(
      Serotype = factor(
        Serotype,
        levels = serotype_order
      )
    )
  
  
  return(
    list(
      results_plot = results_plot,
      counts_all = counts_all,
      forest_data = forest_data,
      star_data = star_data,
      serotype_order = serotype_order
    )
  )
}


# ============================================================
# 16. PREPARE INVASIVE DISEASE PLOT DATA
# ============================================================

disease_plot_data <- prepare_plot_data(
  results = disease_results,
  dataset = data,
  group_name = "disease"
)


# ============================================================
# 17. PREPARE CARRIAGE PLOT DATA
# ============================================================

carriage_plot_data <- prepare_plot_data(
  results = carriage_results,
  dataset = data,
  group_name = "carriage"
)


# ============================================================
# 18. DEFINE PLOT COLOURS
# ============================================================

custom_colors <- c(
  
  "VT_Pre" = "lightgray",
  
  "VT_Post" = "gray30",
  
  "NVT_Pre" = "steelblue",
  
  "NVT_Post" = "darkblue"
)


group_labels <- c(
  
  "VT_Pre" = "VT Pre",
  
  "VT_Post" = "VT Post",
  
  "NVT_Pre" = "NVT Pre",
  
  "NVT_Post" = "NVT Post"
)


forest_colors <- c(
  
  "VT" = "darkgray",
  
  "NVT" = "steelblue"
)


# ============================================================
# 19. FUNCTION TO CREATE THE THREE-PANEL PLOT
#
# PANEL 1:
#   Prevalence
#
# PANEL 2:
#   Corrected OR forest plot
#
# PANEL 3:
#   Significance stars
# ============================================================

make_combined_plot <- function(
    plot_data,
    panel_title
) {
  
  results_plot <- plot_data$results_plot
  
  counts_all <- plot_data$counts_all
  
  forest_data <- plot_data$forest_data
  
  star_data <- plot_data$star_data
  
  serotype_order <- plot_data$serotype_order
  
  
  # ==========================================================
  # 19A. BAR PLOT
  # ==========================================================
  
  bar_plot <- ggplot(
    counts_all,
    
    aes(
      x = Percentage,
      y = Serotype,
      fill = Group
    )
  ) +
    
    geom_col(
      position = "dodge"
    ) +
    
    scale_y_discrete(
      limits = rev(
        serotype_order
      )
    ) +
    
    labs(
      x = "Percentage",
      y = "Serotype",
      title = "Prevalence"
    ) +
    
    scale_fill_manual(
      values = custom_colors,
      labels = group_labels,
      name = "Group"
    ) +
    
    theme_minimal(
      base_size = 12
    ) +
    
    theme(
      
      legend.position = "top",
      
      axis.title.x = element_text(
        size = 16,
        face = "bold",
        color = "black"
      ),
      
      axis.title.y = element_text(
        size = 16,
        face = "bold",
        color = "black"
      ),
      
      axis.text = element_text(
        size = 12,
        face = "bold",
        color = "black"
      ),
      
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        color = "black",
        size = 16
      ),
      
      panel.grid.major.y = element_line(
        linewidth = 0.3
      ),
      
      panel.grid.minor = element_blank()
    )
  
  
  # ==========================================================
  # 19B. FOREST PLOT
  # ==========================================================
  
  forest_plot <- ggplot(
    forest_data,
    
    aes(
      x = corrected_OR,
      y = Serotype
    )
  ) +
    
    geom_vline(
      xintercept = 1,
      linetype = "dashed"
    ) +
    
    # Confidence intervals
    # 6A is excluded because it is the reference
    geom_errorbarh(
      data = forest_data %>%
        filter(
          Serotype != "6A"
        ),
      
      aes(
        xmin = corrected_CI_lower,
        xmax = corrected_CI_upper
      ),
      
      height = 0.3
    ) +
    
    # Corrected OR points
    geom_point(
      aes(
        color = vaccine_type
      ),
      
      size = 3
    ) +
    
    # 6A reference label
    geom_text(
      data = forest_data %>%
        filter(
          Serotype == "6A"
        ),
      
      aes(
        x = 1,
        label = "Reference"
      ),
      
      hjust = -0.15,
      
      color = "black",
      
      size = 4,
      
      fontface = "italic"
    ) +
    
    scale_x_log10(
      
      breaks = c(
        0.1,
        1,
        10,
        100
      ),
      
      labels = c(
        "1e-01",
        "1e+00",
        "1e+01",
        "1e+02"
      )
    ) +
    
    scale_y_discrete(
      limits = rev(
        serotype_order
      )
    ) +
    
    labs(
      
      x = "Odds Ratio (log scale)",
      
      y = "Serotype",
      
      title =
        "Changes Post-PCV10 in Reference to 6A"
    ) +
    
    scale_color_manual(
      
      values = forest_colors,
      
      name = "vaccine type"
    ) +
    
    theme_minimal() +
    
    theme(
      
      legend.position = "top",
      
      axis.title.x = element_text(
        size = 16,
        face = "bold",
        color = "black"
      ),
      
      axis.title.y = element_text(
        size = 16,
        face = "bold",
        color = "black"
      ),
      
      axis.text = element_text(
        size = 14,
        face = "bold",
        color = "black"
      ),
      
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        color = "black",
        size = 16
      ),
      
      plot.margin = margin(
        5,
        0,
        5,
        5
      ),
      
      panel.grid.minor = element_blank()
    )
  
  
  # ==========================================================
  # 19C. STAR PLOT
  # ==========================================================
  
  star_plot <- ggplot(
    star_data,
    
    aes(
      x = 1,
      y = Serotype,
      label = stars
    )
  ) +
    
    geom_text(
      color = "maroon",
      size = 7,
      fontface = "bold"
    ) +
    
    scale_y_discrete(
      limits = rev(
        serotype_order
      )
    ) +
    
    theme_void() +
    
    theme(
      plot.margin = margin(
        35,
        5,
        5,
        0
      )
    )
  
  
  # ==========================================================
  # 19D. COMBINE THREE PANELS
  # ==========================================================
  
  combined_plot <- (
    
    bar_plot +
      
      forest_plot +
      
      star_plot
    
  ) +
    
    plot_layout(
      widths = c(
        1,
        1,
        0.15
      )
    )
  
  
  # ==========================================================
  # 19E. ADD PANEL TITLE
  # ==========================================================
  
  combined_plot <- combined_plot +
    
    plot_annotation(
      
      title = panel_title,
      
      theme = theme(
        plot.title = element_text(
          hjust = 0,
          face = "bold",
          size = 16
        )
      )
    )
  
  
  return(
    combined_plot
  )
}


# ============================================================
# 20. CREATE INVASIVE DISEASE PLOT
# ============================================================

invasive_plot <- make_combined_plot(
  plot_data = disease_plot_data,
  panel_title = "Invasive"
)


# ============================================================
# 21. CREATE CARRIAGE PLOT
# ============================================================

carriage_plot <- make_combined_plot(
  plot_data = carriage_plot_data,
  panel_title = "Carriage"
)


# ============================================================
# 22. DISPLAY INVASIVE PLOT
# ============================================================

print(
  invasive_plot
)


# ============================================================
# 23. DISPLAY CARRIAGE PLOT
# ============================================================

print(
  carriage_plot
)


# ============================================================
# 24. SAVE INVASIVE PLOT
# ============================================================

ggsave(
  "Invasive_6A_reference_corrected_OR_aligned.png",
  
  plot = invasive_plot,
  
  width = 12,
  
  height = 8,
  
  units = "in",
  
  dpi = 300
)


ggsave(
  "Invasive_6A_reference_corrected_OR_aligned.pdf",
  
  plot = invasive_plot,
  
  width = 12,
  
  height = 8,
  
  units = "in"
)


# ============================================================
# 25. SAVE CARRIAGE PLOT
# ============================================================

ggsave(
  "Carriage_6A_reference_corrected_OR_aligned.png",
  
  plot = carriage_plot,
  
  width = 12,
  
  height = 8,
  
  units = "in",
  
  dpi = 300
)


ggsave(
  "Carriage_6A_reference_corrected_OR_aligned.pdf",
  
  plot = carriage_plot,
  
  width = 12,
  
  height = 8,
  
  units = "in"
)


# ============================================================
# 26. FINAL SUMMARY
# ============================================================

cat("\n\n============================================\n")
cat("ANALYSIS COMPLETE\n")
cat("============================================\n\n")

cat(
  "Reference serotype: 6A\n\n"
)

cat(
  "Fisher's exact test:\n"
)

cat(
  "  - Based on raw counts\n"
)

cat(
  "  - BH-adjusted p-values\n\n"
)

cat(
  "Effect size:\n"
)

cat(
  "  - 0.5 Haldane-Anscombe corrected OR\n"
)

cat(
  "  - 95% CI calculated on log OR scale\n\n"
)

cat(
  "Eligibility:\n"
)

cat(
  "  - Invasive: >=5 isolates\n"
)

cat(
  "  - Carriage: >=10 isolates\n\n"
)

cat(
  "Output files:\n\n"
)

cat(
  "Invasive results:\n",
  "  disease_Fisher_vs_6A_raw_and_0.5_corrected.xlsx\n"
)

cat(
  "Carriage results:\n",
  "  carriage_Fisher_vs_6A_raw_and_0.5_corrected.xlsx\n\n"
)

cat(
  "Invasive plots:\n",
  "  Invasive_6A_reference_corrected_OR_aligned.png\n",
  "  Invasive_6A_reference_corrected_OR_aligned.pdf\n\n"
)

cat(
  "Carriage plots:\n",
  "  Carriage_6A_reference_corrected_OR_aligned.png\n",
  "  Carriage_6A_reference_corrected_OR_aligned.pdf\n"
)

cat(
  "\n============================================\n"
)
