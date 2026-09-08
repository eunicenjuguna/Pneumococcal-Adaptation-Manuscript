library(openxlsx)
library(tidyverse)
library(readxl)
library(ggplot2)
library(writexl)
library(patchwork)
library(purrr)

# ============================================================
# 1. READ DATA & DEFINE PERIOD
# ============================================================

data <- read_excel("AMR_prepost.xlsx") %>%
  mutate(
    Period = recode(
      vaccine_period,
      "pre-PCV10"  = "pre",
      "post-PCV10" = "post"
    ),
    Period = factor(Period, levels = c("pre", "post"))
  )


# ============================================================
# 2. PIVOT TO LONG FORMAT
# ============================================================

amr_long <- data %>%
  pivot_longer(
    cols = ends_with("_Res"),
    names_to = "Antibiotic",
    values_to = "Result"
  ) %>%
  mutate(
    Antibiotic = str_remove(Antibiotic, "_Res.*"),
    Antibiotic = str_trim(Antibiotic),
    Resistant  = if_else(Result == "R", 1, 0)
  )


# ============================================================
# 3. DEFINE ANALYSIS GROUPS
# ============================================================

amr_long <- amr_long %>%
  mutate(
    Group = case_when(
      disease_carriage == "disease"  ~ "disease",
      disease_carriage == "carriage" ~ "carriage",
      vaccine_type == "VT"           ~ "VT",
      vaccine_type == "NVT"          ~ "nVT",
      TRUE                           ~ "overall"
    )
  )


# ============================================================
# 4. SUMMARISE COUNTS
# ============================================================

amr_summary <- amr_long %>%
  group_by(Group, Antibiotic, Period) %>%
  summarise(
    Resistant = sum(Resistant),
    Total     = n(),
    .groups   = "drop"
  )

# ============================================================
# 5. FISHER TEST + HALDANE-ANSCOMBE 0.5 CORRECTION
# ============================================================
#
# Fisher's exact test:
#   - calculated from ORIGINAL counts
#   - retained as the exact Fisher p-value
#
# Haldane-Anscombe correction:
#   - +0.5 added to ALL four cells
#   - corrected OR calculated from corrected counts
#   - corrected 95% CI calculated from corrected counts
#   - corrected p-value calculated from corrected log(OR)
#
# FOREST PLOT SIGNIFICANCE:
#   - based on the corrected p-value
#   - BH adjusted corrected p-value
#
# BAR PLOT:
#   - remains based on ORIGINAL observed proportions
# ============================================================

stats_full <- amr_summary %>%
  group_by(Group, Antibiotic) %>%
  filter(n_distinct(Period) == 2) %>%
  summarise(
    
    pre_R =
      Resistant[Period == "pre"],
    
    pre_S =
      Total[Period == "pre"] -
      Resistant[Period == "pre"],
    
    post_R =
      Resistant[Period == "post"],
    
    post_S =
      Total[Period == "post"] -
      Resistant[Period == "post"],
    
    # --------------------------------------------------------
    # Fisher's exact test on ORIGINAL counts
    # --------------------------------------------------------
    
    fisher = list(
      fisher.test(
        matrix(
          c(
            post_R,
            post_S,
            pre_R,
            pre_S
          ),
          nrow = 2,
          byrow = TRUE
        )
      )
    ),
    
    .groups = "drop"
  ) %>%
  mutate(
    
    # ========================================================
    # ORIGINAL FISHER RESULTS
    # ========================================================
    
    OR_raw =
      map_dbl(
        fisher,
        ~ .x$estimate
      ),
    
    CI_low_raw =
      map_dbl(
        fisher,
        ~ .x$conf.int[1]
      ),
    
    CI_high_raw =
      map_dbl(
        fisher,
        ~ .x$conf.int[2]
      ),
    
    fisher_p =
      map_dbl(
        fisher,
        ~ .x$p.value
      ),
    
    # ========================================================
    # HALDANE-ANSCOMBE 0.5 CORRECTION
    # ========================================================
    
    a_corrected =
      post_R + 0.5,
    
    b_corrected =
      post_S + 0.5,
    
    c_corrected =
      pre_R + 0.5,
    
    d_corrected =
      pre_S + 0.5,
    
    # ========================================================
    # CORRECTED ODDS RATIO
    # ========================================================
    
    OR_corrected =
      (
        a_corrected *
          d_corrected
      ) /
      (
        b_corrected *
          c_corrected
      ),
    
    # ========================================================
    # STANDARD ERROR OF CORRECTED LOG OR
    # ========================================================
    
    SE_log_OR =
      sqrt(
        1 / a_corrected +
          1 / b_corrected +
          1 / c_corrected +
          1 / d_corrected
      ),
    
    # ========================================================
    # CORRECTED 95% CI
    # ========================================================
    
    CI_low_corrected =
      exp(
        log(OR_corrected) -
          1.96 * SE_log_OR
      ),
    
    CI_high_corrected =
      exp(
        log(OR_corrected) +
          1.96 * SE_log_OR
      ),
    
    # ========================================================
    # P-VALUE BASED ON THE CORRECTED OR
    # ========================================================
    #
    # This is the p-value that corresponds to the corrected
    # OR and corrected SE.
    #
    # Therefore the significance stars in the forest plot
    # correspond to the 0.5-corrected analysis.
    # ========================================================
    
    corrected_p =
      2 *
      pnorm(
        -abs(
          log(OR_corrected) /
            SE_log_OR
        )
      ),
    
    # ========================================================
    # BH ADJUSTMENT OF CORRECTED P-VALUES
    # ========================================================
    
    p_adj_corrected =
      p.adjust(
        corrected_p,
        method = "BH"
      ),
    
    significant_corrected =
      p_adj_corrected < 0.05,
    
    # ========================================================
    # ORIGINAL OBSERVED PROPORTIONS
    # ========================================================
    
    prop_pre =
      pre_R /
      (pre_R + pre_S),
    
    prop_post =
      post_R /
      (post_R + post_S),
    
    diff =
      prop_post -
      prop_pre,
    
    prop_pre_pct =
      round(
        prop_pre * 100,
        1
      ),
    
    prop_post_pct =
      round(
        prop_post * 100,
        1
      ),
    
    diff_pct =
      round(
        diff * 100,
        1
      ),
    
    # ========================================================
    # REPORTING VALUES
    # ========================================================
    
    OR_increase =
      OR_corrected,
    
    CI_inc_low =
      CI_low_corrected,
    
    CI_inc_high =
      CI_high_corrected,
    
    OR_round =
      round(
        OR_corrected,
        2
      ),
    
    CI_95 =
      paste0(
        "(",
        round(
          CI_low_corrected,
          2
        ),
        "–",
        round(
          CI_high_corrected,
          2
        ),
        ")"
      ),
    
    interpretation =
      case_when(
        
        significant_corrected &
          diff > 0 ~
          "Increased resistance post-PCV10",
        
        significant_corrected &
          diff < 0 ~
          "Decreased resistance post-PCV10",
        
        TRUE ~
          "No significant change"
      )
  ) %>%
  select(
    -fisher
  )


# ============================================================
# 6. ANTIBIOTIC ORDER
# ============================================================

antibiotic_order <- c(
  "TMP", "SMX", "COT", "CHL", "ERY",
  "CLI", "DOX", "TET", "CFX", "PEN"
)

group_order <- c(
  "carriage",
  "disease"
)

antibiotic_display_levels <- rev(
  antibiotic_order
)

plot_font_family <- "Arial"
base_font_size <- 12
plot_title_size <- 14


# ============================================================
# 7. PREPARE STATS FOR PLOTTING
# ============================================================

stats_full <- stats_full %>%
  mutate(
    Antibiotic =
      factor(
        Antibiotic,
        levels = antibiotic_display_levels
      ),
    
    Group =
      factor(
        Group,
        levels = group_order
      )
  ) %>%
  arrange(
    Group,
    Antibiotic
  )


# ============================================================
# 8. BAR PLOT DATA
# ============================================================
#
# IMPORTANT:
# Bar heights use ORIGINAL observed resistance proportions.
# No continuity correction is applied here.
# ============================================================

plot_data <- amr_summary %>%
  mutate(
    Proportion =
      Resistant / Total * 100,
    
    Antibiotic =
      factor(
        Antibiotic,
        levels = antibiotic_display_levels
      ),
    
    Group =
      factor(
        Group,
        levels = group_order
      )
  )


# ============================================================
# 9. BAR PLOT
# ============================================================

bar_plot <- ggplot(
  plot_data,
  aes(
    x = Proportion,
    y = Antibiotic,
    fill = Period
  )
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  facet_wrap(
    ~ Group,
    ncol = 1
  ) +
  scale_fill_manual(
    values = c(
      "pre" = "grey70",
      "post" = "grey30"
    )
  ) +
  scale_y_discrete(
    limits = antibiotic_display_levels,
    drop = FALSE
  ) +
  scale_x_continuous(
    expand = expansion(
      mult = c(0, 0.12)
    )
  ) +
  theme_minimal(
    base_family = plot_font_family,
    base_size = base_font_size
  ) +
  theme(
    text = element_text(
      family = plot_font_family,
      face = "bold",
      color = "black"
    ),
    
    plot.title = element_text(
      size = plot_title_size,
      face = "bold"
    ),
    
    axis.title.x = element_text(
      size = base_font_size,
      face = "bold"
    ),
    
    axis.title.y = element_text(
      size = base_font_size,
      face = "bold"
    ),
    
    axis.text.x = element_text(
      size = base_font_size,
      face = "bold"
    ),
    
    axis.text.y = element_text(
      size = base_font_size,
      face = "bold"
    ),
    
    strip.text = element_text(
      size = base_font_size,
      face = "bold"
    ),
    
    legend.title = element_text(
      size = base_font_size,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = base_font_size,
      face = "bold"
    ),
    
    legend.position = "top"
  ) +
  labs(
    title =
      "Antibiotic Resistance Pre- and Post-PCV10",
    
    x =
      "Proportion Resistant (%)",
    
    y =
      "Antibiotic",
    
    fill =
      "Period"
  )


# ============================================================
# 10. EXPORT OVERALL AMR RESULTS
# ============================================================

write_xlsx(
  list(
    "AMR_Increased_Resistance" =
      stats_full %>%
      select(
        Group,
        Antibiotic,
        
        prop_pre_pct,
        prop_post_pct,
        diff_pct,
        
        # Corrected OR and CI
        OR_round,
        CI_95,
        
        # Raw Fisher and corrected-analysis p-values
        fisher_p,
        corrected_p,
        p_adj_corrected,
        
        interpretation
      )
  ),
  "AMR_prepost_INCREASED_resistance.xlsx"
)

# ============================================================
# 11. FOREST PLOT DATA
# ============================================================
#
# Forest plot uses:
#   - Haldane-Anscombe corrected OR
#   - Haldane-Anscombe corrected 95% CI
#   - BH-adjusted corrected p-value for significance stars
#
# ============================================================

forest_data <- stats_full %>%
  filter(
    is.finite(OR_corrected),
    is.finite(CI_low_corrected),
    is.finite(CI_high_corrected)
  ) %>%
  mutate(
    
    Antibiotic =
      factor(
        Antibiotic,
        levels = antibiotic_display_levels
      ),
    
    Group =
      factor(
        Group,
        levels = group_order
      ),
    
    OR_plot =
      OR_corrected,
    
    CI_low_plot =
      CI_low_corrected,
    
    CI_high_plot =
      CI_high_corrected,
    
    # ========================================================
    # SIGNIFICANCE STARS
    # ========================================================
    #
    # Stars are now based on the BH-adjusted
    # 0.5-corrected p-value.
    # ========================================================
    
    star =
      case_when(
        
        p_adj_corrected < 0.001 ~ "***",
        
        p_adj_corrected < 0.01 ~ "**",
        
        p_adj_corrected < 0.05 ~ "*",
        
        TRUE ~ NA_character_
      )
  )


# ============================================================
# 12. DETERMINE FOREST PLOT LIMITS
# ============================================================
#
# IMPORTANT:
# Do NOT cap the CIs.
#
# All corrected confidence intervals are finite after the
# +0.5 correction, so we allow the x-axis to extend far
# enough to show EVERY confidence interval.
# ============================================================

finite_values <- c(
  forest_data$CI_low_plot,
  forest_data$CI_high_plot,
  forest_data$OR_plot
)

finite_values <-
  finite_values[
    is.finite(finite_values)
  ]

x_min <-
  min(
    finite_values,
    na.rm = TRUE
  )

x_max <-
  max(
    finite_values,
    na.rm = TRUE
  )

# Add space on both sides of the log scale
x_min_plot <-
  x_min / 1.5

x_max_plot <-
  x_max * 1.5


# ============================================================
# 13. FOREST PLOT
# ============================================================

forest_plot <- ggplot(
  forest_data,
  aes(
    x = OR_plot,
    y = Antibiotic
  )
) +
  
  # ----------------------------------------------------------
# NO EFFECT LINE
# ----------------------------------------------------------

geom_vline(
  xintercept = 1,
  linetype = "dashed",
  colour = "grey40",
  linewidth = 0.7
) +
  
  # ----------------------------------------------------------
# CORRECTED 95% CI
# ----------------------------------------------------------

geom_errorbarh(
  aes(
    xmin = CI_low_plot,
    xmax = CI_high_plot
  ),
  height = 0.20,
  linewidth = 0.8
) +
  
  # ----------------------------------------------------------
# CORRECTED OR
# ----------------------------------------------------------

geom_point(
  size = 3
) +
  
  # ----------------------------------------------------------
# SIGNIFICANCE STARS
# ----------------------------------------------------------

geom_text(
  data =
    forest_data %>%
    filter(
      !is.na(star)
    ),
  
  aes(
    x =
      CI_high_plot * 1.12,
    
    y =
      Antibiotic,
    
    label =
      star
  ),
  
  inherit.aes = FALSE,
  
  size = 5,
  
  fontface = "bold",
  
  hjust = 0
) +
  
  # ----------------------------------------------------------
# Y AXIS
# ----------------------------------------------------------

scale_y_discrete(
  limits =
    antibiotic_display_levels,
  
  drop = FALSE
) +
  
  # ----------------------------------------------------------
# LOG ODDS RATIO AXIS
# ----------------------------------------------------------

scale_x_log10(
  
  breaks =
    c(
      0.1,
      0.25,
      0.5,
      1,
      2,
      4,
      10,
      20,
      50,
      100,
      500,
      1000
    ),
  
  limits =
    c(
      x_min_plot,
      x_max_plot
    ),
  
  expand =
    expansion(
      mult =
        c(
          0.02,
          0.02
        )
    )
) +
  
  # ----------------------------------------------------------
# FACETS
# ----------------------------------------------------------

facet_wrap(
  ~ Group,
  ncol = 1,
  scales = "free_y"
) +
  
  # ----------------------------------------------------------
# THEME
# ----------------------------------------------------------

theme_minimal(
  base_family =
    plot_font_family,
  
  base_size =
    base_font_size
) +
  
  theme(
    
    text =
      element_text(
        family =
          plot_font_family,
        
        face =
          "bold",
        
        colour =
          "black"
      ),
    
    axis.title.x =
      element_text(
        size =
          base_font_size,
        
        face =
          "bold",
        
        margin =
          margin(
            t = 10
          )
      ),
    
    axis.title.y =
      element_text(
        size =
          base_font_size,
        
        face =
          "bold"
      ),
    
    axis.text.x =
      element_text(
        size = 10,
        
        face =
          "bold",
        
        colour =
          "black"
      ),
    
    axis.text.y =
      element_text(
        size = 11,
        
        face =
          "bold",
        
        colour =
          "black"
      ),
    
    strip.text =
      element_text(
        size = 12,
        
        face =
          "bold",
        
        colour =
          "black"
      ),
    
    strip.background =
      element_blank(),
    
    panel.grid.major.y =
      element_line(
        colour =
          "grey90",
        
        linewidth =
          0.5
      ),
    
    panel.grid.major.x =
      element_line(
        colour =
          "grey90",
        
        linewidth =
          0.5
      ),
    
    panel.grid.minor.x =
      element_line(
        colour =
          "grey95",
        
        linewidth =
          0.4
      ),
    
    plot.title =
      element_blank(),
    
    plot.margin =
      margin(
        5,
        40,
        10,
        5
      )
  ) +
  
  labs(
    
    x =
      "Odds Ratio (post-PCV10 vs pre-PCV10, log scale)",
    
    y =
      "Antibiotic"
  )


# ============================================================
# 14. REMOVE BAR-PLOT TITLE
# ============================================================

bar_plot <- bar_plot +
  labs(
    title = NULL
  )


# ============================================================
# 15. COMBINE BAR + FOREST PLOT
# ============================================================

combined_plot <-
  bar_plot +
  forest_plot +
  
  plot_layout(
    widths =
      c(
        1.05,
        1.95
      ),
    
    guides =
      "collect"
  ) +
  
  plot_annotation(
    
    title =
      "Change in Antibiotic Resistance After PCV10 Introduction",
    
    theme =
      theme(
        
        plot.title =
          element_text(
            
            family =
              plot_font_family,
            
            face =
              "bold",
            
            size =
              16,
            
            hjust =
              0.5,
            
            margin =
              margin(
                b = 15
              )
          )
      )
  )


# ============================================================
# 16. DISPLAY
# ============================================================

combined_plot

# ============================================================
# 17. IDENTIFY ANTIBIOTICS WITH SIGNIFICANT INCREASE
# ============================================================

sig_antibiotics <- stats_full %>%
  filter(
    significant_corrected,
    diff > 0
  ) %>%
  distinct(
    Group,
    Antibiotic
  )


# ============================================================
# 18. SUBSET DATA TO SIGNIFICANT ANTIBIOTICS
# ============================================================

gpsc_data <- amr_long %>%
  semi_join(
    sig_antibiotics,
    by = c(
      "Group",
      "Antibiotic"
    )
  )


# ============================================================
# 19. SPLIT BY GROUP
# ============================================================

gpsc_data_carriage <-
  gpsc_data %>%
  filter(
    Group == "carriage"
  )

gpsc_data_disease <-
  gpsc_data %>%
  filter(
    Group == "disease"
  )


# ============================================================
# 20. SUMMARISE COUNTS BY GPSC AND ANTIBIOTIC
# ============================================================

gpsc_summary <- function(df) {
  
  df %>%
    group_by(
      GPSC,
      Antibiotic,
      Period
    ) %>%
    summarise(
      Resistant =
        sum(Resistant),
      
      Total =
        n(),
      
      .groups =
        "drop"
    )
}


gpsc_summary_carriage <-
  gpsc_summary(
    gpsc_data_carriage
  )

gpsc_summary_disease <-
  gpsc_summary(
    gpsc_data_disease
  )


# ============================================================
# 21. FISHER + 0.5 CORRECTION PER GPSC
# ============================================================

run_gpsc_fisher <- function(df) {
  
  df %>%
    group_by(
      GPSC,
      Antibiotic
    ) %>%
    
    filter(
      n_distinct(Period) == 2
    ) %>%
    
    summarise(
      
      pre_R =
        Resistant[
          Period == "pre"
        ],
      
      pre_S =
        Total[
          Period == "pre"
        ] -
        Resistant[
          Period == "pre"
        ],
      
      post_R =
        Resistant[
          Period == "post"
        ],
      
      post_S =
        Total[
          Period == "post"
        ] -
        Resistant[
          Period == "post"
        ],
      
      # Fisher exact test on raw counts
      fisher =
        list(
          fisher.test(
            matrix(
              c(
                post_R,
                post_S,
                pre_R,
                pre_S
              ),
              nrow = 2,
              byrow = TRUE
            )
          )
        ),
      
      .groups =
        "drop"
    ) %>%
    
    mutate(
      
      # Raw Fisher OR
      OR_raw =
        map_dbl(
          fisher,
          ~ .x$estimate
        ),
      
      # Raw Fisher exact CI
      CI_low_raw =
        map_dbl(
          fisher,
          ~ .x$conf.int[1]
        ),
      
      CI_high_raw =
        map_dbl(
          fisher,
          ~ .x$conf.int[2]
        ),
      
      # Fisher exact p-value
      p_value =
        map_dbl(
          fisher,
          ~ .x$p.value
        ),
      
      # ------------------------------------------------------
      # Haldane-Anscombe +0.5 correction
      # ------------------------------------------------------
      
      a_corrected =
        post_R + 0.5,
      
      b_corrected =
        post_S + 0.5,
      
      c_corrected =
        pre_R + 0.5,
      
      d_corrected =
        pre_S + 0.5,
      
      # Corrected post/pre OR
      OR_increase =
        (
          a_corrected *
            d_corrected
        ) /
        (
          b_corrected *
            c_corrected
        ),
      
      # Standard error
      SE_log_OR =
        sqrt(
          1 / a_corrected +
            1 / b_corrected +
            1 / c_corrected +
            1 / d_corrected
        ),
      
      # Corrected 95% CI
      CI_low_inc =
        exp(
          log(OR_increase) -
            1.96 * SE_log_OR
        ),
      
      CI_high_inc =
        exp(
          log(OR_increase) +
            1.96 * SE_log_OR
        ),
      
      # BH correction based on Fisher p-values
      p_adj =
        p.adjust(
          p_value,
          method = "BH"
        ),
      
      significant =
        p_adj < 0.05,
      
      OR_round =
        round(
          OR_increase,
          2
        ),
      
      CI_95 =
        paste0(
          "(",
          round(
            CI_low_inc,
            2
          ),
          "–",
          round(
            CI_high_inc,
            2
          ),
          ")"
        )
    ) %>%
    
    select(
      GPSC,
      Antibiotic,
      
      pre_R,
      pre_S,
      post_R,
      post_S,
      
      # Corrected estimates
      OR_round,
      CI_95,
      
      # Fisher p-values
      p_value,
      p_adj,
      significant
    )
}


# ============================================================
# 22. RUN GPSC FISHER ANALYSIS
# ============================================================

gpsc_stats_carriage <-
  run_gpsc_fisher(
    gpsc_summary_carriage
  )

gpsc_stats_disease <-
  run_gpsc_fisher(
    gpsc_summary_disease
  )


# ============================================================
# 23. KEEP ONLY SIGNIFICANT GPSCs
# ============================================================

gpsc_sig_carriage <-
  gpsc_stats_carriage %>%
  filter(
    significant
  )

gpsc_sig_disease <-
  gpsc_stats_disease %>%
  filter(
    significant
  )


# ============================================================
# 24. EXPORT SIGNIFICANT GPSC RESULTS
# ============================================================

write_xlsx(
  list(
    "Carriage_GPSC_sig_increase" =
      gpsc_sig_carriage,
    
    "Disease_GPSC_sig_increase" =
      gpsc_sig_disease
  ),
  "GPSC_sig_increase_disease_carriage.xlsx"
)


# ============================================================
# 25. EXPORT ALL GPSC RESULTS
# ============================================================

write_xlsx(
  list(
    "Carriage_GPSC_All" =
      gpsc_stats_carriage,
    
    "Disease_GPSC_All" =
      gpsc_stats_disease
  ),
  "GPSC_results_disease_carriage_FULL.xlsx"
)


# ============================================================
# 26. FILTER DATA FOR SIGNIFICANT GPSCs
# ============================================================

gpsc_sig_data_carriage <-
  gpsc_data_carriage %>%
  semi_join(
    gpsc_sig_carriage,
    by = c(
      "GPSC",
      "Antibiotic"
    )
  )

gpsc_sig_data_disease <-
  gpsc_data_disease %>%
  semi_join(
    gpsc_sig_disease,
    by = c(
      "GPSC",
      "Antibiotic"
    )
  )


# ============================================================
# 27. SUMMARISE COUNTS BY GPSC, ANTIBIOTIC,
#     SEROTYPE AND PERIOD
# ============================================================

serotype_summary <- function(df) {
  
  df %>%
    group_by(
      GPSC,
      Antibiotic,
      Serotype_merged,
      Period
    ) %>%
    summarise(
      Resistant =
        sum(Resistant),
      
      Total =
        n(),
      
      .groups =
        "drop"
    )
}


serotype_summary_carriage <-
  serotype_summary(
    gpsc_sig_data_carriage
  )

serotype_summary_disease <-
  serotype_summary(
    gpsc_sig_data_disease
  )


# ============================================================
# 28. FISHER + 0.5 CORRECTION PER
#     GPSC × ANTIBIOTIC × SEROTYPE
# ============================================================

run_serotype_fisher <- function(df) {
  
  df %>%
    group_by(
      GPSC,
      Antibiotic,
      Serotype_merged
    ) %>%
    
    filter(
      n_distinct(Period) == 2
    ) %>%
    
    summarise(
      
      pre_R =
        Resistant[
          Period == "pre"
        ],
      
      pre_S =
        Total[
          Period == "pre"
        ] -
        Resistant[
          Period == "pre"
        ],
      
      post_R =
        Resistant[
          Period == "post"
        ],
      
      post_S =
        Total[
          Period == "post"
        ] -
        Resistant[
          Period == "post"
        ],
      
      # Fisher exact test on raw counts
      fisher =
        list(
          fisher.test(
            matrix(
              c(
                post_R,
                post_S,
                pre_R,
                pre_S
              ),
              nrow = 2,
              byrow = TRUE
            )
          )
        ),
      
      .groups =
        "drop"
    ) %>%
    
    mutate(
      
      # Raw Fisher OR
      OR_raw =
        map_dbl(
          fisher,
          ~ .x$estimate
        ),
      
      # Raw Fisher CI
      CI_low_raw =
        map_dbl(
          fisher,
          ~ .x$conf.int[1]
        ),
      
      CI_high_raw =
        map_dbl(
          fisher,
          ~ .x$conf.int[2]
        ),
      
      # Fisher exact p-value
      p_value =
        map_dbl(
          fisher,
          ~ .x$p.value
        ),
      
      # ------------------------------------------------------
      # Haldane-Anscombe +0.5 correction
      # ------------------------------------------------------
      
      a_corrected =
        post_R + 0.5,
      
      b_corrected =
        post_S + 0.5,
      
      c_corrected =
        pre_R + 0.5,
      
      d_corrected =
        pre_S + 0.5,
      
      # Corrected post/pre OR
      OR_increase =
        (
          a_corrected *
            d_corrected
        ) /
        (
          b_corrected *
            c_corrected
        ),
      
      # Standard error
      SE_log_OR =
        sqrt(
          1 / a_corrected +
            1 / b_corrected +
            1 / c_corrected +
            1 / d_corrected
        ),
      
      # Corrected 95% CI
      CI_low_inc =
        exp(
          log(OR_increase) -
            1.96 * SE_log_OR
        ),
      
      CI_high_inc =
        exp(
          log(OR_increase) +
            1.96 * SE_log_OR
        ),
      
      # BH correction based on Fisher p-values
      p_adj =
        p.adjust(
          p_value,
          method = "BH"
        ),
      
      significant =
        p_adj < 0.05,
      
      OR_round =
        round(
          OR_increase,
          2
        ),
      
      CI_95 =
        paste0(
          "(",
          round(
            CI_low_inc,
            2
          ),
          "–",
          round(
            CI_high_inc,
            2
          ),
          ")"
        )
    ) %>%
    
    select(
      GPSC,
      Antibiotic,
      Serotype_merged,
      
      pre_R,
      pre_S,
      post_R,
      post_S,
      
      # Corrected estimates
      OR_round,
      CI_95,
      
      # Fisher p-values
      p_value,
      p_adj,
      significant
    )
}


# ============================================================
# 29. APPLY SEROTYPE FISHER ANALYSIS
# ============================================================

serotype_stats_carriage <-
  run_serotype_fisher(
    serotype_summary_carriage
  )

# No significant GPSCs currently identified for disease,
# so this remains commented as in your original workflow.

# serotype_stats_disease <-
#   run_serotype_fisher(
#     serotype_summary_disease
#   )


# ============================================================
# 30. KEEP ALL SEROTYPES FROM SIGNIFICANT GPSCs
# ============================================================

serotype_results_carriage <-
  serotype_stats_carriage

# serotype_results_disease <-
#   serotype_stats_disease


# ============================================================
# 31. EXPORT ALL SEROTYPE RESULTS FOR SIGNIFICANT GPSCs
# ============================================================

write_xlsx(
  serotype_results_carriage,
  "Serotype_results_in_sig_GPSCs.xlsx"
)
