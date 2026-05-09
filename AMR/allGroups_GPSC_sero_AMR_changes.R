library(openxlsx)
library(tidyverse)
library(readxl)
library(ggplot2)
library(writexl)
library(patchwork)   # 🔴 REQUIRED for combining plots

# -----------------------------
# 1️⃣ Read data & define period
# -----------------------------
data <- read_excel("AMR_prepost.xlsx") %>%
  mutate(
    Period = recode(
      vaccine_period,
      "pre-PCV10"  = "pre",
      "post-PCV10" = "post"
    ),
    Period = factor(Period, levels = c("pre", "post"))
  )

# -----------------------------
# 2️⃣ Pivot to long format
# -----------------------------
amr_long <- data %>%
  pivot_longer(
    cols = ends_with("_Res"),
    names_to  = "Antibiotic",
    values_to = "Result"
  ) %>%
  mutate(
    Antibiotic = str_remove(Antibiotic, "_Res.*"),
    Antibiotic = str_trim(Antibiotic),
    Resistant  = if_else(Result == "R", 1, 0)
  )

# -----------------------------
# 3️⃣ Define analysis groups
# -----------------------------
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

# -----------------------------
# 4️⃣ Summarise counts
# -----------------------------
amr_summary <- amr_long %>%
  group_by(Group, Antibiotic, Period) %>%
  summarise(
    Resistant = sum(Resistant),
    Total     = n(),
    .groups   = "drop"
  )

# -----------------------------
# 5️⃣ Fisher test (raw OR)
# -----------------------------
stats_full <- amr_summary %>%
  group_by(Group, Antibiotic) %>%
  filter(n_distinct(Period) == 2) %>%
  summarise(
    pre_R   = Resistant[Period == "pre"],
    pre_S   = Total[Period == "pre"]  - Resistant[Period == "pre"],
    post_R  = Resistant[Period == "post"],
    post_S  = Total[Period == "post"] - Resistant[Period == "post"],
    fisher  = list(
      fisher.test(matrix(
        c(pre_R, pre_S, post_R, post_S),
        nrow = 2,
        byrow = TRUE
      ))
    ),
    .groups = "drop"
  ) %>%
  mutate(
    OR_raw  = map_dbl(fisher, ~ .x$estimate),
    CI_low_raw  = map_dbl(fisher, ~ .x$conf.int[1]),
    CI_high_raw = map_dbl(fisher, ~ .x$conf.int[2]),
    p_value = map_dbl(fisher, ~ .x$p.value)
  ) %>%
  select(-fisher)

# -----------------------------
# 6️⃣ Convert to INCREASED RESISTANCE OR
# -----------------------------
stats_full <- stats_full %>%
  mutate(
    # proportions
    prop_pre  = pre_R  / (pre_R  + pre_S),
    prop_post = post_R / (post_R + post_S),
    diff      = prop_post - prop_pre,
    
    # BH correction
    p_adj = p.adjust(p_value, method = "BH"),
    significant = p_adj < 0.05,
    
    # percentages
    prop_pre_pct  = round(prop_pre  * 100, 1),
    prop_post_pct = round(prop_post * 100, 1),
    diff_pct      = round(diff * 100, 1),
    
    # FLIPPED OR (↑ resistance)
    OR_increase = 1 / OR_raw,
    CI_inc_low  = 1 / CI_high_raw,
    CI_inc_high = 1 / CI_low_raw,
    
    OR_round = round(OR_increase, 2),
    CI_95 = paste0(
      "(",
      round(CI_inc_low, 2), "–",
      round(CI_inc_high, 2), ")"
    ),
    
    interpretation = case_when(
      p_adj < 0.05 & diff > 0 ~ "Increased resistance post-PCV10",
      p_adj < 0.05 & diff < 0 ~ "Decreased resistance post-PCV10",
      TRUE                   ~ "No significant change"
    )
  )

# -----------------------------
# 7️⃣ Antibiotic order
# -----------------------------
antibiotic_order <- c(
  "TMP", "SMX", "COT", "CHL", "ERY",
  "CLI", "DOX", "TET", "CFX", "PEN"
)

stats_full <- stats_full %>%
  mutate(Antibiotic = factor(Antibiotic, levels = antibiotic_order)) %>%
  arrange(Group, Antibiotic)

# -----------------------------
# 8️⃣ Plot data
# -----------------------------
plot_data <- amr_summary %>%
  mutate(
    Proportion = Resistant / Total * 100,
    Antibiotic = factor(Antibiotic, levels = antibiotic_order)
  )

sig_data <- stats_full %>%
  mutate(
    star = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ NA_character_
    )
  ) %>%
  filter(!is.na(star)) %>%
  left_join(
    plot_data %>%
      group_by(Group, Antibiotic) %>%
      summarise(y_pos = max(Proportion) + 5, .groups = "drop"),
    by = c("Group", "Antibiotic")
  )

# -----------------------------
# 9️⃣ Plot
# -----------------------------
bar_plot <- ggplot(
  plot_data,
  aes(x = Antibiotic, y = Proportion, fill = Period)
) +
  geom_bar(stat = "identity", position = position_dodge(0.8)) +
  facet_wrap(~ Group, ncol = 1) +
  geom_text(
    data = sig_data,
    aes(x = Antibiotic, y = y_pos, label = star),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold"
  ) +
  scale_fill_manual(values = c("pre" = "grey70", "post" = "grey30")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    strip.text  = element_text(face = "bold"),
    legend.position = "top"
  ) +
  labs(
    title = "Antibiotic Resistance Pre- and Post-PCV10",
    x = "Antibiotic",
    y = "Proportion Resistant (%)",
    fill = "Period"
  )

# -----------------------------
# 🔟 Export Excel (INTUITIVE)
# -----------------------------
write_xlsx(
  list(
    "AMR_Increased_Resistance" = stats_full %>%
      select(
        Group, Antibiotic,
        prop_pre_pct, prop_post_pct, diff_pct,
        OR_round, CI_95,
        p_value, p_adj,
        interpretation
      )
  ),
  "AMR_prepost_INCREASED_resistance.xlsx"
)

forest_data <- stats_full %>%
  filter(!is.na(OR_round)) %>%
  mutate(
    Antibiotic = factor(Antibiotic, levels = antibiotic_order)
  )
forest_data <- forest_data %>%
  mutate(
    OR_plot = ifelse(OR_round < 1, 1 / OR_round, OR_round),
    CI_low_plot = ifelse(OR_round < 1, 1 / CI_inc_high, CI_inc_low),
    CI_high_plot = ifelse(OR_round < 1, 1 / CI_inc_low, CI_inc_high),
    
    # stars for significance
    star = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ NA_character_
    )
  )

forest_plot <- ggplot(forest_data,
                      aes(x = Antibiotic, y = OR_plot)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_errorbar(aes(ymin = CI_low_plot, ymax = CI_high_plot), width = 0.2) +
  geom_point(size = 3) +
  geom_text(aes(label = star, y = CI_high_plot * 1.15),
            na.rm = TRUE, size = 4, fontface = "bold") +
  scale_y_log10(
    breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 10, 50),
    limits = c(0.05, 50)
  ) +
  facet_wrap(~ Group, ncol = 1) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold")
  ) +
  labs(
    x = "Antibiotic",
    y = "Odds Ratio (↑ resistance post-PCV10, log scale)",
    title = "Change in Antibiotic Resistance Post-PCV10"
  )


forest_plot
bar_plot


library(tidyverse)
library(readxl)
library(writexl)
library(purrr)

# -----------------------------
# 1️⃣ Load data
# -----------------------------
amr_long <- read_excel("AMR_prepost.xlsx") %>%
  pivot_longer(
    cols = ends_with("_Res"),       # all antibiotic columns
    names_to = "Antibiotic",
    values_to = "Result"
  ) %>%
  mutate(
    Antibiotic = str_remove(Antibiotic, "_Res$"),  # clean names
    Resistant  = if_else(Result == "R", 1, 0),
    Period     = vaccine_period,
    Group      = disease_carriage
  )

# -----------------------------
# 2️⃣ Identify antibiotics with significant increase
# -----------------------------
sig_antibiotics <- stats_full %>%
  filter(p_adj < 0.05, diff > 0) %>%
  distinct(Group, Antibiotic)

# -----------------------------
# 3️⃣ Subset data to only significant antibiotics
# -----------------------------
gpsc_data <- amr_long %>%
  semi_join(sig_antibiotics, by = c("Group", "Antibiotic"))

# -----------------------------
# 4️⃣ Split by Group: Carriage and Disease
# -----------------------------
gpsc_data_carriage <- gpsc_data %>% filter(Group == "carriage")
gpsc_data_disease  <- gpsc_data %>% filter(Group == "disease")

# -----------------------------
# 5️⃣ Summarise counts by GPSC and Antibiotic within each Group
# -----------------------------
gpsc_summary <- function(df){
  df %>%
    group_by(GPSC, Antibiotic, Period) %>%
    summarise(
      Resistant = sum(Resistant),
      Total     = n(),
      .groups = "drop"
    )
}

gpsc_summary_carriage <- gpsc_summary(gpsc_data_carriage)
gpsc_summary_disease  <- gpsc_summary(gpsc_data_disease)

# -----------------------------
# 6️⃣ Function to run Fisher test per GPSC
# -----------------------------
run_gpsc_fisher <- function(df){
  df %>%
    group_by(GPSC, Antibiotic) %>%
    filter(n_distinct(Period) == 2) %>%
    summarise(
      pre_R  = Resistant[Period == "pre"],
      pre_S  = Total[Period == "pre"]  - Resistant[Period == "pre"],
      post_R = Resistant[Period == "post"],
      post_S = Total[Period == "post"] - Resistant[Period == "post"],
      fisher = list(fisher.test(matrix(c(pre_R, pre_S, post_R, post_S), nrow = 2))),
      .groups = "drop"
    ) %>%
    mutate(
      OR_raw      = map_dbl(fisher, ~ .x$estimate),
      CI_low_raw  = map_dbl(fisher, ~ .x$conf.int[1]),
      CI_high_raw = map_dbl(fisher, ~ .x$conf.int[2]),
      p_value     = map_dbl(fisher, ~ .x$p.value),
      OR_increase = ifelse(OR_raw < 1, 1 / OR_raw, OR_raw),
      CI_low_inc  = ifelse(OR_raw < 1, 1 / CI_high_raw, CI_low_raw),
      CI_high_inc = ifelse(OR_raw < 1, 1 / CI_low_raw, CI_high_raw),
      p_adj       = p.adjust(p_value, method = "BH"),
      significant = p_adj < 0.05,
      OR_round    = round(OR_increase, 2),
      CI_95       = paste0("(", round(CI_low_inc, 2), "–", round(CI_high_inc, 2), ")")
    ) %>%
    select(GPSC, Antibiotic, pre_R, pre_S, post_R, post_S, OR_round, CI_95, p_value, p_adj, significant)
}

# -----------------------------
# 7️⃣ Run Fisher test for each Group
# -----------------------------
gpsc_stats_carriage <- run_gpsc_fisher(gpsc_summary_carriage)
gpsc_stats_disease  <- run_gpsc_fisher(gpsc_summary_disease)

# -----------------------------
# 8️⃣ Keep only significant GPSCs
# -----------------------------
gpsc_sig_carriage <- gpsc_stats_carriage %>% filter(significant)
gpsc_sig_disease  <- gpsc_stats_disease  %>% filter(significant)

# -----------------------------
# 9️⃣ Export results to Excel
# -----------------------------
write_xlsx(
  list(
    "Carriage_GPSC_sig_increase" = gpsc_sig_carriage,
    "Disease_GPSC_sig_increase"  = gpsc_sig_disease
  ),
  "GPSC_sig_increase_disease_carriage.xlsx"
)



# -----------------------------
# 0️⃣ Libraries already loaded
# -----------------------------
# library(tidyverse)
# library(readxl)
# library(writexl)
# library(purrr)

# -----------------------------
# 1️⃣ Filter data for only significant GPSCs
# -----------------------------
gpsc_sig_data_carriage <- gpsc_data_carriage %>%
  semi_join(gpsc_sig_carriage, by = c("GPSC", "Antibiotic"))

gpsc_sig_data_disease <- gpsc_data_disease %>%
  semi_join(gpsc_sig_disease, by = c("GPSC", "Antibiotic"))

# -----------------------------
# 2️⃣ Summarise counts by GPSC, Antibiotic, Serotype, Period
# -----------------------------
serotype_summary <- function(df){
  df %>%
    group_by(GPSC, Antibiotic, Serotype_merged, Period) %>%
    summarise(
      Resistant = sum(Resistant),
      Total     = n(),
      .groups = "drop"
    )
}

serotype_summary_carriage <- serotype_summary(gpsc_sig_data_carriage)
serotype_summary_disease  <- serotype_summary(gpsc_sig_data_disease)

# -----------------------------
# 3️⃣ Run Fisher test per GPSC × Antibiotic × Serotype
# -----------------------------
run_serotype_fisher <- function(df){
  df %>%
    group_by(GPSC, Antibiotic, Serotype_merged) %>%
    filter(n_distinct(Period) == 2) %>%  # only serotypes with pre+post
    summarise(
      pre_R  = Resistant[Period == "pre"],
      pre_S  = Total[Period == "pre"] - Resistant[Period == "pre"],
      post_R = Resistant[Period == "post"],
      post_S = Total[Period == "post"] - Resistant[Period == "post"],
      fisher = list(fisher.test(matrix(c(pre_R, pre_S, post_R, post_S), nrow = 2))),
      .groups = "drop"
    ) %>%
    mutate(
      OR_raw      = map_dbl(fisher, ~ .x$estimate),
      CI_low_raw  = map_dbl(fisher, ~ .x$conf.int[1]),
      CI_high_raw = map_dbl(fisher, ~ .x$conf.int[2]),
      p_value     = map_dbl(fisher, ~ .x$p.value),
      
      # Flip OR if <1 to represent ↑ resistance
      OR_increase = ifelse(OR_raw < 1, 1 / OR_raw, OR_raw),
      CI_low_inc  = ifelse(OR_raw < 1, 1 / CI_high_raw, CI_low_raw),
      CI_high_inc = ifelse(OR_raw < 1, 1 / CI_low_raw, CI_high_raw),
      
      p_adj       = p.adjust(p_value, method = "BH"),
      significant = p_adj < 0.05,
      OR_round    = round(OR_increase, 2),
      CI_95       = paste0("(", round(CI_low_inc, 2), "–", round(CI_high_inc, 2), ")")
    ) %>%
    select(GPSC, Antibiotic, Serotype_merged, pre_R, pre_S, post_R, post_S,
           OR_round, CI_95, p_value, p_adj, significant)
}

# -----------------------------
# 4️⃣ Apply function
# -----------------------------
serotype_stats_carriage <- run_serotype_fisher(serotype_summary_carriage)
#serotype_stats_disease  <- run_serotype_fisher(serotype_summary_disease)# no sig gpsvc

# -----------------------------
# 5️⃣ Keep only significant serotypes
# -----------------------------
serotype_sig_carriage <- serotype_stats_carriage %>% filter(significant)
#serotype_sig_disease  <- serotype_stats_disease  %>% filter(significant)

# -----------------------------
# 6️⃣ Export to Excel
# -----------------------------
write_xlsx( serotype_sig_carriage,"Serotype_sig_increase_in_sig_GPSCs.xlsx")


