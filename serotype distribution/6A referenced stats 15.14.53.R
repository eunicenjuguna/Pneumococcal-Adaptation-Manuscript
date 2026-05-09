library(writexl)
library(dplyr)
library(tidyr)
library(tidyverse)
library(readxl)
library(stats)
library(ggplot2)
library(forcats)
library(patchwork)

# Load your data
data <- read_excel("./pre_post_disease_serotypes.xlsx")

# Extract reference counts for 6A (needed for Fisher test)
ref <- data %>% filter(Serotype == "6A")
ref_pre <- ref$Pre
ref_post <- ref$Post

if (is.na(ref_pre) || is.na(ref_post) || ref_pre < 0 || ref_post < 0) {
  stop("Invalid reference serotype counts (non-finite or negative).")
}

# Perform Fisher's exact test comparing each serotype to 6A
results <- data %>%
  filter(Serotype != "6A") %>%
  rowwise() %>%
  mutate(
    fisher_output = list({
      if (any(is.na(c(Pre, Post))) || any(c(Pre, Post) < 0)) {
        NULL
      } else {
        fisher.test(matrix(c(Post, Pre, ref_post, ref_pre), nrow = 2))
      }
    }),
    fisher_p = if (!is.null(fisher_output)) fisher_output$p.value else NA_real_,
    odds_ratio = if (!is.null(fisher_output)) fisher_output$estimate else NA_real_,
    ci_lower = if (!is.null(fisher_output)) fisher_output$conf.int[1] else NA_real_,
    ci_upper = if (!is.null(fisher_output)) fisher_output$conf.int[2] else NA_real_
  ) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(fisher_p, method = "BH")) %>%
  filter(!is.na(odds_ratio)) %>%
  mutate(Significant = p_adj < 0.05,
         stars = case_when(
           p_adj < 0.001 ~ "***",
           p_adj < 0.01  ~ "**",
           p_adj < 0.05  ~ "*",
           TRUE ~ ""
         ))

# Create order of serotypes with 6A first, then the rest ordered by odds ratio
serotype_order <- c("6A", levels(fct_reorder(results$Serotype, results$odds_ratio)))

# Prepare counts data for bar plot including 6A
# 1) Counts for other serotypes from results + pivot
counts <- results %>%
  select(Serotype, vaccine_type, Pre, Post) %>%
  pivot_longer(cols = c("Pre", "Post"), names_to = "Period", values_to = "Count") %>%
  mutate(Group = paste0(vaccine_type, "_", Period))

# 2) Counts for 6A from original data (assigning to NVT groups to match colors)
ref_counts <- data %>%
  filter(Serotype == "6A") %>%
  pivot_longer(cols = c("Pre", "Post"), names_to = "Period", values_to = "Count") %>%
  mutate(vaccine_type = "NVT",    # Assign NVT to 6A for colors
         Group = paste0(vaccine_type, "_", Period))

# Combine counts for bar plot
counts_all <- bind_rows(counts, ref_counts)

# Calculate total counts per Period to get percentages
totals <- counts_all %>%
  group_by(Period) %>%
  summarize(Total = sum(Count), .groups = "drop")

counts_all <- counts_all %>%
  left_join(totals, by = "Period") %>%
  mutate(Percentage = (Count / Total) * 100,
         Serotype = factor(Serotype, levels = serotype_order),
         Group = factor(Group, levels = c("VT_Pre", "VT_Post", "NVT_Pre", "NVT_Post")))

# Define colors for groups
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

# Create bar plot
bar_plot <- ggplot(counts_all, aes(x = Percentage, y = Serotype, fill = Group)) +
  geom_col(position = "dodge") +
  scale_y_discrete(limits = serotype_order) +
  labs(x = "Percentage", y = "Serotype", title = "Prevalence") +
  scale_fill_manual(values = custom_colors, labels = group_labels, name = "Group") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    axis.title.x = element_text(size = 16, face = "bold", color = "black"),
    axis.title.y = element_text(size = 16, face = "bold", color = "black"),
    axis.text = element_text(size = 12, face = "bold", color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold", color = "black", size = 16)
  )

# Prepare results data for forest plot
results <- results %>%
  mutate(Serotype = factor(Serotype, levels = serotype_order))

# Create forest plot
forest_plot <- ggplot(results, aes(x = odds_ratio, y = Serotype)) +
  geom_point(aes(color = vaccine_type), size = 3) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.3) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_log10(
    breaks = c(0.1, 1, 10, 100),
    labels = c("1e-01", "1e+00", "1e+01", "1e+02")
  ) +
  scale_y_discrete(limits = serotype_order) +
  labs(
    x = "Odds Ratio (log scale)",
    y = "Serotype",
    title = "Changes Post-PCV10 in Reference to 6A"
  ) +
  theme_minimal() +
  theme(
    legend.position = "top",
    axis.title.x = element_text(size = 16, face = "bold", color = "black"),
    axis.title.y = element_text(size = 16, face = "bold", color = "black"),
    axis.text = element_text(size = 14, face = "bold", color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold", color = "black", size = 16),
    plot.margin = margin(5, 0, 5, 5)
  ) +
  scale_color_manual(values = c("VT" = "darkgray", "NVT" = "steelblue"), name = "vaccine type")

star_plot <- ggplot(
  results %>%
    filter(stars != "") %>%
    mutate(Serotype = factor(Serotype, levels = serotype_order)),
  aes(x = 1, y = Serotype, label = stars)
) +
  geom_text(color = "maroon", size = 7, fontface = "bold") +
  scale_y_discrete(limits = serotype_order) +
  theme_void() +
  theme(plot.margin = margin(35, 5, 5, 0))

# Combine plots side by side with patchwork
combined_plot <- bar_plot + forest_plot + star_plot + plot_layout(widths = c(1, 1, 0.15))
print(combined_plot)
ggsave("combined_plot.svg", plot = combined_plot, width = 12, height = 8, units = "in", dpi = 300)


# Filter significant serotypes
significant_serotypes <- results %>%
  filter(Significant == TRUE) %>%
  select(Serotype, vaccine_type, Pre, Post, odds_ratio, ci_lower, ci_upper, fisher_p, p_adj, stars)

# Write to Excel
write_xlsx(significant_serotypes, path = "significant_disease_serotypes.xlsx")


# Filter significant serotypes before p-value adjustment
significant_before_correction <- results %>%
  filter(fisher_p < 0.05) %>%
  select(Serotype, vaccine_type, Pre, Post, odds_ratio, ci_lower, ci_upper, fisher_p, stars)

# Write to Excel
write_xlsx(significant_before_correction, path = "significant_disease_serotypes_before_correction.xlsx")





#carriage

# Load data
# ================================
data <- read_excel("./pre_post_carriage_serotypes.xlsx")

# ================================
# Extract reference counts for 6A
# ================================
ref <- data %>% filter(Serotype == "6A")
ref_pre <- ref$Pre
ref_post <- ref$Post

if (is.na(ref_pre) || is.na(ref_post) || ref_pre < 0 || ref_post < 0) {
  stop("Invalid reference serotype counts (non-finite or negative).")
}

# ================================
# Perform Fisher's exact test
# ================================
results <- data %>%
  filter(Serotype != "6A") %>%
  rowwise() %>%
  mutate(
    fisher_output = list({
      if (any(is.na(c(Pre, Post))) || any(c(Pre, Post) < 0)) {
        NULL
      } else {
        # Add a small pseudocount to avoid infinite OR if desired
        # mat <- matrix(c(Post + 0.1, Pre + 0.1, ref_post + 0.1, ref_pre + 0.1), nrow = 2)
        mat <- matrix(c(Post, Pre, ref_post, ref_pre), nrow = 2)
        fisher.test(mat)
      }
    }),
    fisher_p = if (!is.null(fisher_output)) fisher_output$p.value else NA_real_,
    odds_ratio = if (!is.null(fisher_output)) fisher_output$estimate else NA_real_,
    ci_lower = if (!is.null(fisher_output)) fisher_output$conf.int[1] else NA_real_,
    ci_upper = if (!is.null(fisher_output)) fisher_output$conf.int[2] else NA_real_
  ) %>%
  ungroup() %>%
  mutate(
    p_adj = p.adjust(fisher_p, method = "BH"),
    Significant = p_adj < 0.05,
    stars = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE ~ ""
    )
  )

# ================================
# Prepare counts data for bar plot
# ================================
counts <- results %>%
  select(Serotype, vaccine_type, Pre, Post) %>%
  pivot_longer(cols = c("Pre", "Post"), names_to = "Period", values_to = "Count") %>%
  mutate(Group = paste0(vaccine_type, "_", Period))

# 6A counts for reference
ref_counts <- data %>%
  filter(Serotype == "6A") %>%
  pivot_longer(cols = c("Pre", "Post"), names_to = "Period", values_to = "Count") %>%
  mutate(vaccine_type = "NVT",
         Group = paste0(vaccine_type, "_", Period))

# Combine counts
counts_all <- bind_rows(counts, ref_counts)

# Calculate percentages
totals <- counts_all %>%
  group_by(Period) %>%
  summarize(Total = sum(Count), .groups = "drop")

counts_all <- counts_all %>%
  left_join(totals, by = "Period") %>%
  mutate(Percentage = (Count / Total) * 100)

# ================================
# Filter only significant serotypes
# ================================
results_sig <- results %>% filter(Significant == TRUE)
serotype_order_sig <- c("6A", levels(fct_reorder(results_sig$Serotype, results_sig$odds_ratio)))

counts_all_sig <- counts_all %>%
  filter(Serotype %in% results_sig$Serotype | Serotype == "6A") %>%
  mutate(Serotype = factor(Serotype, levels = serotype_order_sig))

# ================================
# Define colors
# ================================
custom_colors <- c(
  "VT_Pre" = "lightgray",
  "VT_Post" = "gray30",
  "NVT_Pre" = "steelblue",
  "NVT_Post" = "darkblue"
)

# ================================
# Create bar plot
# ================================
bar_plot <- ggplot(counts_all_sig, aes(x = Percentage, y = Serotype, fill = Group)) +
  geom_col(position = "dodge") +
  scale_y_discrete(limits = serotype_order_sig) +
  labs(x = "Percentage", y = "Serotype", title = "Prevalence of Significant Serotypes") +
  scale_fill_manual(values = custom_colors, labels = group_labels, name = "Group") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    axis.title.x = element_text(size = 16, face = "bold", color = "black"),
    axis.title.y = element_text(size = 16, face = "bold", color = "black"),
    axis.text = element_text(size = 14, face = "bold", color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold", color = "black", size = 16)
  )

# ================================
# Create forest plot
# ================================
forest_plot <- ggplot(results_sig %>% mutate(Serotype = factor(Serotype, levels = serotype_order_sig)),
                      aes(x = odds_ratio, y = Serotype)) +
  geom_point(aes(color = vaccine_type), size = 3) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.3) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_log10(
    breaks = c(0.1, 1, 10, 100),
    labels = c("1e-01", "1e+00", "1e+01", "1e+02")
  ) +
  scale_y_discrete(limits = serotype_order_sig) +
  labs(x = "Odds Ratio (log scale)", y = "Serotype", title = "Significant Changes Post-PCV10 vs 6A") +
  theme_minimal() +
  theme(
    legend.position = "top",
    axis.title.x = element_text(size = 16, face = "bold", color = "black"),
    axis.title.y = element_text(size = 16, face = "bold", color = "black"),
    axis.text = element_text(size = 14, face = "bold", color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold", color = "black", size = 16),
    plot.margin = margin(5, 0, 5, 5)
  ) +
  scale_color_manual(values = c("VT" = "darkgray", "NVT" = "steelblue"), name = "vaccine type")

star_plot <- ggplot(
  results_sig %>%
    filter(stars != "") %>%
    mutate(Serotype = factor(Serotype, levels = serotype_order_sig)),
  aes(x = 1, y = Serotype, label = stars)
) +
  geom_text(color = "maroon", size = 7, fontface = "bold") +
  scale_y_discrete(limits = serotype_order_sig) +
  theme_void() +
  theme(plot.margin = margin(35, 5, 5, 0))

# ================================
# Combine plots
# ================================
combined_plot <- bar_plot + forest_plot + star_plot + plot_layout(widths = c(1, 1, 0.15))
print(combined_plot)
ggsave("significant_serotypes_plot.svg", plot = combined_plot, width = 12, height = 8, units = "in", dpi = 300)

# ================================
# Export significant serotypes to Excel
# ================================
significant_serotypes <- results_sig %>%
  select(Serotype, vaccine_type, Pre, Post, odds_ratio, ci_lower, ci_upper, fisher_p, p_adj, stars)

write_xlsx(significant_serotypes, path = "significant_carriage_serotypes.xlsx")
