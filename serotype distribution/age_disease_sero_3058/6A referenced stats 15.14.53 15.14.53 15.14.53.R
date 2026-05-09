
library(writexl)
library(dplyr)
library(tidyr)
library(tidyverse)
library(readxl)
library(stats) 
library(ggplot2)
library(forcats)
library(patchwork)

# Step 1: Load your Excel file
# Adjust the file path and sheet name accordingly

data <- read_excel("./disease_under5.xlsx") ####go to line 37

# Get reference counts for serotype 6A
ref <- data %>% filter(Serotype == "6A")
ref_pre <- ref$Pre
ref_post <- ref$Post


# Ensure the reference values are valid
if (is.na(ref_pre) || is.na(ref_post) || ref_pre < 0 || ref_post < 0) {
  stop("Invalid reference serotype counts (non-finite or negative).")
}

# Perform Fisher's test with robust error handling
results <- data %>%
  filter(Serotype != "6A") %>%
  rowwise() %>%
  mutate(
    fisher_output = list({
      # Check for valid numbers
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
  mutate(p_adj = p.adjust(fisher_p, method = "BH"))

# View all results
print(results)
#view(results)
# View significant results (adjusted p < 0.05)
significant_results <- results %>% filter(p_adj < 0.05)
print(significant_results)


write_xlsx(results, path = "./6Aunder5_disease.xlsx")


# Optional: sort by odds ratio or significance
results <- results %>%
  filter(!is.na(odds_ratio)) %>%
  mutate(Serotype = fct_reorder(Serotype, odds_ratio))

# Create forest plot
ggplot(results, aes(x = odds_ratio, y = Serotype)) +
  geom_point(aes(color = p_adj < 0.05), size = 3) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.3) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_log10() +  # log scale for odds ratios
  labs(
    x = "Odds Ratio (log scale)",
    y = "Serotype",
    title = "Forest Plot of Serotype Changes Post-Vaccine",
    subtitle = "Reference: Serotype 6A"
  ) +
  theme_minimal() +
  theme(legend.position = "top") +
  scale_color_manual(values = c("grey40", "red"), name = "Significant (FDR < 0.05)")



results <- results %>%
  filter(!is.na(odds_ratio)) %>%
  mutate(Serotype = fct_reorder(Serotype, odds_ratio))

# Add a column to mark significance
results <- results %>%
  mutate(Significant = p_adj < 0.05)

results <- results %>%
  mutate(Significant = p_adj < 0.05,
         stars = case_when(
           p_adj < 0.001 ~ "***",
           p_adj < 0.01  ~ "**",
           p_adj < 0.05  ~ "*",
           TRUE ~ ""
         ))


forest_plot <- ggplot(results, aes(x = odds_ratio, y = Serotype)) +
  geom_point(aes(color = vaccine_type), size = 3) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.3) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_log10() +
  labs(
    x = "Odds Ratio (log scale)",
    y = "Serotype",
    title = "Changes Post-PCV10 in Reference to 6A"
  ) +
  theme_minimal() +
  theme(legend.position = "top right") +
  scale_color_manual(values = c("VT" = "darkgray", "NVT" = "steelblue")) +
  geom_text(aes(x = ci_upper * 1.2, label = stars), color = "maroon", size = 5)+
  theme(
    legend.position = "top right",
    axis.title.x = element_text(size = 16, face = "bold", color = "black"),
    axis.title.y = element_text(size = 16, face = "bold", color = "black"),
    axis.text = element_text(size = 10,  face = "bold", color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold", color = "black", size = 16),
    plot.subtitle = element_text(hjust = 0.5, face = "plain", color = "black", size = 13)
  )



df <- results %>%
  drop_na(Pre, Post)

# Reshape
counts <- df %>%
  pivot_longer(cols = c("Pre", "Post"), names_to = "Period", values_to = "Count")%>%
  mutate(Group = paste0(vaccine_type, "_", Period))

# Calculate percentages
counts_percent <- counts %>%
  group_by(Period) %>%
  mutate(Total = sum(Count)) %>%
  ungroup() %>%
  mutate(Percentage = (Count / Total) * 100)

custom_colors <- c(
  "VT_Pre" = "lightgray",
  "VT_Post" = "gray30",
  "NVT_Pre" = "steelblue",
  "NVT_Post" = "darkblue"
)

bar_plot <- ggplot(counts_percent, aes(x = Percentage, y = Serotype, fill = Group)) +
  geom_col(position = "dodge") +
  labs(x = "Percentage Proportion", y = "Serotype",title = "Prevalence") +
  scale_fill_manual(values = custom_colors, name = "Key") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top"
  )+
theme(
  legend.position = "top",
  axis.title.x = element_text(size = 16, face = "bold", color = "black"),
  axis.title.y = element_text(size = 16, face = "bold", color = "black"),
  axis.text = element_text(size = 10,  face = "bold", color = "black"),
  plot.title = element_text(hjust = 0.5, face = "bold", color = "black", size = 16),
  plot.subtitle = element_text(hjust = 0.5, face = "plain", color = "black", size = 13)
)


combined_plot <- bar_plot + forest_plot + plot_layout(widths = c(1, 1))
combined_plot

