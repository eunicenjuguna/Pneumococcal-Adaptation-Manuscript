## Disease GPSC analysis and figure generation

## Load packages
library(tidyverse)
library(ggpattern)
library(patchwork)
library(ggnewscale)
library(readxl)
library(writexl)

## Differential disease analysis relative to GPSC 53
comparison_data <- read_excel("./pre_post_disease_GPSCs.xlsx")

reference_counts <- comparison_data %>%
  filter(GPSC == "53")

ref_pre <- reference_counts$Pre
ref_post <- reference_counts$Post

if (is.na(ref_pre) || is.na(ref_post) || ref_pre < 0 || ref_post < 0) {
  stop("Invalid reference GPSC counts detected for GPSC 53.")
}

results <- comparison_data %>%
  filter(GPSC != "53") %>%
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
  mutate(p_adj = p.adjust(fisher_p, method = "BH"))

print(results)

significant_results <- results %>%
  filter(p_adj < 0.05)

print(significant_results)

write_xlsx(results, path = "./53reference_disease.xlsx")

## Read the curated disease results table used for plotting
results <- results %>%
  filter(!is.na(odds_ratio), (Pre + Post) > 1) %>%
  mutate(
    Significant = p_adj < 0.05,
    stars = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

significant_before_correction <- results %>%
  filter(fisher_p < 0.05)

plot_font_family <- "Arial"
forest_point_size <- 10
forest_star_size <- 10
forest_axis_title_size <- 28
forest_axis_text_size <- 26
forest_title_size <- 22
forest_subtitle_size <- 20
disease_axis_title_size <- 32
disease_axis_text_x_size <- 26
disease_axis_text_y_size <- 32
disease_title_size <- 34
disease_strip_text_size <- 30
legend_text_size <- 26
legend_title_size <- 26

write_xlsx(results, path = "significant_diseaseGPSCs.xlsx")
write_xlsx(significant_before_correction, path = "significant_B4correction_diseaseGPSCs.xlsx")

## Load disease isolate metadata and serotype colours
disease_data <- read_excel("./disesea.xlsx")
serotype_colours <- read_csv("serotype_colours.csv")
names(serotype_colours) <- c("Serotype", "Serotype_color")

## Standardise serotype labels
disease_data$Serotype <- sub("^0+", "", disease_data$Serotype)
disease_data$Serotype <- ifelse(disease_data$Serotype %in% c("15B", "15C"), "15B/15C", disease_data$Serotype)
disease_data$Serotype <- ifelse(grepl("19A", disease_data$Serotype), "19A", disease_data$Serotype)
disease_data$Serotype <- ifelse(grepl("23B", disease_data$Serotype), "23B", disease_data$Serotype)
disease_data$Serotype <- ifelse(grepl("6A", disease_data$Serotype), "6A", disease_data$Serotype)

disease_data <- disease_data %>%
  left_join(serotype_colours, by = "Serotype") %>%
  filter(!is.na(Vaccine_Period))

serotype_colors_vector <- serotype_colours %>%
  distinct(Serotype, Serotype_color) %>%
  deframe()

## Restrict the disease plot to GPSCs with more than one isolate
gpscs_above_threshold <- disease_data %>%
  count(GPSC, name = "n") %>%
  filter(n > 1) %>%
  pull(GPSC)

data_subset <- disease_data %>%
  filter(GPSC %in% gpscs_above_threshold)

gpsc_order <- disease_data %>%
  count(GPSC) %>%
  arrange(desc(n)) %>%
  pull(GPSC)

## Generate a complete GPSC-period-serotype table so zero counts are retained
all_combinations <- expand_grid(
  GPSC = unique(data_subset$GPSC),
  Vaccine_Period = unique(data_subset$Vaccine_Period),
  Serotype = unique(data_subset$Serotype)
)

long_data <- all_combinations %>%
  left_join(
    data_subset %>%
      group_by(GPSC, Vaccine_Period, Serotype) %>%
      summarise(count = n(), .groups = "drop"),
    by = c("GPSC", "Vaccine_Period", "Serotype")
  ) %>%
  mutate(count = replace_na(count, 0))

pcv10_serotypes <- c("4", "6B", "14", "18C", "9V", "19F", "23F", "1", "5", "7F")

long_data <- long_data %>%
  mutate(
    Vaccine_Category = if_else(Serotype %in% pcv10_serotypes, "VT", "NVT"),
    GPSC = factor(GPSC, levels = gpsc_order)
  )

VT_data <- filter(long_data, Vaccine_Category == "VT")
NVT_data <- filter(long_data, Vaccine_Category == "NVT")

VT_serotypes <- unique(VT_data$Serotype)
NVT_serotypes <- unique(NVT_data$Serotype)
VT_labels <- setNames(VT_serotypes, VT_serotypes)
NVT_labels <- setNames(NVT_serotypes, NVT_serotypes)

pattern_values <- c("NVT" = "stripe", "VT" = "none")
gpsc_display_order <- c(67, 862, 62, 8, 54, 21, 2, 9, 884, 170, 27, 5, 117, 10, 879, 26, 233, 251, 65, 355, 22, 891, 268, 53)

forest_plot_data <- tibble(GPSC = gpsc_display_order) %>%
  left_join(results, by = "GPSC") %>%
  mutate(GPSC = factor(GPSC, levels = rev(gpsc_display_order)))

finite_forest_values <- c(
  forest_plot_data$odds_ratio,
  forest_plot_data$ci_lower,
  forest_plot_data$ci_upper
)
finite_forest_values <- finite_forest_values[is.finite(finite_forest_values) & finite_forest_values > 0]

forest_x_lower_limit <- min(finite_forest_values, na.rm = TRUE) / 2
forest_ci_upper_cap <- max(finite_forest_values, na.rm = TRUE) * 1.1
star_column_x <- forest_ci_upper_cap * 1.8
forest_x_upper_limit <- star_column_x * 1.35

forest_plot_data <- forest_plot_data %>%
  mutate(
    odds_ratio_plot = if_else(is.finite(odds_ratio) & odds_ratio > 0, odds_ratio, NA_real_),
    ci_lower_plot = if_else(is.finite(ci_lower) & ci_lower > 0, pmax(ci_lower, forest_x_lower_limit), NA_real_),
    ci_upper_plot = if_else(is.finite(ci_upper) & ci_upper > 0, pmin(ci_upper, forest_ci_upper_cap), forest_ci_upper_cap)
  )

significant_star_data <- forest_plot_data %>%
  filter(!is.na(stars), stars != "")

VT_data <- VT_data %>%
  filter(GPSC %in% gpsc_display_order) %>%
  mutate(GPSC = factor(GPSC, levels = gpsc_display_order))

NVT_data <- NVT_data %>%
  filter(GPSC %in% gpsc_display_order) %>%
  mutate(GPSC = factor(GPSC, levels = gpsc_display_order))

## Forest plot: odds ratios for change in disease GPSCs after PCV10 introduction
forest_plot <- ggplot(forest_plot_data, aes(x = odds_ratio_plot, y = GPSC)) +
  geom_errorbarh(
    aes(xmin = ci_lower_plot, xmax = ci_upper_plot),
    height = 0.3,
    na.rm = TRUE
  ) +
  geom_point(aes(color = vaccine_type), size = forest_point_size, na.rm = TRUE) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_text(
    data = significant_star_data,
    aes(x = star_column_x, y = GPSC, label = stars),
    inherit.aes = FALSE,
    hjust = 0,
    color = "maroon4",
    fontface = "bold",
    size = forest_star_size
  ) +
  scale_x_log10(
    limits = c(forest_x_lower_limit, forest_x_upper_limit),
    breaks = breaks_log(n = 4)
  ) +
  scale_color_manual(
    values = c("VT" = "darkgray", "NVT" = "steelblue", "Mixed" = "green"),
    breaks = c("VT", "NVT", "Mixed"),
    labels = c("VT-GPSCs", "NVT-GPSCs", "Mixed VT/NVT-GPSCs"),
    name = NULL,
    guide = guide_legend(order = 1, override.aes = list(size = 6))
  ) +
  labs(
    x = "Odds ratio (log scale)",
    y = "GPSC",
    title = "Odds Ratios for Change in Disease GPSCs After PCV10 Introduction",
    subtitle = "Reference group: GPSC 533"
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_family = plot_font_family) +
  theme(
    legend.position = "right",
    text = element_text(family = plot_font_family, face = "bold", color = "black"),
    axis.title.x = element_text(size = forest_axis_title_size, face = "bold"),
    axis.title.y = element_text(size = forest_axis_title_size, face = "bold"),
    axis.text.x = element_text(size = forest_axis_text_size, face = "bold"),
    axis.text.y = element_text(size = forest_axis_text_size, face = "bold"),
    legend.text = element_text(size = legend_text_size, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = forest_title_size, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = forest_subtitle_size, face = "bold"),
    plot.margin = margin(10, 80, 10, 10)
  )

forest_plot

## Disease plot: serotype composition before and after PCV10 introduction
GPSC_Plot <- ggplot() +
  geom_col_pattern(
    data = VT_data,
    pattern_key_scale_factor = 0.2,
    aes(
      x = count,
      y = factor(Vaccine_Period, levels = c("post", "pre")),
      fill = Serotype,
      pattern = Vaccine_Category
    ),
    color = "black",
    pattern_density = 0.2,
    pattern_spacing = 0.2
  ) +
  facet_grid(GPSC ~ ., switch = "y") +
  scale_fill_manual(
    aesthetics = "fill",
    values = serotype_colors_vector,
    labels = VT_labels,
    breaks = VT_serotypes,
    name = "VTs",
    guide = guide_legend(
      order = 2,
      ncol = 2,
      override.aes = list(
        pattern = "none",
        pattern_fill = NA,
        pattern_colour = NA
      )
    )
  ) +
  scale_pattern_fill_manual(values = pattern_values, guide = "none") +
  new_scale_fill() +
  geom_col_pattern(
    data = NVT_data,
    pattern_key_scale_factor = 0.2,
    aes(
      x = count,
      y = factor(Vaccine_Period, levels = c("post", "pre")),
      fill = Serotype,
      pattern = Vaccine_Category
    ),
    color = "black",
    pattern_density = 0.5,
    pattern_spacing = 0.2
  ) +
  scale_fill_manual(
    aesthetics = "fill",
    values = serotype_colors_vector,
    labels = NVT_labels,
    breaks = NVT_serotypes,
    name = "NVTs",
    guide = guide_legend(
      order = 3,
      ncol = 2,
      override.aes = list(
        pattern = "stripe",
        pattern_fill = "black",
        pattern_colour = "black",
        pattern_density = 0.35,
        pattern_spacing = 0.03
      )
    )
  ) +
  scale_pattern_manual(values = pattern_values, guide = "none") +
  scale_x_continuous(position = "top", expand = c(0, 0)) +
  labs(
    title = "Disease GPSC Composition Before and After PCV10 Introduction",
    x = "Number of samples (n)",
    y = NULL
  ) +
  theme(
    text = element_text(family = plot_font_family, face = "bold", color = "black"),
    axis.text.x = element_text(size = disease_axis_text_x_size, face = "bold"),
    axis.text.y = element_text(size = disease_axis_text_y_size, face = "bold"),
    axis.title.x = element_text(size = disease_axis_title_size, face = "bold"),
    axis.title.y = element_text(size = disease_axis_title_size, face = "bold"),
    legend.text = element_text(size = legend_text_size, face = "bold"),
    legend.title = element_text(size = legend_title_size, face = "bold"),
    plot.title = element_text(size = disease_title_size, hjust = 0, face = "bold"),
    plot.caption = element_text(size = disease_axis_title_size),
    strip.placement = "outside",
    strip.background = element_rect(fill = NA, color = "white"),
    strip.text.y.left = element_text(size = disease_strip_text_size, angle = 0, face = "bold"),
    panel.background = element_rect(fill = NA, color = "white"),
    panel.spacing = unit(0.05, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.key = element_rect(fill = "white", color = NA),
    legend.key.size = unit(0.8, "cm")
  )

## Export the disease-only plot
ggsave("Figure.png", plot = GPSC_Plot, width = 30, height = 20, dpi = 300)

## Combine the disease composition plot and forest plot
large_plot <- GPSC_Plot + forest_plot + guide_area() +
  plot_layout(widths = c(3.8, 2.2, 1.3), guides = "collect") &
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.text = element_text(family = plot_font_family, face = "bold", size = legend_text_size),
    legend.title = element_text(family = plot_font_family, face = "bold", size = legend_title_size)
  )

## Export the final multi-panel disease figure
ggsave("Figure3.png", plot = large_plot, width = 30, height = 20, dpi = 300)
