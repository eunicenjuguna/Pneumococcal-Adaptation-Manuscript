### =========================
## Libraries
## =========================
library(tidyverse)
library(readxl)
library(readr)
library(ggplot2)
library(ggpattern)
library(ggnewscale)

## =========================
## Read data
## =========================
data <- read_excel("Major_lineagePopulation.xlsx")

serotype_colours <- read_csv("serotype_colours.csv",
                             col_names = c("Serotype", "Serotype__color"))

## =========================
## Filter to analysis set FIRST
#repeat in each group
## =========================
data <- data %>%
  filter(
  lineage_class == "mixed",
    lineage_class_vp == "premixed_postNVT",
   !is.na(vaccine_period)
 )

#data <- data %>%
# filter(
#    lineage_class == "mixed",
#   lineage_class_vp == "preVT_postmixed",
#    !is.na(vaccine_period)
# )

#data <- data %>%
# filter(
#   lineage_class == "mixed",
#  lineage_class_vp == "premixed_postmixed",
#   !is.na(vaccine_period)
# )

## =========================
## Serotype cleaning
## =========================
data <- data %>%
  mutate(
    Serotype = sub("^0+", "", Serotype),
    Serotype = case_when(
      Serotype %in% c("15B", "15C") ~ "15B/15C",
      grepl("19A", Serotype) ~ "19A",
      grepl("23B", Serotype) ~ "23B",
      grepl("6A", Serotype) ~ "6A",
      TRUE ~ Serotype
    )
  )

## =========================
## Join colours (SAFE)
## =========================
data <- data %>%
  left_join(serotype_colours, by = "Serotype")

serotype_colors_vector <- data %>%
  select(Serotype, Serotype__color) %>%
  distinct() %>%
  deframe()

## =========================
## Remove rare GPSCs (n ≤ 4)
## =========================
gpsc_counts <- data %>%
  count(GPSC)

aboveFour <- gpsc_counts %>%
  filter(n > 4) %>%
  pull(GPSC)

data_subset <- data %>%
  filter(GPSC %in% aboveFour)

## =========================
## Define GPSC order AFTER filtering  ✅ CRITICAL
## =========================
gpsc_order <- data_subset %>%
  count(GPSC) %>%
  arrange(desc(n)) %>%
  pull(GPSC)

data_subset <- data_subset %>%
  mutate(GPSC = factor(GPSC, levels = gpsc_order))

## =========================
## Expand combinations (avoid missing bars)
## =========================
all_combinations <- expand.grid(
  GPSC = levels(data_subset$GPSC),
  vaccine_period = unique(data_subset$vaccine_period),
  Serotype = unique(data_subset$Serotype)
)

## =========================
## Long format with counts
## =========================
long_data <- all_combinations %>%
  left_join(
    data_subset %>%
      group_by(GPSC, vaccine_period, Serotype) %>%
      summarise(count = n(), .groups = "drop"),
    by = c("GPSC", "vaccine_period", "Serotype")
  ) %>%
  mutate(count = replace_na(count, 0))

## =========================
## VT vs NVT classification (PCV13)
## =========================
VT_serotypes <- c("1", "4", "5", "6B", "7F", "9V", "14", "18C", "19F", "23F")

long_data <- long_data %>%
  mutate(
    Vaccine_Category = ifelse(Serotype %in% VT_serotypes, "VT", "NVT"),
    GPSC = factor(GPSC, levels = gpsc_order)
  )

## =========================
## Plot
## =========================
GPSC_Plot <- ggplot(long_data) +
  geom_col_pattern(
    aes(
      x = factor(vaccine_period, levels = c("pre", "post")),
      y = count,
      fill = Serotype,
      pattern = Vaccine_Category
    ),
    color = "black",
    pattern_density = 0.3,
    pattern_spacing = 0.2
  ) +
  facet_grid(~ GPSC, switch = "x") +
  scale_fill_manual(
    values = serotype_colors_vector,
    name = "Serotype"
  ) +
  scale_pattern_manual(
    values = c(
      VT = "none",
      NVT = "stripe"
    ),
    name = "Vaccine category"
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    x = NULL,
    y = "Number of Isolates (n)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 16, face = "bold"),
    axis.text.y = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    strip.background = element_blank(),
    strip.text = element_text(size = 16, face = "bold"),
    strip.placement = "outside",
    panel.spacing = unit(-0.01, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14)
  )

## =========================
## Save
## =========================
ggsave("premixed_postNVT.png", plot = GPSC_Plot, width = 30, height = 20, dpi = 300)

