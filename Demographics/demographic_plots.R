# Load required libraries
library(tidyverse)
library(readxl)
library(patchwork)

# Import population metadata
meta <- read_excel("KLF_population_ex2011.xlsx")

# Construct structured age groups as an ordered factor
meta <- meta %>%
  mutate(age_group = case_when(
    age_yr < 5 ~ "<5 years",
    age_yr >= 5 & age_yr <= 18 ~ "5-18 years",
    age_yr > 18 ~ ">18 years",
    TRUE ~ NA_character_
  )) %>%
  mutate(age_group = factor(
    age_group,
    levels = c("<5 years", "5-18 years", ">18 years")
  ))

# Programmatically compute the global temporal range
temporal_range <- range(meta$year, na.rm = TRUE)

# Define separate vertical scale constraints
vertical_limit_A <- c(0, 180)  
vertical_limit_B <- c(0, 1500) 

# Define the custom theme for bold, black, larger text
custom_theme <- theme_bw() +
  theme(
    axis.text = element_text(size = 12, face = "bold", color = "black"),
    axis.title = element_text(size = 14, face = "bold", color = "black"),
    strip.text = element_text(size = 12, face = "bold", color = "black") # Updates facet labels too
  )

# ==================== CARRIAGE COHORT SUBPLOTS ====================

carriage_data <- meta %>% filter(disease_carriage == "carriage")

p1 <- carriage_data %>%
  count(age_group, sex) %>%
  ggplot(aes(x = age_group, y = n, fill = sex)) +
  geom_col() +
  scale_y_continuous(limits = vertical_limit_B, expand = expansion(mult = c(0, 0.05))) +
  labs(title = "Carriage Participants", x = NULL, y = "Number of participants", fill = "Sex") +
  custom_theme

p2 <- carriage_data %>%
  count(year, age_group) %>%
  ggplot(aes(x = year, y = n)) +
  geom_line(color = "darkgreen") +
  geom_point(color = "darkgreen") +
  facet_wrap(~age_group) +
  scale_x_continuous(limits = temporal_range, breaks = seq(from = temporal_range[1], to = temporal_range[2], by = 2)) +
  scale_y_continuous(limits = vertical_limit_A) +
  labs(x = "Year", y = "Number of participants") +
  custom_theme

# ==================== DISEASE COHORT SUBPLOTS ====================

disease_data <- meta %>% filter(disease_carriage == "disease")

p3 <- disease_data %>%
  count(age_group, sex) %>%
  ggplot(aes(x = age_group, y = n, fill = sex)) +
  geom_col() +
  scale_y_continuous(limits = vertical_limit_B, expand = expansion(mult = c(0, 0.05))) +
  labs(title = "Disease Participants", x = NULL, y = "Number of participants", fill = "Sex") +
  custom_theme

p4 <- disease_data %>%
  count(year, age_group, source) %>%
  ggplot(aes(x = year, y = n, color = source)) +
  geom_line(linewidth = 1) + 
  geom_point(size = 2) +
  facet_wrap(~age_group) +
  scale_color_manual(values = c("blood" = "red", "cerebrospinal fluid" = "gray", "pleural fluid" = "yellow")) +
  scale_x_continuous(limits = temporal_range, breaks = seq(from = temporal_range[1], to = temporal_range[2], by = 2)) +
  scale_y_continuous(limits = vertical_limit_A) +
  labs(x = "Year", y = "Number of participants", color = "Sample source") +
  custom_theme

# ==================== COHORT ASSEMBLY ====================

final_plotA <- (p2 / p4) + plot_layout(axes = "collect", axis_titles = "collect", guides = "keep")
final_plotB <- (p1 / p3) + plot_layout(axes = "collect", axis_titles = "collect", guides = "keep")

final_plotA
final_plotB



#CONVERT TO A TABLE INSTEAD

carriage_table <- carriage %>%
  count(age_group, sex) %>%
  arrange(age_group, sex)

carriage_time_table <- carriage %>%
  count(year, age_group) %>%
  arrange(age_group, year)

disease_table <- disease %>%
  count(age_group, sex) %>%
  arrange(age_group, sex)

disease_time_table <- disease %>%
  count(year, age_group, source) %>%
  arrange(age_group, year, source)

library(writexl)

write_xlsx(
  list(
    carriage_by_age_sex = carriage_table,
    carriage_over_time = carriage_time_table,
    disease_by_age_sex = disease_table,
    disease_over_time = disease_time_table
  ),
  "KLF_population_summary_tables.xlsx"
)

summary_table <- meta %>%
  group_by(disease_carriage, age_group, sex) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(disease_carriage, age_group, sex)

table1_totals <- meta %>%
  group_by(disease_carriage, age_group, sex) %>%
  summarise(n = n(), .groups = "drop") %>%
  bind_rows(
    meta %>%
      group_by(disease_carriage, age_group) %>%
      summarise(n = n(), sex = "Total", .groups = "drop")
  )

