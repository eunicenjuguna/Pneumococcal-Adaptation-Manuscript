library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)

# 1. Read and Prepare Data
data <- read_excel("QC_fail.xlsx")

# Convert from wide to long format
data_long <- data %>%
  pivot_longer(
    cols = c(Assembly_QC, Taxonomy_QC, Mapping_QC),
    names_to = "QC_type",
    values_to = "QC_result"
  )

# 2. Filter for FAILS and create ordered Age Groups
# Define the order explicitly here
age_levels <- c("Under 5", "5-18", "18+")

fail_data <- data_long %>% 
  filter(QC_result == "FAIL") %>%
  mutate(age_group = case_when(
    age_yr < 5 ~ "Under 5",
    age_yr >= 5 & age_yr < 18 ~ "5-18",
    age_yr >= 18 ~ "18+",
    TRUE ~ NA_character_
  )) %>%
  mutate(age_group = factor(age_group, levels = age_levels))

# 3. Create Custom Theme for bold, black, larger text
custom_theme <- theme_minimal() +
  theme(
    axis.text = element_text(size = 12, face = "bold", color = "black"),
    axis.title = element_text(size = 14, face = "bold", color = "black"),
    strip.text = element_text(size = 12, face = "bold", color = "black"),
    plot.title = element_text(size = 16, face = "bold")
  )

# 4. Generate the Plot
allsamples <- ggplot(fail_data, aes(x = age_group, fill = QC_type)) +
  geom_bar(position = "stack") +
  facet_wrap(~ disease_carriage) +
  labs(
    title = "QC Failures by Age Group and Disease/Carriage",
    x = "Age Group",
    y = "Number of FAILs"
  ) +
  custom_theme

# Display plot
allsamples
