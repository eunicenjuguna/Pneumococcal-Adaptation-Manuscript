library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)

# Read your Excel file
data <- read_excel("QC_fail.xlsx")

# View first rows
head(data)


# Convert from wide to long
data_long <- data %>%
  pivot_longer(
    cols = c(Assembly_QC, Taxonomy_QC, Mapping_QC),
    names_to = "QC_type",
    values_to = "QC_result"
  )

# Keep only FAILs
fail_data <- data_long %>% filter(QC_result == "FAIL")

ggplot(fail_data, aes(x = disease_carriage, fill = QC_type)) +
  geom_bar(position = "stack") +
  facet_wrap(~ source) +
  labs(
    title = "QC Failures by Disease/Carriage and Source",
    x = "Disease vs Carriage",
    y = "Number of FAILs"
  ) +
  theme_minimal()
fail_data <- fail_data %>%
  mutate(age_group = case_when(
    age_yr < 5 ~ "Under 5",
    age_yr >= 5 & age_yr < 18 ~ "5-18",
    age_yr >= 18 ~ "18+",
    TRUE ~ NA_character_
  ))

y<-ggplot(fail_data, aes(x = age_group, fill = QC_type)) +
  geom_bar(position = "stack") +
  facet_wrap(~ disease_carriage) +
  labs(
    title = "QC Failures by Age Group and Disease/Carriage",
    x = "Age Group",
    y = "Number of FAILs"
  ) +
  theme_minimal()

