
# Load libraries
library(tidyverse)
library(readxl)
library(patchwork)


# Read metadata
meta <- read_excel("KLF_population_ex2011.xlsx")

# Create age groups
meta <- meta %>%
  mutate(age_group = case_when(
    age_yr < 5 ~ "<5 years",
    age_yr >= 5 & age_yr <= 18 ~ "5-18 years",
    age_yr > 18 ~ ">18 years",
    TRUE ~ NA_character_
  ))

meta <- meta %>%
  mutate(age_group = factor(age_group,
                            levels = c("<5 years", "5-18 years", ">18 years")))


carriage <- meta %>%
  filter(disease_carriage == "carriage")

p1 <- carriage %>%
  count(age_group, sex) %>%
  ggplot(aes(x = age_group, y = n, fill = sex)) +
  geom_col() +
  labs(title = "Carriage Participants",
       x = NULL,
       y = "Number of participants",
       fill = "Sex") +
  theme_bw()

p2 <- carriage %>%
  count(year, age_group) %>%
  ggplot(aes(x = year, y = n)) +
  geom_line(color = "darkgreen") +
  geom_point(color = "darkgreen") +
  facet_wrap(~age_group, scales = "free_y") +
  labs(x = "Year",
       y = "Number of participants") +
  theme_bw()

disease <- meta %>%
  filter(disease_carriage == "disease")

p3 <- disease %>%
  count(age_group, sex) %>%
  ggplot(aes(x = age_group, y = n, fill = sex)) +
  geom_col() +
  labs(title = "Disease Participants",
       x = NULL,
       y = "Number of participants",
       fill = "Sex") +
  theme_bw()
p4 <- disease %>%
  count(year, age_group, source) %>%
  ggplot(aes(x = year, y = n, color = source)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_wrap(~age_group, scales = "free_y") +
  scale_color_manual(values = c(
    "blood" = "red",
    "cerebrospinal fluid" = "gray",
    "pleural fluid" = "yellow"
  )) +
  labs(x = "Year",
       y = "Number of participants",
       color = "Sample source") +
  theme_bw()



final_plot <- (p1 | p2) / (p3 | p4)

final_plot


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

