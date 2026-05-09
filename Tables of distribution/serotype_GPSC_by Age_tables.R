
# Load libraries
library(tidyverse)
library(readxl)
library(dplyr)
library(tidyr)
library(writexl)
library(readr)
library(ggplot2)



# Read data
data <- read_excel("3035_KLFpopulation_data.xlsx")

data <- data %>%
  mutate(age_yr = case_when(
    age_yr < 5 ~ "<5",
    age_yr >= 5 & age_yr <= 18 ~ "5-18",
    age_yr > 18 ~ ">18"
  ))

data <- data %>%
  mutate(
    vaccine_period = factor(vaccine_period, levels = c("pre", "post")),
    disease_carriage = factor(disease_carriage, levels = c("disease", "carriage"))
  )

data <- data %>%
  mutate(vaccine_type = factor(vaccine_type,
                               levels = c("VT", "NVT")))
serotype_table <- data %>%
  group_by(Serotype_merged, vaccine_type,
           vaccine_period, disease_carriage, age_yr) %>%
  
  summarise(n = n(), .groups = "drop") %>%
  
  unite(Group,
        vaccine_period, disease_carriage, age_yr,
        sep = "_") %>%
  
  pivot_wider(
    names_from = Group,
    values_from = n,
    values_fill = 0
  )

serotype_table <- serotype_table %>%
  mutate(Total = rowSums(across(where(is.numeric)))) %>%
  arrange(vaccine_type, desc(Total))

vt_nvt_summary <- serotype_table %>%
  group_by(vaccine_type) %>%
  summarise(
    total_isolates = sum(Total),
    .groups = "drop"
  ) %>%
  mutate(
    proportion = round(total_isolates / sum(total_isolates), 3)
  )

vt_nvt_summary

write_xlsx(serotype_table,"Serotype_distribution_by_AG.xlsx")
 

# Build lineage table
lineage_table <- data %>%
  mutate(GPSC = as.character(GPSC)) %>%   # IMPORTANT: prevent summing GPSC
  group_by(GPSC, lineage_class,
           vaccine_period, disease_carriage, age_yr) %>%
  
  summarise(n = n(), .groups = "drop") %>%
  
  unite(Group,
        vaccine_period, disease_carriage, age_yr,
        sep = "_") %>%
  
  pivot_wider(
    names_from = Group,
    values_from = n,
    values_fill = 0
  ) %>%
  
  mutate(
    Total = rowSums(across(where(is.numeric)))
  ) %>%
  
  arrange(lineage_class, desc(Total))


# VT / NVT / Mixed summary
vt_nvt_mixed_summary <- lineage_table %>%
  group_by(lineage_class) %>%
  summarise(
    total_isolates = sum(Total),
    .groups = "drop"
  ) %>%
  
  mutate(
    proportion = round(total_isolates / sum(total_isolates), 3)
  )

# View result
vt_nvt_mixed_summary


write_xlsx(lineage_table,"GPSC_distribution_by_AG.xlsx")


# Build lineage table
mixed_data <- data %>%
  filter(lineage_class == "mixed")

lineage_mixed_table <- mixed_data %>%
  mutate(GPSC = as.character(GPSC)) %>%   # IMPORTANT: prevent summing GPSC
  group_by(GPSC, lineage_class_vp,
           vaccine_period, disease_carriage, age_yr) %>%
  
  summarise(n = n(), .groups = "drop") %>%
  
  unite(Group,
        vaccine_period, disease_carriage, age_yr,
        sep = "_") %>%
  
  pivot_wider(
    names_from = Group,
    values_from = n,
    values_fill = 0
  ) %>%
  
  mutate(
    Total = rowSums(across(where(is.numeric)))
  ) %>%
  
  arrange(lineage_class_vp, desc(Total))

write_xlsx(lineage_mixed_table,"mixedGPSC_distribution_by_AG.xlsx")

