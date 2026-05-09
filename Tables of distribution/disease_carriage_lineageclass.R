
library(dplyr)
library(readxl)
library(purrr)
library(writexl)
library(readr)

data <- read_excel("Major_lineagePopulation.xlsx")

data <- data %>%
  mutate(
    lineage_class = factor(lineage_class, levels = c("VT","NVT","mixed")),
    vaccine_period = factor(vaccine_period, levels = c("pre","post")),
    sample_type = factor(disease_carriage, levels = c("disease","carriage"))
  )

count_table <- data %>%
  count(sample_type, lineage_class, vaccine_period) %>%
  pivot_wider(
    names_from = vaccine_period,
    values_from = n,
    values_fill = 0
  )

count_table

stats_results <- count_table %>%
  group_by(sample_type, lineage_class) %>%
  summarise(
    pre = sum(pre),
    post = sum(post),
    .groups = "drop"
  ) %>%
  
  # Build contingency tables
  mutate(
    test = map2(pre, post, ~{
      mat <- matrix(
        c(.x, .y,
          sum(pre) - .x, sum(post) - .y),
        nrow = 2
      )
      
      if(any(mat < 5)){
        fisher.test(mat)
      } else {
        chisq.test(mat)
      }
    }),
    
    p_value = map_dbl(test, ~.x$p.value),
    
    test_used = map_chr(test, ~{
      ifelse(.x$method == "Fisher's Exact Test",
             "Fisher", "Chi-square")
    })
  ) %>%
  
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    significant = ifelse(p_adj < 0.05, "Yes", "No")
  )

stats_results

stats_results <- stats_results %>%
  mutate(
    fold_change = round(post / pre, 2),
    direction = case_when(
      fold_change > 1 ~ "Increase",
      fold_change < 1 ~ "Decrease",
      TRUE ~ "No change"
    )
  )

stats_results
write_xlsx(stats_results,"VT_NVT_Mixed changes.xlsx")



# ----- Plot 1: Lineage class by vaccine period, separated by disease/carriage -----
ggplot(data, aes(x = lineage_class, fill = vaccine_period)) +
  geom_bar(position = position_dodge()) +
  scale_fill_manual(values = c("pre"="lightgrey", "post"="darkgrey")) +
  facet_wrap(~sample_type) +
  labs(
    x = "Lineage class",
    y = "Number of samples",
    title = "Distribution of samples by lineage class (pre vs post vaccine)",
    fill = "Vaccine period"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# ----- Subset mixed lineages -----
mixed_data <- data %>%
  filter(lineage_class == "mixed") %>%
  mutate(GPSC = paste0("GPSC", GPSC)) %>%
  mutate(GPSC = factor(GPSC, levels = c(
    "GPSC5","GPSC9","GPSC10","GPSC21","GPSC22",
    "GPSC54","GPSC61","GPSC62","GPSC65","GPSC67",
    "GPSC92","GPSC170","GPSC184","GPSC251","GPSC862","GPSC879"
  )))

# Define GPSC colours
gpsc_colours <- c(
  "GPSC5"="black", "GPSC9"="green", "GPSC10"="maroon", "GPSC21"="lightgreen",
  "GPSC22"="pink", "GPSC54"="gray", "GPSC61"="blue", "GPSC62"="darkblue",
  "GPSC65"="purple", "GPSC67"="yellow", "GPSC92"="skyblue", "GPSC170"="darkgreen",
  "GPSC184"="lightblue", "GPSC251"="blue", "GPSC862"="orange", "GPSC879"="red"
)

# ----- Plot 2: Mixed lineage GPSCs by vaccine period, separated by disease/carriage -----
ggplot(mixed_data, aes(x = lineage_class_vp, fill = GPSC)) +
  geom_bar() +
  scale_fill_manual(values = gpsc_colours, na.value = "lightgray") +
  facet_wrap(~sample_type) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Mixed lineage (pre vs post)", y = "Number of samples", fill = "GPSC")


###LINEAGE Of interest

serotype_colours <- read_csv("serotype_colours.csv", col_types = cols())

serotype_colours_vec <- setNames(
  serotype_colours$In_silico_serotype__colour,
  serotype_colours$In_silico_serotype
)

premixed_postNVT_pre <- mixed_data %>%
  filter(lineage_class_vp == "premixed_postNVT", vaccine_period=="pre") %>%
  mutate(
    GPSC = factor(GPSC),
    Serotype = factor(Serotype_merged))
premixed_postNVT_post <- mixed_data %>%
  filter(lineage_class_vp == "premixed_postNVT", vaccine_period=="post") %>%
  mutate(
    GPSC = factor(GPSC),
    Serotype = factor(Serotype_merged))

ggplot(premixed_postNVT_pre, aes(x = GPSC, fill = Serotype)) +
  geom_bar() +
  scale_fill_manual(values = serotype_colours_vec, na.value = "lightgray") +
  facet_wrap(~sample_type) +
  labs(
    x = "VT GPSC",
    y = "Number of samples",
    title = "Distribution of VT GPSCs by Serotype",
    fill = "Serotype"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )


preVT_postmixed_pre <- mixed_data %>%
  filter(lineage_class_vp == "preVT_postmixed", vaccine_period=="pre") %>%
  mutate(
    GPSC = factor(GPSC),
    Serotype = factor(Serotype))
preVT_postmixed_post <- mixed_data %>%
  filter(lineage_class_vp == "preVT_postmixed", vaccine_period=="post") %>%
  mutate(
    GPSC = factor(GPSC),
    Serotype = factor(Serotype))

ggplot(preVT_postmixed_post, aes(x = GPSC, fill = Serotype)) +
  geom_bar() +
  scale_fill_manual(values = serotype_colours_vec, na.value = "lightgray") +
  facet_wrap(~sample_type) +
  labs(
    x = "VT GPSC",
    y = "Number of samples",
    title = "Distribution of VT GPSCs by Serotype",
    fill = "Serotype"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

premixed_postmixed_pre <- mixed_data %>%
  filter(lineage_class_vp == "premixed_postmixed", vaccine_period=="pre") %>%
  mutate(
    GPSC = factor(GPSC),
    Serotype = factor(Serotype))
premixed_postmixed_post <- mixed_data %>%
  filter(lineage_class_vp == "premixed_postmixed", vaccine_period=="post") %>%
  mutate(
    GPSC = factor(GPSC),
    Serotype = factor(Serotype))

ggplot(premixed_postmixed_pre, aes(x = GPSC, fill = Serotype)) +
  geom_bar() +
  scale_fill_manual(values = serotype_colours_vec, na.value = "lightgray") +
  facet_wrap(~sample_type) +
  labs(
    x = "VT GPSC",
    y = "Number of samples",
    title = "Distribution of VT GPSCs by Serotype",
    fill = "Serotype"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )
