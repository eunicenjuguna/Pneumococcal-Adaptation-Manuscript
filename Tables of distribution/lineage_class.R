# Load libraries
library(tidyverse)
library(readxl)
library(dplyr)
library(tidyr)
library(writexl)
library(readr)
library(ggplot2)

#Lineages

# Read data
data <- read_excel("Major_lineagePopulation.xlsx")

# Read serotype colours
serotype_colours <- read_csv("serotype_colours.csv", col_types = cols())

serotype_colours_vec <- setNames(
  serotype_colours$In_silico_serotype__colour,
  serotype_colours$In_silico_serotype
)

# Clean data
plot_data <- data %>%
  mutate(
    GPSC = factor(GPSC),
    Serotype = factor(Serotype_merged)
  )

# Plot
ggplot(plot_data,
       aes(x = GPSC, fill = Serotype)) +
  geom_bar() +
  scale_fill_manual(values = serotype_colours_vec, na.value = "lightgray") +
  labs(
    x = "GPSC",
    y = "Number of samples",
    title = "Distribution of GPSCs by Serotype",
    fill = "Serotype"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

# Read data
data <- read_excel("Major_lineagePopulation.xlsx")

# Read serotype colours
serotype_colours <- read_csv("serotype_colours.csv", col_types = cols())

serotype_colours_vec <- setNames(
  serotype_colours$In_silico_serotype__colour,
  serotype_colours$In_silico_serotype
)

# Prepare data
plot_data <- data %>%
  mutate(
    GPSC = factor(GPSC),
    Serotype = factor(Serotype_merged)
  )

# ---- Keep only GPSCs with >1 serotype ----
multi_serotype_gpsc <- plot_data %>%
  group_by(GPSC) %>%
  summarise(n_serotypes = n_distinct(Serotype)) %>%
  filter(n_serotypes > 1)

plot_data_filtered <- plot_data %>%
  filter(GPSC %in% multi_serotype_gpsc$GPSC)

# ---- Plot ----
ggplot(plot_data_filtered,
       aes(x = GPSC, fill = Serotype)) +
  geom_bar() +
  scale_fill_manual(values = serotype_colours_vec, na.value = "lightgray") +
  labs(
    x = "GPSC (only >1 serotype)",
    y = "Number of samples",
    title = "GPSCs with More Than One Serotype",
    fill = "Serotype"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )


#classification plots

# Read major population data

data <- read_excel("Major_lineagePopulation.xlsx")

# Ensure factors are ordered
data$lineage_class <- factor(data$lineage_class, levels = c("VT","NVT","mixed"))
data$vaccine_period <- factor(data$vaccine_period, levels = c("pre", "post"))

ggplot(data, aes(x = lineage_class, fill = vaccine_period)) +
  geom_bar(position = position_dodge()) +   # <-- dodge for side-by-side bars
  scale_fill_manual(values = c(
    "pre" = "lightgrey",
    "post" = "darkgrey"
  )) +
  labs(
    x = "Lineage class",
    y = "Number of samples",
    title = "Distribution of samples by lineage class (pre vs post vaccine)",
    fill = "Vaccine period"
  ) +
  theme_minimal()


#futher classify mixed

# Subset mixed lineages
mixed_data <- data %>%
  filter(lineage_class == "mixed")

# Convert GPSC numbers to character with "GPSC" prefix
mixed_data <- mixed_data %>%
  mutate(GPSC = paste0("GPSC", GPSC))  # e.g., 5 -> GPSC5

# Convert to factor with levels matching the colour vector
mixed_data$GPSC <- factor(mixed_data$GPSC,
                          levels = c("GPSC5","GPSC9","GPSC10","GPSC21","GPSC22",
                                     "GPSC54","GPSC61","GPSC62","GPSC65","GPSC67",
                                     "GPSC92","GPSC170","GPSC184","GPSC251","GPSC862","GPSC879"))

# Define colours
gpsc_colours <- c(
  "GPSC5" = "black",
  "GPSC9" = "green",
  "GPSC10" = "maroon",
  "GPSC21" = "lightgreen",
  "GPSC22" = "pink",
  "GPSC54" = "gray",
  "GPSC61" = "blue",
  "GPSC62" = "darkblue",
  "GPSC65" = "purple",
  "GPSC67" = "yellow",
  "GPSC92" = "skyblue",
  "GPSC170" = "darkgreen",
  "GPSC184" = "lightblue",
  "GPSC251" = "blue",
  "GPSC862" = "orange",
  "GPSC879" = "red"
)

# Plot
ggplot(mixed_data,
       aes(x = lineage_class_vp, fill = GPSC)) +
  geom_bar() +
  scale_fill_manual(values = gpsc_colours, na.value = "lightgray") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#now plot mixed gpsc by their serotypes

library(ggplot2)
library(dplyr)
library(readr)

# Read serotype colours CSV
serotype_colours <- read_csv("serotype_colours.csv", col_types = cols())


# Named vector for colours
serotype_colours_vec <- setNames(
  serotype_colours$In_silico_serotype__colour,
  serotype_colours$In_silico_serotype
)

# Subset mixed lineages and fix columns
mixed_data <- data %>%
  filter(lineage_class == "mixed") %>%
  mutate(
    Serotype = factor(Serotype_merged),
    GPSC = factor(GPSC)  # convert to factor for discrete x-axis
  )

# Plot
ggplot(mixed_data,
       aes(x = GPSC, fill = Serotype)) +
  geom_bar() +
  scale_fill_manual(values = serotype_colours_vec, na.value = "lightgray") +
  labs(
    x = "Mixed GPSC",
    y = "Number of samples",
    title = "Distribution of mixed GPSCs by Serotype"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )


library(ggplot2)
library(dplyr)
library(readr)

# Read serotype colours CSV
serotype_colours <- read_csv("serotype_colours.csv", col_types = cols())


# Named vector for colours
serotype_colours_vec <- setNames(
  serotype_colours$In_silico_serotype__colour,
  serotype_colours$In_silico_serotype
)

# Subset mixed lineages and fix columns
mixed_data <- data %>%
  filter(lineage_class == "NVT") %>%
  mutate(
    Serotype = factor(Serotype_merged),
    GPSC = factor(GPSC)  # convert to factor for discrete x-axis
  )

# Plot
ggplot(mixed_data,
       aes(x = GPSC, fill = Serotype)) +
  geom_bar() +
  scale_fill_manual(values = serotype_colours_vec, na.value = "lightgray") +
  labs(
    x = "Mixed GPSC",
    y = "Number of samples",
    title = "Distribution of mixed GPSCs by Serotype"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

library(ggplot2)
library(dplyr)
library(readr)

# Read serotype colours CSV
serotype_colours <- read_csv("serotype_colours.csv", col_types = cols())


# Named vector for colours
serotype_colours_vec <- setNames(
  serotype_colours$In_silico_serotype__colour,
  serotype_colours$In_silico_serotype
)

# Subset mixed lineages and fix columns
mixed_data <- data %>%
  filter(lineage_class == "VT") %>%
  mutate(
    Serotype = factor(Serotype_merged),
    GPSC = factor(GPSC)  # convert to factor for discrete x-axis
  )

# Plot
ggplot(mixed_data,
       aes(x = GPSC, fill = Serotype)) +
  geom_bar() +
  scale_fill_manual(values = serotype_colours_vec, na.value = "lightgray") +
  labs(
    x = "Mixed GPSC",
    y = "Number of samples",
    title = "Distribution of mixed GPSCs by Serotype"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )



