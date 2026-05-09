## Load packages. Install them first if you haven't used them before with e.g. install.packages("tidyverse")
library(tidyverse)
library(ggpattern)
library(patchwork)
library(ggpubr)
library(ggnewscale)
library(readxl)
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
data <- read_excel("./pre_post_disease_GPSCs.xlsx")


# Get reference counts for serotype 6A
ref <- data %>% filter(GPSC == "53")
ref_pre <- ref$Pre
ref_post <- ref$Post

# Ensure the reference values are valid
if (is.na(ref_pre) || is.na(ref_post) || ref_pre < 0 || ref_post < 0) {
  stop("Invalid reference serotype counts (non-finite or negative).")
}

# Perform Fisher's test with robust error handling
results <- data %>%
  filter(GPSC != "53") %>%
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


write_xlsx(results, path = "./53reference_disease.xlsx")

results <- read_excel("./53reference_disease_edited.xlsx")
# Optional: sort by odds ratio or significance
results <- results %>%
  filter(!is.na(odds_ratio)) %>%
  mutate(GPSC = fct_reorder(as.factor(GPSC), odds_ratio))

# Create forest plot
ggplot(results, aes(x = odds_ratio, y = GPSC)) +
  geom_point(aes(color = p_adj < 0.05), size = 3) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.3) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_log10() +  # log scale for odds ratios
  labs(
    x = "Odds Ratio (log scale)",
    y = "GPSC",
    title = "Forest Plot of GPSC Changes Post-Vaccine",
    subtitle = "Reference: GPSC53"
  ) +
  theme_minimal() +
  theme(legend.position = "top") +
  scale_color_manual(values = c("grey40", "red"), name = "Significant (FDR < 0.05)")



results <- results %>%
  filter(!is.na(odds_ratio)) %>%
  mutate(GPSC = fct_reorder(GPSC, odds_ratio))

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


forest_plot <- ggplot(results, aes(x = odds_ratio, y = GPSC)) +
  geom_point(aes(color = vaccine_type), size = 3) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.3) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_log10() +
  labs(
    x = "Odds Ratio (log scale)",
    y = "GPSC",
    title = "Changes Post-PCV10 in Reference to GPSC53"
  ) +
  theme_minimal() +
  theme(legend.position = "top right") +
  scale_color_manual(values = c("VT" = "darkgray", "NVT" = "steelblue", "Mixed" = "green")) +
  geom_text(aes(x = ci_upper * 1.2, label = stars), color = "maroon", size = 5)+
  theme(
    legend.position = "top right",
    axis.title.x = element_text(size = 16, face = "bold", color = "black"),
    axis.title.y = element_text(size = 16, face = "bold", color = "black"),
    axis.text = element_text(size = 18,  face = "bold", color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold", color = "black", size = 16),
    plot.subtitle = element_text(hjust = 0.5, face = "plain", color = "black", size = 13)
  )

forest_plot2 <- forest_plot + coord_flip()
forest_plot2


# Load in data
#data <- read_csv("metadata_study1.csv")
data <- read_excel("./disesea.xlsx")
serotype_colours <- read_csv("serotype_colours.csv")
colnames(serotype_colours)[1] <- "Serotype"
colnames(serotype_colours)[2] <- "Serotype__color"

# Tidy the serotypes by removing leading zeros and other fun things
data$Serotype <- sub("^0+", "", data$Serotype)

## Some serotype tidying - this will be unique to your data. 
data$Serotype <- ifelse(data$Serotype == "15B" | data$Serotype == "15C", "15B/15C", data$Serotype) # e.g. make all 15B or 15C "15B/15C"
data$Serotype <- ifelse(grepl("19A", data$Serotype), "19A", data$Serotype) # e.g. make all 19A(I), 19A(II) just 19A
data$Serotype <- ifelse(grepl("23B", data$Serotype), "23B", data$Serotype)  # e.g. make all 23B, 23B1 just 23B
data$Serotype <- ifelse(grepl("6A", data$Serotype), "6A", data$Serotype)  # e.g. make all 6A(I), 6A(II) just 6A

# Add the hex code colours for the serotype and GPSC
data <- left_join(data, serotype_colours, by = "Serotype")

# Create a named vector for the colours 
serotype_colors_vector <- setNames(data$Serotype__color, data$Serotype)

# Remove sample with NA vaccination period (just one)
data <- filter(data, !is.na(data$Vaccine_Period))

# Filter data to remove samples belonging to GPSCs with a count below 4 
counts <- table(data$GPSC)
aboveFour <- names(counts[counts > 1])
data_subset <- filter(data, data$GPSC %in% aboveFour)

# Count how many times each GPSC appears in `data`
gpsc_order <- data %>%
  count(GPSC) %>%
  arrange(desc(n))

data <- data %>%
  mutate(GPSC = factor(GPSC, levels = gpsc_order$GPSC))


# Generate all combinations of year and serotype, so there are no missing bars 
all_combinations <- expand.grid(
  GPSC = unique(data_subset$GPSC),
  Vaccine_Period = unique(data_subset$Vaccine_Period),
  Serotype = unique(data_subset$Serotype)
)

# Flip the data into a long format
long_data <- all_combinations %>%
  left_join(data_subset %>%
              group_by(GPSC, Vaccine_Period, Serotype) %>%
              summarise(count = n(), .groups = "drop"),
            by = c("GPSC", "Vaccine_Period", "Serotype")) %>%
  mutate(count = replace_na(count, 0))

# Add in the VT or NVT designation for each serotype. This is for PCV13, but change the serotypes and rename the column for the PCV10 formulation you use
# If you change the name from PCV13, you'll need to change this throughout the script. Control+F and replace :D
long_data$Vaccine_Category <- ifelse(long_data$Serotype == "4" | 
                           long_data$Serotype == "6B" | 
                           long_data$Serotype == "14" | 
                           long_data$Serotype == "18C" | 
                           long_data$Serotype == "9V" | 
                           long_data$Serotype == "19F" | 
                           long_data$Serotype == "23F" | 
                           long_data$Serotype == "1" | 
                           long_data$Serotype == "5" | 
                           long_data$Serotype == "7F" , "VT", "NVT")

## GPSC Plot
long_data <- long_data %>%
  mutate(GPSC = factor(GPSC, levels = gpsc_order$GPSC))

## Filter VT and NVT data - this is to try and make the legend behave.  
VT_data <- filter(long_data , long_data$Vaccine_Category != "NVT")
NVT_data <- filter(long_data , long_data$Vaccine_Category == "NVT")

NVT_labels <- unique(NVT_data$Serotype)
NVT_labels <- setNames(NVT_labels, NVT_labels)

VT_labels <- unique(VT_data$Serotype)
VT_labels <- setNames(VT_labels, VT_labels)

VT_serotypes <- unique(VT_data$Serotype)
NVT_serotypes <- unique(NVT_data$Serotype)

list <- c("NVT" = "stripe", "VT" = "none")
gpsc_order <- c(67,862, 62,8, 54,21,2,9, 884, 170, 27,5, 117,10,879,26,233,251,65,355,22,891,268,53)
long_data <- long_data[long_data$GPSC %in% gpsc_order, ]

# Make sure GPSC is treated as a factor in correct order
long_data$GPSC <- factor(long_data$GPSC, levels = gpsc_order)

# vaccine period serotypes
GPSC_Plot <- ggplot() + 
  geom_col_pattern(data = long_data, pattern_key_scale_factor = 0.2, 
                   aes(x = factor(Vaccine_Period, levels=c("pre", "post")), 
                       y = count, fill = Serotype, pattern = Vaccine_Category), color = "black", 
                   pattern_density = 0.2, pattern_spacing = 0.2)+ 
  facet_grid(~ GPSC,switch = "y") +
  #facet_grid(~ GPSC,switch = "x") + #to make lineages on y axis
  theme(axis.text.x = element_text(angle = 90, 
                                   vjust = 0.1, 
                                   hjust = 1, 
                                   size = 18,
                                   face = "bold",
                                   color = "black"),
        axis.text.y = element_text(size = 25, face = "bold", color = "black"),
        axis.title.x = element_text(size = 25,face = "bold", color = "black"),
        axis.title.y = element_text(size = 25,face = "bold", color = "black"),
        legend.text = element_text(size = 20,face = "bold", color = "black"),
        legend.title = element_text(size = 25,face = "bold", color = "black"),
        plot.caption = element_text(size = 25)) +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.background = element_rect(fill = NA, color = "white"),
        strip.text = element_blank(),
        panel.spacing = unit(-.01, "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.key.size = unit(0.8, "cm")) +
  ylab("Number of Isolates (n)") + 
  xlab(element_blank()) +
  scale_y_continuous(expand = c(0,0)) +
  guides(fill = guide_legend(ncol = 2)) +
  scale_fill_manual(aesthetics = "fill", 
                    values = serotype_colors_vector, 
                    labels = VT_labels,
                    breaks = VT_serotypes, 
                    name = " VT:",
                    guide = guide_legend(order = 1)) +
  scale_pattern_fill_manual(values = list) +
  new_scale_fill() +
  geom_col_pattern(data = long_data, pattern_key_scale_factor = 0.2, aes(x = factor(Vaccine_Period, 
                                                                                    levels=c("pre", "post")), 
                                                                         y = count, fill = Serotype, pattern = Vaccine_Category),
                   color = "black", pattern_density = 0.5, pattern_spacing = 0.2) + 
  scale_fill_manual(aesthetics = "fill", 
                    values = serotype_colors_vector, 
                    labels = NVT_labels,
                    breaks = NVT_serotypes, 
                    name = "NVT:",
                    guide = guide_legend(order = 0)) +
  theme(legend.position = "right", 
        legend.direction = "horizontal",
        legend.key = element_rect(fill = "black")) +
  guides(fill = guide_legend(ncol = 2)) +
  scale_pattern_manual(values = c(
    NVT = "stripe", 
    VT= "none")
    )


#ggsave("Figureall.png", plot = GPSC_Plot, width = 30, height = 20, dpi = 300)
ggsave("Figure.png", plot = GPSC_Plot, width = 30, height = 20, dpi = 300)


large_plot <- ( forest_plot2/GPSC_Plot) + plot_layout(heights = c(2,6)) 

# This will take a good amount of time!
ggsave("Figure3.png", plot = large_plot, width = 30, height = 20, dpi = 300)

