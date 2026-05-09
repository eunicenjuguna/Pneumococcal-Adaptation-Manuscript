############################
# LIBRARIES
############################
library(tidyverse)
library(readxl)
library(writexl)
library(ggplot2)
library(forcats)

############################
# LOAD DATA
############################
data <- read_excel("./disesea2.xlsx")

############################
# SUMMARISE RESISTANCE
############################
AMR_data <- data %>%
  group_by(GPSC, Vaccine_Period) %>%
  summarise(
    TMP = sum(TMP == "R", na.rm = TRUE),
    SMX = sum(SMX == "R", na.rm = TRUE),
    COT = sum(COT == "R", na.rm = TRUE),
    PEN = sum(PEN == "R", na.rm = TRUE),
    CFX = sum(CFX == "R", na.rm = TRUE),
    TET = sum(TET == "R", na.rm = TRUE),
    DOX = sum(DOX == "R", na.rm = TRUE),
    ERY = sum(ERY == "R", na.rm = TRUE),
    CHL = sum(CHL == "R", na.rm = TRUE),
    CLI = sum(CLI == "R", na.rm = TRUE),
    RIF = sum(RIF == "R", na.rm = TRUE),
    Total = n(),
    .groups = "drop"
  )

############################
# REPLACE NA & CALCULATE %
############################
AMR_data <- AMR_data %>%
  mutate(across(c(TMP, SMX, COT, PEN, CFX, TET, DOX, ERY, CHL, CLI, RIF, Total),
                ~replace_na(., 0))) %>%
  mutate(across(c(TMP, SMX, COT, PEN, CFX, TET, DOX, ERY, CHL, CLI, RIF),
                ~ifelse(Total > 0, (. / Total) * 100, 0)))

############################
# RENAME ANTIBIOTICS TO A–J
############################
AMR_data <- AMR_data %>%
  rename(
    A = PEN,
    B = CFX,
    C = TET,
    D = DOX,
    E = CLI,
    F = ERY,
    G = CHL,
    H = COT,
    I = SMX,
    J = TMP
  )

############################
# FILTER: >5 ANTIBIOTICS (≥6)
############################
MDR_5plus <- AMR_data %>%
  mutate(
    n_resistant =
      (A > 0) + (B > 0) + (C > 0) + (D > 0) + (E > 0) +
      (F > 0) + (G > 0) + (H > 0) + (I > 0) + (J > 0)
  ) %>%
  filter(n_resistant > 5) %>%
  select(GPSC, Vaccine_Period)

############################
# LONG FORMAT
############################
long_AMR <- AMR_data %>%
  pivot_longer(
    cols = c(A,B,C,D,E,F,G,H,I,J),
    names_to = "Antibiotic",
    values_to = "Percentage"
  ) %>%
  semi_join(MDR_5plus, by = c("GPSC", "Vaccine_Period"))

############################
# SAFETY CHECK
############################
if (nrow(long_AMR) == 0) {
  stop("No GPSC–Vaccine_Period combinations resistant to >5 antibiotics.")
}

############################
# GPSC ORDER (EDIT IF NEEDED)
############################
gpsc_order <- sort(unique(long_AMR$GPSC))

long_AMR <- long_AMR %>%
  mutate(GPSC = factor(GPSC, levels = gpsc_order),
         Antibiotic = factor(Antibiotic, levels = rev(c("A","B","C","D","E","F","G","H","I","J"))))

############################
# PLOT
############################
AMR_Plot <- ggplot(
  long_AMR,
  aes(
    x = factor(Vaccine_Period, levels = c("pre", "post")),
    y = Antibiotic,
    fill = Percentage
  )
) +
  geom_tile(color = "black") +
  facet_grid(~ GPSC, switch = "x") +
  scale_fill_gradient(
    low = "white",
    high = "red",
    limits = c(0, 100),
    name = "% Resistance"
  ) +
  labs(
    x = "GPSC",
    y = "Antibiotics"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
    axis.text.y = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    panel.spacing = unit(-0.01, "cm")
  )

print(AMR_Plot)
