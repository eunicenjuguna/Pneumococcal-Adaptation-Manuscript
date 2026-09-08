# ============================================================
# GPSC-REFERENCED FISHER'S EXACT TEST


# + 0.5 HALDANE-ANScombe CORRECTED OR / 95% CI
# + PUBLICATION-STYLE GPSC PREVALENCE + FOREST PLOT
#
# REFERENCES:
#   Invasive disease = GPSC 53
#   Carriage         = GPSC 233
#
# STATISTICAL ANALYSIS:
#   Fisher's exact test uses RAW counts
#   BH correction for multiple testing
#   0.5 Haldane-Anscombe correction applied ONLY
#   to OR and 95% CI
#
# PLOTTING:
#   VT serotypes  = solid colours
#   NVT serotypes = striped colours
#
# FOREST PLOT:
#   VT GPSC       = grey
#   NVT GPSC      = blue
#   Mixed GPSC    = green
#
# IMPORTANT:
#   A GPSC is classified as:
#     VT              = contains only VT serotypes
#     NVT             = contains only NVT serotypes
#     Mixed VT/NVT    = contains both
#
# ============================================================


# ============================================================
# 1. LOAD PACKAGES
# ============================================================

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(writexl)
library(scales)

# ggpattern is required for striped NVT bars
if (!requireNamespace("ggpattern", quietly = TRUE)) {
  install.packages("ggpattern")
}

library(ggpattern)


# ============================================================
# 2. LOAD DATA
# ============================================================

data <- read_excel(
  "./KLF_population_ex2011.xlsx"
)


# ============================================================
# 3. CLEAN VARIABLES
# ============================================================

data <- data %>%
  mutate(
    
    # GPSC as character
    GPSC = trimws(
      as.character(GPSC)
    ),
    
    # Remove possible "GPSC" prefix
    GPSC = gsub(
      "^GPSC[[:space:]]*",
      "",
      GPSC,
      ignore.case = TRUE
    ),
    
    # Serotype
    Serotype_merged = trimws(
      as.character(Serotype_merged)
    ),
    
    # Disease/carriage
    disease_carriage = trimws(
      as.character(disease_carriage)
    ),
    
    # Vaccine period
    vaccine_period = trimws(
      as.character(vaccine_period)
    )
  )


# ============================================================
# 4. STANDARDISE DISEASE/CARRIAGE LABELS
# ============================================================

data <- data %>%
  mutate(
    
    disease_carriage = case_when(
      
      tolower(disease_carriage) %in%
        c(
          "disease",
          "invasive",
          "invasive disease"
        ) ~ "disease",
      
      tolower(disease_carriage) %in%
        c(
          "carriage",
          "carrier"
        ) ~ "carriage",
      
      TRUE ~ tolower(disease_carriage)
    )
  )


# ============================================================
# 5. STANDARDISE PRE/POST LABELS
# ============================================================

data <- data %>%
  mutate(
    
    vaccine_period = case_when(
      
      tolower(vaccine_period) %in%
        c(
          "pre",
          "pre-pcv10",
          "pre_pcv10"
        ) ~ "pre",
      
      tolower(vaccine_period) %in%
        c(
          "post",
          "post-pcv10",
          "post_pcv10"
        ) ~ "post",
      
      TRUE ~ tolower(vaccine_period)
    )
  )


# ============================================================
# 6. RAW DATA CHECK
# ============================================================

cat("\n")
cat("============================================================\n")
cat("RAW DATA CHECK\n")
cat("============================================================\n\n")

cat("Dataset dimensions:\n")
print(dim(data))

cat("\nDisease/carriage:\n")
print(
  table(
    data$disease_carriage,
    useNA = "always"
  )
)

cat("\nVaccine period:\n")
print(
  table(
    data$vaccine_period,
    useNA = "always"
  )
)

cat("\nNumber of unique GPSCs:\n")
print(
  length(
    unique(
      na.omit(data$GPSC)
    )
  )
)

cat("\nFirst GPSCs:\n")
print(
  head(
    sort(
      unique(
        na.omit(data$GPSC)
      )
    ),
    30
  )
)


# ============================================================
# 7. DEFINE PCV10 VACCINE-TYPE SEROTYPES
# ============================================================

VT_serotypes <- c(
  "1",
  "4",
  "5",
  "6B",
  "9V",
  "7F",
  "14",
  "18C",
  "19F",
  "23F"
)


# ============================================================
# 8. DEFINE VACCINE TYPE FOR EACH ISOLATE
# ============================================================

data <- data %>%
  mutate(
    
    vaccine_type = case_when(
      
      Serotype_merged %in%
        VT_serotypes ~ "VT",
      
      TRUE ~ "NVT"
    )
  )


# ============================================================
# 9. CHECK REFERENCE GPSCs
# ============================================================

cat("\n")
cat("============================================================\n")
cat("REFERENCE GPSC CHECK\n")
cat("============================================================\n\n")

cat("GPSC 53 — invasive disease:\n")

print(
  data %>%
    filter(
      GPSC == "53",
      disease_carriage == "disease"
    ) %>%
    count(
      vaccine_period
    )
)


cat("\nGPSC 233 — carriage:\n")

print(
  data %>%
    filter(
      GPSC == "233",
      disease_carriage == "carriage"
    ) %>%
    count(
      vaccine_period
    )
)


# ============================================================
# 10. FISHER'S EXACT TEST FUNCTION
# ============================================================

run_raw_fisher <- function(
    pre_gpsc,
    post_gpsc,
    pre_reference,
    post_reference
) {
  
  mat <- matrix(
    c(
      pre_gpsc,
      post_gpsc,
      pre_reference,
      post_reference
    ),
    nrow = 2,
    byrow = TRUE
  )
  
  rownames(mat) <- c(
    "GPSC",
    "Reference"
  )
  
  colnames(mat) <- c(
    "Pre",
    "Post"
  )
  
  ft <- fisher.test(
    mat
  )
  
  tibble(
    
    fisher_OR =
      unname(
        ft$estimate
      ),
    
    fisher_CI_lower =
      ft$conf.int[1],
    
    fisher_CI_upper =
      ft$conf.int[2],
    
    fisher_p =
      ft$p.value
  )
}


# ============================================================
# 11. HALDANE-ANScombe 0.5 CORRECTION
# ============================================================

calculate_0.5_corrected <- function(
    pre_gpsc,
    post_gpsc,
    pre_reference,
    post_reference
) {
  
  # Four cells
  a <- pre_gpsc + 0.5
  b <- post_gpsc + 0.5
  c <- pre_reference + 0.5
  d <- post_reference + 0.5
  
  
  # Corrected OR
  OR <- (
    b * c
  ) / (
    a * d
  )
  
  
  # Standard error
  SE <- sqrt(
    1 / a +
      1 / b +
      1 / c +
      1 / d
  )
  
  
  # Log OR
  log_OR <- log(
    OR
  )
  
  
  # 95% CI
  CI_lower <- exp(
    log_OR -
      1.96 * SE
  )
  
  CI_upper <- exp(
    log_OR +
      1.96 * SE
  )
  
  
  tibble(
    
    corrected_OR =
      OR,
    
    corrected_CI_lower =
      CI_lower,
    
    corrected_CI_upper =
      CI_upper
  )
}


# ============================================================
# 12. GPSC ANALYSIS FUNCTION
# ============================================================

run_GPSC_analysis <- function(
    dataset,
    group_name,
    reference_gpsc,
    minimum_total
) {
  
  
  cat("\n")
  cat("============================================================\n")
  cat(
    toupper(group_name),
    "GPSC ANALYSIS\n"
  )
  cat("============================================================\n\n")
  
  
  cat(
    "Reference GPSC:",
    reference_gpsc,
    "\n"
  )
  
  cat(
    "Eligibility threshold:",
    minimum_total,
    "isolates\n\n"
  )
  
  
  # ----------------------------------------------------------
  # 12A. SUBSET DATA
  # ----------------------------------------------------------
  
  group_data <- dataset %>%
    filter(
      disease_carriage == group_name,
      !is.na(GPSC),
      GPSC != "",
      vaccine_period %in%
        c(
          "pre",
          "post"
        )
    )
  
  
  # ----------------------------------------------------------
  # 12B. TOTAL GPSC COUNTS
  # ----------------------------------------------------------
  
  totals <- group_data %>%
    count(
      GPSC,
      name = "Total"
    )
  
  
  # ----------------------------------------------------------
  # 12C. ELIGIBLE GPSCs
  #
  # Reference is always retained.
  # ----------------------------------------------------------
  
  eligible <- totals %>%
    filter(
      Total >= minimum_total |
        GPSC == reference_gpsc
    )
  
  
  cat("Number of eligible GPSCs:")
  cat(
    nrow(eligible),
    "\n\n"
  )
  
  
  cat("Eligible GPSCs:\n")
  
  print(
    eligible,
    n = Inf
  )
  
  
  # ----------------------------------------------------------
  # 12D. PRE/POST COUNTS
  # ----------------------------------------------------------
  
  counts <- group_data %>%
    filter(
      GPSC %in%
        eligible$GPSC
    ) %>%
    count(
      GPSC,
      vaccine_period
    ) %>%
    pivot_wider(
      names_from =
        vaccine_period,
      values_from =
        n,
      values_fill =
        0
    )
  
  
  # Ensure pre exists
  if (!"pre" %in% names(counts)) {
    counts$pre <- 0
  }
  
  
  # Ensure post exists
  if (!"post" %in% names(counts)) {
    counts$post <- 0
  }
  
  
  # Make sure every eligible GPSC is retained
  counts <- eligible %>%
    select(
      GPSC
    ) %>%
    left_join(
      counts,
      by = "GPSC"
    ) %>%
    mutate(
      
      pre = replace_na(
        pre,
        0
      ),
      
      post = replace_na(
        post,
        0
      )
    )
  
  
  # ----------------------------------------------------------
  # 12E. REFERENCE COUNTS
  # ----------------------------------------------------------
  
  reference_counts <- group_data %>%
    filter(
      GPSC == reference_gpsc
    ) %>%
    count(
      vaccine_period
    ) %>%
    pivot_wider(
      names_from =
        vaccine_period,
      values_from =
        n,
      values_fill =
        0
    )
  
  
  if (!"pre" %in% names(reference_counts)) {
    reference_counts$pre <- 0
  }
  
  
  if (!"post" %in% names(reference_counts)) {
    reference_counts$post <- 0
  }
  
  
  reference_pre <-
    reference_counts$pre[1]
  
  reference_post <-
    reference_counts$post[1]
  
  
  cat("\nReference GPSC counts:\n")
  
  cat(
    "GPSC",
    reference_gpsc,
    "Pre =",
    reference_pre,
    "\n"
  )
  
  cat(
    "GPSC",
    reference_gpsc,
    "Post =",
    reference_post,
    "\n"
  )
  
  
  # ----------------------------------------------------------
  # 12F. RUN FISHER FOR EACH GPSC
  # ----------------------------------------------------------
  
  results_list <- list()
  
  
  test_gpscs <- counts %>%
    filter(
      GPSC != reference_gpsc
    )
  
  
  for (
    i in seq_len(
      nrow(test_gpscs)
    )
  ) {
    
    current_gpsc <-
      test_gpscs$GPSC[i]
    
    current_pre <-
      test_gpscs$pre[i]
    
    current_post <-
      test_gpscs$post[i]
    
    
    fisher_result <-
      run_raw_fisher(
        
        pre_gpsc =
          current_pre,
        
        post_gpsc =
          current_post,
        
        pre_reference =
          reference_pre,
        
        post_reference =
          reference_post
      )
    
    
    corrected_result <-
      calculate_0.5_corrected(
        
        pre_gpsc =
          current_pre,
        
        post_gpsc =
          current_post,
        
        pre_reference =
          reference_pre,
        
        post_reference =
          reference_post
      )
    
    
    results_list[[i]] <-
      tibble(
        
        GPSC =
          current_gpsc,
        
        Pre =
          current_pre,
        
        Post =
          current_post
        
      ) %>%
      
      bind_cols(
        fisher_result
      ) %>%
      
      bind_cols(
        corrected_result
      )
  }
  
  
  results <- bind_rows(
    results_list
  )
  
  
  # ----------------------------------------------------------
  # 12G. ADD TOTAL
  # ----------------------------------------------------------
  
  results <- results %>%
    left_join(
      totals,
      by = "GPSC"
    )
  
  
  # ----------------------------------------------------------
  # 12H. DETERMINE GPSC VT/NVT STATUS
  # ----------------------------------------------------------
  
  gpsc_status <- group_data %>%
    filter(
      GPSC %in%
        eligible$GPSC
    ) %>%
    group_by(
      GPSC
    ) %>%
    summarise(
      
      has_VT =
        any(
          vaccine_type == "VT",
          na.rm = TRUE
        ),
      
      has_NVT =
        any(
          vaccine_type == "NVT",
          na.rm = TRUE
        ),
      
      GPSC_type =
        case_when(
          
          has_VT &
            has_NVT ~
            "Mixed VT/NVT-GPSC",
          
          has_VT ~
            "VT-GPSC",
          
          has_NVT ~
            "NVT-GPSC",
          
          TRUE ~
            NA_character_
        ),
      
      .groups = "drop"
    )
  
  
  # ----------------------------------------------------------
  # 12I. ADD GPSC STATUS
  # ----------------------------------------------------------
  
  results <- results %>%
    left_join(
      gpsc_status,
      by = "GPSC"
    )
  
  
  # ----------------------------------------------------------
  # 12J. BH MULTIPLE TESTING CORRECTION
  # ----------------------------------------------------------
  
  results <- results %>%
    mutate(
      
      p_adj =
        p.adjust(
          fisher_p,
          method = "BH"
        ),
      
      Significant =
        p_adj < 0.05,
      
      stars =
        case_when(
          
          p_adj < 0.001 ~
            "***",
          
          p_adj < 0.01 ~
            "**",
          
          p_adj < 0.05 ~
            "*",
          
          TRUE ~
            ""
        )
    )
  
  
  # ----------------------------------------------------------
  # 12K. CI EXCLUDES 1
  # ----------------------------------------------------------
  
  results <- results %>%
    mutate(
      
      Corrected_CI_excludes_1 =
        corrected_CI_lower > 1 |
        corrected_CI_upper < 1,
      
      Direction =
        case_when(
          
          corrected_OR > 1 ~
            "Increased",
          
          corrected_OR < 1 ~
            "Decreased",
          
          TRUE ~
            "No change"
        )
    )
  
  
  # ----------------------------------------------------------
  # 12L. FINAL COLUMN ORDER
  # ----------------------------------------------------------
  
  results <- results %>%
    select(
      
      GPSC,
      GPSC_type,
      
      Pre,
      Post,
      Total,
      
      fisher_OR,
      fisher_CI_lower,
      fisher_CI_upper,
      fisher_p,
      
      p_adj,
      
      corrected_OR,
      corrected_CI_lower,
      corrected_CI_upper,
      
      Significant,
      stars,
      
      Corrected_CI_excludes_1,
      
      Direction
    ) %>%
    
    arrange(
      fisher_p
    )
  
  
  # ----------------------------------------------------------
  # 12M. RETURN RESULTS
  # ----------------------------------------------------------
  
  return(
    list(
      
      results =
        results,
      
      eligible =
        eligible,
      
      counts =
        counts,
      
      status =
        gpsc_status,
      
      reference_pre =
        reference_pre,
      
      reference_post =
        reference_post
    )
  )
}


# ============================================================
# 13. INVASIVE DISEASE ANALYSIS
#
# GPSC 53 = reference
# >=5 isolates eligible
# ============================================================

disease_analysis <-
  run_GPSC_analysis(
    
    dataset =
      data,
    
    group_name =
      "disease",
    
    reference_gpsc =
      "53",
    
    minimum_total =
      5
  )


disease_results <-
  disease_analysis$results

disease_eligible <-
  disease_analysis$eligible

disease_counts <-
  disease_analysis$counts


# ============================================================
# 14. CARRIAGE ANALYSIS
#
# GPSC 233 = reference
# >=10 isolates eligible
# ============================================================

carriage_analysis <-
  run_GPSC_analysis(
    
    dataset =
      data,
    
    group_name =
      "carriage",
    
    reference_gpsc =
      "233",
    
    minimum_total =
      10
  )


carriage_results <-
  carriage_analysis$results

carriage_eligible <-
  carriage_analysis$eligible

carriage_counts <-
  carriage_analysis$counts


# ============================================================
# 15. PRINT FINAL RESULTS
# ============================================================

cat("\n\n")
cat("============================================================\n")
cat("INVASIVE DISEASE — RESULTS VS GPSC 53\n")
cat("============================================================\n\n")

print(
  disease_results,
  n = Inf
)


cat("\n\n")
cat("============================================================\n")
cat("CARRIAGE — RESULTS VS GPSC 233\n")
cat("============================================================\n\n")

print(
  carriage_results,
  n = Inf
)


# ============================================================
# 16. PRINT SIGNIFICANT RESULTS
# ============================================================

cat("\n\n")
cat("============================================================\n")
cat("SIGNIFICANT INVASIVE DISEASE GPSCs\n")
cat("============================================================\n\n")

print(
  disease_results %>%
    filter(
      Significant == TRUE
    ),
  n = Inf
)


cat("\n\n")
cat("============================================================\n")
cat("SIGNIFICANT CARRIAGE GPSCs\n")
cat("============================================================\n\n")

print(
  carriage_results %>%
    filter(
      Significant == TRUE
    ),
  n = Inf
)


# ============================================================
# 17. SAVE STATISTICAL RESULTS
# ============================================================

write_xlsx(
  list(
    
    Results =
      disease_results,
    
    Eligible_GPSCs =
      disease_eligible,
    
    GPSC_VT_NVT_status =
      disease_analysis$status
    
  ),
  "disease_GPSC_reference_results.xlsx"
)


write_xlsx(
  list(
    
    Results =
      carriage_results,
    
    Eligible_GPSCs =
      carriage_eligible,
    
    GPSC_VT_NVT_status =
      carriage_analysis$status
    
  ),
  "carriage_GPSC_reference_results.xlsx"
)


# ============================================================
# 17B. PRE-PCV10 NVT REFERENCE FREQUENCY TABLES
# ============================================================
#
# These tables document the choice of reference categories:
#   Invasive disease lineage = GPSC 53
#   Carriage lineage         = GPSC 233
#   NVT serotype             = 6A
#
# GPSC classifications use all observations within each population.
# Frequencies themselves use pre-PCV10 observations only. Tables are
# ordered from the lowest to the highest pre-PCV10 frequency, as
# requested; Rank = 1 identifies the most frequent category.

all_gpsc_status <- data %>%
  filter(
    disease_carriage %in% c("disease", "carriage"),
    !is.na(GPSC),
    GPSC != ""
  ) %>%
  group_by(disease_carriage, GPSC) %>%
  summarise(
    has_VT = any(vaccine_type == "VT", na.rm = TRUE),
    has_NVT = any(vaccine_type == "NVT", na.rm = TRUE),
    GPSC_type = case_when(
      has_VT & has_NVT ~ "Mixed VT/NVT-GPSC",
      has_VT ~ "VT-GPSC",
      has_NVT ~ "NVT-GPSC",
      TRUE ~ NA_character_
    ),
    .groups = "drop"
  )

pre_nvt_gpsc_frequency <- data %>%
  filter(
    disease_carriage %in% c("disease", "carriage"),
    vaccine_period == "pre",
    !is.na(GPSC),
    GPSC != ""
  ) %>%
  count(disease_carriage, GPSC, name = "Pre_count") %>%
  left_join(
    all_gpsc_status %>%
      select(disease_carriage, GPSC, GPSC_type),
    by = c("disease_carriage", "GPSC")
  ) %>%
  filter(GPSC_type == "NVT-GPSC") %>%
  group_by(disease_carriage) %>%
  mutate(
    Pre_NVT_GPSC_total = sum(Pre_count),
    Pre_percent = 100 * Pre_count / Pre_NVT_GPSC_total,
    Rank = min_rank(desc(Pre_count)),
    Reference = case_when(
      disease_carriage == "disease" & GPSC == "53" ~ "Yes",
      disease_carriage == "carriage" & GPSC == "233" ~ "Yes",
      TRUE ~ "No"
    )
  ) %>%
  ungroup() %>%
  mutate(
    Population = recode(
      disease_carriage,
      disease = "Invasive disease",
      carriage = "Carriage"
    ),
    Pre_percent = round(Pre_percent, 2)
  ) %>%
  select(
    Population,
    GPSC,
    GPSC_type,
    Pre_count,
    Pre_NVT_GPSC_total,
    Pre_percent,
    Rank,
    Reference
  ) %>%
  arrange(Population, Pre_count, GPSC)

pre_nvt_serotype_frequency <- data %>%
  filter(
    disease_carriage %in% c("disease", "carriage"),
    vaccine_period == "pre",
    vaccine_type == "NVT",
    !is.na(Serotype_merged),
    Serotype_merged != ""
  ) %>%
  count(disease_carriage, Serotype_merged, name = "Pre_count") %>%
  group_by(disease_carriage) %>%
  mutate(
    Pre_NVT_serotype_total = sum(Pre_count),
    Pre_percent = 100 * Pre_count / Pre_NVT_serotype_total,
    Rank = min_rank(desc(Pre_count)),
    Reference = if_else(Serotype_merged == "6A", "Yes", "No")
  ) %>%
  ungroup() %>%
  mutate(
    Population = recode(
      disease_carriage,
      disease = "Invasive disease",
      carriage = "Carriage"
    ),
    Pre_percent = round(Pre_percent, 2)
  ) %>%
  select(
    Population,
    Serotype = Serotype_merged,
    Pre_count,
    Pre_NVT_serotype_total,
    Pre_percent,
    Rank,
    Reference
  ) %>%
  arrange(Population, Pre_count, Serotype)

all_reference_frequencies <- bind_rows(
  pre_nvt_gpsc_frequency %>%
    transmute(
      Population,
      Category = "NVT-only GPSC lineage",
      Category_value = paste0("GPSC ", GPSC),
      Pre_count,
      Denominator = Pre_NVT_GPSC_total,
      Pre_percent,
      Rank,
      Selected_reference = Reference
    ),
  pre_nvt_serotype_frequency %>%
    transmute(
      Population,
      Category = "NVT serotype",
      Category_value = Serotype,
      Pre_count,
      Denominator = Pre_NVT_serotype_total,
      Pre_percent,
      Rank,
      Selected_reference = Reference
    )
) %>%
  arrange(Category, Population, Pre_count, Category_value)

reference_frequency_summary <- all_reference_frequencies %>%
  filter(Selected_reference == "Yes") %>%
  rename(Reference = Category_value) %>%
  select(-Selected_reference)

# The selected references must be the most frequent pre-PCV10
# categories in their respective population and category tables.
if (nrow(reference_frequency_summary) != 4 ||
    any(reference_frequency_summary$Rank != 1)) {
  stop(
    "One or more selected GPSC/serotype references are not ranked first ",
    "in the pre-PCV10 NVT frequency tables."
  )
}

write_xlsx(
  list(
    Pre_PCV10_NVT_frequency = {
      invasive_serotypes <- pre_nvt_serotype_frequency %>%
        filter(Population == "Invasive disease") %>%
        arrange(Pre_count, Serotype)

      carriage_serotypes <- pre_nvt_serotype_frequency %>%
        filter(Population == "Carriage") %>%
        arrange(Pre_count, Serotype)

      invasive_gpscs <- pre_nvt_gpsc_frequency %>%
        filter(Population == "Invasive disease") %>%
        arrange(Pre_count, GPSC)

      carriage_gpscs <- pre_nvt_gpsc_frequency %>%
        filter(Population == "Carriage") %>%
        arrange(Pre_count, GPSC)

      table_rows <- max(
        nrow(invasive_serotypes),
        nrow(carriage_serotypes),
        nrow(invasive_gpscs),
        nrow(carriage_gpscs)
      )

      pad_character <- function(x) {
        c(as.character(x), rep(NA_character_, table_rows - length(x)))
      }

      pad_numeric <- function(x) {
        c(as.numeric(x), rep(NA_real_, table_rows - length(x)))
      }

      tibble(
        Invasive_serotype = pad_character(invasive_serotypes$Serotype),
        Invasive_serotype_n = pad_numeric(invasive_serotypes$Pre_count),
        Carriage_serotype = pad_character(carriage_serotypes$Serotype),
        Carriage_serotype_n = pad_numeric(carriage_serotypes$Pre_count),
        Invasive_GPSC = pad_character(invasive_gpscs$GPSC),
        Invasive_GPSC_n = pad_numeric(invasive_gpscs$Pre_count),
        Carriage_GPSC = pad_character(carriage_gpscs$GPSC),
        Carriage_GPSC_n = pad_numeric(carriage_gpscs$Pre_count)
      )
    }
  ),
  "Pre_PCV10_NVT_reference_frequencies.xlsx"
)

# ============================================================
# 18. GPSC PLOTTING
# ============================================================
#
# REFERENCES:
#   Disease  = GPSC 53
#   Carriage = GPSC 233
#
# BAR PLOTS:
#   VT  = plain/solid colour
#   NVT = striped colour
#
#   Mixed GPSC = both solid and striped segments
#
# FOREST PLOTS:
#   VT-GPSC           = grey
#   NVT-GPSC          = blue
#   Mixed VT/NVT-GPSC = green
#
# ODDS RATIOS:
#   corrected_OR
#   corrected_CI_lower
#   corrected_CI_upper
#
# ============================================================


# ============================================================
# 18. REQUIRED PACKAGES
# ============================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(ggpattern)
library(scales)


# ============================================================
# 19. GET ELIGIBLE GPSCs
# ============================================================
#
# These are taken directly from the statistical analysis.
# Disease: >=5 isolates, plus reference GPSC 53
# Carriage: >=10 isolates, plus reference GPSC 233
#
# ============================================================

disease_gpscs <- disease_eligible %>%
  mutate(
    GPSC = as.character(GPSC)
  ) %>%
  pull(GPSC)


carriage_gpscs <- carriage_eligible %>%
  mutate(
    GPSC = as.character(GPSC)
  ) %>%
  pull(GPSC)


# ============================================================
# 20. CREATE GPSC ORDER
# ============================================================
#
# Order GPSCs by corrected OR.
#
# Reference is deliberately placed at the bottom.
#
# ============================================================

disease_order <- disease_results %>%
  filter(
    GPSC != "53"
  ) %>%
  arrange(
    corrected_OR
  ) %>%
  pull(
    GPSC
  ) %>%
  c(
    "53"
  )


carriage_order <- carriage_results %>%
  filter(
    GPSC != "233"
  ) %>%
  arrange(
    corrected_OR
  ) %>%
  pull(
    GPSC
  ) %>%
  c(
    "233"
  )


# Make sure only eligible GPSCs are included
disease_order <- disease_order[
  disease_order %in% disease_gpscs
]

carriage_order <- carriage_order[
  carriage_order %in% carriage_gpscs
]


# ============================================================
# 21. CREATE NUMERIC GPSC POSITIONS
# ============================================================

disease_y <- tibble(
  
  GPSC = disease_order,
  
  GPSC_position = seq_along(
    disease_order
  )
)


carriage_y <- tibble(
  
  GPSC = carriage_order,
  
  GPSC_position = seq_along(
    carriage_order
  )
)

# ============================================================
# 22. CREATE GPSC VT/NVT STATUS
# ============================================================
#
# VT only           = grey
# NVT only          = blue
# Mixed VT + NVT    = green
#
# IMPORTANT:
# GPSC is explicitly converted to character on BOTH sides
# before joining.
# ============================================================


# -------------------------
# DISEASE
# -------------------------

disease_status <- data %>%
  
  mutate(
    GPSC = as.character(GPSC),
    vaccine_type = as.character(vaccine_type)
  ) %>%
  
  filter(
    disease_carriage == "disease",
    GPSC %in% disease_order
  ) %>%
  
  group_by(
    GPSC
  ) %>%
  
  summarise(
    
    has_VT =
      any(
        vaccine_type == "VT",
        na.rm = TRUE
      ),
    
    has_NVT =
      any(
        vaccine_type == "NVT",
        na.rm = TRUE
      ),
    
    .groups = "drop"
  ) %>%
  
  mutate(
    
    GPSC_type =
      case_when(
        
        has_VT & has_NVT ~
          "Mixed VT/NVT-GPSC",
        
        has_VT ~
          "VT-GPSC",
        
        has_NVT ~
          "NVT-GPSC",
        
        TRUE ~
          NA_character_
      )
  )


# -------------------------
# CARRIAGE
# -------------------------

carriage_status <- data %>%
  
  mutate(
    GPSC = as.character(GPSC),
    vaccine_type = as.character(vaccine_type)
  ) %>%
  
  filter(
    disease_carriage == "carriage",
    GPSC %in% carriage_order
  ) %>%
  
  group_by(
    GPSC
  ) %>%
  
  summarise(
    
    has_VT =
      any(
        vaccine_type == "VT",
        na.rm = TRUE
      ),
    
    has_NVT =
      any(
        vaccine_type == "NVT",
        na.rm = TRUE
      ),
    
    .groups = "drop"
  ) %>%
  
  mutate(
    
    GPSC_type =
      case_when(
        
        has_VT & has_NVT ~
          "Mixed VT/NVT-GPSC",
        
        has_VT ~
          "VT-GPSC",
        
        has_NVT ~
          "NVT-GPSC",
        
        TRUE ~
          NA_character_
      )
  )


# ============================================================
# 23. CHECK GPSC STATUS
# ============================================================

cat("\n================ DISEASE GPSC STATUS ================\n")

print(
  disease_status,
  n = Inf
)


cat("\n================ CARRIAGE GPSC STATUS ================\n")

print(
  carriage_status,
  n = Inf
)


# ============================================================
# 24. CREATE GPSC POSITIONS
# ============================================================

disease_y <- tibble(
  
  GPSC =
    as.character(
      disease_order
    ),
  
  GPSC_position =
    seq_along(
      disease_order
    )
)


carriage_y <- tibble(
  
  GPSC =
    as.character(
      carriage_order
    ),
  
  GPSC_position =
    seq_along(
      carriage_order
    )
)


# ============================================================
# 25. CREATE BAR DATA
# ============================================================

make_bar_data <- function(
    dataset,
    group_name,
    gpsc_order,
    gpsc_positions
) {
  
  
  dataset %>%
    
    mutate(
      GPSC =
        as.character(GPSC),
      
      vaccine_period =
        tolower(
          as.character(
            vaccine_period
          )
        ),
      
      vaccine_type =
        as.character(
          vaccine_type
        ),
      
      Serotype_merged =
        as.character(
          Serotype_merged
        )
    ) %>%
    
    filter(
      
      disease_carriage ==
        group_name,
      
      GPSC %in%
        gpsc_order,
      
      vaccine_period %in%
        c(
          "pre",
          "post"
        ),
      
      !is.na(
        Serotype_merged
      ),
      
      Serotype_merged != ""
    ) %>%
    
    count(
      
      GPSC,
      
      vaccine_period,
      
      Serotype_merged,
      
      vaccine_type,
      
      name =
        "Count"
    ) %>%
    
    left_join(
      
      gpsc_positions,
      
      by =
        "GPSC"
    ) %>%
    
    mutate(
      
      y =
        case_when(
          
          vaccine_period == "pre" ~
            GPSC_position - 0.18,
          
          vaccine_period == "post" ~
            GPSC_position + 0.18
        ),
      
      Pattern =
        case_when(
          
          vaccine_type == "VT" ~
            "none",
          
          vaccine_type == "NVT" ~
            "stripe",
          
          TRUE ~
            "none"
        )
    )
}


disease_bar_data <- make_bar_data(
  
  data,
  
  "disease",
  
  disease_order,
  
  disease_y
)


carriage_bar_data <- make_bar_data(
  
  data,
  
  "carriage",
  
  carriage_order,
  
  carriage_y
)

# ============================================================
# 26. CREATE FOREST DATA
# ============================================================
#
# IMPORTANT:
# disease_results may already contain GPSC_type.
# Therefore we call the newly calculated classification
# GPSC_type_new to avoid .x / .y conflicts.
#
# Forest colours:
#   VT only           = grey
#   NVT only          = blue
#   Mixed VT/NVT      = green
#
# ============================================================


disease_forest_data <- disease_results %>%
  
  mutate(
    GPSC = as.character(GPSC)
  ) %>%
  
  filter(
    GPSC %in% disease_order
  ) %>%
  
  left_join(
    
    disease_status %>%
      mutate(
        GPSC = as.character(GPSC)
      ) %>%
      select(
        GPSC,
        GPSC_type
      ) %>%
      rename(
        GPSC_type_new = GPSC_type
      ),
    
    by = "GPSC"
  ) %>%
  
  left_join(
    
    disease_y,
    
    by = "GPSC"
  ) %>%
  
  mutate(
    
    # Use the newly calculated GPSC classification
    GPSC_type = GPSC_type_new,
    
    # Reference GPSC 53 = OR 1
    corrected_OR =
      ifelse(
        GPSC == "53",
        1,
        corrected_OR
      ),
    
    corrected_CI_lower =
      ifelse(
        GPSC == "53",
        1,
        corrected_CI_lower
      ),
    
    corrected_CI_upper =
      ifelse(
        GPSC == "53",
        1,
        corrected_CI_upper
      )
  ) %>%
  
  select(
    -GPSC_type_new
  )


# ============================================================
# 27. CREATE CARRIAGE FOREST DATA
# ============================================================

carriage_forest_data <- carriage_results %>%
  
  mutate(
    GPSC = as.character(GPSC)
  ) %>%
  
  filter(
    GPSC %in% carriage_order
  ) %>%
  
  left_join(
    
    carriage_status %>%
      mutate(
        GPSC = as.character(GPSC)
      ) %>%
      select(
        GPSC,
        GPSC_type
      ) %>%
      rename(
        GPSC_type_new = GPSC_type
      ),
    
    by = "GPSC"
  ) %>%
  
  left_join(
    
    carriage_y,
    
    by = "GPSC"
  ) %>%
  
  mutate(
    
    # Use the newly calculated GPSC classification
    GPSC_type = GPSC_type_new,
    
    # Reference GPSC 233 = OR 1
    corrected_OR =
      ifelse(
        GPSC == "233",
        1,
        corrected_OR
      ),
    
    corrected_CI_lower =
      ifelse(
        GPSC == "233",
        1,
        corrected_CI_lower
      ),
    
    corrected_CI_upper =
      ifelse(
        GPSC == "233",
        1,
        corrected_CI_upper
      )
  ) %>%
  
  select(
    -GPSC_type_new
  )


# ============================================================
# 28. CHECK FOREST DATA
# ============================================================

cat("\n============================================================\n")
cat("DISEASE FOREST DATA\n")
cat("============================================================\n")

print(
  disease_forest_data %>%
    select(
      GPSC,
      GPSC_position,
      corrected_OR,
      corrected_CI_lower,
      corrected_CI_upper,
      GPSC_type,
      stars
    ),
  n = Inf
)


cat("\n============================================================\n")
cat("CARRIAGE FOREST DATA\n")
cat("============================================================\n")

print(
  carriage_forest_data %>%
    select(
      GPSC,
      GPSC_position,
      corrected_OR,
      corrected_CI_lower,
      corrected_CI_upper,
      GPSC_type,
      stars
    ),
  n = Inf
)


# ============================================================
# 29. CHECK GPSC TYPES
# ============================================================

cat("\n============================================================\n")
cat("DISEASE GPSC TYPES\n")
cat("============================================================\n")

print(
  table(
    disease_forest_data$GPSC_type,
    useNA = "ifany"
  )
)


cat("\n============================================================\n")
cat("CARRIAGE GPSC TYPES\n")
cat("============================================================\n")

print(
  table(
    carriage_forest_data$GPSC_type,
    useNA = "ifany"
  )
)


# ============================================================
# 30. STOP IF GPSC TYPE IS MISSING
# ============================================================

if (
  any(
    is.na(
      disease_forest_data$GPSC_type
    )
  )
) {
  
  cat("\nGPSCs with missing disease type:\n")
  
  print(
    disease_forest_data %>%
      filter(
        is.na(GPSC_type)
      ) %>%
      select(
        GPSC
      )
  )
  
  stop(
    "Disease forest data contains missing GPSC_type."
  )
}


if (
  any(
    is.na(
      carriage_forest_data$GPSC_type
    )
  )
) {
  
  cat("\nGPSCs with missing carriage type:\n")
  
  print(
    carriage_forest_data %>%
      filter(
        is.na(GPSC_type)
      ) %>%
      select(
        GPSC
      )
  )
  
  stop(
    "Carriage forest data contains missing GPSC_type."
  )
}


# ============================================================
# 31. SIGNIFICANCE DATA
# ============================================================

disease_star_data <-
  disease_forest_data %>%
  filter(
    GPSC != "53",
    !is.na(stars),
    stars != ""
  )


carriage_star_data <-
  carriage_forest_data %>%
  filter(
    GPSC != "233",
    !is.na(stars),
    stars != ""
  )


# Fail early if a significant result would be drawn without a star.
# This keeps the figure annotation tied directly to the BH-adjusted
# significance flag used in the exported results tables.
validate_significance_stars <- function(forest_data, panel_name) {
  missing_stars <- forest_data %>%
    filter(
      Significant %in% TRUE,
      is.na(stars) | stars == ""
    )

  unexpected_stars <- forest_data %>%
    filter(
      !(Significant %in% TRUE),
      !is.na(stars),
      stars != ""
    )

  if (nrow(missing_stars) > 0) {
    stop(
      panel_name,
      " significant GPSC(s) missing a star: ",
      paste(missing_stars$GPSC, collapse = ", ")
    )
  }

  if (nrow(unexpected_stars) > 0) {
    stop(
      panel_name,
      " non-significant GPSC(s) unexpectedly marked with a star: ",
      paste(unexpected_stars$GPSC, collapse = ", ")
    )
  }

  message(
    panel_name,
    ": verified stars for all ",
    nrow(missing_stars) + sum(forest_data$Significant %in% TRUE),
    " BH-significant GPSC(s)."
  )
}

validate_significance_stars(disease_forest_data, "Invasive")
validate_significance_stars(carriage_forest_data, "Carriage")


# ============================================================
# 32. COLOURS AND FOREST-PLOT FUNCTION
# ============================================================

# Use the supplied serotype colour document for every serotype in
# either panel. Vaccine types are listed first in the legend;
# non-vaccine types follow and are distinguished by stripes.
observed_serotypes <- sort(
  unique(
    c(
      disease_bar_data$Serotype_merged,
      carriage_bar_data$Serotype_merged
    )
  )
)

observed_serotypes <- observed_serotypes[
  !is.na(observed_serotypes) & observed_serotypes != ""
]

serotype_order <- c(
  VT_serotypes[VT_serotypes %in% observed_serotypes],
  setdiff(observed_serotypes, VT_serotypes)
)

serotype_palette <- utils::read.csv(
  "./serotype_colours.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_palette_columns <- c(
  "In_silico_serotype",
  "In_silico_serotype__colour"
)

if (!all(required_palette_columns %in% names(serotype_palette))) {
  stop(
    "serotype_colours.csv must contain columns: ",
    paste(required_palette_columns, collapse = ", ")
  )
}

serotype_palette <- serotype_palette %>%
  transmute(
    serotype = trimws(as.character(In_silico_serotype)),
    colour = trimws(as.character(In_silico_serotype__colour)),
    match_key = toupper(serotype)
  ) %>%
  filter(serotype != "", colour != "") %>%
  distinct(match_key, .keep_all = TRUE)

palette_match <- match(
  toupper(serotype_order),
  serotype_palette$match_key
)

missing_palette_serotypes <- serotype_order[is.na(palette_match)]

if (length(missing_palette_serotypes) > 0) {
  stop(
    "No colour supplied in serotype_colours.csv for plotted serotype(s): ",
    paste(missing_palette_serotypes, collapse = ", ")
  )
}

serotype_colours <- setNames(
  serotype_palette$colour[palette_match],
  serotype_order
)

# Validate colour strings before ggplot starts drawing the figures.
invalid_palette_colours <- !vapply(
  serotype_colours,
  function(colour) {
    tryCatch(
      {
        grDevices::col2rgb(colour)
        TRUE
      },
      error = function(e) FALSE
    )
  },
  logical(1)
)

if (any(invalid_palette_colours)) {
  stop(
    "Invalid colour value(s) in serotype_colours.csv for: ",
    paste(names(serotype_colours)[invalid_palette_colours], collapse = ", ")
  )
}

gpsc_type_colours <- c(
  "VT-GPSC" = "#8C8C8C",
  "NVT-GPSC" = "#2C7FB8",
  "Mixed VT/NVT-GPSC" = "#18A558"
)


make_forest_plot <- function(
    forest_data,
    star_data,
    gpsc_positions,
    reference_gpsc
) {

  # Match the top-to-bottom ordering used by make_bar_plot().
  y_lookup <- gpsc_positions %>%
    mutate(
      GPSC = as.character(GPSC),
      y_base = rev(seq_along(GPSC))
    ) %>%
    select(GPSC, y_base)

  plot_data <- forest_data %>%
    select(-any_of(c("GPSC_position"))) %>%
    left_join(y_lookup, by = "GPSC")

  # Reference GPSCs are deliberately omitted from the statistical
  # result tables, so add their OR = 1 plotting rows here.
  reference_type <- if (reference_gpsc == "53") {
    disease_status %>%
      filter(GPSC == reference_gpsc) %>%
      pull(GPSC_type)
  } else {
    carriage_status %>%
      filter(GPSC == reference_gpsc) %>%
      pull(GPSC_type)
  }

  reference_row <- y_lookup %>%
    filter(GPSC == reference_gpsc) %>%
    mutate(
      corrected_OR = 1,
      corrected_CI_lower = 1,
      corrected_CI_upper = 1,
      GPSC_type = reference_type[1],
      stars = ""
    )

  plot_data <- plot_data %>%
    filter(GPSC != reference_gpsc) %>%
    bind_rows(reference_row)

  finite_limits <- c(
    plot_data$corrected_CI_lower,
    plot_data$corrected_CI_upper
  )
  finite_limits <- finite_limits[is.finite(finite_limits) & finite_limits > 0]

  ci_limits <- range(c(0.05, 1, finite_limits), na.rm = TRUE)
  lower_limit <- 10^floor(log10(ci_limits[1]))
  upper_ci_power <- 10^ceiling(log10(ci_limits[2]))

  # Put all significance symbols at one shared x coordinate.  The
  # invisible line reserves a consistent column beyond every CI.
  star_column_x <- upper_ci_power * 1.6
  x_limits <- c(lower_limit, star_column_x * 1.45)

  # Use the validated significance table supplied to this function,
  # rather than reconstructing the subset from the plotting rows.
  star_plot_data <- star_data %>%
    mutate(GPSC = as.character(GPSC)) %>%
    filter(Significant %in% TRUE) %>%
    select(GPSC, stars) %>%
    left_join(y_lookup, by = "GPSC") %>%
    mutate(star_x = star_column_x)

  ggplot(
    plot_data,
    aes(
      x = corrected_OR,
      y = y_base,
      colour = GPSC_type
    )
  ) +
    geom_vline(
      xintercept = 1,
      linetype = "dashed",
      colour = "grey55",
      linewidth = 0.45
    ) +
    geom_errorbar(
      aes(
        xmin = corrected_CI_lower,
        xmax = corrected_CI_upper
      ),
      orientation = "y",
      width = 0,
      linewidth = 0.55
    ) +
    geom_point(size = 2.8) +
    geom_vline(
      xintercept = star_column_x,
      colour = "transparent",
      linewidth = 0.1,
      show.legend = FALSE
    ) +
    geom_text(
      data = star_plot_data,
      aes(x = star_x, y = y_base, label = stars),
      inherit.aes = FALSE,
      hjust = 0.5,
      fontface = "bold",
      colour = "black",
      size = 5.2,
      show.legend = FALSE
    ) +
    scale_colour_manual(
      values = gpsc_type_colours,
      breaks = names(gpsc_type_colours),
      name = "GPSC class",
      drop = FALSE
    ) +
    scale_x_log10(
      limits = x_limits,
      breaks = scales::breaks_log(n = 5),
      labels = scales::label_number(accuracy = 0.01),
      expand = expansion(mult = c(0.02, 0.08))
    ) +
    scale_y_continuous(
      breaks = y_lookup$y_base,
      labels = NULL,
      expand = expansion(mult = c(0.015, 0.015))
    ) +
    labs(
      x = "Odds ratio (log scale)",
      y = "GPSC",
      caption = "BH-adjusted significance: * p < 0.05; ** p < 0.01; *** p < 0.001"
    ) +
    coord_cartesian(clip = "off") +
    theme_minimal(base_size = 15) +
    theme(
      axis.title = element_text(face = "bold", size = 16),
      axis.text.x = element_text(face = "bold", size = 12),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3),
      panel.grid.minor.y = element_blank(),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 13),
      legend.text = element_text(size = 11),
      plot.caption = element_text(size = 10, hjust = 0),
      plot.margin = margin(5, 15, 5, 5)
    )
}


# ============================================================
# 30. SIGNIFICANCE DATA
# ============================================================

disease_star_data <-
  disease_forest_data %>%
  filter(
    GPSC != "53",
    !is.na(stars),
    stars != ""
  )


carriage_star_data <-
  carriage_forest_data %>%
  filter(
    GPSC != "233",
    !is.na(stars),
    stars != ""
  )


# ============================================================
# 30. BAR PLOT FUNCTION
# ============================================================

# ============================================================
# CORRECTED GPSC BAR PLOT FUNCTION
# ============================================================

make_bar_plot <- function(
    bar_data,
    gpsc_positions,
    title_text
) {
  
  # ----------------------------------------------------------
  # Create plotting positions
  # ----------------------------------------------------------
  
  # GPSCs are displayed from top to bottom in the order
  # supplied by gpsc_positions.
  
  gpsc_positions_plot <- gpsc_positions %>%
    mutate(
      y_base =
        rev(
          seq_along(GPSC)
        )
    )
  
  
  # ----------------------------------------------------------
  # Add plotting y positions to bar data
  # ----------------------------------------------------------
  
  plot_data <- bar_data %>%
    
    select(
      -any_of(
        c(
          "y",
          "GPSC_position"
        )
      )
    ) %>%
    
    left_join(
      
      gpsc_positions_plot %>%
        select(
          GPSC,
          y_base
        ),
      
      by = "GPSC"
    ) %>%
    
    mutate(
      
      # Two horizontal bars per GPSC
      y =
        case_when(
          
          vaccine_period == "pre" ~
            y_base + 0.18,
          
          vaccine_period == "post" ~
            y_base - 0.18,
          
          TRUE ~
            y_base
        ),
      
      # VT = plain
      # NVT = striped
      Pattern =
        case_when(
          
          vaccine_type == "VT" ~
            "none",
          
          vaccine_type == "NVT" ~
            "stripe",
          
          TRUE ~
            "none"
        )
    )

  panel_serotype_order <- serotype_order[
    serotype_order %in% unique(plot_data$Serotype_merged)
  ]

  max_period_total <- plot_data %>%
    group_by(GPSC, vaccine_period) %>%
    summarise(total = sum(Count), .groups = "drop") %>%
    summarise(maximum = max(total), .groups = "drop") %>%
    pull(maximum)

  period_label_x <- -0.035 * max_period_total
  bar_left_limit <- -0.085 * max_period_total
  
  
  # ----------------------------------------------------------
  # Plot
  # ----------------------------------------------------------
  
  p <- ggplot(
    
    plot_data,
    
    aes(
      x = Count,
      y = y,
      fill = Serotype_merged,
      pattern = Pattern
    )
  ) +
    
    
    # --------------------------------------------------------
  # STACKED HORIZONTAL BARS
  # --------------------------------------------------------
  
  ggpattern::geom_col_pattern(
    
    position = "stack",
    
    orientation = "y",
    
    width = 0.28,
    
    colour = "black",
    
    linewidth = 0.15,
    
    pattern_density = 0.35,
    
    pattern_spacing = 0.025,
    
    pattern_angle = 45,
    
    pattern_colour = "black",
    
    pattern_fill = "black"
  ) +
    
    
    # --------------------------------------------------------
  # PRE / POST LABELS
  # --------------------------------------------------------
  
  geom_text(
    
    data =
      gpsc_positions_plot %>%
      tidyr::expand_grid(
        Period =
          c(
            "pre",
            "post"
          )
      ) %>%
      
      mutate(
        
        y =
          case_when(
            
            Period == "pre" ~
              y_base + 0.18,
            
            Period == "post" ~
              y_base - 0.18
          )
      ),
    
    aes(
      
      x = period_label_x,
      
      y = y,
      
      label = Period
    ),
    
    inherit.aes = FALSE,
    
    hjust = 0.5,
    
    size = 3.5,
    
    fontface = "bold"
  ) +
    
    
    # --------------------------------------------------------
  # SEROTYPE COLOURS
  # --------------------------------------------------------
  
  scale_fill_manual(
    
    values =
      serotype_colours,
    
    breaks =
      panel_serotype_order,
    
    name =
      "Serotype\n(solid = VT; striped = NVT)",
    
    drop =
      TRUE
  ) +
    
    
    # --------------------------------------------------------
  # PATTERN
  # --------------------------------------------------------
  
  scale_pattern_manual(
    
    values =
      c(
        
        "none" =
          "none",
        
        "stripe" =
          "stripe"
      ),
    
    guide =
      "none"
  ) +

    guides(
      fill = guide_legend(
        override.aes = list(
          pattern = ifelse(
            panel_serotype_order %in% VT_serotypes,
            "none",
            "stripe"
          )
        )
      )
    ) +
    
    
    # --------------------------------------------------------
  # Y AXIS
  # --------------------------------------------------------
  
  scale_y_continuous(
    
    breaks =
      gpsc_positions_plot$y_base,
    
    labels =
      gpsc_positions_plot$GPSC,
    
    expand =
      expansion(
        mult =
          c(
            0.015,
            0.015
          )
      )
  ) +
    
    
    # --------------------------------------------------------
  # X AXIS AT TOP
  # --------------------------------------------------------
  
  scale_x_continuous(
    
    position =
      "top",

    limits =
      c(
        bar_left_limit,
        NA
      ),

    breaks = function(limits) {
      breaks <- scales::breaks_pretty(n = 4)(c(0, limits[2]))
      breaks[breaks >= 0]
    },
    
    expand =
      expansion(
        mult =
          c(
            0,
            0.05
          )
      )
  ) +
    
    
    # --------------------------------------------------------
  # LABELS
  # --------------------------------------------------------
  
  labs(
    
    title =
      title_text,
    
    x =
      "Number of samples (n)",
    
    y =
      NULL
  ) +
    
    
    # --------------------------------------------------------
  # THEME
  # --------------------------------------------------------
  
  theme_minimal(
    base_size = 15
  ) +

    coord_cartesian(
      clip = "off"
    ) +
    
    theme(
      
      plot.title =
        element_text(
          size = 22,
          face = "bold",
          hjust = 0
        ),
      
      axis.title.x =
        element_text(
          size = 17,
          face = "bold"
        ),
      
      axis.text.x =
        element_text(
          size = 13,
          face = "bold"
        ),
      
      axis.text.y =
        element_text(
          size = 14,
          face = "bold",
          margin = margin(r = 6)
        ),
      
      axis.ticks =
        element_blank(),
      
      panel.grid.major.y =
        element_line(
          colour = "grey90",
          linewidth = 0.3
        ),
      
      panel.grid.minor =
        element_blank(),
      
      legend.title =
        element_text(
          face = "bold",
          size = 14
        ),
      
      legend.text =
        element_text(
          size = 11
        ),
      
      plot.margin =
        margin(
          5,
          5,
          5,
          25
        )
    )
  
  
  return(p)
}


# ============================================================
# 32. CREATE DISEASE BAR
# ============================================================

disease_bar_plot <- make_bar_plot(

  bar_data =
    disease_bar_data,

  gpsc_positions =
    disease_y,

  title_text =
    "Invasive"
)


# ============================================================
# 33. CREATE DISEASE FOREST
# ============================================================

disease_forest_plot <- make_forest_plot(
  
  forest_data =
    disease_forest_data,
  
  star_data =
    disease_star_data,
  
  gpsc_positions =
    disease_y,
  
  reference_gpsc =
    "53"
)


# ============================================================
# 34. CREATE CARRIAGE BAR
# ============================================================

carriage_bar_plot <- make_bar_plot(
  
  bar_data =
    carriage_bar_data,
  
  gpsc_positions =
    carriage_y,
  
  title_text =
    "Carriage"
)


# ============================================================
# 35. CREATE CARRIAGE FOREST
# ============================================================

carriage_forest_plot <- make_forest_plot(
  
  forest_data =
    carriage_forest_data,
  
  star_data =
    carriage_star_data,
  
  gpsc_positions =
    carriage_y,
  
  reference_gpsc =
    "233"
)


# ============================================================
# 36. COMBINE DISEASE
# ============================================================

disease_final_plot <-
  
  disease_bar_plot +
  
  disease_forest_plot +
  
  plot_layout(
    
    widths =
      c(
        1.15,
        1
      ),
    
    guides =
      "collect"
  ) &
  theme(
    legend.position = "right"
  )


# ============================================================
# 37. COMBINE CARRIAGE
# ============================================================

carriage_final_plot <-
  
  carriage_bar_plot +
  
  carriage_forest_plot +
  
  plot_layout(
    
    widths =
      c(
        1.15,
        1
      ),
    
    guides =
      "collect"
  ) &
  theme(
    legend.position = "right"
  )


# ============================================================
# 38. DISPLAY
# ============================================================

# Display plots when the script is sourced in an interactive R session.
# Avoid creating a slow, unused Rplots.pdf during command-line runs.
if (interactive()) {
  print(
    disease_final_plot
  )

  print(
    carriage_final_plot
  )
}


# ============================================================
# 39. SAVE DISEASE
# ============================================================

ggsave(
  
  filename =
    "Invasive_GPSC_final.png",
  
  plot =
    disease_final_plot,
  
  width =
    15,
  
  height =
    max(
      8,
      length(disease_order) * 0.40
    ),
  
  units =
    "in",
  
  dpi =
    300
)


ggsave(
  
  filename =
    "Invasive_GPSC_final.pdf",
  
  plot =
    disease_final_plot,
  
  width =
    15,
  
  height =
    max(
      8,
      length(disease_order) * 0.40
    ),
  
  units =
    "in"
)


# SVG is the preferred vector format for modern Microsoft Word
# and PowerPoint. It remains sharp when resized on the page.
ggsave(

  filename =
    "Invasive_GPSC_final.svg",

  plot =
    disease_final_plot,

  device =
    svglite::svglite,

  width =
    15,

  height =
    max(
      8,
      length(disease_order) * 0.40
    ),

  units =
    "in"
)


# ============================================================
# 40. SAVE CARRIAGE
# ============================================================

ggsave(
  
  filename =
    "Carriage_GPSC_final.png",
  
  plot =
    carriage_final_plot,
  
  width =
    15,
  
  height =
    max(
      8,
      length(carriage_order) * 0.40
    ),
  
  units =
    "in",
  
  dpi =
    300
)


ggsave(
  
  filename =
    "Carriage_GPSC_final.pdf",
  
  plot =
    carriage_final_plot,
  
  width =
    15,
  
  height =
    max(
      8,
      length(carriage_order) * 0.40
    ),
  
  units =
    "in"
)


ggsave(

  filename =
    "Carriage_GPSC_final.svg",

  plot =
    carriage_final_plot,

  device =
    svglite::svglite,

  width =
    15,

  height =
    max(
      8,
      length(carriage_order) * 0.40
    ),

  units =
    "in"
)


# ============================================================
# 41. SAVE PLOTTING DATA
# ============================================================

write_xlsx(
  
  list(
    
    Results =
      disease_results,
    
    Bar_data =
      disease_bar_data,
    
    Forest_data =
      disease_forest_data,
    
    GPSC_status =
      disease_status
    
  ),
  
  "Invasive_GPSC_plot_data.xlsx"
)


write_xlsx(
  
  list(
    
    Results =
      carriage_results,
    
    Bar_data =
      carriage_bar_data,
    
    Forest_data =
      carriage_forest_data,
    
    GPSC_status =
      carriage_status
    
  ),
  
  "Carriage_GPSC_plot_data.xlsx"
)


# ============================================================
# 42. FINAL CHECKS
# ============================================================

cat("\n")
cat("============================================================\n")
cat("GPSC PLOTTING COMPLETE\n")
cat("============================================================\n")

cat("\nDisease reference: GPSC 53\n")
cat("Carriage reference: GPSC 233\n")

cat("\nDisease GPSCs plotted:\n")
print(disease_order)

cat("\nCarriage GPSCs plotted:\n")
print(carriage_order)

cat("\nForest colours:\n")
cat("  VT-GPSC            = grey\n")
cat("  NVT-GPSC           = blue\n")
cat("  Mixed VT/NVT-GPSC  = green\n")

cat("\nBar patterns:\n")
cat("  VT                 = solid\n")
cat("  NVT                = striped\n")
cat("  Mixed              = solid + striped stacked\n")

cat("\nPlots saved as:\n")
cat("  Invasive_GPSC_final.png\n")
cat("  Invasive_GPSC_final.pdf\n")
cat("  Invasive_GPSC_final.svg\n")
cat("  Carriage_GPSC_final.png\n")
cat("  Carriage_GPSC_final.pdf\n")
cat("  Carriage_GPSC_final.svg\n")
cat("  Pre_PCV10_NVT_reference_frequencies.xlsx\n")
