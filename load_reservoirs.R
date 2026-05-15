library(readr)
library(dplyr)

load_reservoirs <- function(
    path = "coloc-simulation-results.csv.gz",
    trap_quantile = 0.90
) {
  coloc_data <- read_csv(path, show_col_types = FALSE) %>%
    mutate(W = H4 - H3)
  
  # -------------------------
  # 1. Shared signal pool
  # -------------------------
  reservoir_shared <- bind_rows(
    filter(coloc_data, scen == "A_1_B_-1"),
    filter(coloc_data, scen == "A_1_B_1")
  )
  
  shared_h4_q25 <- quantile(reservoir_shared$H4, 0.25, na.rm = TRUE)
  shared_h4_median <- median(reservoir_shared$H4, na.rm = TRUE)
  
  # Weak shared pool:
  # Contains all shared edges below the median.
  # This is retained for reference but may include very weak edges.
  reservoir_shared_weak <- reservoir_shared %>%
    filter(H4 <= shared_h4_median)
  
  # Moderate shared pool:
  # Used for the local bridge trap scenario.
  # This keeps the true clusters relatively weak but avoids near-zero shared edges.
  reservoir_shared_moderate <- reservoir_shared %>%
    filter(H4 >= shared_h4_q25, H4 <= shared_h4_median)
  
  # Strong shared pool:
  # Used for triangle, bowtie, and LD trap scenarios,
  # where the intended true signal structure should be clearly present.
  reservoir_shared_strong <- reservoir_shared %>%
    filter(H4 > shared_h4_median)
  
  # -------------------------
  # 2. Distinct pool from A_0_B_0
  # -------------------------
  a0b0 <- filter(coloc_data, scen == "A_0_B_0")
  
  distinct_h4_median <- median(a0b0$H4, na.rm = TRUE)
  distinct_h4_trap_cut <- quantile(a0b0$H4, trap_quantile, na.rm = TRUE)
  
  # Ordinary background noise
  reservoir_distinct_bg <- a0b0 %>%
    filter(H4 <= distinct_h4_median)
  
  # Dangerous high-H4 distinct tail
  # Used for LD trap / local bridge trap cross-cluster edges
  reservoir_distinct_trap <- a0b0 %>%
    filter(H4 >= distinct_h4_trap_cut)
  
  reservoir_distinct_mid <- a0b0 %>%
    filter(H4 > distinct_h4_median, H4 < distinct_h4_trap_cut)
  
  # -------------------------
  # 3. Ambiguous / mixed pool
  # -------------------------
  reservoir_ambiguous <- filter(coloc_data, scen == "A_1_B_0")
  
  list(
    coloc_data = coloc_data,
    
    reservoir_shared = reservoir_shared,
    reservoir_shared_weak = reservoir_shared_weak,
    reservoir_shared_moderate = reservoir_shared_moderate,
    reservoir_shared_strong = reservoir_shared_strong,
    
    a0b0 = a0b0,
    reservoir_distinct_bg = reservoir_distinct_bg,
    reservoir_distinct_mid = reservoir_distinct_mid,
    reservoir_distinct_trap = reservoir_distinct_trap,
    
    reservoir_ambiguous = reservoir_ambiguous,
    
    shared_h4_q25 = shared_h4_q25,
    shared_h4_median = shared_h4_median,
    distinct_h4_median = distinct_h4_median,
    distinct_h4_trap_cut = distinct_h4_trap_cut
  )
}
