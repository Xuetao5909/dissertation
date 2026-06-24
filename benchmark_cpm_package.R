

library(dplyr)
library(ggplot2)
library(microbenchmark)
library(igraph)

source("load_reservoirs.R")
source("network_generators.R")
source("method.R")

# -----------------------------
# Package check
# -----------------------------
if (!requireNamespace("CliquePercolation", quietly = TRUE)) {
  stop(
    "Package 'CliquePercolation' is not installed. Run:\n",
    "install.packages('CliquePercolation')\n",
    "Then rerun source('benchmark_cpm_package.R')."
  )
}

set.seed(2769)

# -----------------------------
# Settings
# -----------------------------
n_reps <- 50
times_per_network <- 10L
tau_grid <- seq(0.1, 0.7, by = 0.05)
scenarios <- c("triangle", "bowtie", "ld_trap", "local_bridge_trap")

# -----------------------------
# Load reservoirs
# -----------------------------
reservoirs <- load_reservoirs()

# -----------------------------
# Generate one network by scenario
# -----------------------------
generate_network_by_scenario <- function(scenario) {
  if (scenario == "triangle") {
    return(generate_triangle_network(
      reservoir_shared = reservoirs$reservoir_shared_strong,
      reservoir_distinct_bg = reservoirs$reservoir_distinct_bg
    ))
  }
  
  if (scenario == "bowtie") {
    return(generate_bowtie_network(
      reservoir_shared = reservoirs$reservoir_shared_strong,
      reservoir_distinct_bg = reservoirs$reservoir_distinct_bg
    ))
  }
  
  if (scenario == "ld_trap") {
    return(generate_ld_trap_network(
      reservoir_shared = reservoirs$reservoir_shared_strong,
      reservoir_distinct_bg = reservoirs$reservoir_distinct_bg,
      reservoir_distinct_trap = reservoirs$reservoir_distinct_trap
    ))
  }
  
  if (scenario == "local_bridge_trap") {
    return(generate_local_bridge_trap_network(
      reservoir_shared_weak = reservoirs$reservoir_shared_moderate,
      reservoir_distinct_bg = reservoirs$reservoir_distinct_bg,
      reservoir_distinct_trap = reservoirs$reservoir_distinct_trap
    ))
  }
  
  stop("Unknown scenario: ", scenario)
}

# -----------------------------
# Your custom CPM: scans tau_grid internally and chooses best tau by modularity.
# -----------------------------
run_custom_cpm_once <- function(H4_mat) {
  run_cpm_dynamic(
    H4_mat = H4_mat,
    tau_grid = tau_grid,
    clique_size = 3,
    overlap_required = 2
  )
}

# -----------------------------
# Package CPM, FAIR VERSION:
# To mimic your workflow, threshold H4 at each tau first, then run package CPM as unweighted.
# This keeps the rule closer to your method:
#   edge survives iff H4_ij >= tau.
# Note: package output is not used here for scoring; this is runtime comparison only.
# -----------------------------
run_package_cpm_unweighted_grid <- function(H4_mat) {
  out <- vector("list", length(tau_grid))
  
  for (i in seq_along(tau_grid)) {
    tau <- tau_grid[i]
    A_bin <- (H4_mat >= tau) * 1
    diag(A_bin) <- 0
    A_bin <- as.matrix(A_bin)
    storage.mode(A_bin) <- "numeric"
    
    out[[i]] <- CliquePercolation::cpAlgorithm(
      W = A_bin,
      k = 3,
      method = "unweighted"
    )
  }
  
  out
}

# -----------------------------
# Package CPM, PACKAGE-NATIVE WEIGHTED VERSION:
# This uses cpAlgorithm(method = "weighted", I = tau).
# Warning: this is not exactly identical to your edge-threshold workflow,
# because weighted CPM thresholds clique intensity, not each individual edge.
# Keep this optional unless you explicitly want the package-native comparison.
# -----------------------------
run_package_cpm_weighted_grid <- function(H4_mat) {
  out <- vector("list", length(tau_grid))
  
  for (i in seq_along(tau_grid)) {
    tau <- tau_grid[i]
    W <- as.matrix(H4_mat)
    diag(W) <- 0
    storage.mode(W) <- "numeric"
    
    out[[i]] <- CliquePercolation::cpAlgorithm(
      W = W,
      k = 3,
      method = "weighted",
      I = tau
    )
  }
  
  out
}

# Choose what to compare in the main plot.
# Recommended for your PPT: "package_unweighted_grid" because it best matches your thresholded workflow.
include_package_weighted <- TRUE

benchmark_results <- list()
counter <- 1

for (scenario in scenarios) {
  message("Benchmarking scenario: ", scenario)
  
  for (rep in seq_len(n_reps)) {
    net <- generate_network_by_scenario(scenario)
    H4_mat <- net$H4
    
    if (include_package_weighted) {
      bm <- microbenchmark(
        custom = run_custom_cpm_once(H4_mat),
        package_unweighted_grid = run_package_cpm_unweighted_grid(H4_mat),
        package_weighted_grid = run_package_cpm_weighted_grid(H4_mat),
        times = times_per_network
      )
    } else {
      bm <- microbenchmark(
        custom = run_custom_cpm_once(H4_mat),
        package_unweighted_grid = run_package_cpm_unweighted_grid(H4_mat),
        times = times_per_network
      )
    }
    
    benchmark_results[[counter]] <- data.frame(
      scenario = scenario,
      rep = rep,
      implementation = as.character(bm$expr),
      time_ms = bm$time / 1e6,
      stringsAsFactors = FALSE
    )
    
    counter <- counter + 1
  }
}

benchmark_df <- bind_rows(benchmark_results)
write.csv(benchmark_df, "cpm_package_timing_benchmark_raw.csv", row.names = FALSE)

summary_df <- benchmark_df %>%
  group_by(scenario, implementation) %>%
  summarise(
    median_ms = median(time_ms),
    mean_ms = mean(time_ms),
    q25_ms = quantile(time_ms, 0.25),
    q75_ms = quantile(time_ms, 0.75),
    .groups = "drop"
  ) %>%
  group_by(scenario) %>%
  mutate(
    custom_median_ms = median_ms[implementation == "custom"][1],
    relative_to_custom = median_ms / custom_median_ms,
    custom_speedup_vs_this = median_ms / custom_median_ms
  ) %>%
  ungroup()

write.csv(summary_df, "cpm_package_timing_benchmark_summary.csv", row.names = FALSE)
print(summary_df)

# Main PPT plot: custom vs package version that mimics your threshold workflow.
plot_df <- summary_df %>%
  filter(implementation %in% c("custom", "package_unweighted_grid")) %>%
  mutate(
    implementation = recode(
      implementation,
      custom = "Custom CPM",
      package_unweighted_grid = "CliquePercolation package"
    )
  )

p <- ggplot(plot_df, aes(x = scenario, y = median_ms, fill = implementation)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  labs(
    x = "Simulation scenario",
    y = "Median runtime across tau grid (ms)",
    fill = "Implementation",
    title = "CPM timing comparison"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave("cpm_package_timing_comparison.png", p, width = 7, height = 4, dpi = 300)

# Optional plot including package-native weighted version.
p_all <- ggplot(summary_df, aes(x = scenario, y = median_ms, fill = implementation)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  labs(
    x = "Simulation scenario",
    y = "Median runtime across tau grid (ms)",
    fill = "Implementation",
    title = "CPM timing comparison including package-native weighted CPM"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave("cpm_package_timing_comparison_all.png", p_all, width = 8, height = 4.5, dpi = 300)

message("Done. Files written:")
message("  cpm_package_timing_benchmark_raw.csv")
message("  cpm_package_timing_benchmark_summary.csv")
message("  cpm_package_timing_comparison.png")
message("  cpm_package_timing_comparison_all.png")
