
# benchmark_cpm_fixed.R
# Timing comparison for custom k=3 CPM workflow.
# Put this file in the same folder as:
# load_reservoirs.R, network_generators.R, method.R
# and coloc-simulation-results.csv.gz

library(dplyr)
library(tibble)
library(ggplot2)
library(microbenchmark)
library(igraph)

source("load_reservoirs.R")
source("network_generators.R")
source("method.R")

set.seed(2769)

# -----------------------------
# Load reservoirs
# -----------------------------
reservoirs <- load_reservoirs(
  path = "coloc-simulation-results.csv.gz",
  trap_quantile = 0.90
)

# Your actual reservoir object names are:
# reservoirs$reservoir_shared
# reservoirs$reservoir_shared_weak
# reservoirs$reservoir_shared_moderate
# reservoirs$reservoir_shared_strong
# reservoirs$reservoir_distinct_bg
# reservoirs$reservoir_distinct_mid
# reservoirs$reservoir_distinct_trap
# reservoirs$reservoir_ambiguous

# -----------------------------
# Settings
# -----------------------------
n_reps <- 50
times_per_network <- 10L
tau_grid <- seq(0.1, 0.7, by = 0.05)

scenarios <- c(
  "triangle",
  "bowtie",
  "ld_trap",
  "local_bridge_trap"
)

# -----------------------------
# Generate one network by scenario
# -----------------------------
generate_network_by_scenario <- function(scenario) {
  switch(
    scenario,

    "triangle" = generate_triangle_network(
      reservoir_shared = reservoirs$reservoir_shared_strong,
      reservoir_distinct_bg = reservoirs$reservoir_distinct_bg
    ),

    "bowtie" = generate_bowtie_network(
      reservoir_shared = reservoirs$reservoir_shared_strong,
      reservoir_distinct_bg = reservoirs$reservoir_distinct_bg
    ),

    "ld_trap" = generate_ld_trap_network(
      reservoir_shared = reservoirs$reservoir_shared_strong,
      reservoir_distinct_bg = reservoirs$reservoir_distinct_bg,
      reservoir_distinct_trap = reservoirs$reservoir_distinct_trap
    ),

    "local_bridge_trap" = generate_local_bridge_trap_network(
      reservoir_shared_weak = reservoirs$reservoir_shared_moderate,
      reservoir_distinct_bg = reservoirs$reservoir_distinct_bg,
      reservoir_distinct_trap = reservoirs$reservoir_distinct_trap
    ),

    stop(sprintf("Unknown scenario: %s", scenario))
  )
}

# -----------------------------
# Custom CPM: your actual function
# -----------------------------
run_custom_cpm <- function(H4_mat) {
  run_cpm_dynamic(
    H4_mat = H4_mat,
    tau_grid = tau_grid,
    clique_size = 3,
    overlap_required = 2
  )
}

# -----------------------------
# General-purpose CPM-like implementation
# This deliberately keeps the workflow more generic:
# 1. threshold H4
# 2. enumerate k-cliques using igraph::cliques()
# 3. build clique-overlap graph
# 4. score every tau by modularity
#
# It should produce the same type of CPM logic but is less project-specific.
# -----------------------------
run_general_cpm_like <- function(H4_mat,
                                 tau_grid = seq(0.1, 0.7, by = 0.05),
                                 k = 3,
                                 overlap_required = k - 1) {
  best <- NULL
  best_score <- -Inf

  g_weighted <- graph_from_adjacency_matrix(
    H4_mat,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )

  for (tau in tau_grid) {
    A_bin <- (H4_mat >= tau) * 1
    diag(A_bin) <- 0

    g_bin <- graph_from_adjacency_matrix(
      A_bin,
      mode = "undirected",
      diag = FALSE
    )

    clique_list <- cliques(g_bin, min = k, max = k)
    clique_list <- lapply(clique_list, as.integer)

    communities <- build_cpm_communities(
      clique_list = clique_list,
      overlap_required = overlap_required
    )

    temp_mem <- temporary_membership_from_cpm(
      communities = communities,
      n_nodes = nrow(H4_mat)
    )

    score <- modularity(
      g_weighted,
      membership = temp_mem,
      weights = E(g_weighted)$weight
    )

    if (is.finite(score) && score > best_score) {
      best_score <- score
      best <- list(
        tau = tau,
        communities = communities,
        temp_membership = temp_mem,
        objective = score,
        n_cliques = length(clique_list),
        n_communities = length(communities)
      )
    }
  }

  if (is.null(best)) {
    best <- list(
      tau = NA_real_,
      communities = list(),
      temp_membership = seq_len(nrow(H4_mat)),
      objective = NA_real_,
      n_cliques = 0L,
      n_communities = 0L
    )
  }

  best
}

# -----------------------------
# Benchmark loop
# -----------------------------
benchmark_results <- list()
counter <- 1

for (sc in scenarios) {
  message("Benchmarking scenario: ", sc)

  for (rep in seq_len(n_reps)) {
    net <- generate_network_by_scenario(sc)
    H4_mat <- net$H4

    bm <- microbenchmark(
      custom = run_custom_cpm(H4_mat),
      general_like = run_general_cpm_like(
        H4_mat = H4_mat,
        tau_grid = tau_grid,
        k = 3,
        overlap_required = 2
      ),
      times = times_per_network
    )

    benchmark_results[[counter]] <- tibble(
      scenario = sc,
      replicate = rep,
      method = as.character(bm$expr),
      time_ms = as.numeric(bm$time) / 1e6
    )

    counter <- counter + 1
  }
}

benchmark_df <- bind_rows(benchmark_results)

write.csv(
  benchmark_df,
  "cpm_timing_benchmark_raw.csv",
  row.names = FALSE
)

summary_df <- benchmark_df %>%
  group_by(scenario, method) %>%
  summarise(
    median_ms = median(time_ms),
    mean_ms = mean(time_ms),
    q25_ms = quantile(time_ms, 0.25),
    q75_ms = quantile(time_ms, 0.75),
    .groups = "drop"
  ) %>%
  group_by(scenario) %>%
  mutate(
    baseline_median_ms = median_ms[method == "general_like"][1],
    relative_speed = baseline_median_ms / median_ms
  ) %>%
  ungroup()

write.csv(
  summary_df,
  "cpm_timing_benchmark_summary.csv",
  row.names = FALSE
)

print(summary_df)

# -----------------------------
# Plot for PPT
# -----------------------------
p <- ggplot(summary_df, aes(x = scenario, y = median_ms, fill = method)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  labs(
    x = "Simulation scenario",
    y = "Median runtime across tau grid (ms)",
    fill = "Implementation",
    title = "CPM timing comparison"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  "cpm_timing_comparison.png",
  p,
  width = 7,
  height = 4,
  dpi = 300
)

message("Saved:")
message("- cpm_timing_benchmark_raw.csv")
message("- cpm_timing_benchmark_summary.csv")
message("- cpm_timing_comparison.png")
