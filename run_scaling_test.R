library(tibble)
library(dplyr)

# ----------------------------------------
# Expanded LD-trap generator
# This expands the true signal structure itself:
# - group1 and group2 become larger true clusters
# - cross-group trap edges are sampled from the high-H4 A_0_B_0 tail
# - remaining edges are background noise
# ----------------------------------------
generate_scaled_ldtrap_network <- function(
    n = 50,
    group_size = 10,
    reservoir_shared_strong,
    reservoir_distinct_bg,
    reservoir_distinct_trap,
    n_trap_edges = NULL
) {
  if (2 * group_size > n) {
    stop("2 * group_size must be <= n.")
  }
  
  group1 <- seq_len(group_size)
  group2 <- seq(group_size + 1, 2 * group_size)
  
  if (is.null(n_trap_edges)) {
    n_trap_edges <- 2 * group_size
  }
  
  net <- make_empty_network(n)
  
  # 1. Fill all edges with background noise first
  net <- fill_all_remaining_edges(net, reservoir_distinct_bg)
  
  # 2. Add strong within-cluster true shared edges
  net <- fill_edges_from_pairs(net, combn(group1, 2), reservoir_shared_strong)
  net <- fill_edges_from_pairs(net, combn(group2, 2), reservoir_shared_strong)
  
  # 3. Add only a limited number of misleading cross-cluster trap edges
  all_cross_pairs <- expand.grid(group1, group2)
  n_trap_edges <- min(n_trap_edges, nrow(all_cross_pairs))
  selected_idx <- sample(seq_len(nrow(all_cross_pairs)), size = n_trap_edges, replace = FALSE)
  selected_pairs <- t(as.matrix(all_cross_pairs[selected_idx, ]))
  
  net <- fill_edges_from_pairs(net, selected_pairs, reservoir_distinct_trap)
  
  finalise_network(net)
}

# ----------------------------------------
# Run three methods and record runtime only
# ----------------------------------------
run_methods_with_runtime <- function(net, cpm_tau_grid = seq(0.1, 0.7, by = 0.05)) {
  lou_time <- system.time({
    lou <- run_louvain(net$H4)
  })
  
  cpm_time <- system.time({
    cpm <- run_cpm_dynamic(
      H4_mat = net$H4,
      tau_grid = cpm_tau_grid,
      overlap_required = 2
    )
  })
  
  spg_time <- system.time({
    spg <- run_signed_spinglass(net$W)
  })
  
  tibble(
    louvain_time_sec = as.numeric(lou_time["elapsed"]),
    cpm_time_sec = as.numeric(cpm_time["elapsed"]),
    spinglass_time_sec = as.numeric(spg_time["elapsed"]),
    cpm_tau = cpm$tau,
    cpm_n_cliques = cpm$n_cliques,
    cpm_n_communities = cpm$n_communities
  )
}

# ----------------------------------------
# Runtime scaling test
#
# Two scaling modes:
# 1. background_expansion:
#    Keep the original biological signal structure fixed,
#    but increase the total number of background/noise nodes.
#
# 2. expanded_ldtrap:
#    Expand the LD-trap structure itself by increasing the true cluster sizes.
# ----------------------------------------
run_scaling_test <- function(
    n_reps = 20,
    node_sizes = c(10, 20, 30, 50, 100),
    reservoir_shared_strong,
    reservoir_shared_moderate,
    reservoir_distinct_bg,
    reservoir_distinct_trap,
    cpm_tau_grid = seq(0.1, 0.7, by = 0.05)
) {
  results <- list()
  counter <- 1
  
  # -------------------------
  # Version 1:
  # Background-node expansion
  # -------------------------
  for (n in node_sizes) {
    for (rep in seq_len(n_reps)) {
      
      # Triangle with fixed 3-node core
      net <- generate_triangle_network(
        n = n,
        reservoir_shared = reservoir_shared_strong,
        reservoir_distinct_bg = reservoir_distinct_bg,
        core = c(1, 2, 3)
      )
      
      timing <- run_methods_with_runtime(net, cpm_tau_grid)
      
      results[[counter]] <- bind_cols(
        tibble(
          scaling_type = "background_expansion",
          scenario = "triangle",
          n_nodes = n,
          true_cluster_size = 3,
          replicate = rep
        ),
        timing
      )
      counter <- counter + 1
      
      # LD trap with fixed two 3-node clusters
      net <- generate_ld_trap_network(
        n = n,
        reservoir_shared = reservoir_shared_strong,
        reservoir_distinct_bg = reservoir_distinct_bg,
        reservoir_distinct_trap = reservoir_distinct_trap,
        group1 = c(1, 2, 3),
        group2 = c(4, 5, 6)
      )
      
      timing <- run_methods_with_runtime(net, cpm_tau_grid)
      
      results[[counter]] <- bind_cols(
        tibble(
          scaling_type = "background_expansion",
          scenario = "ld_trap",
          n_nodes = n,
          true_cluster_size = 3,
          replicate = rep
        ),
        timing
      )
      counter <- counter + 1
      
      # Local bridge trap with fixed weak 4-node clusters
      net <- generate_local_bridge_trap_network(
        n = n,
        reservoir_shared_weak = reservoir_shared_moderate,
        reservoir_distinct_bg = reservoir_distinct_bg,
        reservoir_distinct_trap = reservoir_distinct_trap,
        group1 = c(1, 2, 3, 4),
        group2 = c(5, 6, 7, 8)
      )
      
      timing <- run_methods_with_runtime(net, cpm_tau_grid)
      
      results[[counter]] <- bind_cols(
        tibble(
          scaling_type = "background_expansion",
          scenario = "local_bridge_trap",
          n_nodes = n,
          true_cluster_size = 4,
          replicate = rep
        ),
        timing
      )
      counter <- counter + 1
    }
  }
  
  # -------------------------
  # Version 2:
  # Expanded LD-trap structure
  # -------------------------
  ldtrap_group_sizes <- tibble(
    n_nodes = c(10, 20, 30, 50, 100),
    true_cluster_size = c(3, 5, 7, 10, 20)
  ) %>%
    filter(n_nodes %in% node_sizes)
  
  for (i in seq_len(nrow(ldtrap_group_sizes))) {
    n <- ldtrap_group_sizes$n_nodes[i]
    group_size <- ldtrap_group_sizes$true_cluster_size[i]
    
    for (rep in seq_len(n_reps)) {
      net <- generate_scaled_ldtrap_network(
        n = n,
        group_size = group_size,
        reservoir_shared_strong = reservoir_shared_strong,
        reservoir_distinct_bg = reservoir_distinct_bg,
        reservoir_distinct_trap = reservoir_distinct_trap
      )
      
      timing <- run_methods_with_runtime(net, cpm_tau_grid)
      
      results[[counter]] <- bind_cols(
        tibble(
          scaling_type = "expanded_ldtrap",
          scenario = "ld_trap",
          n_nodes = n,
          true_cluster_size = group_size,
          replicate = rep
        ),
        timing
      )
      counter <- counter + 1
    }
  }
  
  bind_rows(results)
}

# ----------------------------------------
# Summarise scaling runtime
# ----------------------------------------
summarise_scaling_runtime <- function(scaling_res) {
  scaling_res %>%
    group_by(scaling_type, scenario, n_nodes, true_cluster_size) %>%
    summarise(
      louvain_mean_time_sec = mean(louvain_time_sec, na.rm = TRUE),
      louvain_median_time_sec = median(louvain_time_sec, na.rm = TRUE),
      louvain_total_time_sec = sum(louvain_time_sec, na.rm = TRUE),
      
      cpm_mean_time_sec = mean(cpm_time_sec, na.rm = TRUE),
      cpm_median_time_sec = median(cpm_time_sec, na.rm = TRUE),
      cpm_total_time_sec = sum(cpm_time_sec, na.rm = TRUE),
      
      spinglass_mean_time_sec = mean(spinglass_time_sec, na.rm = TRUE),
      spinglass_median_time_sec = median(spinglass_time_sec, na.rm = TRUE),
      spinglass_total_time_sec = sum(spinglass_time_sec, na.rm = TRUE),
      
      cpm_mean_n_cliques = mean(cpm_n_cliques, na.rm = TRUE),
      cpm_mean_tau = mean(cpm_tau, na.rm = TRUE),
      .groups = "drop"
    )
}

