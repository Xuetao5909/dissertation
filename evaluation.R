library(dplyr)

# ----------------------------------------
# Adjusted Rand Index
# Implemented directly to avoid an additional package dependency
# ----------------------------------------
compute_ari <- function(pred, truth) {
  if (length(pred) != length(truth)) {
    stop("pred and truth must have the same length.")
  }
  
  tab <- table(pred, truth)
  n <- sum(tab)
  
  if (n <= 1) return(NA_real_)
  
  choose2 <- function(x) x * (x - 1) / 2
  
  sum_comb_c <- sum(choose2(tab))
  sum_comb_rows <- sum(choose2(rowSums(tab)))
  sum_comb_cols <- sum(choose2(colSums(tab)))
  total_comb <- choose2(n)
  
  expected_index <- sum_comb_rows * sum_comb_cols / total_comb
  max_index <- 0.5 * (sum_comb_rows + sum_comb_cols)
  
  denom <- max_index - expected_index
  
  if (denom == 0) {
    return(NA_real_)
  }
  
  (sum_comb_c - expected_index) / denom
}

# ----------------------------------------
# Ground-truth labels for non-overlapping scenarios
# Bowtie is excluded from ARI because the true structure is overlapping
# ----------------------------------------
get_truth_labels <- function(scenario) {
  switch(
    scenario,
    "triangle" = c(1, 1, 1, 2, 2, 2, 2, 2, 2, 2),
    "ld_trap" = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 3),
    "local_bridge_trap" = c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3),
    "bowtie" = NULL,
    stop(sprintf("Unknown scenario for truth labels: %s", scenario))
  )
}

# ----------------------------------------
# Safe extraction helper for evaluation outputs
# ----------------------------------------
get_metric <- function(x, name, default = NA) {
  if (!is.null(x[[name]])) x[[name]] else default
}

# ----------------------------------------
# Generic core-cluster success check
# Used for triangle-like single-core scenarios
# ----------------------------------------
check_core_cluster_success <- function(mem, core = c(1, 2, 3)) {
  core_same <- length(unique(mem[core])) == 1
  
  list(
    success = core_same,
    false_merge = FALSE,
    over_split = !core_same,
    noise_merge = NA,
    forced_assignment = FALSE,
    left_recovered = NA,
    right_recovered = NA,
    core_separated = NA,
    core_recovered = core_same
  )
}

# ----------------------------------------
# Generic success criterion for hard-clustering methods
# Applicable to Louvain and Spinglass under two-cluster scenarios
# ----------------------------------------
check_two_cluster_success <- function(mem, group1 = c(1, 2, 3), group2 = c(4, 5, 6)) {
  g1_same <- length(unique(mem[group1])) == 1
  g2_same <- length(unique(mem[group2])) == 1
  separate <- unique(mem[group1])[1] != unique(mem[group2])[1]
  
  g1_cluster_id <- unique(mem[group1])[1]
  g2_cluster_id <- unique(mem[group2])[1]
  
  g1_has_noise <- sum(mem == g1_cluster_id) > length(group1)
  g2_has_noise <- sum(mem == g2_cluster_id) > length(group2)
  noise_merge <- g1_has_noise || g2_has_noise
  
  core_recovered <- g1_same && g2_same && separate
  strict_success <- core_recovered && !noise_merge
  
  list(
    success = strict_success,
    false_merge = length(unique(mem[c(group1, group2)])) == 1,
    over_split = !(g1_same && g2_same),
    noise_merge = noise_merge,
    forced_assignment = FALSE,
    left_recovered = g1_same,
    right_recovered = g2_same,
    core_separated = separate,
    core_recovered = core_recovered
  )
}

# ----------------------------------------
# Hard-clustering reference criterion for the bowtie scenario
# Hard-clustering methods cannot explicitly represent overlapping membership
# ----------------------------------------
check_hard_bowtie <- function(mem, left = c(1, 2, 3), right = c(3, 4, 5), hub = 3) {
  left_same <- length(unique(mem[left])) == 1
  right_same <- length(unique(mem[right])) == 1
  
  false_merge <- length(unique(mem[c(left, right)])) == 1
  over_split <- !(left_same || right_same) && !false_merge
  forced_assignment <- (left_same || right_same) && !false_merge
  
  list(
    success = FALSE,
    false_merge = false_merge,
    over_split = over_split,
    noise_merge = NA,
    forced_assignment = forced_assignment,
    left_recovered = left_same,
    right_recovered = right_same,
    core_separated = !false_merge,
    core_recovered = forced_assignment
  )
}

# ----------------------------------------
# CPM: triangle scenario
# ----------------------------------------
check_cpm_triangle <- function(comms, core = c(1, 2, 3)) {
  # For a 3-node core, recovery of at least 2 nodes within a CPM community
  # is treated as successful localisation of the core signal
  has_core <- any(sapply(comms, function(x) sum(x %in% core) >= 2))
  
  list(
    success = has_core,
    false_merge = FALSE,
    left_recovered = NA,
    right_recovered = NA,
    core_separated = NA,
    core_recovered = has_core
  )
}

# ----------------------------------------
# CPM: bowtie scenario
# ----------------------------------------
check_cpm_bowtie <- function(comms, left = c(1, 2, 3), right = c(3, 4, 5), hub = 3) {
  if (length(comms) == 0) {
    return(list(
      success = FALSE,
      false_merge = FALSE,
      hub_overlap = FALSE,
      recovered_left = 0,
      recovered_right = 0,
      left_recovered = FALSE,
      right_recovered = FALSE,
      core_separated = FALSE,
      core_recovered = FALSE
    ))
  }
  
  max_left <- max(sapply(comms, function(x) sum(left %in% x)))
  max_right <- max(sapply(comms, function(x) sum(right %in% x)))
  hub_overlap <- sum(sapply(comms, function(x) hub %in% x)) >= 2
  
  non_hub_left <- setdiff(left, hub)
  non_hub_right <- setdiff(right, hub)
  
  false_merge <- any(sapply(comms, function(x) {
    any(non_hub_left %in% x) && any(non_hub_right %in% x)
  }))
  
  left_recovered <- max_left >= 2
  right_recovered <- max_right >= 2
  core_separated <- !false_merge
  core_recovered <- left_recovered && right_recovered && core_separated
  
  list(
    success = core_recovered && hub_overlap,
    false_merge = false_merge,
    hub_overlap = hub_overlap,
    recovered_left = max_left,
    recovered_right = max_right,
    left_recovered = left_recovered,
    right_recovered = right_recovered,
    core_separated = core_separated,
    core_recovered = core_recovered
  )
}

# ----------------------------------------
# CPM: LD trap scenario
# ----------------------------------------
check_cpm_ldtrap <- function(comms, group1 = c(1, 2, 3), group2 = c(4, 5, 6)) {
  # For each 3-node group, recovery of at least 2 nodes is treated
  # as successful capture of the underlying signal
  has_g1 <- any(sapply(comms, function(x) sum(x %in% group1) >= 2))
  has_g2 <- any(sapply(comms, function(x) sum(x %in% group2) >= 2))
  
  false_merge <- any(sapply(comms, function(x) {
    any(group1 %in% x) && any(group2 %in% x)
  }))
  
  core_separated <- !false_merge
  core_recovered <- has_g1 && has_g2 && core_separated
  
  list(
    success = core_recovered,
    false_merge = false_merge,
    left_recovered = has_g1,
    right_recovered = has_g2,
    core_separated = core_separated,
    core_recovered = core_recovered
  )
}

# ----------------------------------------
# CPM: local bridge trap / scenario 4
# ----------------------------------------
check_cpm_local_bridge_trap <- function(comms, group1 = c(1, 2, 3, 4), group2 = c(5, 6, 7, 8)) {
  # For each 4-node weak-signal group, recovery of at least 3 nodes
  # is treated as successful recovery of the local structure
  has_g1 <- any(sapply(comms, function(x) sum(x %in% group1) >= 3))
  has_g2 <- any(sapply(comms, function(x) sum(x %in% group2) >= 3))
  
  false_merge <- any(sapply(comms, function(x) {
    any(group1 %in% x) && any(group2 %in% x)
  }))
  
  core_separated <- !false_merge
  core_recovered <- has_g1 && has_g2 && core_separated
  
  list(
    success = core_recovered,
    false_merge = false_merge,
    left_recovered = has_g1,
    right_recovered = has_g2,
    core_separated = core_separated,
    core_recovered = core_recovered
  )
}

# ----------------------------------------
# Summary table
# Four methods: Louvain, Infomap, CPM, and Signed Spinglass
# ----------------------------------------
summarise_results <- function(res) {
  res %>%
    group_by(scenario) %>%
    summarise(
      # -------------------------
      # Louvain
      # -------------------------
      louvain_success_rate = mean(louvain_success, na.rm = TRUE),
      louvain_core_recovered_rate = mean(louvain_core_recovered, na.rm = TRUE),
      louvain_left_recovered_rate = mean(louvain_left_recovered, na.rm = TRUE),
      louvain_right_recovered_rate = mean(louvain_right_recovered, na.rm = TRUE),
      louvain_core_separated_rate = mean(louvain_core_separated, na.rm = TRUE),
      louvain_false_merge_rate = mean(louvain_false_merge, na.rm = TRUE),
      louvain_over_split_rate = mean(louvain_over_split, na.rm = TRUE),
      louvain_noise_merge_rate = mean(louvain_noise_merge, na.rm = TRUE),
      louvain_forced_assignment_rate = mean(louvain_forced_assignment, na.rm = TRUE),
      louvain_mean_ari = mean(louvain_ari, na.rm = TRUE),
      louvain_mean_time_sec = mean(louvain_time_sec, na.rm = TRUE),
      
      # -------------------------
      # Infomap
      # -------------------------
      infomap_success_rate = mean(infomap_success, na.rm = TRUE),
      infomap_core_recovered_rate = mean(infomap_core_recovered, na.rm = TRUE),
      infomap_left_recovered_rate = mean(infomap_left_recovered, na.rm = TRUE),
      infomap_right_recovered_rate = mean(infomap_right_recovered, na.rm = TRUE),
      infomap_core_separated_rate = mean(infomap_core_separated, na.rm = TRUE),
      infomap_false_merge_rate = mean(infomap_false_merge, na.rm = TRUE),
      infomap_over_split_rate = mean(infomap_over_split, na.rm = TRUE),
      infomap_noise_merge_rate = mean(infomap_noise_merge, na.rm = TRUE),
      infomap_forced_assignment_rate = mean(infomap_forced_assignment, na.rm = TRUE),
      infomap_mean_ari = mean(infomap_ari, na.rm = TRUE),
      infomap_mean_time_sec = mean(infomap_time_sec, na.rm = TRUE),
      
      # -------------------------
      # CPM
      # -------------------------
      cpm_success_rate = mean(cpm_success, na.rm = TRUE),
      cpm_core_recovered_rate = mean(cpm_core_recovered, na.rm = TRUE),
      cpm_left_recovered_rate = mean(cpm_left_recovered, na.rm = TRUE),
      cpm_right_recovered_rate = mean(cpm_right_recovered, na.rm = TRUE),
      cpm_core_separated_rate = mean(cpm_core_separated, na.rm = TRUE),
      cpm_false_merge_rate = mean(cpm_false_merge, na.rm = TRUE),
      cpm_mean_tau = mean(cpm_tau, na.rm = TRUE),
      cpm_mean_ari = mean(cpm_ari, na.rm = TRUE),
      cpm_mean_time_sec = mean(cpm_time_sec, na.rm = TRUE),
      
      # -------------------------
      # Signed Spinglass
      # -------------------------
      spinglass_success_rate = mean(spinglass_success, na.rm = TRUE),
      spinglass_core_recovered_rate = mean(spinglass_core_recovered, na.rm = TRUE),
      spinglass_left_recovered_rate = mean(spinglass_left_recovered, na.rm = TRUE),
      spinglass_right_recovered_rate = mean(spinglass_right_recovered, na.rm = TRUE),
      spinglass_core_separated_rate = mean(spinglass_core_separated, na.rm = TRUE),
      spinglass_false_merge_rate = mean(spinglass_false_merge, na.rm = TRUE),
      spinglass_over_split_rate = mean(spinglass_over_split, na.rm = TRUE),
      spinglass_noise_merge_rate = mean(spinglass_noise_merge, na.rm = TRUE),
      spinglass_forced_assignment_rate = mean(spinglass_forced_assignment, na.rm = TRUE),
      spinglass_mean_ari = mean(spinglass_ari, na.rm = TRUE),
      spinglass_mean_time_sec = mean(spinglass_time_sec, na.rm = TRUE),
      
      .groups = "drop"
    )
}
