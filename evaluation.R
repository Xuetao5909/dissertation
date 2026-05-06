library(dplyr)

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
  
  strict_success <- g1_same && g2_same && separate && !noise_merge
  
  list(
    success = strict_success,
    false_merge = length(unique(mem[c(group1, group2)])) == 1,
    over_split = !(g1_same && g2_same),
    noise_merge = noise_merge,           # Indicates whether background/noise nodes have been absorbed into either target cluster
    forced_assignment = FALSE            # Included for completeness; set to FALSE outside the bowtie scenario
  )
}

# ----------------------------------------
# Hard-clustering reference criterion for the bowtie scenario
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
    noise_merge = NA,                    # A clean definition of noise merge is not straightforward in the bowtie setting; recorded as NA
    forced_assignment = forced_assignment # Indicates whether the hub structure has been reduced to a single hard assignment
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
    false_merge = FALSE
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
      recovered_right = 0
    ))
  }
  
  max_left <- max(sapply(comms, function(x) sum(left %in% x)))
  max_right <- max(sapply(comms, function(x) sum(right %in% x)))
  hub_overlap <- sum(sapply(comms, function(x) hub %in% x)) >= 2
  
  non_hub_left <- setdiff(left, hub)
  non_hub_right <- setdiff(right, hub)
  
  false_merge <- any(sapply(comms, function(x) {
    sum(non_hub_left %in% x) >= 2 && sum(non_hub_right %in% x) >= 2
  }))
  
  list(
    success = (max_left >= 2) && (max_right >= 2) && hub_overlap && !false_merge,
    false_merge = false_merge,
    hub_overlap = hub_overlap,
    recovered_left = max_left,
    recovered_right = max_right
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
  
  list(
    success = has_g1 && has_g2 && !false_merge,
    false_merge = false_merge
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
  
  list(
    success = has_g1 && has_g2 && !false_merge,
    false_merge = false_merge
  )
}

# ----------------------------------------
# Summary table
# The third method is Spinglass rather than SBM
# ----------------------------------------
summarise_results <- function(res) {
  res %>%
    group_by(scenario) %>%
    summarise(
      louvain_success_rate = mean(louvain_success, na.rm = TRUE),
      louvain_false_merge_rate = mean(louvain_false_merge, na.rm = TRUE),
      louvain_over_split_rate = mean(louvain_over_split, na.rm = TRUE),
      louvain_noise_merge_rate = mean(louvain_noise_merge, na.rm = TRUE),
      louvain_forced_assignment_rate = mean(louvain_forced_assignment, na.rm = TRUE),
      
      cpm_success_rate = mean(cpm_success, na.rm = TRUE),
      cpm_false_merge_rate = mean(cpm_false_merge, na.rm = TRUE),
      cpm_mean_tau = mean(cpm_tau, na.rm = TRUE),
      
      spinglass_success_rate = mean(spinglass_success, na.rm = TRUE),
      spinglass_false_merge_rate = mean(spinglass_false_merge, na.rm = TRUE),
      spinglass_over_split_rate = mean(spinglass_over_split, na.rm = TRUE),
      
      # Corrected: these fields now refer to the Spinglass outputs,
      # rather than incorrectly reusing the Louvain variables
      spinglass_noise_merge_rate = mean(spinglass_noise_merge, na.rm = TRUE),
      spinglass_forced_assignment_rate = mean(spinglass_forced_assignment, na.rm = TRUE),
      .groups = "drop"
    )
}