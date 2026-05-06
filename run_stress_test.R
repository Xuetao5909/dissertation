library(tibble)

run_stress_test <- function(
    n_reps = 100,
    scenarios = c("triangle", "bowtie", "ld_trap", "local_bridge_trap"),
    reservoir_shared,
    reservoir_shared_weak,
    reservoir_distinct_bg,
    reservoir_distinct_trap,
    reservoir_ambiguous = NULL,
    cpm_tau_grid = seq(0.1, 0.7, by = 0.05)
) {
  results <- list()
  counter <- 1
  
  for (sc in scenarios) {
    for (rep in seq_len(n_reps)) {
      
      net <- switch(
        sc,
        "triangle" = generate_triangle_network(
          reservoir_shared = reservoir_shared,
          reservoir_distinct_bg = reservoir_distinct_bg
        ),
        
        "bowtie" = generate_bowtie_network(
          reservoir_shared = reservoir_shared,
          reservoir_distinct_bg = reservoir_distinct_bg
        ),
        
        "ld_trap" = generate_ld_trap_network(
          reservoir_shared = reservoir_shared,
          reservoir_distinct_bg = reservoir_distinct_bg,
          reservoir_distinct_trap = reservoir_distinct_trap
        ),
        
        "local_bridge_trap" = generate_local_bridge_trap_network(
          reservoir_shared_weak = reservoir_shared_weak,
          reservoir_distinct_bg = reservoir_distinct_bg,
          reservoir_distinct_trap = reservoir_distinct_trap
        ),
        
        stop(sprintf("Unknown scenario: %s", sc))
      )
      
      # -------------------------
      # Run methods
      # -------------------------
      lou <- run_louvain(net$H4)
      
      cpm <- run_cpm_dynamic(
        H4_mat = net$H4,
        tau_grid = cpm_tau_grid,
        overlap_required = 2
      )
      
      spg <- run_signed_spinglass(net$W)
      
      # -------------------------
      # Evaluate by scenario
      # -------------------------
      if (sc == "triangle") {
        core_lou <- unique(lou$membership[c(1, 2, 3)])[1]
        noise_lou <- sum(lou$membership == core_lou) > 3
        eval_lou <- list(
          success = length(unique(lou$membership[c(1, 2, 3)])) == 1 && !noise_lou,
          false_merge = FALSE,
          over_split = length(unique(lou$membership[c(1, 2, 3)])) != 1,
          noise_merge = noise_lou,
          forced_assignment = FALSE
        )
        
        core_spg <- unique(spg$membership[c(1, 2, 3)])[1]
        noise_spg <- sum(spg$membership == core_spg) > 3
        eval_spg <- list(
          success = length(unique(spg$membership[c(1, 2, 3)])) == 1 && !noise_spg,
          false_merge = FALSE,
          over_split = length(unique(spg$membership[c(1, 2, 3)])) != 1,
          noise_merge = noise_spg,
          forced_assignment = FALSE
        )
        
        eval_cpm <- check_cpm_triangle(
          comms = cpm$communities,
          core = c(1, 2, 3)
        )
      }
      
      if (sc == "bowtie") {
        eval_lou <- check_hard_bowtie(lou$membership)
        eval_spg <- check_hard_bowtie(spg$membership)
        eval_cpm <- check_cpm_bowtie(cpm$communities)
      }
      
      if (sc == "ld_trap") {
        eval_lou <- check_two_cluster_success(
          mem = lou$membership,
          group1 = c(1, 2, 3),
          group2 = c(4, 5, 6)
        )
        
        eval_spg <- check_two_cluster_success(
          mem = spg$membership,
          group1 = c(1, 2, 3),
          group2 = c(4, 5, 6)
        )
        
        eval_cpm <- check_cpm_ldtrap(
          comms = cpm$communities,
          group1 = c(1, 2, 3),
          group2 = c(4, 5, 6)
        )
      }
      
      if (sc == "local_bridge_trap") {
        eval_lou <- check_two_cluster_success(
          mem = lou$membership,
          group1 = c(1, 2, 3, 4),
          group2 = c(5, 6, 7, 8)
        )
        
        eval_spg <- check_two_cluster_success(
          mem = spg$membership,
          group1 = c(1, 2, 3, 4),
          group2 = c(5, 6, 7, 8)
        )
        
        eval_cpm <- check_cpm_local_bridge_trap(
          comms = cpm$communities,
          group1 = c(1, 2, 3, 4),
          group2 = c(5, 6, 7, 8)
        )
      }
      
      # -------------------------
      # Save one replicate result
      # -------------------------
      results[[counter]] <- tibble(
        scenario = sc,
        replicate = rep,
        
        louvain_success = eval_lou$success,
        louvain_false_merge = eval_lou$false_merge,
        louvain_over_split = eval_lou$over_split,
        louvain_noise_merge = eval_lou$noise_merge,
        louvain_forced_assignment = eval_lou$forced_assignment,
        louvain_objective = lou$objective,
        
        cpm_success = eval_cpm$success,
        cpm_false_merge = ifelse(is.null(eval_cpm$false_merge), NA, eval_cpm$false_merge),
        cpm_tau = cpm$tau,
        cpm_objective = cpm$objective,
        cpm_n_cliques = cpm$n_cliques,
        cpm_n_communities = cpm$n_communities,
        cpm_hub_overlap = ifelse(is.null(eval_cpm$hub_overlap), NA, eval_cpm$hub_overlap),
        
        spinglass_success = eval_spg$success,
        spinglass_false_merge = eval_spg$false_merge,
        spinglass_over_split = eval_spg$over_split,
        spinglass_noise_merge = eval_spg$noise_merge,
        spinglass_forced_assignment = eval_spg$forced_assignment,
        spinglass_objective = spg$objective
      )
      
      counter <- counter + 1
    }
  }
  
  dplyr::bind_rows(results)
}