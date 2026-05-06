sample_edge_from_reservoir <- function(df, n = 1) {
  idx <- sample(seq_len(nrow(df)), size = n, replace = TRUE)
  df[idx, c("H3", "H4", "W")]
}

make_empty_network <- function(n = 10) {
  list(
    H3 = matrix(0, n, n),
    H4 = matrix(0, n, n),
    W  = matrix(0, n, n)
  )
}

fill_symmetric_edge <- function(net, i, j, H3, H4) {
  net$H3[i, j] <- net$H3[j, i] <- H3
  net$H4[i, j] <- net$H4[j, i] <- H4
  net$W[i, j]  <- net$W[j, i]  <- H4 - H3
  net
}

fill_edges_from_pairs <- function(net, pairs, reservoir) {
  for (k in seq_len(ncol(pairs))) {
    i <- pairs[1, k]
    j <- pairs[2, k]
    s <- sample_edge_from_reservoir(reservoir, 1)
    net <- fill_symmetric_edge(net, i, j, s$H3, s$H4)
  }
  net
}

fill_all_remaining_edges <- function(net, reservoir) {
  n <- nrow(net$H4)
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      if (net$H4[i, j] == 0 && net$H3[i, j] == 0) {
        s <- sample_edge_from_reservoir(reservoir, 1)
        net <- fill_symmetric_edge(net, i, j, s$H3, s$H4)
      }
    }
  }
  net
}

fill_cross_group_edges <- function(net, group1, group2, reservoir) {
  for (i in group1) {
    for (j in group2) {
      s <- sample_edge_from_reservoir(reservoir, 1)
      net <- fill_symmetric_edge(net, i, j, s$H3, s$H4)
    }
  }
  net
}

finalise_network <- function(net) {
  diag(net$H3) <- 0
  diag(net$H4) <- 0
  diag(net$W)  <- 0
  net
}

generate_triangle_network <- function(
    n = 10,
    reservoir_shared,
    reservoir_distinct_bg,
    core = c(1, 2, 3)
) {
  net <- make_empty_network(n)
  
  # fill all edges with low/background distinct noise first
  net <- fill_all_remaining_edges(net, reservoir_distinct_bg)
  
  # overwrite core triangle with true shared edges
  core_pairs <- combn(core, 2)
  net <- fill_edges_from_pairs(net, core_pairs, reservoir_shared)
  
  finalise_network(net)
}

generate_bowtie_network <- function(
    n = 10,
    reservoir_shared,
    reservoir_distinct_bg,
    left = c(1, 2, 3),
    right = c(3, 4, 5)
) {
  net <- make_empty_network(n)
  
  # fill all edges with low/background distinct noise first
  net <- fill_all_remaining_edges(net, reservoir_distinct_bg)
  
  # overwrite left and right triangles with true shared edges
  net <- fill_edges_from_pairs(net, combn(left, 2), reservoir_shared)
  net <- fill_edges_from_pairs(net, combn(right, 2), reservoir_shared)
  
  finalise_network(net)
}

generate_ld_trap_network <- function(
    n = 10,
    reservoir_shared,
    reservoir_distinct_bg,
    reservoir_distinct_trap,
    group1 = c(1, 2, 3),
    group2 = c(4, 5, 6)
) {
  net <- make_empty_network(n)
  
  # 1. first fill everything with ordinary distinct/background edges
  net <- fill_all_remaining_edges(net, reservoir_distinct_bg)
  
  # 2. overwrite within-group edges with true shared edges
  net <- fill_edges_from_pairs(net, combn(group1, 2), reservoir_shared)
  net <- fill_edges_from_pairs(net, combn(group2, 2), reservoir_shared)
  
  # 3. overwrite cross-group edges with dangerous high-H4 distinct tail edges
  net <- fill_cross_group_edges(net, group1, group2, reservoir_distinct_trap)
  
  finalise_network(net)
}
generate_local_bridge_trap_network <- function(
    n = 10,
    reservoir_shared_weak,
    reservoir_distinct_bg,
    reservoir_distinct_trap,
    group1 = c(1, 2, 3, 4),
    group2 = c(5, 6, 7, 8),
    bridge_pairs = matrix(c(
      4, 5,
      4, 6,
      2, 6
    ), nrow = 2)
) {
  net <- make_empty_network(n)
  
  
  net <- fill_all_remaining_edges(net, reservoir_distinct_bg)
  
  
  net <- fill_edges_from_pairs(net, combn(group1, 2), reservoir_shared_weak)
  net <- fill_edges_from_pairs(net, combn(group2, 2), reservoir_shared_weak)
  
  
  net <- fill_edges_from_pairs(net, bridge_pairs, reservoir_distinct_trap)
  
  finalise_network(net)
}