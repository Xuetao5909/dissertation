library(igraph)

# ----------------------------------------
# Louvain on continuous H4 matrix
# ----------------------------------------
run_louvain <- function(H4_mat) {
  g <- graph_from_adjacency_matrix(
    H4_mat,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )
  
  clu <- cluster_louvain(g, weights = E(g)$weight)
  
  list(
    membership = membership(clu),
    objective = modularity(clu)
  )
}

# ----------------------------------------
# Build CPM communities from clique list
# overlap_required = 2 corresponds to strict CPM for k=3
# ----------------------------------------
build_cpm_communities <- function(clique_list, overlap_required = 2) {
  n_clq <- length(clique_list)
  
  # no clique
  if (n_clq == 0) {
    return(list())
  }
  
  # single clique
  if (n_clq == 1) {
    return(list(sort(unique(as.integer(clique_list[[1]])))))
  }
  
  adj <- matrix(0, n_clq, n_clq)
  
  for (i in seq_len(n_clq - 1)) {
    for (j in seq.int(i + 1, n_clq)) {
      ov <- length(intersect(clique_list[[i]], clique_list[[j]]))
      if (ov >= overlap_required) {
        adj[i, j] <- 1
        adj[j, i] <- 1
      }
    }
  }
  
  g_clq <- graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
  comps <- components(g_clq)
  
  communities <- vector("list", comps$no)
  for (k in seq_len(comps$no)) {
    idx <- which(comps$membership == k)
    communities[[k]] <- sort(unique(unlist(clique_list[idx])))
  }
  
  communities
}

# ----------------------------------------
# Temporary mapping:
# - unique-membership node -> corresponding community label
# - overlap node -> assigned a new temporary label
# - noise node -> assigned its own temporary label
# Used only for scoring / modularity optimisation
# ----------------------------------------
temporary_membership_from_cpm <- function(communities, n_nodes) {
  node_to_comms <- vector("list", n_nodes)
  
  for (k in seq_along(communities)) {
    for (v in communities[[k]]) {
      node_to_comms[[v]] <- c(node_to_comms[[v]], k)
    }
  }
  
  temp_mem <- integer(n_nodes)
  next_label <- length(communities) + 1
  
  for (v in seq_len(n_nodes)) {
    if (length(node_to_comms[[v]]) == 0) {
      # noise node
      temp_mem[v] <- next_label
      next_label <- next_label + 1
    } else if (length(node_to_comms[[v]]) == 1) {
      # unique community membership
      temp_mem[v] <- node_to_comms[[v]][1]
    } else {
      # overlap node: assign temporarily to its first CPM community
      # for modularity scoring only
      temp_mem[v] <- node_to_comms[[v]][1]
    }
  }
  
  temp_mem
}

# ----------------------------------------
# Dynamic CPM:
# 1. threshold H4 by tau
# 2. find true k-cliques
# 3. build CPM communities
# 4. use temporary membership to score by modularity
# 5. choose tau with best modularity
# ----------------------------------------
run_cpm_dynamic <- function(
    H4_mat,
    tau_grid = seq(0.1, 0.7, by = 0.05),
    clique_size = 3,
    overlap_required = 2
) {
  best <- NULL
  best_score <- -Inf
  
  score_trace <- data.frame(
    tau = numeric(0),
    objective = numeric(0),
    n_cliques = integer(0),
    n_communities = integer(0),
    stringsAsFactors = FALSE
  )
  
  for (tau in tau_grid) {
    # 1. threshold to binary graph
    A_bin <- (H4_mat >= tau) * 1
    diag(A_bin) <- 0
    
    g_bin <- graph_from_adjacency_matrix(
      A_bin,
      mode = "undirected",
      diag = FALSE
    )
    
    # 2. true cliques of size = clique_size
    clq <- cliques(g_bin, min = clique_size, max = clique_size)
    clq <- lapply(clq, as.integer)
    
    # 3. build communities
    communities <- build_cpm_communities(
      clique_list = clq,
      overlap_required = overlap_required
    )
    
    # 4. temporary membership for scoring
    temp_mem <- temporary_membership_from_cpm(
      communities = communities,
      n_nodes = nrow(H4_mat)
    )
    
    # 5. score on continuous H4 graph
    g_weighted <- graph_from_adjacency_matrix(
      H4_mat,
      mode = "undirected",
      weighted = TRUE,
      diag = FALSE
    )
    
    score <- modularity(
      g_weighted,
      membership = temp_mem,
      weights = E(g_weighted)$weight
    )
    
    score_trace <- rbind(
      score_trace,
      data.frame(
        tau = tau,
        objective = score,
        n_cliques = length(clq),
        n_communities = length(communities),
        stringsAsFactors = FALSE
      )
    )
    
    if (is.finite(score) && score > best_score) {
      best_score <- score
      best <- list(
        tau = tau,
        communities = communities,
        temp_membership = temp_mem,
        objective = score,
        n_cliques = length(clq),
        n_communities = length(communities),
        A_bin = A_bin
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
      n_communities = 0L,
      A_bin = matrix(0, nrow(H4_mat), ncol(H4_mat))
    )
  }
  
  best$score_trace <- score_trace
  best
}

# ----------------------------------------
# Signed Spinglass baseline on continuous signed matrix
# W = H4 - H3
#
# NOTE:
# - This is a signed graph clustering baseline
# - No thresholding is used here
# ----------------------------------------
run_signed_spinglass <- function(W_mat, seed = 42) {
  g_signed <- graph_from_adjacency_matrix(
    W_mat,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )
  
  set.seed(seed)
  res <- cluster_spinglass(g_signed, implementation = "neg")
  
  list(
    membership = membership(res),
    objective = modularity(res)
  )
}

